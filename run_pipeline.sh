#!/usr/bin/env bash
# 一键顺序跑：3 Tier / 13 Step 漏斗 pipeline
# 用法：
#   bash run_pipeline.sh [CONFIG]             # 正常跑
#   CLEAN=1 bash run_pipeline.sh [CONFIG]     # 先清理旧产物再跑
#
# 模式判断：
#   tier1.enabled: true  → 3 Tier 流程（新）
#   tier1.enabled: false → 旧 4 Phase 流程（merge_and_select.py）

set -euo pipefail

CFG="${1:-configs/config.yaml}"
TS="$(date +%Y%m%d_%H%M%S)"
LOGDIR="logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/pipeline_${TS}.log"

# —— MPNN 参数（按需修改）——
MPNN_ARGS=(
  --num-per-pdb 20000
  --shards 5
  --temps 0.10,0.15,0.20,0.25,0.30
  --seed 4242
  --build-his-seeds
)

run(){ echo -e "\n[RUN] $*" | tee -a "$LOG"; "$@" 2>&1 | tee -a "$LOG"; }
note(){ echo -e "\n[NOTE] $*" | tee -a "$LOG"; }
ok(){ echo -e "[OK] $*" | tee -a "$LOG"; }

# 0) 可选清理
if [[ "${CLEAN:-0}" == "1" ]]; then
  note "清理旧产物..."
  rm -rf mpnn_outputs/* his_seeds/* esm_scores/* \
         foldx/repaired/* foldx/batches/* \
         results/screening/* results/final_top10k.csv results/merged_all.csv \
         results/tier1_candidates.csv results/tier2_candidates.csv results/final_candidates.csv \
         tier2/structures/* tier2/pka/* tier2/rosetta/* tier2/rmsd/* \
         phase_c/structures/* phase_c/pka/* phase_c/rosetta/* phase_c/rmsd/* \
         logs/processed_batches.json \
         configs/mpnn/chain_id.jsonl configs/mpnn/fixed_positions.jsonl || true
  ok "清理完成。"
fi

# 基础检查
[[ -f "$CFG" ]] || { echo "[ERR] 找不到配置：$CFG" | tee -a "$LOG"; exit 1; }
mkdir -p results/screening results/audit foldx/repaired foldx/batches \
         esm_scores mpnn_outputs his_seeds configs/mpnn \
         tier2/structures tier2/pka tier2/rosetta tier2/rmsd

# 读取模式
TIER1_ENABLED=$(python3 -c "
import yaml; c=yaml.safe_load(open('${CFG}'))
print(str(c.get('tier1',{}).get('enabled',False)).lower())
")
TIER2_ENABLED=$(python3 -c "
import yaml; c=yaml.safe_load(open('${CFG}'))
print(str(c.get('tier2',{}).get('enabled',False)).lower())
")

# ═══════════════════════════════════════════════════════════════════════════════
# Tier 1: 生成 + 高通量筛选（Steps 1-8）
# ═══════════════════════════════════════════════════════════════════════════════

note "Step 1: 接口热点扫描"
run python scripts/scan_interface.py "$CFG"
ok "his_hotspots.csv 就绪。"

note "Step 2: MPNN 设计"
run python scripts/run_mpnn_design.py "$CFG" "${MPNN_ARGS[@]}"
ok "MPNN 输出就绪。"

note "Step 3: ESM 评分"
run python scripts/run_esm_chunk.py "$CFG"
ok "for_foldx.csv 就绪。"

note "Step 4: FoldX RepairPDB"
run python scripts/repair_pdbs.py "$CFG"
ok "Repaired PDB 就绪。"

note "Step 5: 预计算 WT 相互作用能"
run python scripts/precompute_wt_ac.py "$CFG"
ok "WT_ac.csv 就绪。"

note "Step 6: 生成 FoldX 批次"
run python scripts/make_mutlist_chunk.py "$CFG"
ok "FoldX 批次就绪。"

note "Step 7: FoldX 评估"
run python scripts/run_foldx_batch.py "$CFG"
ok "foldx_summary.csv 就绪。"

# FoldX 中间文件清理（保留 summary CSV 和 batch_seqs.csv）
note "清理 FoldX 中间文件..."
find foldx/batches -name "*.pdb" -not -name "*_Repair.pdb" -delete 2>/dev/null || true
find foldx/batches -name "*.fxout" -delete 2>/dev/null || true
ok "FoldX 中间文件已清理。"

if [[ "$TIER1_ENABLED" == "true" ]]; then
    # ── 3 Tier 模式 ─────────────────────────────────────────────────────────
    note "Step 8: Tier 1 自适应筛选"
    run python scripts/tier1_filter.py "$CFG"
    ok "tier1_candidates.csv 就绪。"

    if [[ "$TIER2_ENABLED" == "true" ]]; then
        # ═════════════════════════════════════════════════════════════════════
        # Tier 2: 结构评估 + 筛选（Steps 9-12）
        # ═════════════════════════════════════════════════════════════════════

        # —— PyRosetta 线 ——
        note "Step 9a: PyRosetta 突变体建模"
        run conda run -n pyrosetta python scripts/build_structures.py "$CFG"
        ok "突变体结构就绪。"

        note "Step 9b: pKa 预测"
        run python scripts/run_pka.py "$CFG"
        ok "pKa 预测就绪。"

        note "Step 9c: Rosetta 评分 (dddG_elec + pH-score)"
        run conda run -n pyrosetta python scripts/run_rosetta_eval.py "$CFG"
        ok "Rosetta 评分就绪。"

        # —— SimpleFold 3x 线 ——
        note "Step 10: SimpleFold 3x 采样"
        run conda run -n simplefold python scripts/run_simplefold_3x.py "$CFG"
        ok "SimpleFold 3x 就绪。"

        note "Step 11: CDR RMSD"
        run python scripts/run_rmsd.py "$CFG"
        ok "CDR RMSD 就绪。"

        # —— Tier 2 筛选 ——
        note "Step 12: Tier 2 筛选 (pKa 相对排序 + RMSD 门槛)"
        run python scripts/tier2_filter.py "$CFG"
        ok "tier2_candidates.csv 就绪。"

        # ═════════════════════════════════════════════════════════════════════
        # Tier 3: 精排（Step 13）
        # ═════════════════════════════════════════════════════════════════════

        note "Step 13: Tier 3 精排 (dddG_elec 排序 + 软标记)"
        run python scripts/merge_and_rank.py "$CFG"
        ok "final_candidates.csv 就绪。"
    else
        note "Tier 2 已禁用，跳过 Steps 9-13。"
        note "Tier 1 结果: results/tier1_candidates.csv"
    fi
else
    # ── 旧模式 (4 Phase / 12 Step) ──────────────────────────────────────────
    PHASE_C=$(python3 -c "
    import yaml; c=yaml.safe_load(open('${CFG}'))
    print(str(c.get('phase_c',{}).get('enabled',False)).lower())
    ")

    if [[ "$PHASE_C" == "true" ]]; then
        note "Step 7.5: 预合并（生成 Phase C 候选列表）"
        run python scripts/merge_and_select.py "$CFG"
        ok "Phase C 候选列表就绪。"

        note "Step 8: Phase C - 结构生成"
        run conda run -n pyrosetta python scripts/build_structures.py "$CFG"
        ok "突变体结构就绪。"

        note "Step 9: Phase C - pKa 预测"
        run python scripts/run_pka.py "$CFG"
        ok "pKa 预测就绪。"

        note "Step 10: Phase C - Rosetta 评分"
        run conda run -n pyrosetta python scripts/run_rosetta_eval.py "$CFG"
        ok "Rosetta 评分就绪。"

        note "Step 11: Phase C - CDR RMSD"
        run python scripts/run_rmsd.py "$CFG"
        ok "CDR RMSD 就绪。"
    else
        note "Phase C 已禁用，跳过 Steps 8-11。"
    fi

    note "Step 12: 合并与筛选"
    run python scripts/merge_and_select.py "$CFG"
    ok "最终结果 -> results/final_top10k.csv"
fi

note "全部完成。详细日志：$LOG"
