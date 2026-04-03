#!/usr/bin/env bash
# 一键顺序跑：MPNN → ESM-1b → FoldX（修复 / WT / 批次 / 评估）→ 合并筛选
# 用法：
#   bash run_pipeline.sh [CONFIG]             # 正常跑
#   CLEAN=1 bash run_pipeline.sh [CONFIG]     # 先清理旧产物再跑
#
# 说明：
# - 默认使用 configs/config.yaml
# - 需要修改 MPNN_ARGS 里参数时，直接改下面那一行

set -euo pipefail

CFG="${1:-configs/config.yaml}"
TS="$(date +%Y%m%d_%H%M%S)"
LOGDIR="logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/pipeline_${TS}.log"

# —— 你可以在这里改 MPNN 的参数（只是占位，正式开跑前手动改）——
MPNN_ARGS=(
  --num-per-pdb 20000         # 每个PDB的设计数（示例）
  --shards 5                # 分片数（示例）
  --temps 0.10,0.15,0.20,0.25,0.30   # 采样温度（示例）
  --seed 4242               # 随机种子
  --build-his-seeds         # 同步生成 His 种子库
)

run(){ echo -e "\n[RUN] $*" | tee -a "$LOG"; "$@" 2>&1 | tee -a "$LOG"; }
note(){ echo -e "\n[NOTE] $*" | tee -a "$LOG"; }
ok(){ echo -e "[OK] $*" | tee -a "$LOG"; }

# 0) 可选清理
if [[ "${CLEAN:-0}" == "1" ]]; then
  note "清理旧产物（保留 experiments/1E62/data 和 experiments/1E62/data）…"
  rm -rf mpnn_outputs/* his_seeds/* esm_scores/* \
         foldx/repaired/* foldx/batches/* \
         results/screening/* results/final_top10k.csv \
         logs/processed_batches.json \
         configs/mpnn/chain_id.jsonl configs/mpnn/fixed_positions.jsonl || true
  ok "清理完成。"
fi

# 基础检查 & 目录
[[ -f "$CFG" ]] || { echo "[ERR] 找不到配置：$CFG" | tee -a "$LOG"; exit 1; }
mkdir -p results/screening foldx/repaired foldx/batches esm_scores mpnn_outputs his_seeds configs/mpnn

# 1) 接口热点扫描（生成 results/screening/his_hotspots.csv）
note "Step 1: 接口热点扫描（His 偏向位点识别）"
run python scripts/scan_interface.py "$CFG"
ok "his_hotspots.csv 就绪。"

# 2) MPNN 设计（按 config.design.region；并生成 his seeds）
note "Step 2: MPNN 设计（生成设计库 + His 种子）"
run python scripts/run_mpnn_design.py "$CFG" "${MPNN_ARGS[@]}"
ok "MPNN 输出与 His seeds 就绪。"

# 3) ESM-1b 评分与合并（自动读取 mpnn + his_seed，生成 results/screening/for_foldx.csv）
note "Step 3: ESM-1b 评分与合并（for_foldx.csv）"
run python scripts/run_esm_chunk.py "$CFG"
ok "for_foldx.csv 就绪。"

# 4) FoldX：RepairPDB（每个PDB一次）
note "Step 4: FoldX RepairPDB"
run python scripts/repair_pdbs.py "$CFG"
ok "foldx/repaired/*_Repair.pdb 就绪。"

# 5) 预计算 WT 的 AC（pH=7.4/6.0），缓存到 foldx/repaired/WT_ac.csv
note "Step 5: 预计算 WT 相互作用能"
run python scripts/precompute_wt_ac.py "$CFG"
ok "WT_ac.csv 就绪。"

# 6) 生成 FoldX 批次（individual_list.txt + batch_seqs.csv）
note "Step 6: 生成 FoldX 批次（individual_list.txt）"
run python scripts/make_mutlist_chunk.py "$CFG"
ok "foldx/batches/<pdb>/batch_XXXXX 就绪。"

# 7) FoldX 评估（BuildModel + AnalyseComplex，两点 pH）
note "Step 7: FoldX 评估（两点 pH）"
run python scripts/run_foldx_batch.py "$CFG"
ok "各批次 foldx_summary.csv 就绪。"

# 8) 合并 & 自适应筛选（输出 results/final_top10k.csv，并包含 WT 对比）
note "Step 8: 合并与自适应筛选（含 WT 对比）"
run python scripts/merge_and_select.py "$CFG"
ok "最终结果 -> results/final_top10k.csv"

note "全部完成 ✅ 详细日志：$LOG"