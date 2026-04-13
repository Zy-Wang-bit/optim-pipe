#!/bin/bash
# 批量轨迹分析: 30 个 MD 轨迹 (15 变体 × 2 pH)
# 在 Slurm 计算节点上运行 (需 optim-pipe 环境)
set -euo pipefail

REPO="/public/home/ziyang/code/optim-pipe"
MD_DIR="${REPO}/experiments/sdab_variants/md"
CONFIG="${REPO}/third_party/molecular_dynamics/configs/md_config_sdab.yaml"
MD_CODE="${REPO}/third_party/molecular_dynamics"

cd "$MD_CODE"

# ---- Phase 1: 轨迹分析 ----
echo "=== Phase 1: Trajectory Analysis (30 sims) ==="
for id in $(seq -w 0 14); do
    variant="sdab_v2_00${id}"
    for ph in 7.4 6.0; do
        traj_dir="${MD_DIR}/${variant}/pH_${ph}"
        if [[ ! -f "${traj_dir}/production.xtc" ]]; then
            echo "SKIP: ${variant}/pH_${ph} — no trajectory"
            continue
        fi
        if [[ -f "${traj_dir}/analysis/summary.json" ]]; then
            echo "SKIP: ${variant}/pH_${ph} — already analyzed"
            continue
        fi
        echo "[$(date '+%H:%M:%S')] Analyzing ${variant} pH=${ph}..."
        python analyze_trajectory.py \
            --traj "$traj_dir" \
            --config "$CONFIG" 2>&1 | tail -3
    done
done

# ---- Phase 2: 双 pH 对比 ----
echo ""
echo "=== Phase 2: pH Comparison (15 variants) ==="
for id in $(seq -w 0 14); do
    variant="sdab_v2_00${id}"
    variant_dir="${MD_DIR}/${variant}"
    if [[ -f "${variant_dir}/ph_comparison.json" ]]; then
        echo "SKIP: ${variant} — already compared"
        continue
    fi
    echo "[$(date '+%H:%M:%S')] Comparing pH for ${variant}..."
    python compare_ph.py \
        --variant-dir "$variant_dir" \
        --base-ph 7.4 --target-ph 6.0 2>&1 | tail -2
done

echo ""
echo "=== All analysis complete at $(date) ==="
