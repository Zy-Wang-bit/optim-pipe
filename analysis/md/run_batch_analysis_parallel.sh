#!/bin/bash
# 批量轨迹分析 — 30 路并行 (computer1, 128 cores)
set -uo pipefail

REPO="/public/home/ziyang/code/optim-pipe"
MD_DIR="${REPO}/experiments/sdab_variants/md"
CONFIG="${REPO}/third_party/molecular_dynamics/configs/md_config_sdab.yaml"
MD_CODE="${REPO}/third_party/molecular_dynamics"
LOG_DIR="${MD_DIR}/analysis_logs"
mkdir -p "$LOG_DIR"

cd "$MD_CODE"

echo "=== Phase 1: Parallel Trajectory Analysis (30 sims, $(nproc) cores) ==="
echo "Start: $(date)"

PIDS=()
for id in $(seq -w 0 14); do
    variant="sdab_v2_00${id}"
    for ph in 7.4 6.0; do
        traj_dir="${MD_DIR}/${variant}/pH_${ph}"
        log="${LOG_DIR}/${variant}_pH${ph}.log"

        if [[ ! -f "${traj_dir}/production.xtc" ]]; then
            echo "SKIP: ${variant}/pH_${ph} — no trajectory"
            continue
        fi
        if [[ -f "${traj_dir}/analysis/summary.json" ]]; then
            echo "SKIP: ${variant}/pH_${ph} — already done"
            continue
        fi

        echo "LAUNCH: ${variant} pH=${ph}"
        python analyze_trajectory.py \
            --traj "$traj_dir" \
            --config "$CONFIG" > "$log" 2>&1 &
        PIDS+=($!)
    done
done

echo ""
echo "Launched ${#PIDS[@]} analysis jobs. Waiting..."
echo ""

# Wait and report
FAIL=0
for pid in "${PIDS[@]}"; do
    if ! wait "$pid"; then
        FAIL=$((FAIL+1))
    fi
done

echo "=== Phase 1 complete: $((${#PIDS[@]}-FAIL))/${#PIDS[@]} succeeded ==="
if (( FAIL > 0 )); then
    echo "Failed logs:"
    grep -l "Error\|Traceback" "$LOG_DIR"/*.log 2>/dev/null
fi

# ---- Phase 2: pH comparison (sequential, fast) ----
echo ""
echo "=== Phase 2: pH Comparison ==="
for id in $(seq -w 0 14); do
    variant="sdab_v2_00${id}"
    variant_dir="${MD_DIR}/${variant}"
    if [[ -f "${variant_dir}/ph_comparison.json" ]]; then
        echo "SKIP: ${variant}"
        continue
    fi
    echo "Compare: ${variant}"
    python compare_ph.py \
        --variant-dir "$variant_dir" \
        --base-ph 7.4 --target-ph 6.0 2>&1 | tail -1
done

echo ""
echo "=== ALL DONE at $(date) ==="
echo "Summaries: find ${MD_DIR} -name summary.json | wc -l"
echo "Comparisons: find ${MD_DIR} -name ph_comparison.json | wc -l"
