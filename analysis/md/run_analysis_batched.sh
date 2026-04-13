#!/bin/bash
# 分批分析: 5 路并行, 避免 IO 瓶颈
set -uo pipefail

CD="/public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics"
CONFIG="/public/home/ziyang/code/optim-pipe/third_party/molecular_dynamics/configs/md_config_sdab.yaml"
MD_DIR="/public/home/ziyang/code/optim-pipe/experiments/sdab_variants/md"
PARALLEL=5

cd "$CD"

echo "=== Batched analysis ($PARALLEL parallel) ==="
echo "Start: $(date)"

# Build job list
JOBS=()
for id in $(seq -w 0 14); do
    for ph in 7.4 6.0; do
        traj="${MD_DIR}/sdab_v2_00${id}/pH_${ph}"
        if [[ -f "${traj}/analysis/summary.json" ]]; then
            continue  # skip done
        fi
        JOBS+=("$traj")
    done
done

echo "Jobs: ${#JOBS[@]}"

# Process in batches
i=0
while (( i < ${#JOBS[@]} )); do
    PIDS=()
    for (( j=0; j<PARALLEL && i<${#JOBS[@]}; j++, i++ )); do
        traj="${JOBS[$i]}"
        echo "[$(date '+%H:%M:%S')] START $(basename $(dirname $traj))/$(basename $traj)"
        python analyze_trajectory.py --traj "$traj" --config "$CONFIG" > /dev/null 2>&1 &
        PIDS+=($!)
    done
    for pid in "${PIDS[@]}"; do
        wait "$pid"
    done
    echo "[$(date '+%H:%M:%S')] Batch done ($(find ${MD_DIR} -name summary.json | wc -l)/30)"
done

echo ""
echo "=== Analysis done at $(date) ==="
echo "summary.json: $(find ${MD_DIR} -name summary.json | wc -l)/30"

# compare_ph
echo "=== compare_ph ==="
for id in $(seq -w 0 14); do
    python compare_ph.py --variant-dir "${MD_DIR}/sdab_v2_00${id}" --base-ph 7.4 --target-ph 6.0 > /dev/null 2>&1
done
echo "ph_comparison: $(find ${MD_DIR} -name ph_comparison.json | wc -l)/15"

echo "=== ALL DONE at $(date) ==="
