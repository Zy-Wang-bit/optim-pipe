#!/bin/bash
# Migrate a running GROMACS production MD to a different GPU using checkpoint restart.
# Intended use: during sdab R4 multi-sim run, slow sims on one GPU can be relocated
# to a freed GPU after other sims on that GPU finish.
#
# Usage:
#   migrate_sim.sh <variant> <ph> <new_gpu>
# Example:
#   migrate_sim.sh r4_06 7.4 6
#
# What it does:
#   1. Finds the gmx mdrun PID running for <variant>/pH_<ph>
#   2. Sends SIGINT -> gmx flushes production.cpt and exits cleanly
#   3. Waits for child and parent python run_md.py to exit
#   4. Relaunches gmx mdrun with -cpi production.cpt -append on <new_gpu>,
#      detached via nohup, appending to the same per-sim log
set -euo pipefail

VARIANT=${1:?need variant, e.g. r4_06}
PH=${2:?need pH, e.g. 7.4}
NEW_GPU=${3:?need target GPU index 0-7}

ROOT=/root/sdab_r4_md
WORKDIR=$ROOT/runs/$VARIANT/pH_$PH
LOG=$ROOT/logs/${VARIANT}_pH${PH}.log

if [ ! -d "$WORKDIR" ]; then
    echo "[err] workdir not found: $WORKDIR" >&2
    exit 1
fi
if [ ! -f "$WORKDIR/production.cpt" ]; then
    echo "[err] no production.cpt yet in $WORKDIR — sim not in production phase" >&2
    exit 1
fi

cd "$WORKDIR"

MDRUN_PID=$(pgrep -f "gmx mdrun.*production" | while read p; do
    if readlink -f /proc/$p/cwd 2>/dev/null | grep -q "$WORKDIR$"; then
        echo $p; break
    fi
done)

if [ -z "$MDRUN_PID" ]; then
    echo "[warn] no running gmx mdrun for $VARIANT pH $PH — assuming already stopped"
else
    echo "[migrate] SIGINT -> gmx mdrun pid=$MDRUN_PID"
    kill -INT "$MDRUN_PID"
    for i in $(seq 1 60); do
        if ! kill -0 "$MDRUN_PID" 2>/dev/null; then break; fi
        sleep 2
    done
    if kill -0 "$MDRUN_PID" 2>/dev/null; then
        echo "[err] gmx mdrun still alive after 120s" >&2
        exit 2
    fi
    echo "[migrate] gmx mdrun exited cleanly"
fi

PYPID=$(pgrep -f "run_md.py.*${VARIANT}\.pdb.*--ph ${PH}" || true)
if [ -n "$PYPID" ]; then
    echo "[migrate] also killing parent python run_md.py pid=$PYPID"
    kill -TERM $PYPID 2>/dev/null || true
    sleep 2
fi

echo "[migrate] restarting on GPU $NEW_GPU ..."
source /etc/profile.d/modules.sh 2>/dev/null || true
module load gromacs/2024.2 openmpi/5.0.3 2>/dev/null || true

{
    echo "# === migrated to GPU $NEW_GPU at $(date -Iseconds) ==="
} >> "$LOG"

CUDA_VISIBLE_DEVICES=$NEW_GPU nohup gmx mdrun \
    -s production.tpr \
    -deffnm production \
    -cpi production.cpt \
    -append \
    -ntmpi 1 -ntomp 5 \
    >> "$LOG" 2>&1 &

NEW_PID=$!
disown $NEW_PID 2>/dev/null || true
echo "[migrate] new gmx mdrun pid=$NEW_PID on GPU $NEW_GPU"
echo "[migrate] log: $LOG"
