# AF3 1E62 CeS Complex Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate five AlphaFold 3 structures for the 1E62 heavy/light chain antibody complexed with the HBsAg CeS sequence supplied by the user.

**Architecture:** Run the installed AlphaFold 3 Apptainer image directly on the already available `gpu6` compute node. Split AF3 into a CPU data-pipeline step and a GPU inference step, both launched as direct shell commands from `gpu6`; do not submit Slurm jobs from this plan.

**Tech Stack:** AlphaFold 3 `3.0.1` via Apptainer, direct execution on `gpu6`, one visible NVIDIA A40 GPU, Python standard library for JSON/result validation.

---

## File Structure

- Create `experiments/1E62/R4/af3_ces_wt/input/1e62_wt_hbsag_ces.json`: AlphaFold 3 input JSON with the exact heavy chain, light chain, and CeS antigen sequences.
- Create `experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_data_pipeline_gpu6.sh`: direct-run script for the CPU-only AF3 data pipeline on `gpu6`.
- Create `experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_infer_gpu6.sh`: direct-run script for AF3 GPU inference on `gpu6`.
- Produce `experiments/1E62/R4/af3_ces_wt/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json`: AF3 data-pipeline output containing MSA/template features.
- Produce `experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces/seed-*_sample-0/model.cif`: five AF3 complex structures, one per randomized seed.
- Produce `experiments/1E62/R4/af3_ces_wt/summary/confidence_summary.tsv`: compact confidence summary for deciding whether the structures are usable.
- Produce `experiments/1E62/R4/af3_ces_wt/summary/ranked_structures.tsv`: manifest mapping rank, seed, sample, score, source path, and copied path.
- Copy the five ranked structures to `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank1.cif` through `AF3-1E62-CeS-WT-rank5.cif`.
- Copy rank 1 to `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT.cif` as the top-ranked convenience alias.

## Assumptions and Guardrails

- Do not run `sbatch`, `srun`, or any Slurm submission command for this task.
- Run all compute commands only from `gpu6`, because the user has confirmed this node has enough resources.
- The current direct-run resource check is `hostname=gpu6` and one visible `NVIDIA A40` with `46068 MiB` total memory.
- The AF3 install remains at `/public/apps/alphafold3/alphafold3`.
- The AF3 model directory remains `/public/shared/alphafold3/models`.
- The fast local AF3 database path remains `/ssd/alphafold3`; the fallback AF3 database path remains `/public/shared/alphafold3`.
- The job uses five randomized seeds, `[75801088, 635954951, 495296019, 1129861070, 2026110668]`, with one diffusion sample per seed. This returns exactly five structures while still sampling five different random initializations.
- This first run predicts the full-length supplied HBsAg CeS sequence. If confidence shows the full-length membrane regions dominate or destabilize the complex, create a separate antigenic-loop/truncation run rather than editing this run in place.

---

### Task 1: Create AF3 Workspace and Input JSON

**Files:**
- Create: `experiments/1E62/R4/af3_ces_wt/input/1e62_wt_hbsag_ces.json`

- [ ] **Step 1: Create the AF3 workspace directories**

Run:

```bash
mkdir -p experiments/1E62/R4/af3_ces_wt/input
mkdir -p experiments/1E62/R4/af3_ces_wt/scripts
mkdir -p experiments/1E62/R4/af3_ces_wt/logs
mkdir -p experiments/1E62/R4/af3_ces_wt/summary
```

Expected: all four directories exist and `mkdir` prints no output.

- [ ] **Step 2: Write the AlphaFold 3 input JSON**

Create `experiments/1E62/R4/af3_ces_wt/input/1e62_wt_hbsag_ces.json` with this exact content:

```json
{
  "name": "1e62_wt_hbsag_ces",
  "modelSeeds": [75801088, 635954951, 495296019, 1129861070, 2026110668],
  "sequences": [
    {
      "protein": {
        "id": "H",
        "sequence": "EMQLVESGGGLVQPGGSLRLSCAASGFTFSDYWMNWVRQSPGKGLEWVAEVKLKSNNYATHYAESVKGRFTISRDDSKNSVYLQMNSLKTEDTAVYYCASGFDYWGQGTLVTVSS"
      }
    },
    {
      "protein": {
        "id": "L",
        "sequence": "DIVMSQSPDSLAVSLGERATITCKSSQSLLYSSDQNNYLAWSQQKPGQPPKLLIYWASIRDSGVPDRFSGSGTDFTLTISSLQAEDVAVYYCHQYYTYPFTFGQGTKLEIK"
      }
    },
    {
      "protein": {
        "id": "A",
        "sequence": "MENTTSGFLGPLLVLQAGFFLLTRILTIPQSLDSWWTSLNFLGGAPTCPGQNSQSPTSNHSPTSCPPICPGYRWMCLRRFIIFLFILLLCLIFLLVLLDYQGMLPVCPLLPGTSTTSTGPCKTCTIPAQGTSMFPSCCCTKPSDGNCTCIPIPSSWAFARFLWEWASVRFSWLSLLVPFVQWFVGLSPTVWLSVIWMMWYWGPSLYNILSPFLPLLPIFFCLWVYI"
      }
    }
  ],
  "dialect": "alphafold3",
  "version": 1
}
```

- [ ] **Step 3: Validate JSON syntax**

Run:

```bash
python -m json.tool experiments/1E62/R4/af3_ces_wt/input/1e62_wt_hbsag_ces.json >/tmp/1e62_wt_hbsag_ces.validated.json
```

Expected: command exits with status `0` and prints no error.

- [ ] **Step 4: Validate chain IDs, sequence lengths, and seeds**

Run:

```bash
python - <<'PY'
import json
from pathlib import Path

path = Path("experiments/1E62/R4/af3_ces_wt/input/1e62_wt_hbsag_ces.json")
data = json.loads(path.read_text())
chains = []
for item in data["sequences"]:
    protein = item["protein"]
    chains.append((protein["id"], len(protein["sequence"])))
print(chains)
print(data["modelSeeds"])
assert chains == [("H", 115), ("L", 111), ("A", 226)]
assert data["modelSeeds"] == [75801088, 635954951, 495296019, 1129861070, 2026110668]
assert data["dialect"] == "alphafold3"
assert data["version"] == 1
PY
```

Expected output:

```text
[('H', 115), ('L', 111), ('A', 226)]
[75801088, 635954951, 495296019, 1129861070, 2026110668]
```

- [ ] **Step 5: Record local worktree status**

Run:

```bash
git status --short
```

Expected: runtime files under `experiments/` may appear only if they are not ignored locally; do not stage runtime inputs or outputs unless the user explicitly asks.

---

### Task 2: Verify Direct Compute Context on gpu6

**Files:**
- No files created or modified.

- [ ] **Step 1: Confirm this shell is on gpu6**

Run:

```bash
hostname
```

Expected output:

```text
gpu6
```

- [ ] **Step 2: Confirm the visible GPU is idle enough for AF3 inference**

Run:

```bash
nvidia-smi --query-gpu=index,name,memory.total,memory.used,utilization.gpu --format=csv,noheader
```

Expected output should show one visible `NVIDIA A40`; the pre-run check observed:

```text
0, NVIDIA A40, 46068 MiB, 3 MiB, 0 %
```

- [ ] **Step 3: Confirm AF3 paths are visible from gpu6**

Run:

```bash
test -s /public/apps/alphafold3/alphafold3/alphafold3.0.1.sif
test -s /public/apps/alphafold3/alphafold3/run_alphafold.py
test -s /public/shared/alphafold3/models/af3.bin
test -d /ssd/alphafold3
test -d /public/shared/alphafold3
```

Expected: all five `test` commands exit with status `0` and print no output.

---

### Task 3: Run AF3 Data Pipeline Directly on gpu6

**Files:**
- Create: `experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_data_pipeline_gpu6.sh`
- Produce: `experiments/1E62/R4/af3_ces_wt/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json`
- Produce: `experiments/1E62/R4/af3_ces_wt/logs/af3_ces_data_pipeline_gpu6.log`

- [ ] **Step 1: Write the direct data-pipeline script**

Create `experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_data_pipeline_gpu6.sh` with this exact content:

```bash
#!/bin/bash
set -euo pipefail

if [ "$(hostname)" != "gpu6" ]; then
  echo "ERROR: this script must run directly on gpu6; current host is $(hostname)" >&2
  exit 2
fi

projectDir=/public/home/ziyang/code/optim-pipe/experiments/1E62/R4/af3_ces_wt
af3Bin=/public/apps/alphafold3/alphafold3
af3Database=/public/shared/alphafold3
af3Ssd=/ssd/alphafold3
logFile="${projectDir}/logs/af3_ces_data_pipeline_gpu6.log"

if [ -d "${projectDir}/msa/1e62_wt_hbsag_ces" ]; then
  echo "ERROR: output directory already exists: ${projectDir}/msa/1e62_wt_hbsag_ces" >&2
  echo "Move or remove it before rerunning to avoid timestamped AF3 output directories." >&2
  exit 3
fi

mkdir -p "${projectDir}/msa" "${projectDir}/logs"
cd "${af3Bin}"

exec > >(tee "${logFile}") 2>&1

apptainer exec \
  --bind "${projectDir}":/root/af_input \
  --bind "${projectDir}":/root/af_output \
  --bind "${af3Database}/models":/root/models \
  --bind "${af3Ssd}":/root/public_databases \
  --bind "${af3Database}":/root/public_databases_fallback \
  alphafold3.0.1.sif \
  python run_alphafold.py \
    --json_path=/root/af_input/input/1e62_wt_hbsag_ces.json \
    --model_dir=/root/models \
    --db_dir=/root/public_databases \
    --db_dir=/root/public_databases_fallback \
    --output_dir=/root/af_output/msa \
    --run_inference=False
```

- [ ] **Step 2: Check script syntax**

Run:

```bash
bash -n experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_data_pipeline_gpu6.sh
```

Expected: command exits with status `0` and prints no output.

- [ ] **Step 3: Run the data pipeline directly on gpu6**

Run:

```bash
bash experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_data_pipeline_gpu6.sh
```

Expected: command runs in the current `gpu6` shell. It does not call `sbatch` or Slurm. The log includes `Skipping model inference`.

- [ ] **Step 4: Verify data-pipeline JSON exists**

Run:

```bash
test -s experiments/1E62/R4/af3_ces_wt/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json
```

Expected: command exits with status `0` and prints no output.

- [ ] **Step 5: Verify the data JSON still has the expected chains and seeds**

Run:

```bash
python - <<'PY'
import json
from pathlib import Path

path = Path("experiments/1E62/R4/af3_ces_wt/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json")
data = json.loads(path.read_text())
chains = []
for item in data["sequences"]:
    protein = item["protein"]
    chains.append((protein["id"], len(protein["sequence"])))
print(chains)
print(data["modelSeeds"])
assert chains == [("H", 115), ("L", 111), ("A", 226)]
assert data["name"] == "1e62_wt_hbsag_ces"
assert data["modelSeeds"] == [75801088, 635954951, 495296019, 1129861070, 2026110668]
PY
```

Expected output:

```text
[('H', 115), ('L', 111), ('A', 226)]
[75801088, 635954951, 495296019, 1129861070, 2026110668]
```

---

### Task 4: Run AF3 Inference Directly on gpu6

**Files:**
- Create: `experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_infer_gpu6.sh`
- Produce: `experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_model.cif`
- Produce: `experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_summary_confidences.json`
- Produce: `experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces/ranking_scores.csv`
- Produce: `experiments/1E62/R4/af3_ces_wt/logs/af3_ces_infer_gpu6.log`

- [ ] **Step 1: Confirm the data-pipeline JSON is present**

Run:

```bash
test -s experiments/1E62/R4/af3_ces_wt/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json
```

Expected: command exits with status `0` and prints no output.

- [ ] **Step 2: Write the direct inference script**

Create `experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_infer_gpu6.sh` with this exact content:

```bash
#!/bin/bash
set -euo pipefail

if [ "$(hostname)" != "gpu6" ]; then
  echo "ERROR: this script must run directly on gpu6; current host is $(hostname)" >&2
  exit 2
fi

projectDir=/public/home/ziyang/code/optim-pipe/experiments/1E62/R4/af3_ces_wt
af3Bin=/public/apps/alphafold3/alphafold3
af3Database=/public/shared/alphafold3
af3Ssd=/ssd/alphafold3
logFile="${projectDir}/logs/af3_ces_infer_gpu6.log"

if [ ! -s "${projectDir}/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json" ]; then
  echo "ERROR: missing data-pipeline JSON: ${projectDir}/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json" >&2
  exit 3
fi

if [ -d "${projectDir}/inference/1e62_wt_hbsag_ces" ]; then
  echo "ERROR: output directory already exists: ${projectDir}/inference/1e62_wt_hbsag_ces" >&2
  echo "Move or remove it before rerunning to avoid timestamped AF3 output directories." >&2
  exit 4
fi

mkdir -p "${projectDir}/inference" "${projectDir}/logs"
cd "${af3Bin}"

export CUDA_VISIBLE_DEVICES=0
exec > >(tee "${logFile}") 2>&1

apptainer exec --nv \
  --bind "${projectDir}":/root/af_input \
  --bind "${projectDir}":/root/af_output \
  --bind "${af3Database}/models":/root/models \
  --bind "${af3Ssd}":/root/public_databases \
  --bind "${af3Database}":/root/public_databases_fallback \
  alphafold3.0.1.sif \
  python run_alphafold.py \
    --json_path=/root/af_input/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json \
    --model_dir=/root/models \
    --db_dir=/root/public_databases \
    --db_dir=/root/public_databases_fallback \
    --output_dir=/root/af_output/inference \
    --gpu_device=0 \
    --num_diffusion_samples=1 \
    --run_data_pipeline=False
```

- [ ] **Step 3: Check script syntax**

Run:

```bash
bash -n experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_infer_gpu6.sh
```

Expected: command exits with status `0` and prints no output.

- [ ] **Step 4: Run inference directly on gpu6**

Run:

```bash
bash experiments/1E62/R4/af3_ces_wt/scripts/af3_ces_infer_gpu6.sh
```

Expected: command runs in the current `gpu6` shell. It does not call `sbatch` or Slurm. The log includes `Skipping data pipeline`, `Writing outputs`, and `Done running 1 fold jobs`.

- [ ] **Step 5: Verify exactly five inference outputs exist**

Run:

```bash
test -s experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_model.cif
test -s experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_summary_confidences.json
test -s experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces/ranking_scores.csv
python - <<'PY'
import csv
from pathlib import Path

job = Path("experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces")
with (job / "ranking_scores.csv").open() as f:
    rows = list(csv.DictReader(f))
print(len(rows))
assert len(rows) == 5
for row in rows:
    path = job / f"seed-{row['seed']}_sample-{row['sample']}" / "model.cif"
    assert path.exists() and path.stat().st_size > 0, path
PY
```

Expected: the three `test` commands exit with status `0`; the Python check prints `5`.

---

### Task 5: Summarize Confidence and Promote the Five Ranked Structures

**Files:**
- Create: `experiments/1E62/R4/af3_ces_wt/summary/confidence_summary.tsv`
- Create: `experiments/1E62/R4/af3_ces_wt/summary/ranked_structures.tsv`
- Create: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank1.cif`
- Create: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank2.cif`
- Create: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank3.cif`
- Create: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank4.cif`
- Create: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank5.cif`
- Create: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT.cif`

- [ ] **Step 1: Extract ranking, copy five structures, and write summaries**

Run:

```bash
python - <<'PY'
import csv
import json
import shutil
from pathlib import Path

job = Path("experiments/1E62/R4/af3_ces_wt/inference/1e62_wt_hbsag_ces")
summary_path = job / "1e62_wt_hbsag_ces_summary_confidences.json"
ranking_path = job / "ranking_scores.csv"
summary_dir = Path("experiments/1E62/R4/af3_ces_wt/summary")
data_pdb_dir = Path("experiments/1E62/R4/data/pdb")
confidence_path = summary_dir / "confidence_summary.tsv"
manifest_path = summary_dir / "ranked_structures.tsv"

summary = json.loads(summary_path.read_text())
with ranking_path.open() as f:
    rows = list(csv.DictReader(f))

ranked = sorted(rows, key=lambda row: float(row["ranking_score"]), reverse=True)
assert len(ranked) == 5, f"Expected 5 ranked structures, found {len(ranked)}"

summary_dir.mkdir(parents=True, exist_ok=True)
data_pdb_dir.mkdir(parents=True, exist_ok=True)

manifest_rows = []
for rank, row in enumerate(ranked, start=1):
    src = job / f"seed-{row['seed']}_sample-{row['sample']}" / "model.cif"
    dst = data_pdb_dir / f"AF3-1E62-CeS-WT-rank{rank}.cif"
    assert src.exists() and src.stat().st_size > 0, src
    shutil.copy2(src, dst)
    manifest_rows.append([
        rank,
        row["seed"],
        row["sample"],
        row["ranking_score"],
        str(src),
        str(dst),
    ])

top_alias = data_pdb_dir / "AF3-1E62-CeS-WT.cif"
shutil.copy2(data_pdb_dir / "AF3-1E62-CeS-WT-rank1.cif", top_alias)

best = ranked[0]
chain_order = ["H", "L", "A"]
pair_iptm = summary.get("chain_pair_iptm", [])
pair_pae = summary.get("chain_pair_pae_min", [])

records = [
    ("job", "1e62_wt_hbsag_ces"),
    ("num_returned_structures", "5"),
    ("randomized_model_seeds", "75801088,635954951,495296019,1129861070,2026110668"),
    ("num_diffusion_samples_per_seed", "1"),
    ("rank1_seed", best["seed"]),
    ("rank1_sample", best["sample"]),
    ("top_ranking_score", best["ranking_score"]),
    ("ptm", str(summary.get("ptm", "NA"))),
    ("iptm", str(summary.get("iptm", "NA"))),
    ("fraction_disordered", str(summary.get("fraction_disordered", "NA"))),
    ("has_clash", str(summary.get("has_clash", "NA"))),
]

for i, left in enumerate(chain_order):
    for j, right in enumerate(chain_order):
        if i < len(pair_iptm) and j < len(pair_iptm[i]):
            records.append((f"chain_pair_iptm_{left}_{right}", str(pair_iptm[i][j])))
        if i < len(pair_pae) and j < len(pair_pae[i]):
            records.append((f"chain_pair_pae_min_{left}_{right}", str(pair_pae[i][j])))

with confidence_path.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerows(records)

with manifest_path.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["rank", "seed", "sample", "ranking_score", "source_cif", "copied_cif"])
    writer.writerows(manifest_rows)

print(confidence_path)
print(manifest_path)
PY
```

Expected output:

```text
experiments/1E62/R4/af3_ces_wt/summary/confidence_summary.tsv
experiments/1E62/R4/af3_ces_wt/summary/ranked_structures.tsv
```

- [ ] **Step 2: Display the confidence and ranked-structure summaries**

Run:

```bash
column -t -s $'\t' experiments/1E62/R4/af3_ces_wt/summary/confidence_summary.tsv
column -t -s $'\t' experiments/1E62/R4/af3_ces_wt/summary/ranked_structures.tsv
```

Expected: the confidence table includes `num_returned_structures`, `rank1_seed`, `top_ranking_score`, `ptm`, `iptm`, `chain_pair_iptm_H_A`, `chain_pair_iptm_L_A`, `chain_pair_pae_min_H_A`, and `chain_pair_pae_min_L_A` when those AF3 metrics are present. The ranked-structure table has exactly five data rows.

- [ ] **Step 3: Verify the five promoted ranked structures exist**

Run:

```bash
for rank in 1 2 3 4 5; do
  test -s "experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank${rank}.cif"
done
test -s experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT.cif
ls -lh experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank*.cif experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT.cif
```

Expected: `ls` prints five non-empty ranked files plus the top-ranked alias `AF3-1E62-CeS-WT.cif`.

- [ ] **Step 4: Check whether each promoted AF3 output contains all chain IDs**

Run:

```bash
python - <<'PY'
from pathlib import Path
import re

for rank in range(1, 6):
    path = Path(f"experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank{rank}.cif")
    text = path.read_text(errors="replace")
    present = {chain: bool(re.search(rf"\b{chain}\b", text)) for chain in ["H", "L", "A"]}
    print(rank, present)
    assert all(present.values()), path
PY
```

Expected output:

```text
1 {'H': True, 'L': True, 'A': True}
2 {'H': True, 'L': True, 'A': True}
3 {'H': True, 'L': True, 'A': True}
4 {'H': True, 'L': True, 'A': True}
5 {'H': True, 'L': True, 'A': True}
```

---

### Task 6: Completion Report

**Files:**
- Read: `experiments/1E62/R4/af3_ces_wt/summary/confidence_summary.tsv`
- Read: `experiments/1E62/R4/af3_ces_wt/summary/ranked_structures.tsv`
- Read: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT.cif`
- Read: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank1.cif`
- Read: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank2.cif`
- Read: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank3.cif`
- Read: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank4.cif`
- Read: `experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank5.cif`
- Read: `experiments/1E62/R4/af3_ces_wt/logs/af3_ces_data_pipeline_gpu6.log`
- Read: `experiments/1E62/R4/af3_ces_wt/logs/af3_ces_infer_gpu6.log`

- [ ] **Step 1: Collect final paths and confidence highlights**

Run:

```bash
python - <<'PY'
from pathlib import Path

structure = Path("experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT.cif")
ranked = sorted(Path("experiments/1E62/R4/data/pdb").glob("AF3-1E62-CeS-WT-rank*.cif"))
summary = Path("experiments/1E62/R4/af3_ces_wt/summary/confidence_summary.tsv")
manifest = Path("experiments/1E62/R4/af3_ces_wt/summary/ranked_structures.tsv")
print(f"structure={structure} exists={structure.exists()} size={structure.stat().st_size if structure.exists() else 0}")
for path in ranked:
    print(f"ranked_structure={path} exists={path.exists()} size={path.stat().st_size if path.exists() else 0}")
print(summary.read_text())
print(manifest.read_text())
PY
```

Expected: output includes one top-ranked alias and five ranked structures with `exists=True`, positive file sizes, the confidence table, and the ranked-structure manifest.

- [ ] **Step 2: State interpretation boundary**

Report these points to the user:

```text
AF3 produced the requested full-length 1E62 heavy/light + HBsAg CeS complex structures.
The five ranked structure paths are experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank1.cif through experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT-rank5.cif.
The top-ranked convenience alias is experiments/1E62/R4/data/pdb/AF3-1E62-CeS-WT.cif.
Use chain-pair metrics H-A and L-A to judge antibody-antigen interface confidence.
Because HBsAg is membrane-associated and the run used the full-length supplied sequence, low interface confidence or high disorder should trigger a separate antigenic-loop/truncation run rather than silent reuse of this structure.
```

- [ ] **Step 3: Do not commit runtime outputs**

Run:

```bash
git status --short
```

Expected: do not stage generated files under `experiments/` unless the user explicitly asks to preserve them in git. The structures can be used from the shared workspace paths without committing them.

---

## Self-Review

**Spec coverage:** The user supplied heavy chain, light chain, and HBsAg CeS sequences and asked to obtain the complex structure. Task 1 creates the AF3 input with those exact sequences and randomized seeds. Task 2 verifies that the current shell is directly on `gpu6` with visible A40 resources. Task 3 runs CPU data-pipeline work directly on `gpu6`. Task 4 runs GPU inference directly on `gpu6` with one diffusion sample per seed, producing exactly five structures. Task 5 promotes all five ranked structures and summarizes confidence. Task 6 reports the five final structure paths and interpretation boundary.

**Placeholder scan:** The plan contains exact paths, exact sequences, exact direct-run scripts, exact commands, and exact expected checks. It does not require unspecified script names, unspecified sequences, or unspecified directories.

**Type consistency:** The AF3 job name is consistently `1e62_wt_hbsag_ces`. Chain IDs are consistently `H`, `L`, and `A`. The data-pipeline path used by inference is consistently `experiments/1E62/R4/af3_ces_wt/msa/1e62_wt_hbsag_ces/1e62_wt_hbsag_ces_data.json`. Final ranked files use `AF3-1E62-CeS-WT-rank1.cif` through `AF3-1E62-CeS-WT-rank5.cif`.
