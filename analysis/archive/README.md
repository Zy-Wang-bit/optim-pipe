# Analysis Archive

This directory keeps historical analysis scripts that are no longer active
entrypoints.

## `initial_design_generation_scripts_20260703.zip`

Archived tracked snapshot of `analysis/initial_design_generation/`.

- Created: 2026-07-03
- Contents: 61 files with their original relative paths, including 52 `.py`
  scripts plus the historical config and schema files.
- The live tree no longer keeps a partial copy of this workflow.

Restore when needed:

```bash
python3 - <<'PY'
from zipfile import ZipFile
with ZipFile("analysis/archive/initial_design_generation_scripts_20260703.zip") as zf:
    zf.extractall(".")
PY
```

After restore, run from the repository root and expect historical environment
paths/results under `results/initial_design_generation/`.

## `per_target_15k_candidate_pools_20260705.zip`

Archived copy of the two per-target 15K delivery candidate pools from
`results/initial_design_generation/per_target_15k_candidate_pools/`.

- Created: 2026-07-05
- Contents: 9 files with their original relative paths.
- Main CSVs:
  - `final_candidate_pool_1E62_15k_draft.csv`: 15,000 rows, 78 columns.
  - `final_candidate_pool_sdAb_15k_draft.csv`: 15,000 rows, 78 columns.
- Includes the audit reports, input manifest, hard-failure table, control
  anchor plan, and per-target comparison note.
- SHA256:
  `af3bd3cc67a020d24b094766702b888a44051f62121ef219cff9a2bb0a6d125b`

Restore when needed:

```bash
python3 - <<'PY'
from zipfile import ZipFile
with ZipFile("analysis/archive/per_target_15k_candidate_pools_20260705.zip") as zf:
    zf.extractall(".")
PY
```

## 2026-07-03 cleanup index

The current tree keeps the runnable 3-tier pipeline and removes historical
analysis surfaces. The last full pre-cleanup commit is `07ec8b9`
(`chore: prepare initial design deliverables`).

Restore any removed path with:

```bash
git checkout 07ec8b9 -- <path>
```

Removed from the live tree:

- `analysis/initial_design_generation/` historical workflow snapshot
- `analysis/window_selection/` completed window-selection workflow
- `analysis/visualize/` one-off plotting and slide scripts
- `analysis/md/` legacy MD post-analysis scripts
- standalone Rosetta analysis CLIs not used by `scripts/run_rosetta_eval.py`
- historical batch/node config matrix under `configs/`

Kept in the live tree:

- `run_pipeline.sh`
- `scripts/`
- `analysis/rosetta/calc_ph_score_antibody.py`
- `analysis/rosetta/batch_calc_dddg_elec.py`
- `analysis/pka/`
- `analysis/r4/`
- `analysis/naming/`
- `third_party/molecular_dynamics/`

## Retired 8-step / Phase C pipeline

The old `merge_and_select.py` / `phase_c` entrypoint was removed from current
HEAD. To reproduce those historical outputs, checkout commit `d7726f1`
(`feat: extend pipeline to 4 Phase / 12 Step with optional Phase C evaluation`)
or `64cfb02` (`feat: restructure pipeline from 4 Phase / 12 Step to 3 Tier /
13 Step`).
