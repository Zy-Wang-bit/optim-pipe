# optim-pipe 项目重组 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganize the optim-pipe project: clean up 28 root items to ~11 visible, consolidate external dependencies, restructure experiments, establish cross-session task management, and write CLAUDE.md work protocol.

**Architecture:** File moves + path reference updates. No logic changes. All config-driven paths flow through YAML files; hardcoded absolute paths in experiment scripts also need updating. SimpleFold weights move from `artifacts/` into `third_party/ml-simplefold/weights/`, with CLI default and wrapper default updated.

**Tech Stack:** Bash (file ops), Python YAML configs, shell scripts

**Spec:** `docs/superpowers/specs/2026-04-03-project-reorganization-design.md`

---

### Task 1: Create target directory structure

**Files:**
- Create: `docs/design/`, `docs/reports/`, `.tasks/active/`, `.tasks/done/sdab-pka/`, `experiments/1E62/data/`, `third_party/foldx/`, `third_party/molecular_dynamics/`, `third_party/ml-simplefold/weights/`

- [ ] **Step 1: Create all target directories**

```bash
cd /public/home/ziyang/code/optim-pipe
mkdir -p docs/design docs/reports
mkdir -p .tasks/active .tasks/done/sdab-pka
mkdir -p experiments/1E62/data
mkdir -p third_party/foldx
mkdir -p third_party/molecular_dynamics
mkdir -p third_party/ml-simplefold/weights
```

- [ ] **Step 2: Verify directories exist**

Run: `ls -d docs/design docs/reports .tasks/active .tasks/done/sdab-pka experiments/1E62/data third_party/foldx third_party/molecular_dynamics third_party/ml-simplefold/weights`
Expected: All 8 paths listed without error.

- [ ] **Step 3: Commit**

```bash
git add docs/.gitkeep .tasks/.gitkeep
git commit -m "chore: create target directory structure for project reorganization"
```

Note: empty directories won't be tracked by git unless they contain files. The subsequent tasks will populate them.

---

### Task 2: Move documentation files to docs/

**Files:**
- Move: `ROADMAP.md` → `docs/ROADMAP.md`
- Move: `Phase1-combination_design.md` → `docs/design/Phase1-combination_design.md`
- Move: `Phase1-binary_classification.md` → `docs/design/Phase1-binary_classification.md`
- Move: `analysis/correlation_report.md` → `docs/reports/correlation_report.md`
- Move: `analysis/deep_metric_analysis.md` → `docs/reports/deep_metric_analysis.md`

- [ ] **Step 1: Move files**

```bash
cd /public/home/ziyang/code/optim-pipe
mv ROADMAP.md docs/ROADMAP.md
mv Phase1-combination_design.md docs/design/Phase1-combination_design.md
mv Phase1-binary_classification.md docs/design/Phase1-binary_classification.md
mv analysis/correlation_report.md docs/reports/correlation_report.md
mv analysis/deep_metric_analysis.md docs/reports/deep_metric_analysis.md
```

- [ ] **Step 2: Verify moves**

Run: `ls docs/ROADMAP.md docs/design/Phase1-*.md docs/reports/*.md`
Expected: All 5 files listed.

- [ ] **Step 3: Commit**

```bash
git add docs/ analysis/
git add -u ROADMAP.md Phase1-combination_design.md Phase1-binary_classification.md
git commit -m "chore: move documentation to docs/"
```

---

### Task 3: Move task management files to .tasks/

**Files:**
- Move: `task_plan.md` → `.tasks/done/sdab-pka/plan.md`
- Move: `progress.md` → `.tasks/done/sdab-pka/progress.md`
- Move: `findings.md` → `.tasks/done/sdab-pka/findings.md`

- [ ] **Step 1: Move files**

```bash
cd /public/home/ziyang/code/optim-pipe
mv task_plan.md .tasks/done/sdab-pka/plan.md
mv progress.md .tasks/done/sdab-pka/progress.md
mv findings.md .tasks/done/sdab-pka/findings.md
```

- [ ] **Step 2: Verify**

Run: `ls .tasks/done/sdab-pka/`
Expected: `findings.md  plan.md  progress.md`

- [ ] **Step 3: Commit**

```bash
git add .tasks/
git add -u task_plan.md progress.md findings.md
git commit -m "chore: move task files to .tasks/done/sdab-pka/"
```

---

### Task 4: Restructure experiments/1E62 (move R1/R2/R3)

**Files:**
- Move: `experiments/1E62_R1/` → `experiments/1E62/R1/`
- Move: `experiments/1E62_R2/` → `experiments/1E62/R2/`
- Move: `experiments/1E62_R3/` → `experiments/1E62/R3/`

- [ ] **Step 1: Move experiment directories**

```bash
cd /public/home/ziyang/code/optim-pipe
mv experiments/1E62_R1 experiments/1E62/R1
mv experiments/1E62_R2 experiments/1E62/R2
mv experiments/1E62_R3 experiments/1E62/R3
```

- [ ] **Step 2: Verify structure**

Run: `ls experiments/1E62/`
Expected: `R1  R2  R3  data`

- [ ] **Step 3: Commit**

```bash
git add experiments/1E62/
git commit -m "chore: restructure 1E62 experiments into experiments/1E62/R1-R3"
```

---

### Task 5: Migrate data/ to experiments/1E62/data/

**Files:**
- Move: `data/pdbs/n{1-4}-0_Repair.pdb` → `experiments/1E62/data/AF3-abag-n{1-4}.pdb` (rename)
- Move: `data/wt/*` → `experiments/1E62/data/`
- Move: `data/molecular_dynamics/` → `third_party/molecular_dynamics/`
- Delete: `data/pdb_AF3/`, `data/pdbs/n{1-4}-0.pdb` (non-Repair), `data/` (empty shell after)

- [ ] **Step 1: Move and rename PDB Repair files**

```bash
cd /public/home/ziyang/code/optim-pipe
cp data/pdbs/n1-0_Repair.pdb experiments/1E62/data/AF3-abag-n1.pdb
cp data/pdbs/n2-0_Repair.pdb experiments/1E62/data/AF3-abag-n2.pdb
cp data/pdbs/n3-0_Repair.pdb experiments/1E62/data/AF3-abag-n3.pdb
cp data/pdbs/n4-0_Repair.pdb experiments/1E62/data/AF3-abag-n4.pdb
```

- [ ] **Step 2: Move WT files**

```bash
cp data/wt/ab_wt.pdb experiments/1E62/data/ab_wt.pdb
cp data/wt/heavy.fasta experiments/1E62/data/heavy.fasta
cp data/wt/light.fasta experiments/1E62/data/light.fasta
cp data/wt/antigen.fasta experiments/1E62/data/antigen.fasta
cp data/wt/original.fasta experiments/1E62/data/original.fasta
cp data/wt/wt_binding_energy.csv experiments/1E62/data/wt_binding_energy.csv
```

- [ ] **Step 3: Move molecular dynamics scripts**

```bash
cp -r data/molecular_dynamics/* third_party/molecular_dynamics/
```

- [ ] **Step 4: Verify all copied files**

Run: `ls -la experiments/1E62/data/ && ls third_party/molecular_dynamics/`
Expected: 10 files in data/ (4 PDB + 6 from wt), MD scripts in molecular_dynamics/.

- [ ] **Step 5: Delete source data/**

```bash
rm -rf data/pdb_AF3
rm data/pdbs/n1-0.pdb data/pdbs/n2-0.pdb data/pdbs/n3-0.pdb data/pdbs/n4-0.pdb
rm data/pdbs/n1-0_Repair.pdb data/pdbs/n2-0_Repair.pdb data/pdbs/n3-0_Repair.pdb data/pdbs/n4-0_Repair.pdb
rm data/wt/ab_wt.pdb data/wt/heavy.fasta data/wt/light.fasta data/wt/antigen.fasta data/wt/original.fasta data/wt/wt_binding_energy.csv
rm -rf data/molecular_dynamics
rmdir data/pdbs data/wt data
```

- [ ] **Step 6: Verify data/ is gone**

Run: `test -d data && echo "STILL EXISTS" || echo "REMOVED"`
Expected: `REMOVED`

- [ ] **Step 7: Commit**

```bash
git add experiments/1E62/data/ third_party/molecular_dynamics/
git add -u data/
git commit -m "chore: migrate data/ into experiments/1E62/data/ and third_party/"
```

---

### Task 6: Consolidate third_party/ (merge tools/ and artifacts/)

**Files:**
- Move: `tools/foldx/foldx` → `third_party/foldx/foldx`
- Move: `tools/foldx/` other files → `third_party/foldx/`
- Move: `artifacts/*.ckpt` → `third_party/ml-simplefold/weights/`
- Delete: `tools/`, `artifacts/`

- [ ] **Step 1: Move FoldX**

```bash
cd /public/home/ziyang/code/optim-pipe
cp -r tools/foldx/* third_party/foldx/
```

- [ ] **Step 2: Move SimpleFold weights**

```bash
cp artifacts/simplefold_3B.ckpt third_party/ml-simplefold/weights/
cp artifacts/simplefold_1.6B.ckpt third_party/ml-simplefold/weights/
cp artifacts/plddt_module_1.6B.ckpt third_party/ml-simplefold/weights/
cp artifacts/plddt.ckpt third_party/ml-simplefold/weights/
```

- [ ] **Step 3: Verify copies**

Run: `ls -lh third_party/foldx/foldx && ls -lh third_party/ml-simplefold/weights/`
Expected: FoldX binary (~83MB), 4 checkpoint files (11G, 6G, 442M, 442M).

- [ ] **Step 4: Delete source directories**

```bash
rm -rf tools
rm -rf artifacts
```

- [ ] **Step 5: Commit**

```bash
git add third_party/foldx/ third_party/ml-simplefold/weights/
git add -u tools/ artifacts/
git commit -m "chore: consolidate tools/ and artifacts/ into third_party/"
```

---

### Task 7: Delete obsolete directories

**Files:**
- Delete: `results/` (153M historical data)
- Delete: `nano01/` (obsolete nanobody template)
- Delete: `molecules/` (empty directory)

- [ ] **Step 1: Delete obsolete items**

```bash
cd /public/home/ziyang/code/optim-pipe
rm -rf results
rm -rf nano01
rmdir molecules
```

- [ ] **Step 2: Add results/ to .gitignore for future pipeline runs**

The main pipeline (`run_pipeline.sh`) creates `results/` at runtime. Add it to `.gitignore` so it doesn't re-clutter the repo.

Create/update `.gitignore`:
```
# Pipeline runtime output (recreated by run_pipeline.sh)
results/
foldx/
mpnn_outputs/
his_seeds/
esm_scores/
logs/
```

- [ ] **Step 3: Verify deletions**

Run: `ls -d results nano01 molecules 2>&1`
Expected: All three report "No such file or directory"

- [ ] **Step 4: Commit**

```bash
git add .gitignore
git add -u results/ nano01/ molecules/
git commit -m "chore: remove obsolete results/, nano01/, molecules/"
```

---

### Task 8: Update main config.yaml paths

**Files:**
- Modify: `configs/config.yaml`

The main config currently has:
- `foldx_bin: "/public/home/ziyang/code/optim-pipe/foldX/foldx"` (line 3, note: already points to `foldX/` not `tools/`)
- `pdb_dir: "data/pdbs"` (line 4)
- `wt_dir: "data/wt"` (line 5)

- [ ] **Step 1: Update config.yaml**

Edit `configs/config.yaml`:

```yaml
paths:
  foldx_bin: "third_party/foldx/foldx"
  pdb_dir: "experiments/1E62/data"
  wt_dir: "experiments/1E62/data"
  wt_files:
    A: "heavy.fasta"
    B: "light.fasta"
    C: "antigen.fasta"
  mpnn_out_dir: "mpnn_outputs"
  his_seed_dir: "his_seeds"
  esm_out_dir: "esm_scores"
  foldx_dir: "foldx"
  results_dir: "results"
  mpnn_runner: "third_party/ProteinMPNN/protein_mpnn_run.py"
```

Note: `mpnn_runner` was `ProteinMPNN/protein_mpnn_run.py`, now needs `third_party/` prefix. Verify this is loaded relative to project root.

- [ ] **Step 2: Validate YAML syntax**

Run: `python -c "import yaml; yaml.safe_load(open('configs/config.yaml')); print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add configs/config.yaml
git commit -m "chore: update config.yaml paths for new directory structure"
```

---

### Task 9: Update all PDB-specific config files

**Files:**
- Modify: `configs/config_foldx_n1-0.yaml`
- Modify: `configs/config_foldx_n2-0.yaml`
- Modify: `configs/config_foldx_n3-0.yaml`
- Modify: `configs/config_foldx_n4-0.yaml`
- Modify: `configs/config_foldx_n1-0_chainA.yaml`
- Modify: `configs/config_foldx_n1-0_chainB.yaml`
- Modify: `configs/config_foldx_n1-0_ext_HL.yaml`
- Modify: `configs/config_esm_n1-0.yaml`
- Modify: `configs/config_esm_n2-0.yaml`
- Modify: `configs/config_esm_n3-0.yaml`
- Modify: `configs/config_esm_n4-0.yaml`
- Modify: `configs/config_his.yaml`

All these configs share the same pattern. In each file, update:

| Old | New |
|-----|-----|
| `foldx_bin: tools/foldx/foldx` | `foldx_bin: third_party/foldx/foldx` |
| `pdb_dir: data/pdbs` | `pdb_dir: experiments/1E62/data` |
| `wt_dir: data/wt` | `wt_dir: experiments/1E62/data` |
| `mpnn_runner: third_party/ProteinMPNN/protein_mpnn_run.py` | (no change, already correct) |

- [ ] **Step 1: Batch update all config files with sed**

```bash
cd /public/home/ziyang/code/optim-pipe
for f in configs/config_foldx_*.yaml configs/config_esm_*.yaml configs/config_his.yaml; do
  sed -i 's|foldx_bin: tools/foldx/foldx|foldx_bin: third_party/foldx/foldx|g' "$f"
  sed -i 's|pdb_dir: data/pdbs|pdb_dir: experiments/1E62/data|g' "$f"
  sed -i 's|pdb_dir: "data/pdbs"|pdb_dir: "experiments/1E62/data"|g' "$f"
  sed -i 's|wt_dir: data/wt|wt_dir: experiments/1E62/data|g' "$f"
  sed -i 's|wt_dir: "data/wt"|wt_dir: "experiments/1E62/data"|g' "$f"
done
```

- [ ] **Step 2: Verify no old paths remain in configs**

Run: `grep -r "data/pdbs\|data/wt\|tools/foldx" configs/`
Expected: No output (all references updated).

- [ ] **Step 3: Validate all YAML files**

```bash
for f in configs/*.yaml; do
  python -c "import yaml; yaml.safe_load(open('$f')); print('OK: $f')"
done
```
Expected: All print `OK`.

- [ ] **Step 4: Commit**

```bash
git add configs/
git commit -m "chore: update all config YAML files for new directory structure"
```

---

### Task 10: Update experiment config files

**Files:**
- Modify: `experiments/1E62/R3/config_r3_eval.yaml`
- Modify: `experiments/1E62/R2/config_r2_pka.yaml`

- [ ] **Step 1: Update R3 eval config**

Edit `experiments/1E62/R3/config_r3_eval.yaml`:

| Old | New |
|-----|-----|
| `foldx_bin: tools/foldx/foldx` | `foldx_bin: third_party/foldx/foldx` |
| `wt_pdb: data/wt/ab_wt.pdb` | `wt_pdb: experiments/1E62/data/ab_wt.pdb` |
| `template_pdb: data/pdbs/n1-0_Repair.pdb` | `template_pdb: experiments/1E62/data/AF3-abag-n1.pdb` |
| `r3_csv: experiments/1E62_R3/design/r3_combinations.csv` | `r3_csv: experiments/1E62/R3/design/r3_combinations.csv` |
| `output_dir: experiments/1E62_R3/results` | `output_dir: experiments/1E62/R3/results` |

- [ ] **Step 2: Update R2 pKa config**

Edit `experiments/1E62/R2/config_r2_pka.yaml`:

| Old | New |
|-----|-----|
| `template_pdb: data/pdbs/n1-0_Repair.pdb` | `template_pdb: experiments/1E62/data/AF3-abag-n1.pdb` |
| `wt_pdb: data/wt/ab_wt.pdb` | `wt_pdb: experiments/1E62/data/ab_wt.pdb` |
| `elisa_csv: experiments/1E62_R2/wet_lab/elisa_summary.csv` | `elisa_csv: experiments/1E62/R2/wet_lab/elisa_summary.csv` |
| `output_dir: experiments/1E62_R2/results` | `output_dir: experiments/1E62/R2/results` |

- [ ] **Step 3: Validate both configs**

```bash
python -c "import yaml; yaml.safe_load(open('experiments/1E62/R3/config_r3_eval.yaml')); print('R3 OK')"
python -c "import yaml; yaml.safe_load(open('experiments/1E62/R2/config_r2_pka.yaml')); print('R2 OK')"
```

- [ ] **Step 4: Commit**

```bash
git add experiments/1E62/R3/config_r3_eval.yaml experiments/1E62/R2/config_r2_pka.yaml
git commit -m "chore: update experiment config files for new paths"
```

---

### Task 11: Update experiment Python scripts with hardcoded paths

**Files:**
- Modify: `experiments/1E62/R1/run_r1_pka.py` (lines 24-26)
- Modify: `experiments/1E62/R1/migrate_naming.py` (line 18)
- Modify: `experiments/1E62/R2/scripts/run_foldx_r2_binding.py` (lines 12-14)
- Modify: `experiments/1E62/R2/scripts/run_rosetta_interface_ddg.py` (lines 10-12)
- Modify: `experiments/1E62/R2/scripts/analyze_simplefold_3x_rmsd.py` (lines 18-20)
- Modify: `experiments/1E62/R2/migrate_naming.py` (line 20)
- Modify: `experiments/1E62/R2/verify_mapping.py`
- Modify: `experiments/1E62/R3/scripts/calibrate_weights.py` (lines 365-369)
- Modify: `experiments/1E62/R3/scripts/build_slides.py` (line 25)
- Modify: `experiments/1E62/R3/scripts/generate_report_figures.py` (lines 16-26)
- Modify: `experiments/1E62/R3/scripts/build_r3_structures.py`
- Modify: `experiments/1E62/R3/scripts/analyze_simplefold_3x_rmsd.py` (lines 17-20)
- Modify: `experiments/1E62/R3/migrate_naming.py` (line 19)
- Modify: `experiments/1E62/R3/eval_r3.sh`
- Modify: `analysis/elisa_vs_computation.py` (lines 26-27)
- Modify: `scripts/build_mutants_from_csv.py` (line 75)

All these files contain hardcoded references to the old directory structure. The changes are mechanical string replacements.

- [ ] **Step 1: Batch replace old paths in all experiment scripts**

Run these sed replacements across all Python and shell files under experiments/1E62/:

```bash
cd /public/home/ziyang/code/optim-pipe

# Replace 1E62_R1 → 1E62/R1, 1E62_R2 → 1E62/R2, 1E62_R3 → 1E62/R3
find experiments/1E62/ -type f \( -name "*.py" -o -name "*.sh" -o -name "*.md" \) -exec \
  sed -i 's|experiments/1E62_R1|experiments/1E62/R1|g; s|experiments/1E62_R2|experiments/1E62/R2|g; s|experiments/1E62_R3|experiments/1E62/R3|g' {} +

# Replace data/pdbs/n1-0_Repair.pdb → experiments/1E62/data/AF3-abag-n1.pdb
find experiments/1E62/ -type f \( -name "*.py" -o -name "*.sh" \) -exec \
  sed -i 's|data/pdbs/n1-0_Repair\.pdb|experiments/1E62/data/AF3-abag-n1.pdb|g; s|data/pdbs/n2-0_Repair\.pdb|experiments/1E62/data/AF3-abag-n2.pdb|g; s|data/pdbs/n3-0_Repair\.pdb|experiments/1E62/data/AF3-abag-n3.pdb|g; s|data/pdbs/n4-0_Repair\.pdb|experiments/1E62/data/AF3-abag-n4.pdb|g' {} +

# Replace data/pdbs/n1-0.pdb (non-Repair) → experiments/1E62/data/AF3-abag-n1.pdb
find experiments/1E62/ -type f -name "*.py" -exec \
  sed -i 's|data/pdbs/n1-0\.pdb|experiments/1E62/data/AF3-abag-n1.pdb|g; s|data/pdbs/n2-0\.pdb|experiments/1E62/data/AF3-abag-n2.pdb|g' {} +

# Replace data/wt paths
find experiments/1E62/ -type f \( -name "*.py" -o -name "*.sh" \) -exec \
  sed -i 's|data/wt/ab_wt\.pdb|experiments/1E62/data/ab_wt.pdb|g; s|data/wt|experiments/1E62/data|g' {} +

# Replace tools/foldx references
find experiments/1E62/ -type f \( -name "*.py" -o -name "*.sh" \) -exec \
  sed -i 's|tools/foldx/foldx|third_party/foldx/foldx|g' {} +
```

- [ ] **Step 2: Update analysis/elisa_vs_computation.py**

```bash
sed -i 's|experiments/1E62_R2|experiments/1E62/R2|g' analysis/elisa_vs_computation.py
```

- [ ] **Step 3: Update scripts/build_mutants_from_csv.py**

```bash
sed -i 's|data/pdbs/n1-0\.pdb|experiments/1E62/data/AF3-abag-n1.pdb|g' scripts/build_mutants_from_csv.py
sed -i 's|data/1E62_R2|experiments/1E62/R2|g' scripts/build_mutants_from_csv.py
```

- [ ] **Step 4: Verify no old path references remain in active code**

```bash
grep -rn "data/pdbs\|data/wt\|tools/foldx\|1E62_R[123]" scripts/ analysis/ experiments/1E62/ configs/ --include="*.py" --include="*.sh" --include="*.yaml" | grep -v "archive/" | grep -v ".md"
```
Expected: No output (all references updated). Markdown files may still have old references in prose — that's OK.

- [ ] **Step 5: Commit**

```bash
git add experiments/1E62/ analysis/ scripts/
git commit -m "chore: update hardcoded paths in experiment and analysis scripts"
```

---

### Task 12: Update SimpleFold default checkpoint path

**Files:**
- Modify: `third_party/ml-simplefold/src/simplefold/cli.py:20`
- Modify: `third_party/ml-simplefold/src/simplefold/wrapper.py:55`
- Modify: `third_party/ml-simplefold/configs/paths/default.yaml:15`

- [ ] **Step 1: Update cli.py default**

In `third_party/ml-simplefold/src/simplefold/cli.py`, line 20:

Old:
```python
parser.add_argument("--ckpt_dir", type=str, default="artifacts", help="Directory to save the checkpoint.")
```

New:
```python
parser.add_argument("--ckpt_dir", type=str, default="third_party/ml-simplefold/weights", help="Directory to save the checkpoint.")
```

- [ ] **Step 2: Update wrapper.py default**

In `third_party/ml-simplefold/src/simplefold/wrapper.py`, line 55:

Old:
```python
ckpt_dir="./artifacts",
```

New:
```python
ckpt_dir="./third_party/ml-simplefold/weights",
```

- [ ] **Step 3: Update default.yaml**

In `third_party/ml-simplefold/configs/paths/default.yaml`, line 15:

Old:
```yaml
output_dir: ${paths.root_dir}artifacts/
```

New:
```yaml
output_dir: ${paths.root_dir}third_party/ml-simplefold/weights/
```

- [ ] **Step 4: Verify SimpleFold weight files are accessible**

```bash
ls -lh third_party/ml-simplefold/weights/*.ckpt
```
Expected: 4 files listed with correct sizes.

- [ ] **Step 5: Commit**

```bash
git add third_party/ml-simplefold/src/simplefold/cli.py third_party/ml-simplefold/src/simplefold/wrapper.py third_party/ml-simplefold/configs/paths/default.yaml
git commit -m "chore: update SimpleFold weight paths from artifacts/ to third_party/ml-simplefold/weights/"
```

---

### Task 13: Update run_pipeline.sh

**Files:**
- Modify: `run_pipeline.sh`

- [ ] **Step 1: Update path references in run_pipeline.sh**

Line 34 comment: change `data/pdbs` and `data/wt` references to new locations.

```bash
sed -i 's|data/pdbs|experiments/1E62/data|g; s|data/wt|experiments/1E62/data|g' run_pipeline.sh
```

No other changes needed — `run_pipeline.sh` uses `results/` as a runtime working directory (covered by .gitignore), and reads paths from config.yaml.

- [ ] **Step 2: Verify**

Run: `grep -n "data/pdbs\|data/wt\|tools/foldx" run_pipeline.sh`
Expected: No output.

- [ ] **Step 3: Commit**

```bash
git add run_pipeline.sh
git commit -m "chore: update run_pipeline.sh path references"
```

---

### Task 14: Update .claude/settings files

**Files:**
- Modify: `.claude/settings.json`
- Modify: `.claude/settings.local.json`

- [ ] **Step 1: Update settings.json**

Replace `experiments/1E62_R3` with `experiments/1E62/R3` in `.claude/settings.json`.

```bash
sed -i 's|experiments/1E62_R3|experiments/1E62/R3|g' .claude/settings.json
```

- [ ] **Step 2: Update settings.local.json**

Replace `tools/foldx/foldx` with `third_party/foldx/foldx` in `.claude/settings.local.json`.

```bash
sed -i 's|tools/foldx/foldx|third_party/foldx/foldx|g' .claude/settings.local.json
```

- [ ] **Step 3: Commit**

```bash
git add .claude/settings.json .claude/settings.local.json
git commit -m "chore: update .claude settings for new paths"
```

---

### Task 15: Write CLAUDE.md

**Files:**
- Create: `CLAUDE.md` (root)

- [ ] **Step 1: Write CLAUDE.md with project overview and work protocol**

```markdown
# optim-pipe

pH 依赖性抗体亲和力优化的自动化计算流水线。

## 环境

- **主环境**: `conda activate optim-pipe` (Python 3.11, PyTorch 2.8, fair-esm 2.0, BioPython 1.85)
- **ProteinMPNN**: `conda activate proteinmpnn`
- **SimpleFold**: `conda activate simplefold`
- **PyRosetta**: `conda activate pyrosetta`

## 目录结构

| 目录 | 职责 |
|------|------|
| `scripts/` | 核心 pipeline 脚本（8 步） |
| `configs/` | YAML 配置文件 |
| `analysis/` | 后分析代码（pKa、Rosetta、RMSD、ELISA 关联） |
| `experiments/` | 按抗体系统/轮次组织的实验数据 |
| `docs/` | 路线图、设计文档、分析报告 |
| `third_party/` | 外部工具（FoldX、ProteinMPNN、SimpleFold、MD） |
| `archive/` | 历史归档 |
| `.tasks/` | 跨对话任务管理 |

## 跨对话任务协议

### 新对话启动
1. 读取 `.tasks/active/` 下所有活跃任务的 `plan.md` 和 `progress.md`
2. 用"五问重启检查"验证对上下文的理解

### 接受新任务
1. 在 `.tasks/active/<task-name>/` 创建三个文件：`plan.md`、`findings.md`、`progress.md`
2. 与用户确认 `plan.md` 后开始执行

### 执行中
- 每完成一个阶段更新 `progress.md`
- 重要发现记入 `findings.md`
- 做重大决策前重新读取 `plan.md`

### 任务完成
- 将任务目录从 `active/` 移到 `done/`

### 文件格式

**plan.md** 包含：目标、阶段清单（含检查项和状态）、关键问题、已做决策。

**progress.md** 包含：每个会话的操作记录、创建/修改的文件、"五问重启检查"表。

**findings.md** 包含：需求、研究发现、技术决策、资源列表。

### 五问重启检查

| 问题 | 说明 |
|------|------|
| 我在哪里？ | 当前阶段 |
| 我要去哪里？ | 下一阶段 |
| 目标是什么？ | 任务目标 |
| 我学到了什么？ | 见 findings.md |
| 我做了什么？ | 见 progress.md |

## 编码规范

- 配置路径统一使用 `configs/config.yaml`，不硬编码绝对路径
- 实验数据放 `experiments/<系统>/<轮次>/`，分析代码放 `analysis/`
- 新增外部工具放 `third_party/`
- Pipeline 运行时输出（results/, foldx/, logs/ 等）由 .gitignore 管理，不提交
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "chore: write CLAUDE.md with project overview and cross-session work protocol"
```

---

### Task 16: Update README.md

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Update README.md directory structure and path references**

Rewrite the "目录结构" section to reflect new layout. Update all path references:

- Remove `results/`, `nano01/`, `molecules/`, `data/`, `tools/`, `artifacts/`, `Phase1-*.md`, `ROADMAP.md` from listing
- Add `docs/`, `.tasks/`, `third_party/foldx/`, `third_party/molecular_dynamics/`
- Update `experiments/` to show `1E62/R1/R2/R3` structure
- Update FoldX path reference (line 51): `tools/foldx/foldx` → `third_party/foldx/foldx`
- Update ROADMAP link: `[ROADMAP.md](ROADMAP.md)` → `[ROADMAP.md](docs/ROADMAP.md)`
- Remove `results/` directory section (lines 133-140)
- Update `artifacts/` reference to `third_party/ml-simplefold/weights/`

- [ ] **Step 2: Verify no stale paths in README**

```bash
grep -n "data/pdbs\|data/wt\|tools/foldx\|artifacts/\|nano01\|molecules\|1E62_R[123]" README.md
```
Expected: No output.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "chore: update README.md for new directory structure"
```

---

### Task 17: Final verification

- [ ] **Step 1: Check root directory item count**

```bash
ls /public/home/ziyang/code/optim-pipe/ | wc -l
```
Expected: ~11 visible items (README.md, CLAUDE.md, run_pipeline.sh, configs, scripts, analysis, experiments, docs, third_party, archive, workspace).

- [ ] **Step 2: Scan for any remaining broken path references**

```bash
grep -rn "data/pdbs\|data/wt\|tools/foldx\|/artifacts/" scripts/ configs/ analysis/ run_pipeline.sh experiments/1E62/ --include="*.py" --include="*.sh" --include="*.yaml" 2>/dev/null
```
Expected: No output.

- [ ] **Step 3: Verify FoldX binary is accessible**

```bash
third_party/foldx/foldx --help 2>&1 | head -3
```
Expected: FoldX help output.

- [ ] **Step 4: Verify SimpleFold weights**

```bash
ls -lh third_party/ml-simplefold/weights/*.ckpt | awk '{print $5, $NF}'
```
Expected: 4 files with correct sizes (11G, 6.0G, 442M, 442M).

- [ ] **Step 5: Validate all config files**

```bash
for f in configs/*.yaml; do
  python -c "import yaml; yaml.safe_load(open('$f'))" && echo "OK: $f" || echo "FAIL: $f"
done
```
Expected: All print `OK`.

- [ ] **Step 6: Verify task management structure**

```bash
ls .tasks/done/sdab-pka/
ls .tasks/active/
```
Expected: `plan.md findings.md progress.md` in done/sdab-pka/; active/ is empty.

- [ ] **Step 7: Verify CLAUDE.md is populated**

```bash
wc -l CLAUDE.md
```
Expected: ~80+ lines (non-empty).
