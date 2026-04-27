# optim-pipe Agent Guide

This file is the Codex entry point for the repository. The project was originally maintained with Claude Code; keep `CLAUDE.md` and this file aligned when changing shared agent workflow rules.

## Project Context

`optim-pipe` is an automated computational pipeline for pH-dependent antibody affinity optimization. The main workflow is a 3-tier funnel:

- Tier 1: high-throughput generation and filtering with ProteinMPNN, ESM, and FoldX.
- Tier 2: structure evaluation with PyRosetta and SimpleFold 3x.
- Tier 3: final ranking with `dddG_elec`.

The current biological focus is pH-sensitive engineering of anti-HBsAg antibody `1E62`, with active work around R2/R3/R4-style design and wet-lab feedback.

## Worktree Guardrails

- This repository often has local experimental data and unfinished analysis changes. Run `git status --short` before editing, and never revert or overwrite unrelated user changes.
- Keep runtime outputs out of commits. The ignored runtime directories include `results/`, `foldx/`, `phase_c/`, `tier2/`, `mpnn_outputs/`, `his_seeds/`, `esm_scores/`, `logs/`, `workspace/`, and `experiments/`.
- `third_party/` is generally external and ignored, except tracked code under `third_party/molecular_dynamics/`.
- Treat `.claude/`, `.tasks/`, and `.agents/` as local agent state unless the user explicitly asks to inspect or migrate them.

## Environment

Use conda environments rather than the empty `.venv/` at the repo root.

| Component | Environment | Use |
| --- | --- | --- |
| Main pipeline | `conda activate optim-pipe` | Core Python scripts, ESM, FoldX orchestration |
| ProteinMPNN | `conda activate proteinmpnn` | Sequence design |
| SimpleFold | `conda activate simplefold` | Tier 2 structure prediction |
| PyRosetta | `conda activate pyrosetta` | Mutant modeling and Rosetta scoring |

GROMACS for MD work:

```bash
module load gromacs/2024.2 openmpi/5.0.3
```

Slurm cluster notes from the Claude workflow:

- Jump/login host: `121.192.179.122`; internal nodes use `10.10.10.x`.
- Available GPU partitions: A40 `gpu[5-9]`, RTX3090 `gpu[14-17]`.
- Available CPU partition: AMD `computer[1-4]`.
- Do not run compute jobs on the login/management node. Use it only for light file operations, Git, `scp`, and similar control-plane work.
- If the user has already reserved a compute node, enter it directly with `ssh <node>`; the main cluster filesystem is shared across these nodes.
- Do not assume a compute node remains available. The user must explicitly tell you which node is available and how many cores can be used.
- CPU-heavy jobs such as FoldX and PyRosetta should have at least as many parallel tasks as available cores. For FoldX-like batches, estimate `ceil(candidate_count / chunk_size) * pH_count` and tune `chunk_size` so task count is near or slightly above the core count.

H100 standalone node notes from Claude memory:

- Connect with `ssh h100`; if the network is unstable, retry after `proxyon`.
- Hardware: 8 x H100 80GB HBM3.
- It is not part of the shared main-cluster filesystem. Sync code/data with `rsync -avz` or `scp`; its working path is `/home/ziyang/`, not `/public/home/ziyang/`.
- It has no Slurm/PBS scheduler. Run processes directly and assign GPUs with `CUDA_VISIBLE_DEVICES` when parallelizing.
- Conda environments are not guaranteed to match the main cluster; check with `conda env list` before using it.

## Common Commands

Run from the repository root unless a script explicitly says otherwise.

```bash
conda activate optim-pipe
bash run_pipeline.sh
bash run_pipeline.sh configs/config_foldx_n1-0.yaml
CLEAN=1 bash run_pipeline.sh
```

MD module:

```bash
cd third_party/molecular_dynamics
python run_md.py --pdb input.pdb --ph 7.4 6.0 --output-dir experiments/1E62/R3/md/
python analyze_trajectory.py --traj experiments/1E62/R3/md/HE1H/pH_7.4/
python compare_ph.py --variant-dir experiments/1E62/R3/md/HE1H/ --base-ph 7.4 --target-ph 6.0
```

When touching MD code, the focused test suite is:

```bash
python -m pytest third_party/molecular_dynamics/tests
```

## Task Management

Use `.tasks/active/<task-name>/` only for large work that is likely to span multiple conversations, has several implementation phases, or needs durable findings/progress notes.

Do not create task files for small one-session fixes, direct Q&A, or one-off config/doc edits. If a small task grows unexpectedly, add task files then.

For a new Codex session:

1. If the user says to continue a named task, read only that task's `plan.md`, `progress.md`, and relevant `findings.md`.
2. If no task is named, answer the user's request directly and do not proactively scan all of `.tasks/active/`.
3. When resuming a task, reconstruct: current phase, next phase, goal, findings, and completed work.

Within one session, keep one active task context. If the user introduces an unrelated large task, create or use a separate task directory instead of merging contexts.

## Agent Team Workflow

When the user explicitly authorizes sub-agents, delegation, or parallel agent work, prefer an agent team shape for work that is independent and parallelizable:

- The main session owns scope, decisions, integration, and the final user-facing answer.
- Parallel workers handle disjoint outputs such as separate plots, pages, modules, files, scripts, or data-cleaning batches.
- A reviewer agent can be used for independent checking when the output is large, cross-file, or user-facing.
- Do not use agent teams for tightly coupled, serial, or judgment-heavy work where the main session must inspect each step directly.

This is adapted from Claude memory's "Agent Team Default" rule, but Codex only delegates when the user has authorized agent usage.

## Coding Conventions

- Default configuration path is `configs/config.yaml`; avoid hard-coded absolute paths.
- Store experiment-round data under `experiments/<system>/<round>/` and analysis code under `analysis/`.
- Add new external tools under `third_party/`; keep project-owned wrappers or analysis modules outside ignored vendor payloads unless there is a clear reason.
- All pipeline scripts are expected to run from the repository root with relative paths such as `mpnn_outputs/` and `foldx/`.
- ProteinMPNN and SimpleFold are invoked as external subprocess tools and require their dedicated conda environments.
- Tier 2 Step 9a/9c require `pyrosetta`; Step 10 requires `simplefold`; other pipeline steps normally use `optim-pipe`.

Mutation naming appears in several historical formats, such as R1 `A_E1H`, R3 `H1E>H`, and FoldX `SA40A`. Final data files should use the normalized form such as `HE1H`; conversion logic lives in `analysis/naming/convert.py`.

## Documentation

- `README.md`: project overview, environment, directory map, and pipeline summary.
- `CLAUDE.md`: legacy Claude Code workflow notes; use it as source history when behavior diverges.
- `docs/ROADMAP.md`: higher-level project direction.
- `docs/reports/`: analysis reports.
- `docs/wet_lab_delivery/`: wet-lab delivery materials.
- `docs/slides/`: presentation materials.
