# sdAb active learning

This package runs a local, replayable active-learning loop for sdAb pH-switch design. It expects wet-lab labels and precomputed ESMC embeddings as input CSV files. It does not run ESMC, submit wet-lab nominations, or write generated experiment outputs into git.

The model target is the pH-switch utility:

```text
utility_recomputed = A7 * (1 - A6)
```

`run_round.py` recomputes this value from the current config. It ignores any stale `utility_input` column.

## Training rows

The main model trains only on rows where:

```text
experiment_status == "measured"
label_is_direct_measurement == True
usable_for_training == True
```

Legacy summaries, pseudo labels, and unknown rows can stay in the table for provenance or diagnostics, but they do not enter production round training. The warm-start adapter writes legacy R4 summary rows with `experiment_status=excluded`, `label_source_type=legacy_training_summary`, and `usable_for_training=False`.

## Embeddings

Embedding CSV files must contain:

```text
sequence_id
sequence_sha256
emb_0
emb_1
...
```

The join checks both `sequence_id` and `sequence_sha256`. Candidates without embeddings are marked `missing_embedding` during `run_round.py` and are not scored.

## Round artifacts

One round writes:

```text
training_data_round_<round>.csv
candidate_pool_round_<round>.csv
model_config_round_<round>.json
round_provenance_round_<round>.json
prediction_table_round_<round>.csv
selected_candidates_round_<round>.csv
```

`model_config_round_<round>.json` stores only model settings and the scikit-learn version. Input hashes, WT sequence hash, embedding metadata, objective config, acquisition config, random seed, timestamp, and git commit go into `round_provenance_round_<round>.json`.

`selected_candidates_round_<round>.csv` is a reproducible dry-run selection artifact. Treat it as a wet-lab nomination list only after it was generated from confirmed labels, production ESMC embeddings, and explicit approval.

Structural metrics are review context only in this implementation. They are not model features and are not hard filters.

## Example commands

Build a warm-start training table:

```bash
python -m active_learning build-sdab-warm-start \
  --config-yaml active_learning/example_config.yaml \
  --wt-fasta experiments/sdab/R2/data/sdab.fasta \
  --legacy-r4-training experiments/sdab/R4/data/training_data.csv \
  --out-csv experiments/sdab/R5/active_learning/warm_start_training_rows_R5.csv
```

Add `--direct-wet-labels <confirmed_direct_wet_labels.csv>` when confirmed pH-specific labels are available. A legacy-only warm-start file is useful for provenance, but it will fail production selection because legacy summary rows are excluded from training.

Build a full single-mutant candidate pool:

```bash
python -m active_learning build-candidate-pool \
  --config-yaml active_learning/example_config.yaml \
  --wt-fasta experiments/sdab/R2/data/sdab.fasta \
  --system-id sdab \
  --round-id R5 \
  --candidate-mode full_single \
  --out-csv experiments/sdab/R5/active_learning/candidate_pool_input_R5.csv
```

Run one dry-run selection round:

```bash
python -m active_learning run-round \
  --config-yaml active_learning/example_config.yaml \
  --training-csv experiments/sdab/R5/active_learning/warm_start_training_rows_R5.csv \
  --candidate-csv experiments/sdab/R5/active_learning/candidate_pool_input_R5.csv \
  --embedding-csv <precomputed_esmc_embeddings.csv> \
  --round-id R5 \
  --top-k 24 \
  --out-dir experiments/sdab/R5/active_learning
```

Optional status inputs can be added with `--failed-synthesis-csv`, `--invalid-candidates-csv`, and `--previously-selected-csv`. When present, `run_round.py` uses them with the training table to compute pre-embedding availability before it checks for missing embeddings.

Run a small retrospective replay from Python. Use a table with direct measured rows from at least two rounds; legacy-only warm-start rows are excluded and will not form a heldout set.

```python
import pandas as pd

from active_learning.tools import load_embedding_table
from active_learning.replay import leave_round_out_enrichment

training = pd.read_csv("experiments/sdab/replay/direct_training_rows_by_round.csv")
embeddings = load_embedding_table("experiments/sdab/R5/active_learning/esmc_embeddings.csv")

result = leave_round_out_enrichment(
    training,
    embeddings,
    heldout_round="R4",
    top_k=24,
)
print(result)
```

Run the active-learning tests:

```bash
conda run -p /data/ziyang/mamba/envs/optim-pipe python -m pytest active_learning/tests -v
```
