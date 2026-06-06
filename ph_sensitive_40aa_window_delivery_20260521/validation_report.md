# Validation Report

- PASS: 42
- WARN: 0
- FAIL: 0

| Category | Check | Status | Detail |
|---|---|---:|---|
| blocker_dry_run | 1E62:capacity_failure | PASS | legal_windows=150, min_capacity=10000, expected_target_blocker=False |
| blocker_dry_run | 1E62:insufficient_AF3_models | PASS | qualified_AF3_models=49, threshold=15 |
| blocker_dry_run | sdAb:capacity_failure | PASS | legal_windows=83, min_capacity=10000, expected_target_blocker=False |
| blocker_dry_run | sdAb:insufficient_AF3_models | PASS | qualified_AF3_models=73, threshold=15 |
| blocker_dry_run | window_no_prior_is_not_missing_wet_data | PASS | no window-level unknown prior triggers missing_wet_data |
| context_only | context_constraints_are_non_scoring | PASS | context_only constraints do not request scoring use |
| context_only | legacy_md_inventory_1e62.csv | PASS | allowed_usage_values=['context_only', 'mechanism_context_only'] |
| context_only | legacy_md_inventory_sdab.csv | PASS | allowed_usage_values=['context_only', 'mechanism_context_only'] |
| context_only | score_fields_ignore_context_only | PASS | context_only/legacy terms absent from scoring fields |
| fasta_length | 1E62:AeS | PASS | length=226, expected=226 |
| fasta_length | 1E62:heavy | PASS | length=115, expected=115 |
| fasta_length | 1E62:light | PASS | length=113, expected=113 |
| fasta_length | sdAb:AeS | PASS | length=226, expected=226 |
| fasta_length | sdAb:sdab | PASS | length=122, expected=122 |
| schema | candidate_windows.csv | PASS | all required columns present |
| schema | hbsag_aes_regions.csv | PASS | all required columns present |
| schema | prior_constraints_table.csv | PASS | all required columns present |
| schema | reference_sequence_map.csv | PASS | all required columns present |
| schema | wet_observation_table_1e62.csv | PASS | all required columns present |
| schema | wet_observation_table_sdab.csv | PASS | all required columns present |
| schema | window_design_capacity_table.csv | PASS | all required columns present |
| schema | window_mutation_mask.csv | PASS | all required columns present |
| schema | window_scores.csv | PASS | all required columns present |
| source_policy | 1E62:experiments/1E62/R1/1e62_R1_kd.csv | PASS | allowed scoring source exists |
| source_policy | 1E62:experiments/1E62/R2/wet_lab/elisa_summary.csv | PASS | allowed scoring source exists |
| source_policy | 1E62:experiments/1E62/R2/wet_lab/expression.csv | PASS | allowed scoring source exists |
| source_policy | 1E62:experiments/1E62/R2/wet_lab/variant_sequences.csv | PASS | allowed scoring source exists |
| source_policy | prior_constraints_non_AeS_hard_filters | PASS | all hard filters are AeS-applicable |
| source_policy | prior_constraints_table.csv | PASS | no excluded sources in scoring inputs |
| source_policy | sdAb:experiments/sdab/R2/data/hs32-92_8ng_ml_Elisa_results.csv | PASS | allowed scoring source exists |
| source_policy | sdAb:experiments/sdab/R2/data/sdab_elisa_raw.csv | PASS | allowed scoring source exists |
| source_policy | sdAb:experiments/sdab/R2/data/sdab_variants.csv | PASS | allowed scoring source exists |
| source_policy | sdAb:experiments/sdab/R4/data/training_data.csv | PASS | allowed scoring source exists |
| source_policy | wet_observation_table_1e62.csv | PASS | no excluded sources in scoring inputs |
| source_policy | wet_observation_table_sdab.csv | PASS | no excluded sources in scoring inputs |
| wet_source_rows | 1E62 R1 KD | PASS | rows=52, expected=52 |
| wet_source_rows | 1E62 R2 ELISA | PASS | rows=20, expected=20 |
| wet_source_rows | 1E62 R2 expression | PASS | rows=20, expected=20 |
| wet_source_rows | sdAb R2 combo | PASS | rows=61, expected=61 |
| wet_source_rows | sdAb R2 raw ELISA | PASS | rows=576, expected=576 |
| wet_source_rows | sdAb R2 single | PASS | rows=31, expected=31 |
| wet_source_rows | sdAb R4 training | PASS | rows=93, expected=93 |
