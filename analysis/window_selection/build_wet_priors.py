#!/usr/bin/env python
from __future__ import annotations

import pandas as pd

from analysis.window_selection.common import (
    diff_mutations,
    load_yaml,
    make_arg_parser,
    mutation_positions,
    read_fasta,
    repo_path,
    split_list,
    write_csv,
)


def row_base(target: str, variant_id: str, assay_type: str, endpoint: str, value, source_file: str) -> dict:
    return {
        "target": target,
        "variant_id": variant_id,
        "assay_type": assay_type,
        "endpoint": endpoint,
        "observation_value": value,
        "derived_endpoint": "",
        "normalization_status": "raw",
        "direction": "higher_is_better",
        "chain": "",
        "mutations": "",
        "positions": "",
        "source_file": source_file,
        "source_antigen": "AeS",
        "use_scope": "window_scoring",
    }


def build_1e62() -> pd.DataFrame:
    heavy = read_fasta("experiments/1E62/data/heavy.fasta")[2]
    light = read_fasta("experiments/1E62/data/light.fasta")[2]
    rows = []

    kd_path = "experiments/1E62/R1/1e62_R1_kd.csv"
    kd = pd.read_csv(repo_path(kd_path))
    for _, rec in kd.iterrows():
        h_muts = diff_mutations(heavy, rec["heavy_sequence"], "H")
        l_muts = diff_mutations(light, rec["light_sequence"], "L")
        muts = h_muts + l_muts
        chains = sorted({m[0] for m in muts})
        row = row_base("1E62", rec["variant_id"], "KD", "KD_ratio", rec["KD_ratio"], kd_path)
        row.update(
            {
                "normalization_status": "standardized",
                "direction": "lower_is_better",
                "chain": "|".join(chains),
                "mutations": "|".join(muts),
                "positions": "|".join(str(x) for x in mutation_positions(muts)),
            }
        )
        rows.append(row)

    elisa_path = "experiments/1E62/R2/wet_lab/elisa_summary.csv"
    elisa = pd.read_csv(repo_path(elisa_path))
    for _, rec in elisa.iterrows():
        muts = split_list(rec.get("mutations"))
        for endpoint in ["Ae_pH6.0", "Ae_pH7.4", "Ae_ratio"]:
            row = row_base("1E62", rec["variant"], "ELISA", endpoint, rec.get(endpoint), elisa_path)
            row.update(
                {
                    "derived_endpoint": "acid_retention" if endpoint == "Ae_ratio" else "",
                    "normalization_status": "derived" if endpoint == "Ae_ratio" else "raw",
                    "direction": "ratio_acid_over_neutral" if endpoint == "Ae_ratio" else "higher_is_better",
                    "chain": "|".join(sorted({m[0] for m in muts if m})),
                    "mutations": "|".join(muts),
                    "positions": "|".join(str(x) for x in mutation_positions(muts)),
                }
            )
            rows.append(row)

    exp_path = "experiments/1E62/R2/wet_lab/expression.csv"
    expr = pd.read_csv(repo_path(exp_path))
    for _, rec in expr.iterrows():
        row = row_base("1E62", rec["编号"], "expression", "concentration_ug_ul", rec["浓度ug/ul"], exp_path)
        rows.append(row)

    return pd.DataFrame(rows)


def build_sdab() -> pd.DataFrame:
    rows = []
    variants_path = "experiments/sdab/R2/data/sdab_variants.csv"
    variants = pd.read_csv(repo_path(variants_path))
    name_to_mut = {}
    for _, rec in variants.iterrows():
        mut = rec["mutant_id"]
        name_to_mut[rec["elisa_name"]] = mut
        row = row_base("sdAb", rec["elisa_name"], "single_variant_metadata", "expressed", rec.get("expressed"), variants_path)
        row.update(
            {
                "normalization_status": "raw",
                "chain": "A",
                "mutations": mut,
                "positions": str(int(rec["position"])),
            }
        )
        rows.append(row)

    raw_path = "experiments/sdab/R2/data/sdab_elisa_raw.csv"
    raw = pd.read_csv(repo_path(raw_path))
    for _, rec in raw.iterrows():
        mut = name_to_mut.get(rec["sample"], "")
        for endpoint, value in [("od_ph74", rec["od_ph74"]), ("od_ph60", rec["od_ph60"])]:
            row = row_base("sdAb", rec["sample"], "ELISA_raw", endpoint, value, raw_path)
            row.update(
                {
                    "chain": "A" if mut else "",
                    "mutations": mut,
                    "positions": "|".join(str(x) for x in mutation_positions([mut])),
                    "derived_endpoint": f"concentration_{rec['concentration_ng_ml']}_replicate_{rec['replicate']}",
                }
            )
            rows.append(row)

    combo_path = "experiments/sdab/R2/data/hs32-92_8ng_ml_Elisa_results.csv"
    combo = pd.read_csv(repo_path(combo_path))
    for _, rec in combo.iterrows():
        muts = split_list(rec.get("Mutations"))
        for endpoint in ["C-pH60_Avg", "C-pH74_avg", "ratio", "D-pH60_Avg", "D-pH74_avg", "比值＞1.5"]:
            row = row_base("sdAb", rec["ID"], "ELISA_combo_8ng_ml", endpoint, rec.get(endpoint), combo_path)
            row.update(
                {
                    "derived_endpoint": "combo_ratio" if endpoint in {"ratio", "比值＞1.5"} else "",
                    "normalization_status": "derived" if endpoint in {"ratio", "比值＞1.5"} else "raw",
                    "direction": "ratio_acid_over_neutral" if endpoint in {"ratio", "比值＞1.5"} else "higher_is_better",
                    "chain": "A",
                    "mutations": "|".join(muts),
                    "positions": "|".join(str(x) for x in mutation_positions(muts)),
                    "source_antigen": "AeS_or_table_label_unresolved",
                }
            )
            rows.append(row)

    r4_path = "experiments/sdab/R4/data/training_data.csv"
    r4 = pd.read_csv(repo_path(r4_path))
    for _, rec in r4.iterrows():
        muts = split_list(rec.get("mutations_unified"))
        for endpoint in ["log_pH74", "log_pH6", "log_ratio"]:
            row = row_base("sdAb", rec["name"], "R4_training", endpoint, rec.get(endpoint), r4_path)
            row.update(
                {
                    "derived_endpoint": endpoint,
                    "normalization_status": "standardized",
                    "direction": "ratio_acid_over_neutral" if endpoint == "log_ratio" else "higher_is_better",
                    "chain": "A",
                    "mutations": "|".join(muts),
                    "positions": "|".join(str(x) for x in mutation_positions(muts)),
                }
            )
            rows.append(row)

    return pd.DataFrame(rows)


def build_manual_constraints() -> pd.DataFrame:
    data = load_yaml("analysis/window_selection/config/manual_prior_constraints.yaml")
    rows = data.get("constraints", [])
    return pd.DataFrame(rows)


def main() -> None:
    make_arg_parser("Build wet observation tables and prior constraints.").parse_args()
    one = build_1e62()
    sdab = build_sdab()
    constraints = build_manual_constraints()
    write_csv(one, "wet_observation_table_1e62.csv")
    write_csv(sdab, "wet_observation_table_sdab.csv")
    write_csv(constraints, "prior_constraints_table.csv")
    print(f"Wrote 1E62 wet rows={len(one)}, sdAb wet rows={len(sdab)}, constraints={len(constraints)}.")


if __name__ == "__main__":
    main()
