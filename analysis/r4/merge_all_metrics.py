#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Combine validation_metrics.csv (Step 4) + md_metrics.csv (MD deltas) into
all_metrics.csv so Step 5 has model + FoldX + Rosetta + MD signals side-by-side.

Inputs:
- experiments/sdab/R4/structures/validation_metrics.csv  (16 rows, from Step 4)
- experiments/sdab/R4/md_metrics.csv                     (16 rows, from analyze_and_compare_r4)

Output:
- experiments/sdab/R4/structures/all_metrics.csv
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

VAL = REPO_ROOT / "experiments/sdab/R4/structures/validation_metrics.csv"
MD = REPO_ROOT / "experiments/sdab/R4/md_metrics.csv"
OUT = REPO_ROOT / "experiments/sdab/R4/structures/all_metrics.csv"

# MD delta columns (from compare_ph.py output). Signs encoded for scoring:
# - delta_rmsd          : pH6 − pH7.4 of antibody RMSD → positive = pH6 more destabilized (good)
# - delta_rmsf_h1/h2/h3 : positive = pH6 more flexible in CDR (good for pH-selective dissociation)
# - delta_n_hbond       : negative = pH6 loses interface H-bonds (good)
# - delta_buried_sasa   : negative = pH6 has smaller buried SASA (good)
# - delta_n_contacts    : negative = pH6 has fewer interface contacts (good)
MD_DELTA_COLS = [
    "delta_rmsd",
    "delta_rmsf_global",
    "delta_rmsf_h1",
    "delta_rmsf_h2",
    "delta_rmsf_h3",
    "delta_n_hbond",
    "delta_buried_sasa",
    "delta_n_contacts",
]


def main():
    if not VAL.exists():
        sys.exit(f"[err] missing {VAL}")
    if not MD.exists():
        sys.exit(f"[err] missing {MD} (run analyze_and_compare_r4.py first)")

    val = pd.read_csv(VAL)
    md = pd.read_csv(MD)

    md_keep = ["name"] + [c for c in MD_DELTA_COLS if c in md.columns]
    if "status" in md.columns:
        md_keep.append("status")
        md = md.rename(columns={"status": "md_status"})
        md_keep = ["md_status" if c == "status" else c for c in md_keep]

    out = val.merge(md[md_keep + (["md_status"] if "md_status" in md.columns else [])],
                    on="name", how="left")
    out.to_csv(OUT, index=False)

    print(f"[OK] wrote {OUT} ({len(out)} rows)")
    # Cross-method rank agreement on pH selectivity signals
    cols = ["model_score", "delta_delta", "dddG_elec", "delta_rmsd",
            "delta_n_hbond", "delta_buried_sasa", "delta_n_contacts"]
    cols = [c for c in cols if c in out.columns]
    ok = out.dropna(subset=cols)
    if len(ok) >= 3:
        print("\nSpearman rank correlation (all rows with complete metrics):")
        ref = "delta_delta"
        for c in cols:
            if c == ref:
                continue
            rho, p = spearmanr(ok[ref], ok[c])
            direction = "same-direction" if rho > 0 else "opposite"
            print(f"  {ref} vs {c:22s}: rho={rho:+.3f} (p={p:.3g})  [{direction}]")

    print()
    preview_cols = ["name", "track", "mutations_unified", "model_score",
                    "delta_delta", "dddG_elec", "delta_rmsd", "delta_n_hbond",
                    "delta_n_contacts"]
    preview_cols = [c for c in preview_cols if c in out.columns]
    print(out.sort_values("delta_delta", ascending=False)[preview_cols].to_string(index=False))


if __name__ == "__main__":
    main()
