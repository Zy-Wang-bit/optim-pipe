#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Merge Step 4 FoldX + Rosetta metrics with the R4 candidate mutlist.

Inputs:
- experiments/sdab/R4/structures/mutlist.csv            (16 rows)
- experiments/sdab/R4/foldx/batches/sdab/batch_000/foldx_summary.csv
- experiments/sdab/R4/tier2/rosetta/dddg_elec.csv       (16 rows)
- experiments/sdab/R4/tier2/rosetta/ph_score.csv        (16 rows)

FoldX summary rows are keyed by mutant PDB filename (e.g. sdab_Repair_1.pdb),
which maps to mutlist row N via the `Dif_build_pH*` ordering. We use the same
ordering the BuildModel writes: line i of individual_list.txt ↔ sdab_Repair_i.pdb.

Rosetta rows are keyed by variant_id = PDB stem (we will rename the PDBs to
`r4_NN` before running Rosetta, so variant_id == mutlist.name).

Output: experiments/sdab/R4/structures/validation_metrics.csv (16 rows).
"""

import re
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

MUTLIST = REPO_ROOT / "experiments/sdab/R4/structures/mutlist.csv"
FOLDX_SUMMARY = REPO_ROOT / "experiments/sdab/R4/foldx/batches/sdab/batch_000/foldx_summary.csv"
DDDG_CSV = REPO_ROOT / "experiments/sdab/R4/tier2/rosetta/dddg_elec.csv"
PHSCORE_CSV = REPO_ROOT / "experiments/sdab/R4/tier2/rosetta/ph_scores.csv"
OUT_CSV = REPO_ROOT / "experiments/sdab/R4/structures/validation_metrics.csv"

MPDB_RE = re.compile(r"^sdab_Repair_(\d+)\.pdb$")


def load_foldx(mutlist: pd.DataFrame) -> pd.DataFrame:
    fx = pd.read_csv(FOLDX_SUMMARY)
    # mpdb column example: sdab_Repair_1.pdb → index 1 (1-based line in individual_list.txt)
    idx = fx["mpdb"].apply(lambda s: int(MPDB_RE.match(s).group(1)))
    fx = fx.assign(mutlist_idx=idx - 1)  # 0-based into mutlist
    fx = fx.sort_values("mutlist_idx").reset_index(drop=True)
    fx["name"] = mutlist.iloc[fx["mutlist_idx"].values]["name"].values
    keep = ["name", "dG_pH7_4", "dG_pH6_0", "ddG_pH7_4", "ddG_pH6_0",
            "delta", "delta_wt", "delta_delta"]
    out = fx[keep].copy()
    out["foldx_status"] = out["delta_delta"].apply(
        lambda v: "ok" if pd.notna(v) and v != "" else "missing")
    return out


def load_rosetta() -> tuple[pd.DataFrame, pd.DataFrame]:
    dddg = pd.read_csv(DDDG_CSV)
    phs = pd.read_csv(PHSCORE_CSV)
    dddg = dddg.rename(columns={"variant_id": "name"})
    phs = phs.rename(columns={"variant_id": "name"})
    dddg = dddg[["name", "dddG_elec", "ddG_elec_pH7", "ddG_elec_pH5",
                 "n_his_binder", "status"]]
    dddg = dddg.rename(columns={"status": "rosetta_dddg_status"})
    phs = phs[["name", "ph_score", "status"]]
    phs = phs.rename(columns={"status": "rosetta_ph_status"})
    return dddg, phs


def main():
    mutlist = pd.read_csv(MUTLIST)
    fx = load_foldx(mutlist)
    dddg, phs = load_rosetta()

    df = mutlist.merge(fx, on="name", how="left") \
                .merge(dddg, on="name", how="left") \
                .merge(phs, on="name", how="left")

    df["rosetta_status"] = df.apply(
        lambda r: "ok" if r.get("rosetta_dddg_status") == "success" and
                  r.get("rosetta_ph_status") == "success" else
                  f"dddg={r.get('rosetta_dddg_status')};ph={r.get('rosetta_ph_status')}",
        axis=1,
    )

    cols = [
        "name", "track", "mutations_unified", "order", "contains_D110H", "confidence",
        "log_pH74_pred", "log_pH6_pred", "score",
        "ddG_pH7_4", "ddG_pH6_0", "delta", "delta_wt", "delta_delta",
        "dddG_elec", "ddG_elec_pH7", "ddG_elec_pH5", "ph_score", "n_his_binder",
        "foldx_status", "rosetta_status",
    ]
    df = df.rename(columns={"score": "model_score"})
    cols = [c if c != "score" else "model_score" for c in cols]
    cols = ["name", "track", "mutations_unified", "order", "contains_D110H", "confidence",
            "log_pH74_pred", "log_pH6_pred", "model_score",
            "ddG_pH7_4", "ddG_pH6_0", "delta", "delta_wt", "delta_delta",
            "dddG_elec", "ddG_elec_pH7", "ddG_elec_pH5", "ph_score", "n_his_binder",
            "foldx_status", "rosetta_status"]
    out = df[cols].sort_values("delta_delta", ascending=False).reset_index(drop=True)
    out.to_csv(OUT_CSV, index=False)

    print(f"[OK] wrote {OUT_CSV} ({len(out)} rows)")
    print(f"  foldx_status: {out['foldx_status'].value_counts().to_dict()}")
    print(f"  rosetta_status ok: {(out['rosetta_status']=='ok').sum()}/{len(out)}")
    print()
    print(out[["name", "track", "mutations_unified", "model_score",
               "delta_delta", "dddG_elec", "ph_score"]].to_string(index=False))

    # Rank consistency diagnostics (larger-is-better for all four)
    ok = out.dropna(subset=["delta_delta", "dddG_elec", "ph_score", "model_score"])
    if len(ok) >= 3:
        print("\nSpearman rank correlation across scoring sources:")
        for a, b in [("model_score", "delta_delta"),
                     ("model_score", "dddG_elec"),
                     ("delta_delta", "dddG_elec"),
                     ("delta_delta", "ph_score"),
                     ("dddG_elec", "ph_score")]:
            rho, p = spearmanr(ok[a], ok[b])
            print(f"  {a} vs {b}: rho={rho:+.3f} (p={p:.3g})")


if __name__ == "__main__":
    main()
