#!/usr/bin/env python3
"""综合排名: 整合 FoldX + MD + 序列相似度"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
MD_DIR = REPO / "experiments/sdab_variants/md"
CAND_CSV = REPO / "workspace/sdab_r2_v2/results/final_candidates.csv"
PLOTS = MD_DIR / "analysis_plots"

variants = [f"sdab_v2_{i:04d}" for i in range(15)]


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(CAND_CSV)

    rows = []
    for i, vid in enumerate(variants):
        r = cand.iloc[i]

        s74_path = MD_DIR / vid / "pH_7.4" / "analysis" / "summary.json"
        s60_path = MD_DIR / vid / "pH_6.0" / "analysis" / "summary.json"
        comp_path = MD_DIR / vid / "ph_comparison.json"

        s74 = json.load(open(s74_path)) if s74_path.exists() else {}
        s60 = json.load(open(s60_path)) if s60_path.exists() else {}
        comp = json.load(open(comp_path)) if comp_path.exists() else {}
        deltas = comp.get("deltas", {})

        row = {
            "variant_id": vid,
            "source": r.get("source", ""),
            "delta_foldx": r.get("delta", None),
            "dG_pH7_4": r.get("dG_pH7_4", None),
            "dddG_elec": r.get("dddG_elec", None),
            # pH 7.4 stability
            "rmsd_mean_74": s74.get("rmsd_rmsd_mean", None),
            "rmsd_converge_74": s74.get("rmsd_converge_ns", None),
            # pH-switch deltas (pH6.0 - pH7.4)
            "d_contacts": deltas.get("delta_n_contacts", None),
            "d_hbond": deltas.get("delta_n_hbond", None),
            "d_rmsd": deltas.get("delta_rmsd", None),
            "d_rmsf": deltas.get("delta_rmsf_global", None),
            # pairwise His occupancy pH 7.4
            "His99_occ_74": s74.get("pairwise_dist_His99-THR123_occupancy_4.5A", None),
            "His110_occ_74": s74.get("pairwise_dist_His110-LYS122_occupancy_4.5A", None),
            "His61_occ_74": s74.get("pairwise_dist_His61-THR118_occupancy_4.5A", None),
            "ARG72_His99_74": s74.get("pairwise_dist_ARG72-His99_occupancy_4.5A", None),
            # pairwise His occupancy pH 6.0
            "His99_occ_60": s60.get("pairwise_dist_His99-THR123_occupancy_4.5A", None),
            "His110_occ_60": s60.get("pairwise_dist_His110-LYS122_occupancy_4.5A", None),
            "His61_occ_60": s60.get("pairwise_dist_His61-THR118_occupancy_4.5A", None),
            # absolute values
            "contacts_74": s74.get("contacts_n_contacts_mean", None),
            "contacts_60": s60.get("contacts_n_contacts_mean", None),
            "hbond_74": s74.get("hbond_n_hbond_mean", None),
            "hbond_60": s60.get("hbond_n_hbond_mean", None),
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    # Merge sequence identity
    id_df = pd.read_csv(PLOTS / "sequence_identity.csv")
    df = df.merge(
        id_df[["variant_id", "n_mutations", "identity_total", "identity_CDR", "identity_FR"]],
        on="variant_id", how="left",
    )

    # ---- Composite ranking ----
    # d_contacts < 0 = pH6 接触减少 = 好的 pH-switch → rank ascending
    # d_hbond < 0 = pH6 氢键减少 → rank ascending
    # rmsd_mean_74 low = pH7.4 结构稳定 → rank ascending
    rank_cols = []
    for col, ascending in [("d_contacts", True), ("d_hbond", True), ("rmsd_mean_74", True)]:
        rname = f"rank_{col}"
        if df[col].notna().any():
            df[rname] = df[col].rank(ascending=ascending)
            rank_cols.append(rname)

    df["avg_rank"] = df[rank_cols].mean(axis=1)
    df["final_rank"] = df["avg_rank"].rank(ascending=True).astype(int)
    df = df.sort_values("final_rank")

    df.to_csv(PLOTS / "final_ranking.csv", index=False)

    # ---- Print summary ----
    sep = "=" * 140
    print(sep)
    hdr = (f"{'Rank':>4}  {'Variant':<14}  {'Src':<13}  {'dG74':>6}  {'delta':>6}  "
           f"{'RMSD74':>7}  {'dCont':>7}  {'dHbond':>7}  "
           f"{'H99@74':>6}  {'H110@74':>6}  {'H61@74':>6}  "
           f"{'ID':>5}  {'nMut':>4}")
    print(hdr)
    print(sep)
    for _, r in df.iterrows():
        def fmt(v, f_str):
            return f_str % v if pd.notna(v) and v is not None else "   N/A"
        h61 = fmt(r.His61_occ_74, "%.3f")
        print(
            f"{int(r.final_rank):4d}  {r.variant_id:<14}  {str(r.source):<13}  "
            f"{r.dG_pH7_4:6.2f}  {r.delta_foldx:6.3f}  "
            f"{fmt(r.rmsd_mean_74, '%7.2f')}  "
            f"{fmt(r.d_contacts, '%+7.1f')}  {fmt(r.d_hbond, '%+7.1f')}  "
            f"{fmt(r.His99_occ_74, '%6.3f')}  {fmt(r.His110_occ_74, '%6.3f')}  "
            f"{h61}  "
            f"{r.identity_total:5.1f}  {int(r.n_mutations):4d}"
        )
    print(sep)

    # pH-switch direction
    good = df[df["d_contacts"] < 0]
    bad = df[df["d_contacts"] >= 0]
    print()
    print(f"pH-switch 方向 (d_contacts < 0 = pH6 接触减少 = 好):")
    print(f"  正确: {len(good)}/15")
    print(f"  异常: {len(bad)}/15")

    # His occupancy change
    print()
    print("His-抗原距离占有率变化 (pH7.4 → pH6.0):")
    for _, r in df.iterrows():
        d99 = (r.His99_occ_60 or 0) - (r.His99_occ_74 or 0)
        d110 = (r.His110_occ_60 or 0) - (r.His110_occ_74 or 0)
        print(f"  {r.variant_id}  His99: {d99:+.3f}  His110: {d110:+.3f}")

    print(f"\n结果: {PLOTS / 'final_ranking.csv'}")


if __name__ == "__main__":
    main()
