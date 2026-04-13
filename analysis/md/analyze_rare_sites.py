#!/usr/bin/env python3
"""非共有位点功能分析: 低频突变 (<5/15) 对单个变体的影响"""

import re
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
CANDIDATES_CSV = REPO / "workspace/sdab_r2_v2/results/final_candidates.csv"
MD_DIR = REPO / "experiments/sdab_variants/md"
OUTPUT_DIR = REPO / "experiments/sdab_variants/md/analysis_plots"

CDR = {
    "FR1": (1, 25), "CDR1": (26, 35), "FR2": (36, 49), "CDR2": (50, 65),
    "FR3": (66, 94), "CDR3": (95, 112), "FR4": (113, 122),
}


def get_region(pos):
    for name, (s, e) in CDR.items():
        if s <= pos <= e:
            return name
    return "?"


def parse_mutations(mut_str):
    muts = {}
    for m in str(mut_str).split(","):
        match = re.match(r"([A-Z])A(\d+)([A-Z])", m.strip())
        if match:
            muts[int(match.group(2))] = (match.group(1), match.group(3))
    return muts


def load_rmsf(variant_id, ph):
    path = MD_DIR / variant_id / f"pH_{ph}" / "analysis" / "rmsf_per_residue.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    df_a = df[df["chain"].str.contains("chain_A", na=False)]
    return dict(zip(df_a["resid"].astype(int), df_a["rmsf"]))


def load_contacts(variant_id, ph):
    path = MD_DIR / variant_id / f"pH_{ph}" / "analysis" / "contacts_frequency.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    ab_contacts = df[df["ab_chain"].str.contains("chain_A", na=False)].groupby("ab_resid")["frequency"].sum()
    return ab_contacts.to_dict()


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cand = pd.read_csv(CANDIDATES_CSV)
    variants = [f"sdab_v2_{i:04d}" for i in range(15)]

    # Count mutation frequency across all variants
    all_positions = Counter()
    variant_mut_map = {}
    for i, (_, row) in enumerate(cand.iterrows()):
        vid = variants[i]
        muts = parse_mutations(row["mutations"])
        variant_mut_map[vid] = muts
        for pos in muts:
            all_positions[pos] += 1

    # Rare positions: frequency < 5/15
    rare_positions = {pos: cnt for pos, cnt in all_positions.items() if cnt < 5}

    # ---- Per-variant rare site analysis ----
    rows = []
    for vid in variants:
        rmsf_74 = load_rmsf(vid, "7.4")
        rmsf_60 = load_rmsf(vid, "6.0")
        contacts_74 = load_contacts(vid, "7.4")
        contacts_60 = load_contacts(vid, "6.0")

        for pos, (wt, mut) in variant_mut_map[vid].items():
            freq = all_positions[pos]
            if freq >= 5:
                continue  # only rare sites

            region = get_region(pos)
            is_his = mut == "H"

            # ΔRMSF at site
            d_rmsf = None
            if pos in rmsf_74 and pos in rmsf_60:
                d_rmsf = rmsf_60[pos] - rmsf_74[pos]

            # ΔRMSF in ±3 window
            d_rmsf_window = []
            for p in range(pos - 3, pos + 4):
                if p in rmsf_74 and p in rmsf_60:
                    d_rmsf_window.append(rmsf_60[p] - rmsf_74[p])
            d_rmsf_window_mean = np.mean(d_rmsf_window) if d_rmsf_window else None

            # Contact change
            c74 = contacts_74.get(pos, 0)
            c60 = contacts_60.get(pos, 0)
            d_contact = c60 - c74 if (c74 > 0 or c60 > 0) else None
            is_interface = c74 > 0.1 or c60 > 0.1

            rows.append({
                "variant_id": vid,
                "position": pos,
                "wt_aa": wt,
                "mut_aa": mut,
                "mutation": f"{wt}{pos}{mut}",
                "region": region,
                "frequency": f"{freq}/15",
                "is_his": is_his,
                "is_interface": is_interface,
                "delta_rmsf": round(d_rmsf, 3) if d_rmsf is not None else None,
                "delta_rmsf_window3": round(d_rmsf_window_mean, 3) if d_rmsf_window_mean is not None else None,
                "contact_74": round(c74, 3),
                "contact_60": round(c60, 3),
                "delta_contact": round(d_contact, 3) if d_contact is not None else None,
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_DIR / "rare_site_analysis.csv", index=False)
    print(f"Rare site analysis: {len(df)} entries saved")

    # ---- Case analysis: most significant rare mutations ----
    print("\n=== Notable Rare Mutations ===")

    # By |ΔRMSF|
    if df["delta_rmsf"].notna().any():
        df_rmsf = df.dropna(subset=["delta_rmsf"]).sort_values("delta_rmsf", key=abs, ascending=False)
        print("\nTop 10 by |ΔRMSF|:")
        for _, r in df_rmsf.head(10).iterrows():
            print(f"  {r['variant_id']}  {r['mutation']:>6s}  {r['region']}  freq={r['frequency']}  "
                  f"ΔRMSF={r['delta_rmsf']:+.2f}  interface={'Y' if r['is_interface'] else 'N'}  "
                  f"His={'Y' if r['is_his'] else 'N'}")

    # By |Δcontact|
    if df["delta_contact"].notna().any():
        df_contact = df.dropna(subset=["delta_contact"]).sort_values("delta_contact", key=abs, ascending=False)
        print("\nTop 10 by |Δcontact|:")
        for _, r in df_contact.head(10).iterrows():
            print(f"  {r['variant_id']}  {r['mutation']:>6s}  {r['region']}  freq={r['frequency']}  "
                  f"Δcontact={r['delta_contact']:+.3f}  c7.4={r['contact_74']:.2f}  c6.0={r['contact_60']:.2f}")

    # Special attention: rare His mutations
    his_rare = df[df["is_his"]]
    if not his_rare.empty:
        print("\n=== Rare His Mutations ===")
        for _, r in his_rare.iterrows():
            print(f"  {r['variant_id']}  {r['mutation']:>6s}  {r['region']}  freq={r['frequency']}  "
                  f"ΔRMSF={r['delta_rmsf']}  Δcontact={r['delta_contact']}  interface={r['is_interface']}")

    print(f"\nAll outputs in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
