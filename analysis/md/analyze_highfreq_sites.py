#!/usr/bin/env python3
"""高频突变位点 pH 敏感性分析 + 序列相似度"""

import json
import re
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

# WT sequence (from WT_Repair.pdb, chain A, 122 residues)
WT_SEQ = None  # Will be loaded from PDB


def get_region(pos):
    for name, (s, e) in CDR.items():
        if s <= pos <= e:
            return name
    return "?"


def load_wt_seq():
    """Load WT sequence from PDB."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("wt", str(REPO / "experiments/sdab/repaired/WT_Repair.pdb"))
    from Bio.Data.IUPACData import protein_letters_3to1
    def three_to_one(name):
        return protein_letters_3to1.get(name.lower(), "X")
    seq = {}
    for res in s[0]["A"].get_residues():
        if res.id[0] != " ":
            continue
        try:
            seq[res.id[1]] = three_to_one(res.resname)
        except KeyError:
            seq[res.id[1]] = "X"
    return seq


def parse_mutations(mut_str):
    """Parse 'EA1S,QA3T,...' → dict {pos: (wt, mut)}"""
    muts = {}
    for m in str(mut_str).split(","):
        match = re.match(r"([A-Z])A(\d+)([A-Z])", m.strip())
        if match:
            muts[int(match.group(2))] = (match.group(1), match.group(3))
    return muts


def compute_identity(variant_muts, wt_seq):
    """Compute sequence identity with WT, total + per-region."""
    total = len(wt_seq)
    mutated = len(variant_muts)
    identity_total = (total - mutated) / total * 100

    region_identity = {}
    for region, (s, e) in CDR.items():
        region_len = sum(1 for p in wt_seq if s <= p <= e)
        region_muts = sum(1 for p in variant_muts if s <= p <= e)
        if region_len > 0:
            region_identity[region] = (region_len - region_muts) / region_len * 100
        else:
            region_identity[region] = 100.0

    cdr_len = sum(1 for p in wt_seq if any(s <= p <= e for name, (s, e) in CDR.items() if "CDR" in name))
    cdr_muts = sum(1 for p in variant_muts if any(s <= p <= e for name, (s, e) in CDR.items() if "CDR" in name))
    fr_len = total - cdr_len
    fr_muts = mutated - cdr_muts

    region_identity["CDR_all"] = (cdr_len - cdr_muts) / max(cdr_len, 1) * 100
    region_identity["FR_all"] = (fr_len - fr_muts) / max(fr_len, 1) * 100

    return identity_total, region_identity


def load_rmsf(variant_id, ph):
    path = MD_DIR / variant_id / f"pH_{ph}" / "analysis" / "rmsf_per_residue.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    df_a = df[df["chain"].str.contains("chain_A", na=False)]
    return dict(zip(df_a["resid"].astype(int), df_a["rmsf"]))


def load_contacts(variant_id, ph):
    """Load interface contact frequency for antibody residues."""
    path = MD_DIR / variant_id / f"pH_{ph}" / "analysis" / "contacts_frequency.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    # Sum contact frequency per antibody residue
    ab_contacts = df[df["ab_chain"].str.contains("chain_A", na=False)].groupby("ab_resid")["frequency"].sum()
    return ab_contacts.to_dict()


def load_summary(variant_id, ph):
    path = MD_DIR / variant_id / f"pH_{ph}" / "analysis" / "summary.json"
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cand = pd.read_csv(CANDIDATES_CSV)
    variants = [f"sdab_v2_{i:04d}" for i in range(15)]
    wt_seq = load_wt_seq()

    # ---- Mutation frequency analysis ----
    all_positions = Counter()
    variant_mut_map = {}

    for i, (_, row) in enumerate(cand.iterrows()):
        vid = variants[i]
        muts = parse_mutations(row["mutations"])
        variant_mut_map[vid] = muts
        for pos in muts:
            all_positions[pos] += 1

    # High-frequency positions (≥10/15)
    high_freq = {pos: cnt for pos, cnt in all_positions.items() if cnt >= 10}
    sorted_hf = sorted(high_freq.keys())

    # ---- Sequence identity ----
    identity_rows = []
    for i, vid in enumerate(variants):
        muts = variant_mut_map[vid]
        id_total, id_regions = compute_identity(muts, wt_seq)
        row = {
            "variant_id": vid,
            "n_mutations": len(muts),
            "identity_total": round(id_total, 1),
            "identity_CDR": round(id_regions["CDR_all"], 1),
            "identity_FR": round(id_regions["FR_all"], 1),
        }
        for region in CDR:
            row[f"identity_{region}"] = round(id_regions[region], 1)
        identity_rows.append(row)

    identity_df = pd.DataFrame(identity_rows)
    identity_df.to_csv(OUTPUT_DIR / "sequence_identity.csv", index=False)
    print("Sequence identity saved.")

    # ---- High-freq site ΔRMSF ----
    hf_rmsf_data = []

    for vid in variants:
        rmsf_74 = load_rmsf(vid, "7.4")
        rmsf_60 = load_rmsf(vid, "6.0")
        contacts_74 = load_contacts(vid, "7.4")
        contacts_60 = load_contacts(vid, "6.0")

        for pos in sorted_hf:
            has_mut = pos in variant_mut_map[vid]
            mut_type = f"{variant_mut_map[vid][pos][0]}{pos}{variant_mut_map[vid][pos][1]}" if has_mut else "-"

            d_rmsf = None
            if pos in rmsf_74 and pos in rmsf_60:
                d_rmsf = rmsf_60[pos] - rmsf_74[pos]

            d_contact = None
            if pos in contacts_74 or pos in contacts_60:
                c74 = contacts_74.get(pos, 0)
                c60 = contacts_60.get(pos, 0)
                d_contact = c60 - c74

            hf_rmsf_data.append({
                "variant_id": vid,
                "position": pos,
                "region": get_region(pos),
                "mutation": mut_type,
                "frequency": high_freq[pos],
                "has_mutation": has_mut,
                "delta_rmsf": d_rmsf,
                "delta_contact_freq": d_contact,
            })

    hf_df = pd.DataFrame(hf_rmsf_data)
    hf_df.to_csv(OUTPUT_DIR / "highfreq_site_analysis.csv", index=False)

    # ---- Summary table: per-position median ----
    summary_rows = []
    for pos in sorted_hf:
        sub = hf_df[(hf_df["position"] == pos) & (hf_df["has_mutation"])]
        summary_rows.append({
            "position": pos,
            "region": get_region(pos),
            "frequency": f"{high_freq[pos]}/15",
            "common_mutation": sub["mutation"].mode().values[0] if len(sub) > 0 else "-",
            "delta_rmsf_median": round(sub["delta_rmsf"].median(), 3) if sub["delta_rmsf"].notna().any() else None,
            "delta_rmsf_iqr": round(sub["delta_rmsf"].quantile(0.75) - sub["delta_rmsf"].quantile(0.25), 3) if sub["delta_rmsf"].notna().sum() > 2 else None,
            "delta_contact_median": round(sub["delta_contact_freq"].median(), 3) if sub["delta_contact_freq"].notna().any() else None,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTPUT_DIR / "highfreq_site_summary.csv", index=False)
    print("High-frequency site summary saved.")

    # ---- Boxplot: ΔRMSF per high-freq position ----
    # Only positions with enough data
    plot_positions = [p for p in sorted_hf if hf_df[(hf_df["position"] == p) & hf_df["delta_rmsf"].notna()].shape[0] >= 5]

    if plot_positions:
        fig, ax = plt.subplots(figsize=(max(14, len(plot_positions) * 0.5), 6))
        box_data = []
        labels = []
        colors = []
        for pos in plot_positions:
            vals = hf_df[(hf_df["position"] == pos) & hf_df["delta_rmsf"].notna()]["delta_rmsf"].values
            box_data.append(vals)
            region = get_region(pos)
            is_his = pos in [52, 55, 61, 99, 108, 109, 110]
            labels.append(f"{pos}\n{'*' if is_his else ''}{region}")
            colors.append("salmon" if is_his else "lightblue" if "CDR" in region else "lightgray")

        bp = ax.boxplot(box_data, patch_artist=True, labels=labels)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)

        ax.axhline(0, color="black", linestyle="--", linewidth=0.5)
        ax.set_ylabel("ΔRMSF (pH6.0 − pH7.4) [Å]", fontsize=11)
        ax.set_xlabel("Mutation Position", fontsize=11)
        ax.set_title("High-Frequency Sites: ΔRMSF Distribution (red=His, blue=CDR, gray=FR)", fontsize=11)
        plt.xticks(fontsize=7, rotation=90)
        plt.tight_layout()
        fig.savefig(OUTPUT_DIR / "highfreq_delta_rmsf_boxplot.png", dpi=200)
        plt.close()
        print("Boxplot saved.")

    # ---- Correlation analysis: presence of RA72K vs global metrics ----
    # Check if RA72K (10/15) correlates with better pH-switch
    print("\n=== RA72K (10/15) correlation ===")
    for vid in variants:
        s74 = load_summary(vid, "7.4")
        s60 = load_summary(vid, "6.0")
        has_r72k = 72 in variant_mut_map[vid]
        n_hb_74 = s74.get("hbond_n_hbond_mean", None)
        n_hb_60 = s60.get("hbond_n_hbond_mean", None)
        if n_hb_74 is not None and n_hb_60 is not None:
            d_hb = n_hb_60 - n_hb_74
            print(f"  {vid}: R72K={'Y' if has_r72k else 'N'}  ΔHbond={d_hb:+.1f}")

    print(f"\nAll outputs in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
