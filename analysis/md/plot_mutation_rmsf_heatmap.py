#!/usr/bin/env python3
"""突变位点 ΔRMSF 热图: X=突变位点, Y=15 变体, 颜色=ΔRMSF(pH6.0 - pH7.4)"""

import re
import sys
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


def get_region(pos):
    for name, (s, e) in CDR.items():
        if s <= pos <= e:
            return name
    return "?"


def parse_mutations(mut_str):
    """Parse 'EA1S,QA3T,...' → list of (pos, wt, mut)"""
    muts = []
    for m in str(mut_str).split(","):
        match = re.match(r"([A-Z])A(\d+)([A-Z])", m.strip())
        if match:
            muts.append((int(match.group(2)), match.group(1), match.group(3)))
    return muts


def load_rmsf(variant_id, ph):
    """Load rmsf_per_residue.csv for a variant/pH, return {resid: rmsf}"""
    path = MD_DIR / variant_id / f"pH_{ph}" / "analysis" / "rmsf_per_residue.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    # Only chain A (antibody) — handle GROMACS segid format
    df_a = df[df["chain"].str.contains("chain_A", na=False)]
    return dict(zip(df_a["resid"].astype(int), df_a["rmsf"]))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load candidates
    cand = pd.read_csv(CANDIDATES_CSV)
    variants = [f"sdab_v2_{i:04d}" for i in range(15)]

    # Parse all mutations, collect all positions
    variant_muts = {}  # variant_id → set of positions
    all_positions = set()
    for i, (_, row) in enumerate(cand.iterrows()):
        vid = variants[i]
        muts = parse_mutations(row["mutations"])
        positions = {pos for pos, wt, mut in muts}
        variant_muts[vid] = positions
        all_positions |= positions

    sorted_positions = sorted(all_positions)
    n_pos = len(sorted_positions)
    n_var = len(variants)
    pos_to_idx = {p: i for i, p in enumerate(sorted_positions)}

    # Build ΔRMSF matrix
    data = np.full((n_var, n_pos), np.nan)
    for vi, vid in enumerate(variants):
        rmsf_74 = load_rmsf(vid, "7.4")
        rmsf_60 = load_rmsf(vid, "6.0")
        if not rmsf_74 or not rmsf_60:
            print(f"  WARN: {vid} missing RMSF data")
            continue
        for pos in variant_muts[vid]:
            if pos in rmsf_74 and pos in rmsf_60:
                pi = pos_to_idx[pos]
                data[vi, pi] = rmsf_60[pos] - rmsf_74[pos]

    # Save CSV
    csv_df = pd.DataFrame(data, index=variants, columns=[str(p) for p in sorted_positions])
    csv_df.to_csv(OUTPUT_DIR / "mutation_delta_rmsf.csv")

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(max(20, n_pos * 0.35), 8))

    # Mask NaN for display
    masked = np.ma.masked_invalid(data)
    vmax = max(abs(np.nanmin(data)), abs(np.nanmax(data)), 1.0)
    im = ax.imshow(masked, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   interpolation="nearest")

    # X ticks: position numbers with region annotation
    x_labels = []
    for p in sorted_positions:
        region = get_region(p)
        x_labels.append(f"{p}\n{region}")
    ax.set_xticks(range(n_pos))
    ax.set_xticklabels(x_labels, fontsize=6, rotation=90)
    ax.set_yticks(range(n_var))
    ax.set_yticklabels(variants, fontsize=8)

    # CDR region color bands at top
    for name, (s, e) in CDR.items():
        idxs = [pos_to_idx[p] for p in sorted_positions if s <= p <= e]
        if idxs and "CDR" in name:
            ax.axvspan(min(idxs) - 0.5, max(idxs) + 0.5, alpha=0.08, color="orange")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("ΔRMSF (pH6.0 − pH7.4) [Å]", fontsize=10)

    ax.set_xlabel("Mutation Position (A chain)", fontsize=11)
    ax.set_ylabel("Variant", fontsize=11)
    ax.set_title("Mutation Site ΔRMSF Heatmap (blue=pH6 more stable, red=pH6 more flexible)", fontsize=12)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "mutation_delta_rmsf_heatmap.png", dpi=200)
    plt.close()

    print(f"Heatmap saved: {OUTPUT_DIR / 'mutation_delta_rmsf_heatmap.png'}")
    print(f"Data saved: {OUTPUT_DIR / 'mutation_delta_rmsf.csv'}")


if __name__ == "__main__":
    main()
