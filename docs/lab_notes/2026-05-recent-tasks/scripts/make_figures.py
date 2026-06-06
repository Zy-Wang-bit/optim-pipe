#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parents[1]
FIG = OUT / "figures"

COLORS = {
    "blue": "#5B7C99",
    "sage": "#7F8F84",
    "clay": "#B7A99A",
    "rose": "#B88C8C",
    "stone": "#8A9199",
    "paper": "#F7F5F1",
    "grid": "#D8D1C7",
}


def setup():
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 220,
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#5D6268",
            "axes.labelcolor": "#333333",
            "xtick.color": "#333333",
            "ytick.color": "#333333",
            "axes.grid": True,
            "grid.color": COLORS["grid"],
            "grid.linewidth": 0.6,
            "grid.alpha": 0.55,
        }
    )


def save(fig, name):
    path = FIG / name
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def plot_r2_heatmap():
    src = ROOT / ".tasks/active/hbsag-genotype-generalization/artifacts/r2_per_type_phenotype_table.csv"
    df = pd.read_csv(src)
    keep = ["com4", "com8", "com9", "com16", "com17", "com18", "com20"]
    df = df[df["variant"].isin(keep)].set_index("variant").loc[keep]
    values = df[["Ae_pH7.4", "B_pH7.4", "D1_pH7.4"]].to_numpy(float)

    fig, ax = plt.subplots(figsize=(5.8, 3.2))
    im = ax.imshow(values, cmap="YlGnBu", vmin=0, vmax=4.1, aspect="auto")
    ax.set_xticks(range(3), ["Ae", "B", "D1"])
    ax.set_yticks(range(len(keep)), keep)
    ax.set_title("R2 pH 7.4 ELISA signal")
    ax.set_xlabel("HBsAg genotype")
    ax.set_ylabel("Variant")
    ax.grid(False)

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center", fontsize=8, color="#222222")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("OD-like signal")
    return save(fig, "fig1_r2_elisa_key_variants.png")


def plot_antigen_map():
    positions = [110, 113, 114, 122, 126, 131, 134, 143, 159, 160, 161, 168]
    labels = {
        110: "Ce",
        113: "Ce",
        114: "B/Ce/D1",
        122: "D1",
        126: "Ce",
        131: "B/Ce/D1",
        134: "D1",
        143: "Ce/D1",
        159: "D1",
        160: "Ce",
        161: "Ce/D1",
        168: "D1",
    }
    regions = [
        (117, 128, "117-128"),
        (124, 147, "124-147"),
        (158, 168, "158-168"),
    ]

    label_offsets = {
        110: (-0.5, 0.00),
        113: (-1.1, 0.10),
        114: (1.1, -0.02),
        122: (0.0, 0.00),
        126: (0.0, 0.00),
        131: (0.0, 0.00),
        134: (0.0, 0.00),
        143: (0.0, 0.00),
        159: (-1.0, 0.10),
        160: (0.0, -0.02),
        161: (1.25, 0.10),
        168: (0.0, 0.00),
    }

    fig, ax = plt.subplots(figsize=(8.0, 2.8))
    ax.hlines(0, 99, 169, color=COLORS["stone"], linewidth=3)
    for start, end, name in regions:
        ax.axvspan(start, end, color=COLORS["clay"], alpha=0.22)
        ax.text((start + end) / 2, 0.78, name, ha="center", va="center", fontsize=9)
    for idx, pos in enumerate(positions):
        height = 0.25 + 0.12 * (idx % 3)
        color = COLORS["rose"] if "D1" in labels[pos] else COLORS["blue"]
        ax.vlines(pos, 0, height, color=color, linewidth=1.8)
        ax.scatter([pos], [height], s=42, color=color, zorder=3)
        dx, dy = label_offsets[pos]
        ax.text(pos + dx, height + 0.08 + dy, f"{pos}\n{labels[pos]}", ha="center", va="bottom", fontsize=7)

    ax.set_xlim(99, 169)
    ax.set_ylim(-0.25, 1.05)
    ax.set_yticks([])
    ax.set_xlabel("HBsAg S-protein position")
    ax.set_title("Type-sensitive positions around the antigenic loop")
    ax.grid(axis="x", alpha=0.35)
    return save(fig, "fig2_hbsag_type_sensitive_map.png")


def plot_foldx_auc():
    auc = pd.DataFrame(
        {
            "Genotype": ["Ae", "B", "D1"],
            "AUC": [0.552, 0.521, 0.594],
        }
    )
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    ax.bar(auc["Genotype"], auc["AUC"], color=[COLORS["blue"], COLORS["clay"], COLORS["rose"]], width=0.58)
    ax.axhline(0.5, color="#5D6268", linestyle="--", linewidth=1.0)
    ax.set_ylim(0.45, 0.65)
    ax.set_ylabel("ROC AUC")
    ax.set_title("FoldX genotype validation stayed weak")
    for i, row in auc.iterrows():
        ax.text(i, row["AUC"] + 0.008, f"{row['AUC']:.3f}", ha="center", fontsize=9)
    return save(fig, "fig3_foldx_auc_boundary.png")


def plot_af3_ces_ranks():
    src = ROOT / "experiments/1E62/R4/af3_ces_wt_corrected_light/summary/ranked_structures.tsv"
    df = pd.read_csv(src, sep="\t")
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    ax.plot(df["rank"], df["ranking_score"], marker="o", color=COLORS["sage"], linewidth=2)
    ax.scatter(df["rank"], df["ranking_score"], s=44, color=COLORS["sage"])
    ax.set_xticks(df["rank"])
    ax.set_xlabel("AF3 rank")
    ax.set_ylabel("ranking score")
    ax.set_ylim(0.42, 0.74)
    ax.set_title("Corrected-light CeS AF3 templates")
    for _, row in df.iterrows():
        ax.text(row["rank"], row["ranking_score"] + 0.012, f"{row['ranking_score']:.3f}", ha="center", fontsize=8)
    return save(fig, "fig4_af3_ces_rank_scores.png")


def plot_r3_md():
    data = [
        ("base", 1.3, -39),
        ("T1", 10.9, 78),
        ("42E", 1.9, -165),
        ("42A", 1.7, -97),
        ("com16diff", 3.9, -154),
        ("com17diff", 23.3, 45),
        ("cross1", 1.6, -131),
        ("cross2", 1.9, 40),
        ("cross3", 1.3, -92),
        ("com17clean", 8.3, 27),
    ]
    df = pd.DataFrame(data, columns=["variant", "rmsd74", "delta_contacts"])
    colors = [COLORS["sage"] if x < 0 and r < 4 else COLORS["rose"] for x, r in zip(df["delta_contacts"], df["rmsd74"])]

    fig, ax = plt.subplots(figsize=(6.8, 3.4))
    ax.bar(df["variant"], df["delta_contacts"], color=colors, width=0.72)
    ax.axhline(0, color="#5D6268", linewidth=0.9)
    ax.set_ylabel("Delta contacts, pH6.0 - pH7.4")
    ax.set_title("1E62 R3 MD: acid-state contact loss")
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    for i, row in df.iterrows():
        y = row["delta_contacts"]
        va = "bottom" if y >= 0 else "top"
        offset = 6 if y >= 0 else -6
        ax.text(i, y + offset, f"{int(y)}", ha="center", va=va, fontsize=8)
    return save(fig, "fig5_1e62_r3_md_delta_contacts.png")


def plot_sdab_overview():
    sections = pd.DataFrame(
        {
            "Section": ["MPNN full/min/Y37", "ELISA model", "FR-His"],
            "Variants": [12, 10, 10],
            "Directional": [3, 9, 7],
        }
    )
    x = np.arange(len(sections))
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    ax.bar(x - 0.18, sections["Variants"], width=0.36, label="delivered", color=COLORS["blue"])
    ax.bar(x + 0.18, sections["Directional"], width=0.36, label="MD direction OK", color=COLORS["sage"])
    ax.set_xticks(x, sections["Section"], rotation=15, ha="right")
    ax.set_ylabel("variant count")
    ax.set_title("sdab delivery package: 32 variants")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    for i, row in sections.iterrows():
        ax.text(i - 0.18, row["Variants"] + 0.35, str(row["Variants"]), ha="center", fontsize=8)
        ax.text(i + 0.18, row["Directional"] + 0.35, str(row["Directional"]), ha="center", fontsize=8)
    return save(fig, "fig6_sdab_32_variant_overview.png")


def write_source_index(generated):
    rows = [
        ("fig1_r2_elisa_key_variants.png", ".tasks/active/hbsag-genotype-generalization/artifacts/r2_per_type_phenotype_table.csv"),
        ("fig2_hbsag_type_sensitive_map.png", ".tasks/active/hbsag-genotype-generalization/artifacts/antigen_difference_table.md"),
        ("fig3_foldx_auc_boundary.png", ".tasks/active/hbsag-genotype-generalization/artifacts/foldx_validation_legacy/findings.md"),
        ("fig4_af3_ces_rank_scores.png", "experiments/1E62/R4/af3_ces_wt_corrected_light/summary/ranked_structures.tsv"),
        ("fig5_1e62_r3_md_delta_contacts.png", "docs/wet_lab_delivery/1e62_r3_design_summary.md"),
        ("fig6_sdab_32_variant_overview.png", "docs/wet_lab_delivery/sdab_design_summary.md"),
    ]
    with (OUT / "figure_sources.tsv").open("w", encoding="utf-8") as f:
        f.write("figure\tsource\n")
        for row in rows:
            f.write("\t".join(row) + "\n")


def main():
    setup()
    generated = [
        plot_r2_heatmap(),
        plot_antigen_map(),
        plot_foldx_auc(),
        plot_af3_ces_ranks(),
        plot_r3_md(),
        plot_sdab_overview(),
    ]
    write_source_index(generated)
    for path in generated:
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
