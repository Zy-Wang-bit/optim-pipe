#!/usr/bin/env python3
"""Build mutation-position distribution plots and PyMOL hotspot coloring inputs.

This script uses only the Python standard library because the project runtime
does not guarantee matplotlib, numpy, or PyMOL Python bindings.
"""

import csv
import html
import math
import re
from collections import Counter, defaultdict, OrderedDict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "results/initial_design_generation/mutation_hotspot_visualization"

TARGETS = {
    "Ab_1E62": {
        "short": "1E62",
        "candidate_csv": ROOT
        / "results/initial_design_generation/per_target_15k_candidate_pools/final_candidate_pool_1E62_15k_draft.csv",
        "pdb": ROOT
        / "results/initial_design_generation/p0_mpnn_backbones/af3_1E62_seed7_model5_rank1.pdb",
        "design_chain": "L",
        "antigen_chain": "C",
        "mutation_chain_prefix": "L",
        "window_start": 1,
        "window_end": 40,
        "window_label": "VL:1-40",
    },
    "Ab_sdAb": {
        "short": "sdAb",
        "candidate_csv": ROOT
        / "results/initial_design_generation/per_target_15k_candidate_pools/final_candidate_pool_sdAb_15k_draft.csv",
        "pdb": ROOT
        / "results/initial_design_generation/p0_mpnn_backbones/af3_sdAb_seed16_model5_rank1.pdb",
        "design_chain": "A",
        "antigen_chain": "B",
        "mutation_chain_prefix": "A",
        "window_start": 72,
        "window_end": 111,
        "window_label": "VHH:72-111",
    },
}

MUTATION_RE = re.compile(r"^([A-Za-z])([A-Za-z])(\d+)([A-Za-z])$")

AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_mutation_list(text):
    if text is None:
        return []
    text = str(text).strip()
    if not text or text.lower() in {"nan", "none", "na"}:
        return []
    parsed = []
    for token in text.split(";"):
        token = token.strip()
        if not token:
            continue
        match = MUTATION_RE.match(token)
        if not match:
            continue
        chain_prefix, wt, pos, new = match.groups()
        parsed.append(
            {
                "token": token,
                "chain_prefix": chain_prefix.upper(),
                "wildtype_aa": wt.upper(),
                "position": int(pos),
                "mutant_aa": new.upper(),
            }
        )
    return parsed


def read_chain_residues(pdb_path, chain_id):
    residues = OrderedDict()
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[21] != chain_id:
                continue
            resseq = line[22:26].strip()
            icode = line[26].strip()
            if icode:
                continue
            if not resseq.lstrip("-").isdigit():
                continue
            position = int(resseq)
            residues.setdefault(position, line[17:20].strip())
    return residues


def collect_target_frequency(target_name, cfg):
    position_to_mutated_rows = Counter()
    position_to_his_rows = Counter()
    position_to_mutants = defaultdict(Counter)
    position_to_tokens = defaultdict(Counter)
    parse_skipped_tokens = Counter()
    row_count = 0
    row_with_mutations = 0

    with cfg["candidate_csv"].open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row_count += 1
            mutations = parse_mutation_list(row.get("mutation_list", ""))
            if mutations:
                row_with_mutations += 1
            seen_positions = set()
            for item in mutations:
                if item["chain_prefix"] != cfg["mutation_chain_prefix"]:
                    parse_skipped_tokens[item["token"]] += 1
                    continue
                position = item["position"]
                if position not in seen_positions:
                    position_to_mutated_rows[position] += 1
                    seen_positions.add(position)
                if item["mutant_aa"] == "H":
                    position_to_his_rows[position] += 1
                position_to_mutants[position][item["mutant_aa"]] += 1
                position_to_tokens[position][item["token"]] += 1

    pdb_residues = read_chain_residues(cfg["pdb"], cfg["design_chain"])
    rows = []
    for position, residue3 in pdb_residues.items():
        if not (cfg["window_start"] <= position <= cfg["window_end"]):
            continue
        count = position_to_mutated_rows[position]
        mutant_counts = position_to_mutants[position]
        token_counts = position_to_tokens[position]
        rows.append(
            {
                "target": target_name,
                "target_short": cfg["short"],
                "pdb_path": str(cfg["pdb"].relative_to(ROOT)),
                "design_chain": cfg["design_chain"],
                "position": position,
                "pdb_residue_name": residue3,
                "pdb_wildtype_aa": AA3_TO_AA1.get(residue3, "X"),
                "mutation_row_count": count,
                "mutation_frequency": count / float(row_count) if row_count else 0.0,
                "mutation_frequency_percent": 100.0 * count / float(row_count) if row_count else 0.0,
                "his_mutation_count": position_to_his_rows[position],
                "his_mutation_frequency_percent": 100.0
                * position_to_his_rows[position]
                / float(row_count)
                if row_count
                else 0.0,
                "mutant_aa_counts": format_counter(mutant_counts),
                "mutation_token_counts": format_counter(token_counts),
                "window_label": cfg["window_label"],
            }
        )

    return {
        "target": target_name,
        "short": cfg["short"],
        "row_count": row_count,
        "row_with_mutations": row_with_mutations,
        "rows": rows,
        "parse_skipped_tokens": parse_skipped_tokens,
    }


def format_counter(counter):
    if not counter:
        return ""
    return ";".join("{}:{}".format(k, counter[k]) for k in sorted(counter, key=lambda x: (-counter[x], x)))


def write_frequency_csv(target_results):
    all_rows = []
    for result in target_results:
        all_rows.extend(result["rows"])
        target_file = OUTDIR / "mutation_position_frequency_{}.csv".format(result["short"])
        with target_file.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(result["rows"][0].keys()))
            writer.writeheader()
            writer.writerows(result["rows"])

    all_file = OUTDIR / "mutation_position_frequency_all_targets.csv"
    with all_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)


def hotspot_color(freq_percent):
    if freq_percent >= 25:
        return "#b2182b"
    if freq_percent >= 15:
        return "#ef8a62"
    if freq_percent >= 5:
        return "#fddbc7"
    if freq_percent > 0:
        return "#fee090"
    return "#d9d9d9"


def write_svg(result):
    rows = result["rows"]
    width = 1280
    height = 520
    margin_left = 70
    margin_right = 30
    margin_top = 70
    margin_bottom = 95
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    max_freq = max(row["mutation_frequency_percent"] for row in rows) or 1.0
    y_max = max(5.0, math.ceil(max_freq / 5.0) * 5.0)
    bar_gap = 1.0
    bar_width = max(2.0, plot_width / len(rows) - bar_gap)

    top_rows = sorted(rows, key=lambda r: (-r["mutation_frequency_percent"], r["position"]))[:12]
    top_positions = {row["position"] for row in top_rows}

    elements = []
    elements.append(
        '<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">'.format(
            w=width, h=height
        )
    )
    elements.append('<rect width="100%" height="100%" fill="white"/>')
    elements.append(
        '<text x="{x}" y="32" font-family="Arial" font-size="22" font-weight="700">{title}</text>'.format(
            x=margin_left,
            title=html.escape(
                "{} mutation position distribution, {} 15K draft".format(
                    result["short"], TARGETS[result["target"]]["window_label"]
                )
            ),
        )
    )
    elements.append(
        '<text x="{x}" y="54" font-family="Arial" font-size="13" fill="#555">40 aa design window only. Bar height = fraction of candidates mutated at this residue. Full 15K draft including controls.</text>'.format(
            x=margin_left
        )
    )

    for tick in range(0, int(y_max) + 1, 5):
        y = margin_top + plot_height - (tick / y_max) * plot_height
        elements.append(
            '<line x1="{x1}" y1="{y:.1f}" x2="{x2}" y2="{y:.1f}" stroke="#e5e5e5" stroke-width="1"/>'.format(
                x1=margin_left, x2=width - margin_right, y=y
            )
        )
        elements.append(
            '<text x="{x}" y="{y:.1f}" text-anchor="end" dominant-baseline="middle" font-family="Arial" font-size="11" fill="#555">{tick}%</text>'.format(
                x=margin_left - 8, y=y, tick=tick
            )
        )

    elements.append(
        '<line x1="{x}" y1="{y}" x2="{x2}" y2="{y}" stroke="#333" stroke-width="1.2"/>'.format(
            x=margin_left, x2=width - margin_right, y=margin_top + plot_height
        )
    )
    elements.append(
        '<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2}" stroke="#333" stroke-width="1.2"/>'.format(
            x=margin_left, y1=margin_top, y2=margin_top + plot_height
        )
    )

    for idx, row in enumerate(rows):
        x = margin_left + idx * (plot_width / len(rows))
        bar_h = (row["mutation_frequency_percent"] / y_max) * plot_height
        y = margin_top + plot_height - bar_h
        elements.append(
            '<rect x="{x:.2f}" y="{y:.2f}" width="{bw:.2f}" height="{bh:.2f}" fill="{fill}"/>'.format(
                x=x,
                y=y,
                bw=bar_width,
                bh=bar_h,
                fill=hotspot_color(row["mutation_frequency_percent"]),
            )
        )
        if row["position"] in top_positions:
            elements.append(
                '<text x="{x:.1f}" y="{y:.1f}" transform="rotate(-55 {x:.1f} {y:.1f})" font-family="Arial" font-size="10" fill="#111">{label}</text>'.format(
                    x=x + bar_width / 2,
                    y=max(18, y - 5),
                    label=html.escape("{}{}".format(row["pdb_wildtype_aa"], row["position"])),
                )
            )

    tick_step = 10
    for position in range(rows[0]["position"], rows[-1]["position"] + 1, tick_step):
        idx = position - rows[0]["position"]
        if 0 <= idx < len(rows):
            x = margin_left + idx * (plot_width / len(rows)) + bar_width / 2
            elements.append(
                '<line x1="{x:.1f}" y1="{y1}" x2="{x:.1f}" y2="{y2}" stroke="#333" stroke-width="1"/>'.format(
                    x=x, y1=margin_top + plot_height, y2=margin_top + plot_height + 4
                )
            )
            elements.append(
                '<text x="{x:.1f}" y="{y}" text-anchor="middle" font-family="Arial" font-size="10" fill="#555">{pos}</text>'.format(
                    x=x, y=margin_top + plot_height + 18, pos=position
                )
            )

    elements.append(
        '<text x="{x}" y="{y}" text-anchor="middle" font-family="Arial" font-size="13" fill="#333">Residue position on design chain</text>'.format(
            x=margin_left + plot_width / 2, y=height - 36
        )
    )
    elements.append(
        '<text x="18" y="{y}" text-anchor="middle" transform="rotate(-90 18 {y})" font-family="Arial" font-size="13" fill="#333">Mutation frequency</text>'.format(
            y=margin_top + plot_height / 2
        )
    )

    legend_x = width - 420
    legend_y = 28
    legend = [
        ("#b2182b", ">=25%"),
        ("#ef8a62", "15-25%"),
        ("#fddbc7", "5-15%"),
        ("#fee090", ">0-5%"),
        ("#d9d9d9", "0%"),
    ]
    for i, (fill, label) in enumerate(legend):
        x = legend_x + i * 78
        elements.append('<rect x="{x}" y="{y}" width="18" height="12" fill="{fill}"/>'.format(x=x, y=legend_y, fill=fill))
        elements.append(
            '<text x="{x}" y="{y}" font-family="Arial" font-size="11" fill="#444">{label}</text>'.format(
                x=x + 23, y=legend_y + 10, label=label
            )
        )

    elements.append("</svg>")
    (OUTDIR / "mutation_position_distribution_{}.svg".format(result["short"])).write_text(
        "\n".join(elements), encoding="utf-8"
    )


def write_combined_svg(target_results):
    # Keep this simple: stack the two per-target SVG charts by embedding them as
    # foreign-object-free repeated panels from the same drawing helper would be noisy.
    lines = [
        "<!doctype html>",
        '<html lang="en"><head><meta charset="utf-8"><title>Mutation hotspot distributions</title>',
        "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222} img{max-width:100%;border:1px solid #ddd;margin:12px 0 28px} table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:6px 9px}</style>",
        "</head><body>",
        "<h1>Mutation Hotspot Distribution</h1>",
        "<p>Each bar shows how often a residue is mutated inside the fixed 40 aa design window. Controls are included; they are 20 rows per target and have negligible impact on frequency.</p>",
    ]
    for result in target_results:
        short = result["short"]
        lines.append("<h2>{}</h2>".format(html.escape(short)))
        lines.append(
            '<img src="mutation_position_distribution_{}.svg" alt="{} mutation position distribution">'.format(
                html.escape(short), html.escape(short)
            )
        )
    lines.append("</body></html>")
    (OUTDIR / "mutation_position_distribution_overview.html").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def update_pdb_bfactor(target_name, cfg, rows):
    freq_by_pos = {row["position"]: row["mutation_frequency_percent"] for row in rows}
    out_pdb = OUTDIR / "{}_mutation_frequency_bfactor.pdb".format(cfg["short"])
    with cfg["pdb"].open() as source, out_pdb.open("w") as dest:
        for line in source:
            if line.startswith(("ATOM  ", "HETATM")) and line[21] == cfg["design_chain"]:
                resseq = line[22:26].strip()
                icode = line[26].strip()
                if not icode and resseq.lstrip("-").isdigit():
                    freq = freq_by_pos.get(int(resseq), 0.0)
                    line = line[:60] + "{:6.2f}".format(freq) + line[66:]
            dest.write(line)
    return out_pdb


def pymol_selection_for_positions(object_name, chain, positions):
    if not positions:
        return "none"
    chunks = []
    sorted_positions = sorted(positions)
    for start in range(0, len(sorted_positions), 40):
        chunk = "+".join(str(pos) for pos in sorted_positions[start : start + 40])
        chunks.append("({} and chain {} and resi {})".format(object_name, chain, chunk))
    return " or ".join(chunks)


def write_pymol_script(target_name, cfg, rows, bfactor_pdb):
    object_name = "{}_mutation_hotspots".format(cfg["short"])
    pml_path = OUTDIR / "{}_mutation_frequency_hotspots.pml".format(cfg["short"])
    bin_defs = [
        ("mut_hot_very_high", [r["position"] for r in rows if r["mutation_frequency_percent"] >= 25], "[0.70, 0.05, 0.10]"),
        ("mut_hot_high", [r["position"] for r in rows if 15 <= r["mutation_frequency_percent"] < 25], "[0.95, 0.35, 0.10]"),
        ("mut_hot_medium", [r["position"] for r in rows if 5 <= r["mutation_frequency_percent"] < 15], "[1.00, 0.75, 0.20]"),
        ("mut_hot_low", [r["position"] for r in rows if 0 < r["mutation_frequency_percent"] < 5], "[1.00, 0.92, 0.55]"),
    ]
    top_rows = sorted(rows, key=lambda r: (-r["mutation_frequency_percent"], r["position"]))[:12]
    window_positions = [row["position"] for row in rows]

    lines = [
        "# PyMOL mutation hotspot coloring for {}".format(cfg["short"]),
        "# Run from the repository root with: pymol {}".format(str(pml_path.relative_to(ROOT))),
        "load {}, {}".format(str(bfactor_pdb.relative_to(ROOT)), object_name),
        "hide everything, {}".format(object_name),
        "show cartoon, {}".format(object_name),
        "bg_color white",
        "color gray80, {}".format(object_name),
        "set cartoon_transparency, 0.55, {} and not chain {}".format(object_name, cfg["design_chain"]),
        "set cartoon_transparency, 0.15, {} and chain {}".format(object_name, cfg["design_chain"]),
        "select design_chain, {} and chain {}".format(object_name, cfg["design_chain"]),
        "select antigen_chain, {} and chain {}".format(object_name, cfg["antigen_chain"]),
        "color slate, antigen_chain",
        "set_color mut_hot_very_high, [0.70, 0.05, 0.10]",
        "set_color mut_hot_high, [0.95, 0.35, 0.10]",
        "set_color mut_hot_medium, [1.00, 0.75, 0.20]",
        "set_color mut_hot_low, [1.00, 0.92, 0.55]",
    ]
    for color_name, positions, _rgb in bin_defs:
        if not positions:
            continue
        selection = pymol_selection_for_positions(object_name, cfg["design_chain"], positions)
        lines.append("color {}, {}".format(color_name, selection))
    mutated_positions = [r["position"] for r in rows if r["mutation_frequency_percent"] > 0]
    top_positions = [r["position"] for r in top_rows]
    lines.extend(
        [
            "show sticks, {}".format(pymol_selection_for_positions(object_name, cfg["design_chain"], mutated_positions)),
            "set stick_radius, 0.18",
            "select design_window, {}".format(pymol_selection_for_positions(object_name, cfg["design_chain"], window_positions)),
            "show spheres, {}".format(pymol_selection_for_positions(object_name, cfg["design_chain"], top_positions)),
            "set sphere_scale, 0.36",
            "# B-factor also stores mutation frequency percent on the design chain.",
            "# Alternative gradient coloring:",
            "# spectrum b, blue_white_red, design_chain, minimum=0, maximum=35",
        ]
    )
    for row in top_rows:
        lines.append(
            'label {} and chain {} and resi {} and name CA, "{}{} {:.1f}%"'.format(
                object_name,
                cfg["design_chain"],
                row["position"],
                row["pdb_wildtype_aa"],
                row["position"],
                row["mutation_frequency_percent"],
            )
        )
    lines.extend(
        [
            "set label_size, 14",
            "set label_color, black",
            "orient design_chain",
            "zoom design_chain, 8",
        ]
    )
    pml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(target_results):
    lines = [
        "# Mutation Hotspot Visualization Report",
        "",
        "This report summarizes mutation-position frequency within the fixed 40 aa design windows in the per-target 15K draft pools and lists the PyMOL inputs generated for spatial hotspot visualization.",
        "",
        "Scope: fixed 40 aa design window only, using the full 15K draft pool per target, including the 20 control / anchor rows. Controls are retained because they are part of the current draft and have negligible effect on per-position frequency.",
        "",
        "PyMOL note: local PyMOL was not available in this environment, so this stage generated `.pml` scripts and B-factor-annotated PDB files for direct use on a machine with PyMOL.",
        "",
        "## Summary",
        "",
        "| target | design_window | candidate_rows | rows_with_mutations | mutated_positions | top_position | top_frequency | top_12_positions |",
        "| --- | --- | ---: | ---: | ---: | --- | ---: | --- |",
    ]
    for result in target_results:
        rows = result["rows"]
        mutated = [row for row in rows if row["mutation_row_count"] > 0]
        top = sorted(rows, key=lambda r: (-r["mutation_frequency_percent"], r["position"]))[:12]
        top_position = "{}{}".format(top[0]["pdb_wildtype_aa"], top[0]["position"]) if top else ""
        top_freq = "{:.2f}%".format(top[0]["mutation_frequency_percent"]) if top else "0%"
        top_12 = "; ".join(
            "{}{}={:.1f}%".format(row["pdb_wildtype_aa"], row["position"], row["mutation_frequency_percent"])
            for row in top
        )
        lines.append(
            "| {target} | {window} | {rows_n} | {mut_rows} | {pos_n} | {top_pos} | {top_freq} | {top12} |".format(
                target=result["short"],
                window=TARGETS[result["target"]]["window_label"],
                rows_n=result["row_count"],
                mut_rows=result["row_with_mutations"],
                pos_n=len(mutated),
                top_pos=top_position,
                top_freq=top_freq,
                top12=top_12,
            )
        )

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `mutation_position_frequency_1E62.csv`",
            "- `mutation_position_frequency_sdAb.csv`",
            "- `mutation_position_frequency_all_targets.csv`",
            "- `mutation_position_distribution_1E62.svg`",
            "- `mutation_position_distribution_sdAb.svg`",
            "- `mutation_position_distribution_overview.html`",
            "- `1E62_mutation_frequency_bfactor.pdb`",
            "- `sdAb_mutation_frequency_bfactor.pdb`",
            "- `1E62_mutation_frequency_hotspots.pml`",
            "- `sdAb_mutation_frequency_hotspots.pml`",
            "",
            "## PyMOL Usage",
            "",
            "From the project root:",
            "",
            "```bash",
            "pymol results/initial_design_generation/mutation_hotspot_visualization/1E62_mutation_frequency_hotspots.pml",
            "pymol results/initial_design_generation/mutation_hotspot_visualization/sdAb_mutation_frequency_hotspots.pml",
            "```",
            "",
            "Color bins in the scripts:",
            "",
            "- deep red: mutation frequency >= 25%",
            "- orange: 15-25%",
            "- gold: 5-15%",
            "- pale yellow: >0-5%",
            "- gray: no mutation in the current 15K draft",
            "",
            "The design-chain B-factor field also stores mutation frequency as percent, so PyMOL can alternatively use `spectrum b, blue_white_red, design_chain, minimum=0, maximum=35`.",
        ]
    )
    (OUTDIR / "mutation_hotspot_visualization_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    target_results = []
    for target_name, cfg in TARGETS.items():
        result = collect_target_frequency(target_name, cfg)
        target_results.append(result)

    write_frequency_csv(target_results)
    for result in target_results:
        cfg = TARGETS[result["target"]]
        write_svg(result)
        bfactor_pdb = update_pdb_bfactor(result["target"], cfg, result["rows"])
        write_pymol_script(result["target"], cfg, result["rows"], bfactor_pdb)
    write_combined_svg(target_results)
    write_report(target_results)

    print("Wrote mutation hotspot visualization outputs to {}".format(OUTDIR.relative_to(ROOT)))


if __name__ == "__main__":
    main()
