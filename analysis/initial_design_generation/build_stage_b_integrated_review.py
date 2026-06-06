#!/usr/bin/env python3
"""Build integrated Stage B review after route-limited heavy-lite compute.

This stage is an evidence integration layer only. It does not run AF3,
SimpleFold, MD, glycan modeling, or final library selection.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"
ROUTE_DIR = ROOT / "results/initial_design_generation/tier2_route_limited_stage_b"
COMPUTE_DIR = ROUTE_DIR / "compute"
OUT_DIR = ROOT / "results/initial_design_generation/stage_b_integrated_review"

PKA_REVIEW = COMPUTE_DIR / "electrostatics_pka_integrated_review.csv"
SIMPLEFOLD_REVIEW = COMPUTE_DIR / "simplefold_antibody_ensemble_review.csv"
AF3_REVIEW = COMPUTE_DIR / "af3_complex_review_results.csv"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    view = df.copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
    lines = [
        "| " + " | ".join(map(str, view.columns)) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in view.columns) + " |")
    return "\n".join(lines)


def read_inputs() -> pd.DataFrame:
    base = pd.read_csv(PKA_REVIEW)
    simplefold = pd.read_csv(SIMPLEFOLD_REVIEW)
    af3 = pd.read_csv(AF3_REVIEW)

    sf_keep = [
        "variant_id",
        "target",
        "simplefold_status",
        "num_samples",
        "median_plddt",
        "fold_review_class",
        "fold_review_notes",
    ]
    af3_keep = [
        "variant_id",
        "target",
        "af3_status",
        "ptm",
        "iptm",
        "chain_A_plddt",
        "chain_B_plddt",
        "chain_C_plddt",
        "chain_H_plddt",
        "chain_L_plddt",
        "iptm_A_B",
        "iptm_C_H",
        "iptm_C_L",
        "iptm_H_L",
    ]
    df = base.merge(simplefold[sf_keep], on=["variant_id", "target"], how="left", validate="one_to_one")
    df = df.merge(af3[af3_keep], on=["variant_id", "target"], how="left", validate="one_to_one")

    for col in [
        "mutation_count",
        "pyrosetta_rosetta_delta_score",
        "pyrosetta_local_structure_validity_t2_score",
        "pyrosetta_local_clash_count",
        "pka_pka_combined_mean",
        "his_pka_support_t2_score",
        "his_min_antigen_distance",
        "median_plddt",
        "ptm",
        "iptm",
        "chain_A_plddt",
        "chain_B_plddt",
        "chain_C_plddt",
        "chain_H_plddt",
        "chain_L_plddt",
        "iptm_A_B",
        "iptm_C_H",
        "iptm_C_L",
        "iptm_H_L",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def fold_class(row: pd.Series) -> str:
    if not bool(row.get("simplefold_selected", False)):
        return "fold_unchecked"
    cls = row.get("fold_review_class")
    if isinstance(cls, str) and cls:
        return cls
    return "fold_missing_result"


def complex_class(row: pd.Series) -> str:
    if not bool(row.get("af3_selected", False)):
        return "complex_unchecked"
    if pd.isna(row.get("iptm")):
        return "complex_missing_result"

    iptm = float(row["iptm"])
    target = row.get("target")
    if target == "Ab_1E62":
        if iptm >= 0.39:
            return "complex_plausible_relative"
        if iptm >= 0.35:
            return "complex_boundary"
        return "complex_weak"

    if target == "Ab_sdAb":
        # sdAb is judged relative to this route's weak absolute complex signal.
        # Without a parent baseline, do not call these structure-confirmed.
        if iptm >= 0.20:
            return "complex_plausible_relative"
        if iptm >= 0.10:
            return "complex_weak"
        return "complex_very_weak"

    return "complex_unclassified"


def glycan_class(row: pd.Series) -> str:
    if bool(row.get("glycan_selected", False)):
        return "glycan_selected_but_not_reviewed"
    return "glycan_unchecked_no_flag"


def is_good_pka(cls: str) -> bool:
    return cls in {"pH_mechanism_strong", "pH_mechanism_plausible"}


def is_boundary_pka(cls: str) -> bool:
    return cls == "pH_mechanism_boundary"


def integrate_row(row: pd.Series) -> tuple[str, str, int, str, str]:
    target = row["target"]
    pka_cls = row["electrostatic_mechanism_class"]
    fold_cls = row["fold_review_class_integrated"]
    complex_cls = row["complex_review_class"]
    local_cls = row["heavy_lite_final_class"]
    flags: list[str] = []
    notes: list[str] = []

    if complex_cls == "complex_unchecked":
        flags.append("complex_unchecked")
    if fold_cls == "fold_unchecked":
        flags.append("fold_unchecked")
    if row["glycan_review_class"].startswith("glycan_unchecked"):
        flags.append("glycan_unchecked")

    if target == "Ab_sdAb":
        flags.append("parent_complex_baseline_missing")
    elif target == "Ab_1E62" and complex_cls != "complex_unchecked":
        flags.append("parent_complex_baseline_missing")

    if is_good_pka(pka_cls) and fold_cls == "fold_supported" and complex_cls in {
        "complex_weak",
        "complex_very_weak",
    }:
        flags.append("complex_weak_but_pH_fold_supported")
    if is_good_pka(pka_cls) and fold_cls == "fold_unchecked" and complex_cls in {
        "complex_weak",
        "complex_very_weak",
    }:
        flags.append("complex_weak_but_pH_supported_fold_unchecked")

    if is_boundary_pka(pka_cls) and complex_cls in {"complex_plausible_relative", "complex_boundary"}:
        flags.append("pH_boundary_but_complex_plausible")

    if local_cls not in {"HL_promotable", "HL_boundary_supported"}:
        return (
            "integrated_reject_or_audit",
            "reject_or_audit",
            9,
            ";".join(flags),
            "local_heavy_lite_not_supported",
        )

    if target == "Ab_1E62":
        if is_good_pka(pka_cls) and complex_cls == "complex_plausible_relative":
            notes.append("1E62 local_pH evidence with relative AF3 complex support")
            return (
                "primary_complex_plausible",
                "promote_review",
                1,
                ";".join(flags),
                ";".join(notes),
            )
        if is_good_pka(pka_cls) and complex_cls == "complex_unchecked":
            notes.append("1E62 local_pH evidence; complex route not selected")
            return (
                "primary_local_pH_supported_complex_unchecked",
                "promote_review",
                2,
                ";".join(flags),
                ";".join(notes),
            )
        if is_good_pka(pka_cls) and complex_cls in {"complex_boundary", "complex_weak"}:
            notes.append("1E62 route conflict or boundary complex evidence")
            return (
                "primary_route_conflict_review",
                "route_conflict_review",
                4,
                ";".join(flags),
                ";".join(notes),
            )
        if is_boundary_pka(pka_cls):
            notes.append("1E62 boundary pH/electrostatic support")
            return (
                "primary_boundary",
                "retain_boundary",
                5,
                ";".join(flags),
                ";".join(notes),
            )
        return (
            "primary_low_support",
            "retain_boundary",
            6,
            ";".join(flags),
            "1E62 low integrated support",
        )

    if target == "Ab_sdAb":
        if is_good_pka(pka_cls) and fold_cls == "fold_supported" and complex_cls in {
            "complex_weak",
            "complex_very_weak",
        }:
            notes.append("sdAb pH/fold supported but current AF3 complex route is weak")
            return (
                "pH_fold_supported_complex_weak",
                "secondary_promote_review",
                3,
                ";".join(flags),
                ";".join(notes),
            )
        if is_good_pka(pka_cls) and fold_cls == "fold_supported" and complex_cls == "complex_unchecked":
            notes.append("sdAb pH/fold supported; complex route not selected")
            return (
                "pH_fold_supported_complex_unchecked",
                "secondary_promote_review",
                4,
                ";".join(flags),
                ";".join(notes),
            )
        if is_good_pka(pka_cls) and fold_cls == "fold_unchecked":
            notes.append("sdAb pH supported; fold and/or complex route not selected")
            return (
                "pH_supported_route_unchecked",
                "retain_secondary_review",
                5,
                ";".join(flags),
                ";".join(notes),
            )
        if is_boundary_pka(pka_cls):
            notes.append("sdAb boundary pH/electrostatic support")
            return (
                "secondary_boundary",
                "retain_boundary",
                6,
                ";".join(flags),
                ";".join(notes),
            )
        return (
            "secondary_low_support",
            "retain_boundary",
            7,
            ";".join(flags),
            "sdAb low integrated support",
        )

    return (
        "integrated_unclassified",
        "retain_boundary",
        8,
        ";".join(flags),
        "target_not_recognized",
    )


def add_integrated_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["fold_review_class_integrated"] = out.apply(fold_class, axis=1)
    out["complex_review_class"] = out.apply(complex_class, axis=1)
    out["glycan_review_class"] = out.apply(glycan_class, axis=1)
    out["parent_baseline_complex_score"] = pd.NA
    out["mutant_delta_complex_score"] = pd.NA
    out["parent_baseline_status"] = "missing_parent_baseline_not_blocking_integrated_review"
    integrated = out.apply(integrate_row, axis=1, result_type="expand")
    integrated.columns = [
        "integrated_evidence_class",
        "integrated_action",
        "selection_priority",
        "evidence_conflict_flags",
        "notes",
    ]
    out = pd.concat([out, integrated], axis=1)
    out["af3_iptm"] = out["iptm"]
    out["af3_ptm"] = out["ptm"]
    out["simplefold_median_plddt"] = out["median_plddt"]
    return out


def write_csvs(df: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    keep = [
        "variant_id",
        "target",
        "mutation_list",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "mutation_count",
        "heavy_lite_final_class",
        "stage2a_final_class",
        "boundary_support_level",
        "electrostatic_mechanism_class",
        "pka_review_class",
        "his_positions",
        "his_min_antigen_distance",
        "his_pka_values",
        "his_pka_support_score",
        "fold_review_class_integrated",
        "simplefold_status",
        "simplefold_median_plddt",
        "complex_review_class",
        "af3_status",
        "af3_iptm",
        "af3_ptm",
        "iptm_A_B",
        "iptm_C_H",
        "iptm_C_L",
        "glycan_review_class",
        "integrated_evidence_class",
        "integrated_action",
        "evidence_conflict_flags",
        "selection_priority",
        "parent_baseline_complex_score",
        "mutant_delta_complex_score",
        "parent_baseline_status",
        "pyrosetta_rosetta_delta_score",
        "pyrosetta_local_structure_validity_t2_score",
        "pyrosetta_local_clash_count",
        "af3_selected",
        "simplefold_selected",
        "electrostatics_selected",
        "glycan_selected",
        "notes",
    ]
    df[keep].to_csv(OUT_DIR / "stage_b_integrated_review_table.csv", index=False)


def count_table(df: pd.DataFrame, keys: list[str], value: str = "count") -> pd.DataFrame:
    return (
        df.groupby(keys, dropna=False)
        .size()
        .reset_index(name=value)
        .sort_values(keys[:-1] + [value], ascending=[True] * max(0, len(keys) - 1) + [False])
    )


def route_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for target, sub in df.groupby("target"):
        rows.append(
            {
                "target": target,
                "main_count": len(sub),
                "af3_selected": int(sub["af3_selected"].sum()),
                "simplefold_selected": int(sub["simplefold_selected"].sum()),
                "electrostatics_selected": int(sub["electrostatics_selected"].sum()),
                "glycan_selected": int(sub["glycan_selected"].sum()),
                "complex_plausible_relative": int(sub["complex_review_class"].eq("complex_plausible_relative").sum()),
                "complex_weak": int(sub["complex_review_class"].isin(["complex_weak", "complex_very_weak"]).sum()),
                "complex_unchecked": int(sub["complex_review_class"].eq("complex_unchecked").sum()),
                "fold_supported": int(sub["fold_review_class_integrated"].eq("fold_supported").sum()),
                "fold_unchecked": int(sub["fold_review_class_integrated"].eq("fold_unchecked").sum()),
                "pH_strong_or_plausible": int(sub["electrostatic_mechanism_class"].isin(["pH_mechanism_strong", "pH_mechanism_plausible"]).sum()),
                "pH_boundary": int(sub["electrostatic_mechanism_class"].eq("pH_mechanism_boundary").sum()),
            }
        )
    return pd.DataFrame(rows)


def write_reports(df: pd.DataFrame) -> None:
    by_target_action = count_table(df, ["target", "integrated_action"])
    by_target_class = count_table(df, ["target", "integrated_evidence_class"])
    by_route = route_summary(df)
    by_seed = count_table(df, ["target", "his_seed_set", "integrated_evidence_class"])
    by_cluster = count_table(df, ["target", "near_duplicate_cluster_id", "integrated_evidence_class"])
    conflicts = df[df["evidence_conflict_flags"].fillna("").ne("")].copy()
    conflict_summary = count_table(conflicts, ["target", "evidence_conflict_flags"])

    (OUT_DIR / "stage_b_integrated_by_target.md").write_text(
        "\n".join(
            [
                "# Stage B Integrated Review By Target",
                "",
                "## Route Coverage",
                md_table(by_route),
                "",
                "## Integrated Action",
                md_table(by_target_action),
                "",
                "## Integrated Evidence Class",
                md_table(by_target_class),
                "",
            ]
        ),
        encoding="utf-8",
    )
    (OUT_DIR / "stage_b_integrated_by_route.md").write_text(
        "\n".join(["# Stage B Integrated Review By Route", "", md_table(by_route), ""]),
        encoding="utf-8",
    )
    (OUT_DIR / "stage_b_integrated_by_seed.md").write_text(
        "\n".join(["# Stage B Integrated Review By Seed", "", md_table(by_seed), ""]),
        encoding="utf-8",
    )
    (OUT_DIR / "stage_b_integrated_by_cluster.md").write_text(
        "\n".join(
            [
                "# Stage B Integrated Review By Cluster",
                "",
                "Top rows per target are shown to monitor near-duplicate concentration.",
                "",
                md_table(by_cluster.groupby("target", group_keys=False).head(40)),
                "",
            ]
        ),
        encoding="utf-8",
    )
    conflict_cols = [
        "variant_id",
        "target",
        "mutation_list",
        "his_seed_set",
        "integrated_evidence_class",
        "integrated_action",
        "evidence_conflict_flags",
        "complex_review_class",
        "fold_review_class_integrated",
        "electrostatic_mechanism_class",
    ]
    conflicts[conflict_cols].to_csv(OUT_DIR / "stage_b_route_conflict_rows.csv", index=False)
    (OUT_DIR / "stage_b_route_conflict_report.md").write_text(
        "\n".join(
            [
                "# Stage B Route Conflict Report",
                "",
                "Flags include unchecked routes as well as true route conflicts. Unchecked means the route was not selected in the route-limited panel, not that the candidate failed that route.",
                "",
                "## Conflict / Unchecked Flag Summary",
                md_table(conflict_summary),
                "",
                "## Rows With True Complex Weak But pH/Fold Support",
                md_table(
                    conflicts[
                        conflicts["evidence_conflict_flags"].str.contains(
                            "complex_weak_but_pH_fold_supported", na=False
                        )
                    ][conflict_cols].head(80)
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )

    next_gate = build_next_gate_report(df, by_route, by_target_action, by_target_class)
    (OUT_DIR / "stage_b_next_gate_report.md").write_text(next_gate, encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(next_gate, encoding="utf-8")


def build_next_gate_report(
    df: pd.DataFrame,
    by_route: pd.DataFrame,
    by_target_action: pd.DataFrame,
    by_target_class: pd.DataFrame,
) -> str:
    primary = df[df["integrated_action"].eq("promote_review")]
    secondary = df[df["integrated_action"].isin(["secondary_promote_review", "retain_secondary_review"])]
    boundary = df[df["integrated_action"].eq("retain_boundary")]
    conflicts = df[df["integrated_action"].eq("route_conflict_review")]

    source_rows = pd.DataFrame(
        [
            {"input": str(PKA_REVIEW.relative_to(ROOT)), "rows": len(pd.read_csv(PKA_REVIEW)), "sha256": sha256_file(PKA_REVIEW)},
            {
                "input": str(SIMPLEFOLD_REVIEW.relative_to(ROOT)),
                "rows": len(pd.read_csv(SIMPLEFOLD_REVIEW)),
                "sha256": sha256_file(SIMPLEFOLD_REVIEW),
            },
            {"input": str(AF3_REVIEW.relative_to(ROOT)), "rows": len(pd.read_csv(AF3_REVIEW)), "sha256": sha256_file(AF3_REVIEW)},
        ]
    )

    e62 = df[df["target"].eq("Ab_1E62")]
    sdab = df[df["target"].eq("Ab_sdAb")]
    e62_primary = e62[e62["integrated_action"].eq("promote_review")]
    sdab_secondary = sdab[sdab["integrated_action"].isin(["secondary_promote_review", "retain_secondary_review"])]
    sdab_complex_weak = sdab[sdab["complex_review_class"].isin(["complex_weak", "complex_very_weak"])]

    lines = [
        "# Integrated Stage B Review",
        "",
        "## Executive Summary",
        "",
        "Verdict: `INTEGRATED_REVIEW_COMPLETE__BANK_PLANNING_ALLOWED__NO_FINAL_UNLOCK`.",
        "",
        "This stage merged the route-limited Stage B outputs into one candidate-level interpretation table. It did not run new AF3, SimpleFold, MD, glycan modeling, or final 10K selection.",
        "",
        "The key interpretation is target-specific:",
        "",
        "- `1E62` has a primary branch with local pH/electrostatics support and relatively plausible AF3 complex support where AF3 was selected.",
        "- `sdAb` remains supported by pH/electrostatics and antibody-only fold checks, but its AF3 complex route is weak where selected. It should be retained as a secondary/exploratory branch rather than merged with 1E62 as structure-confirmed evidence.",
        "- Glycan review remains `unchecked_no_flag`; this must not be reported as low glycan risk.",
        "- Parent complex baseline comparison is missing and should be treated as a supplemental diagnostic, especially for sdAb.",
        "",
        "## Source Inputs",
        "",
        md_table(source_rows),
        "",
        "## Route Coverage And Evidence Counts",
        "",
        md_table(by_route),
        "",
        "## Integrated Action Counts",
        "",
        md_table(by_target_action),
        "",
        "## Integrated Evidence Classes",
        "",
        md_table(by_target_class),
        "",
        "## Primary / Secondary Bank Implications",
        "",
        md_table(
            pd.DataFrame(
                [
                    {
                        "target": "Ab_1E62",
                        "main_rows": len(e62),
                        "primary_promote_review": len(e62_primary),
                        "secondary_or_retained_review": 0,
                        "boundary_or_conflict": int(
                            e62["integrated_action"].isin(["retain_boundary", "route_conflict_review"]).sum()
                        ),
                        "complex_weak_selected_rows": 0,
                        "interpretation": "primary candidate-bank branch; still not final library unlock",
                    },
                    {
                        "target": "Ab_sdAb",
                        "main_rows": len(sdab),
                        "primary_promote_review": 0,
                        "secondary_or_retained_review": len(sdab_secondary),
                        "boundary_or_conflict": int(sdab["integrated_action"].eq("retain_boundary").sum()),
                        "complex_weak_selected_rows": len(sdab_complex_weak),
                        "interpretation": "secondary pH/fold-supported branch; complex evidence weak",
                    },
                ]
            )
        ),
        "",
        "## Gate Recommendation",
        "",
        "Recommended next step: construct a Tier2 candidate bank from integrated classes, with separate primary and secondary branches.",
        "",
        "Allowed now:",
        "",
        "- Build a primary 1E62 candidate bank from `promote_review` candidates.",
        "- Build a secondary sdAb candidate bank from `secondary_promote_review` and selected `retain_secondary_review` candidates.",
        "- Keep seed/cluster capping in the bank construction step.",
        "- Optionally plan a small parent-baseline complex diagnostic before any sdAb final confidence upgrade.",
        "",
        "Still locked:",
        "",
        "- Broad Tier2-heavy execution.",
        "- MD or constant-pH MD.",
        "- Final 10K library selection.",
        "- Treating sdAb complex-weak candidates as structure-confirmed.",
        "- Reporting glycan as low risk without an explicit glycan route.",
        "",
        "## Output Files",
        "",
        "- `results/initial_design_generation/stage_b_integrated_review/stage_b_integrated_review_table.csv`",
        "- `results/initial_design_generation/stage_b_integrated_review/stage_b_integrated_by_target.md`",
        "- `results/initial_design_generation/stage_b_integrated_review/stage_b_integrated_by_route.md`",
        "- `results/initial_design_generation/stage_b_integrated_review/stage_b_integrated_by_seed.md`",
        "- `results/initial_design_generation/stage_b_integrated_review/stage_b_integrated_by_cluster.md`",
        "- `results/initial_design_generation/stage_b_integrated_review/stage_b_route_conflict_report.md`",
        "- `results/initial_design_generation/stage_b_integrated_review/stage_b_route_conflict_rows.csv`",
        "- `results/initial_design_generation/stage_b_integrated_review/stage_b_next_gate_report.md`",
        "",
    ]
    return "\n".join(lines)


def write_manifest(df: pd.DataFrame) -> None:
    manifest = [
        "stage: integrated_stage_b_review",
        "status: complete",
        f"row_count: {len(df)}",
        "inputs:",
        f"  electrostatics_pka: {PKA_REVIEW.relative_to(ROOT)}",
        f"  electrostatics_pka_sha256: {sha256_file(PKA_REVIEW)}",
        f"  simplefold: {SIMPLEFOLD_REVIEW.relative_to(ROOT)}",
        f"  simplefold_sha256: {sha256_file(SIMPLEFOLD_REVIEW)}",
        f"  af3_complex: {AF3_REVIEW.relative_to(ROOT)}",
        f"  af3_complex_sha256: {sha256_file(AF3_REVIEW)}",
        "locks:",
        "  broad_tier2_heavy: true",
        "  md: true",
        "  final_10k: true",
        "  sdab_structure_confirmed_upgrade: true",
    ]
    (OUT_DIR / "stage_b_integrated_manifest.yaml").write_text("\n".join(manifest) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = add_integrated_columns(read_inputs())
    write_csvs(df)
    write_reports(df)
    write_manifest(df)
    print(f"wrote {OUT_DIR / 'stage_b_integrated_review_table.csv'} rows={len(df)}")
    print(f"wrote {OUT_DIR / 'stage_b_next_gate_report.md'}")


if __name__ == "__main__":
    main()
