#!/usr/bin/env python3
"""Build Tier2 candidate bank from the Integrated Stage B review.

The bank is not a final 10K library. It organizes evidence-backed variants
into primary / secondary branches and emits expansion templates for later
backfill or constrained local expansion.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"
SOURCE = ROOT / "results/initial_design_generation/stage_b_integrated_review/stage_b_integrated_review_table.csv"
OUT_DIR = ROOT / "results/initial_design_generation/tier2_candidate_bank"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def short_hash(text: str, n: int = 10) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:n]


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


def read_source() -> pd.DataFrame:
    df = pd.read_csv(SOURCE)
    if df.duplicated(["variant_id", "target"]).any():
        dupes = df[df.duplicated(["variant_id", "target"], keep=False)][["variant_id", "target"]]
        raise ValueError(f"duplicate variant keys in integrated review:\n{dupes.head(20)}")
    return df


def split_mutations(mutation_list: object) -> tuple[list[str], list[str]]:
    if not isinstance(mutation_list, str) or not mutation_list:
        return [], []
    muts = [m.strip() for m in mutation_list.split(";") if m.strip()]
    his = [m for m in muts if m.endswith("H")]
    non_his = [m for m in muts if not m.endswith("H")]
    return his, non_his


def bank_role(row: pd.Series) -> tuple[str, str, str, int]:
    cls = row["integrated_evidence_class"]
    action = row["integrated_action"]

    if row["target"] == "Ab_1E62":
        branch = "primary_1E62"
        if cls == "primary_complex_plausible":
            return branch, "primary_complex_plausible_template", "high_priority_template", 1
        if cls == "primary_local_pH_supported_complex_unchecked":
            return branch, "primary_local_pH_supported_route_unchecked", "conditional_template", 2
        if action == "retain_boundary":
            return branch, "primary_boundary_representative", "boundary_representative_only", 5
        return branch, "primary_manual_review", "manual_review_only", 6

    if row["target"] == "Ab_sdAb":
        branch = "secondary_sdAb"
        if cls == "pH_fold_supported_complex_weak":
            return branch, "secondary_pH_fold_supported_complex_weak", "secondary_template_complex_diagnostic_required", 3
        if cls == "pH_fold_supported_complex_unchecked":
            return branch, "secondary_pH_fold_supported_complex_unchecked", "secondary_template_complex_diagnostic_required", 4
        if cls == "pH_supported_route_unchecked":
            return branch, "secondary_pH_supported_route_unchecked", "low_priority_secondary_template", 5
        if action == "retain_boundary":
            return branch, "secondary_boundary_representative", "boundary_representative_only", 7
        return branch, "secondary_manual_review", "manual_review_only", 8

    return "unknown", "manual_review", "manual_review_only", 9


def expansion_policy(row: pd.Series) -> dict[str, object]:
    cls = row["integrated_evidence_class"]
    target = row["target"]

    policy = {
        "allowed_backfill": False,
        "allowed_local_expansion": False,
        "max_family_size": 0,
        "seed_cap": 0,
        "cluster_cap": 0,
        "priority_level": int(row["bank_priority"]),
        "template_status": "not_expansion_template",
        "complex_diagnostic_required": False,
    }

    if target == "Ab_1E62" and cls == "primary_complex_plausible":
        policy.update(
            {
                "allowed_backfill": True,
                "allowed_local_expansion": True,
                "max_family_size": 80,
                "seed_cap": 16,
                "cluster_cap": 5,
                "template_status": "active_primary_template",
            }
        )
    elif target == "Ab_1E62" and cls == "primary_local_pH_supported_complex_unchecked":
        policy.update(
            {
                "allowed_backfill": True,
                "allowed_local_expansion": False,
                "max_family_size": 40,
                "seed_cap": 10,
                "cluster_cap": 4,
                "template_status": "active_primary_backfill_template",
            }
        )
    elif target == "Ab_1E62" and cls == "primary_boundary":
        policy.update(
            {
                "allowed_backfill": True,
                "allowed_local_expansion": False,
                "max_family_size": 10,
                "seed_cap": 4,
                "cluster_cap": 2,
                "template_status": "limited_boundary_template",
            }
        )
    elif target == "Ab_sdAb" and cls == "pH_fold_supported_complex_weak":
        policy.update(
            {
                "allowed_backfill": True,
                "allowed_local_expansion": False,
                "max_family_size": 30,
                "seed_cap": 10,
                "cluster_cap": 3,
                "template_status": "secondary_template_complex_weak",
                "complex_diagnostic_required": True,
            }
        )
    elif target == "Ab_sdAb" and cls == "pH_fold_supported_complex_unchecked":
        policy.update(
            {
                "allowed_backfill": True,
                "allowed_local_expansion": False,
                "max_family_size": 20,
                "seed_cap": 8,
                "cluster_cap": 3,
                "template_status": "secondary_template_complex_unchecked",
                "complex_diagnostic_required": True,
            }
        )
    elif target == "Ab_sdAb" and cls == "pH_supported_route_unchecked":
        policy.update(
            {
                "allowed_backfill": True,
                "allowed_local_expansion": False,
                "max_family_size": 10,
                "seed_cap": 5,
                "cluster_cap": 2,
                "template_status": "low_priority_secondary_template",
                "complex_diagnostic_required": True,
            }
        )
    elif target == "Ab_sdAb" and cls == "secondary_boundary":
        policy.update(
            {
                "allowed_backfill": False,
                "allowed_local_expansion": False,
                "max_family_size": 0,
                "seed_cap": 0,
                "cluster_cap": 0,
                "template_status": "boundary_representative_only",
                "complex_diagnostic_required": True,
            }
        )
    return policy


def add_bank_columns(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        branch, role, expansion_role, priority = bank_role(row)
        his, non_his = split_mutations(row.get("mutation_list"))
        item = row.to_dict()
        item["primary_or_secondary_branch"] = branch
        item["complex_evidence_level"] = row["complex_review_class"]
        item["fold_evidence_level"] = row["fold_review_class_integrated"]
        item["pH_mechanism_level"] = row["electrostatic_mechanism_class"]
        item["glycan_status"] = row["glycan_review_class"]
        item["recommended_bank_role"] = role
        item["recommended_expansion_role"] = expansion_role
        item["bank_priority"] = priority
        item["his_mutations"] = ";".join(his)
        item["rescue_pattern"] = ";".join(non_his) if non_his else "none"
        item["seed_pattern"] = row.get("his_seed_set", "")
        item["mutation_pattern"] = row.get("mutation_list", "")
        item["bank_notes"] = bank_note(pd.Series(item))
        rows.append(item)

    out = pd.DataFrame(rows)
    policies = out.apply(expansion_policy, axis=1, result_type="expand")
    return pd.concat([out, policies], axis=1)


def bank_note(row: pd.Series) -> str:
    notes: list[str] = []
    if row["primary_or_secondary_branch"] == "primary_1E62":
        notes.append("1E62 primary branch")
    if row["primary_or_secondary_branch"] == "secondary_sdAb":
        notes.append("sdAb secondary branch")
    if row["complex_evidence_level"] in {"complex_weak", "complex_very_weak"}:
        notes.append("complex evidence weak; not structure-confirmed")
    if row["complex_evidence_level"] == "complex_unchecked":
        notes.append("complex route unchecked")
    if row["fold_evidence_level"] == "fold_unchecked":
        notes.append("fold route unchecked")
    if row["glycan_status"] == "glycan_unchecked_no_flag":
        notes.append("glycan unchecked; do not report low risk")
    if row["recommended_expansion_role"].endswith("diagnostic_required"):
        notes.append("parent/complex diagnostic recommended before confidence upgrade")
    return "; ".join(notes)


def template_id(row: pd.Series) -> str:
    target_short = "1E62" if row["target"] == "Ab_1E62" else "sdAb"
    digest = short_hash(f"{row['target']}|{row['variant_id']}|{row['mutation_pattern']}")
    return f"tpl_{target_short}_{digest}"


def write_csvs(bank: pd.DataFrame) -> pd.DataFrame:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    bank_cols = [
        "variant_id",
        "target",
        "mutation_list",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "mutation_count",
        "integrated_evidence_class",
        "integrated_action",
        "primary_or_secondary_branch",
        "complex_evidence_level",
        "fold_evidence_level",
        "pH_mechanism_level",
        "glycan_status",
        "recommended_bank_role",
        "recommended_expansion_role",
        "bank_priority",
        "allowed_backfill",
        "allowed_local_expansion",
        "max_family_size",
        "seed_cap",
        "cluster_cap",
        "complex_diagnostic_required",
        "template_status",
        "selection_priority",
        "pyrosetta_rosetta_delta_score",
        "pyrosetta_local_structure_validity_t2_score",
        "pyrosetta_local_clash_count",
        "his_min_antigen_distance",
        "his_pka_support_score",
        "simplefold_median_plddt",
        "af3_iptm",
        "af3_ptm",
        "evidence_conflict_flags",
        "parent_baseline_status",
        "bank_notes",
    ]
    bank[bank_cols].to_csv(OUT_DIR / "tier2_candidate_bank.csv", index=False)

    templates = bank[bank["template_status"].ne("not_expansion_template")].copy()
    templates["template_id"] = templates.apply(template_id, axis=1)
    template_cols = [
        "template_id",
        "target",
        "variant_id",
        "seed_pattern",
        "mutation_pattern",
        "rescue_pattern",
        "integrated_evidence_class",
        "recommended_expansion_role",
        "allowed_backfill",
        "allowed_local_expansion",
        "max_family_size",
        "seed_cap",
        "cluster_cap",
        "priority_level",
        "template_status",
        "complex_diagnostic_required",
        "glycan_status",
        "bank_notes",
    ]
    templates[template_cols].to_csv(OUT_DIR / "tier2_candidate_bank_expansion_templates.csv", index=False)
    return templates


def count_table(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    return (
        df.groupby(keys, dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(keys[:-1] + ["count"], ascending=[True] * max(0, len(keys) - 1) + [False])
    )


def write_reports(bank: pd.DataFrame, templates: pd.DataFrame) -> None:
    by_target = count_table(bank, ["target", "primary_or_secondary_branch", "recommended_bank_role"])
    by_evidence = count_table(bank, ["target", "integrated_evidence_class", "recommended_expansion_role"])
    by_seed = count_table(bank, ["target", "his_seed_set", "recommended_bank_role"])
    by_cluster = count_table(bank, ["target", "near_duplicate_cluster_id", "recommended_bank_role"])
    conflicts = bank[
        bank["evidence_conflict_flags"].fillna("").ne("")
        | bank["complex_diagnostic_required"].astype(bool)
        | bank["glycan_status"].eq("glycan_unchecked_no_flag")
    ].copy()

    (OUT_DIR / "tier2_candidate_bank_by_target.md").write_text(
        "\n".join(["# Tier2 Candidate Bank By Target", "", md_table(by_target), ""]),
        encoding="utf-8",
    )
    (OUT_DIR / "tier2_candidate_bank_by_evidence_class.md").write_text(
        "\n".join(["# Tier2 Candidate Bank By Evidence Class", "", md_table(by_evidence), ""]),
        encoding="utf-8",
    )
    (OUT_DIR / "tier2_candidate_bank_by_seed.md").write_text(
        "\n".join(["# Tier2 Candidate Bank By Seed", "", md_table(by_seed), ""]),
        encoding="utf-8",
    )
    (OUT_DIR / "tier2_candidate_bank_by_cluster.md").write_text(
        "\n".join(
            [
                "# Tier2 Candidate Bank By Cluster",
                "",
                md_table(by_cluster.groupby("target", group_keys=False).head(60)),
                "",
            ]
        ),
        encoding="utf-8",
    )

    conflict_cols = [
        "variant_id",
        "target",
        "mutation_list",
        "integrated_evidence_class",
        "recommended_bank_role",
        "recommended_expansion_role",
        "complex_evidence_level",
        "fold_evidence_level",
        "glycan_status",
        "complex_diagnostic_required",
        "evidence_conflict_flags",
        "bank_notes",
    ]
    conflicts[conflict_cols].to_csv(OUT_DIR / "tier2_candidate_bank_conflict_rows.csv", index=False)
    (OUT_DIR / "tier2_candidate_bank_conflict_report.md").write_text(
        "\n".join(
            [
                "# Tier2 Candidate Bank Conflict / Diagnostic Report",
                "",
                "This report keeps route-unchecked, glycan-unchecked, and sdAb complex-diagnostic requirements explicit. These are not automatic rejects.",
                "",
                "## Conflict / Diagnostic Rows By Target",
                md_table(count_table(conflicts, ["target", "recommended_bank_role"])),
                "",
                "## Diagnostic Flags",
                md_table(count_table(conflicts, ["target", "complex_evidence_level", "glycan_status"])),
                "",
            ]
        ),
        encoding="utf-8",
    )

    report = build_next_gate_report(bank, templates, by_target, by_evidence)
    (OUT_DIR / "tier2_candidate_bank_next_gate_report.md").write_text(report, encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(report, encoding="utf-8")
    (OUT_DIR / "tier2_candidate_bank_manifest.yaml").write_text(
        "\n".join(
            [
                "stage: tier2_candidate_bank_construction",
                "status: complete",
                f"source: {SOURCE.relative_to(ROOT)}",
                f"source_sha256: {sha256_file(SOURCE)}",
                f"bank_rows: {len(bank)}",
                f"template_rows: {len(templates)}",
                "locks:",
                "  final_10k: true",
                "  broad_tier2_heavy: true",
                "  md: true",
                "  sdab_structure_confirmed_upgrade: true",
                "  glycan_low_risk_claim: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def build_next_gate_report(
    bank: pd.DataFrame,
    templates: pd.DataFrame,
    by_target: pd.DataFrame,
    by_evidence: pd.DataFrame,
) -> str:
    branch_summary = (
        bank.groupby(["target", "primary_or_secondary_branch"], dropna=False)
        .agg(
            bank_rows=("variant_id", "size"),
            backfill_allowed=("allowed_backfill", "sum"),
            local_expansion_allowed=("allowed_local_expansion", "sum"),
            complex_diagnostic_required=("complex_diagnostic_required", "sum"),
        )
        .reset_index()
    )
    template_summary = count_table(templates, ["target", "template_status"])
    source_hash = sha256_file(SOURCE)

    return "\n".join(
        [
            "# Tier2 Candidate Bank Construction",
            "",
            "## Executive Summary",
            "",
            "Verdict: `TIER2_CANDIDATE_BANK_BUILT__EXPANSION_TEMPLATE_PLANNING_ALLOWED__FINAL_10K_LOCKED`.",
            "",
            "This stage converts the Integrated Stage B Review into a candidate bank and expansion-template table. It does not run new AF3, SimpleFold, MD, glycan modeling, broad Tier2-heavy, or final 10K selection.",
            "",
            "Key interpretation:",
            "",
            "- `1E62` is the primary branch. Its highest-confidence templates come from `primary_complex_plausible`, followed by local pH-supported but complex-unchecked candidates.",
            "- `sdAb` is retained as a secondary branch. Its candidates remain pH/fold-supported or pH-supported, but complex evidence is weak or unchecked, so they are not structure-confirmed templates.",
            "- Glycan status remains `glycan_unchecked_no_flag`; this is not a low-risk call.",
            "- Parent complex baseline comparison is still missing and should be handled as a small supplemental diagnostic before upgrading sdAb confidence.",
            "",
            "## Source",
            "",
            md_table(
                pd.DataFrame(
                    [
                        {
                            "source": str(SOURCE.relative_to(ROOT)),
                            "rows": len(pd.read_csv(SOURCE)),
                            "sha256": source_hash,
                        }
                    ]
                )
            ),
            "",
            "## Branch Summary",
            "",
            md_table(branch_summary),
            "",
            "## Bank Role Counts",
            "",
            md_table(by_target),
            "",
            "## Evidence / Expansion Role Counts",
            "",
            md_table(by_evidence),
            "",
            "## Expansion Template Summary",
            "",
            md_table(template_summary),
            "",
            "## Next Gate",
            "",
            "Allowed now:",
            "",
            "- Use `tier2_candidate_bank.csv` as the organized candidate mother bank.",
            "- Use `tier2_candidate_bank_expansion_templates.csv` to define evidence-guided backfill and constrained local expansion.",
            "- Keep 1E62 and sdAb in separate primary / secondary branches.",
            "- Prepare a small parent-baseline complex diagnostic, especially for sdAb, before any confidence upgrade.",
            "- Prepare backfill rules against the production / reserve pools using these template families and caps.",
            "",
            "Still locked:",
            "",
            "- Final 10K library selection.",
            "- Broad Tier2-heavy or broad AF3/SimpleFold review.",
            "- MD or constant-pH MD.",
            "- Treating sdAb as structure-confirmed.",
            "- Reporting glycan low risk.",
            "",
            "## Output Files",
            "",
            "- `results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank.csv`",
            "- `results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank_expansion_templates.csv`",
            "- `results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank_by_target.md`",
            "- `results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank_by_evidence_class.md`",
            "- `results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank_by_seed.md`",
            "- `results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank_by_cluster.md`",
            "- `results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank_conflict_report.md`",
            "- `results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank_next_gate_report.md`",
            "",
        ]
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    bank = add_bank_columns(read_source())
    templates = write_csvs(bank)
    write_reports(bank, templates)
    print(f"wrote {OUT_DIR / 'tier2_candidate_bank.csv'} rows={len(bank)}")
    print(f"wrote {OUT_DIR / 'tier2_candidate_bank_expansion_templates.csv'} rows={len(templates)}")
    print(f"wrote {OUT_DIR / 'tier2_candidate_bank_next_gate_report.md'}")


if __name__ == "__main__":
    main()
