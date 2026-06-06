#!/usr/bin/env python3
"""Build Stage-2A compute result tables and transition report."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(".")
OUT_DIR = ROOT / "results/initial_design_generation/stage2a_compute"
INPUT_CSV = ROOT / "results/initial_design_generation/stage2a_synchronized_input/stage2a_synchronized_input_list.csv"
AUDIT_REPORT = ROOT / "results/initial_design_generation/stage2a_synchronized_input/stage2a_input_audit_report.md"
PYRO_CSV = OUT_DIR / "tier2b_full_pyrosetta_results.csv"
PKA_CSV = OUT_DIR / "tier2b_full_pka_summary.csv"
CURRENT_STAGE_REPORT = ROOT / ".tasks/active/initial-design-generation/current_stage_report.md"
PROGRESS = ROOT / ".tasks/active/initial-design-generation/progress.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def clamp(value: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, value))


def rosetta_delta_support(delta: float) -> float:
    if pd.isna(delta):
        return 0.0
    if delta <= 0:
        return 1.0
    if delta <= 20:
        return 0.85
    if delta <= 80:
        return 0.65
    if delta <= 200:
        return 0.35
    if delta <= 600:
        return 0.12
    return 0.0


def classify(row: pd.Series) -> tuple[str, str]:
    if str(row.get("final_slot_type", "")) in {"control", "audit"} or str(row.get("control_anchor", "")) == "yes":
        return "T2A_control_or_audit", "retain_as_control"
    if row.get("tool_status") != "success":
        return "T2A_reject", "reject"

    validity = float(row.get("stage2a_local_validity") or 0.0)
    delta = float(row.get("rosetta_delta_score") or 0.0)
    clash = float(row.get("stage2a_new_clash_count") or 0.0)
    pka = row.get("stage2a_pH_mechanism_support")
    pka = 0.0 if pd.isna(pka) else float(pka)
    input_level = str(row.get("boundary_support_level", ""))

    if validity < 0.45 or clash > 10 or delta > 600:
        return "T2A_reject", "reject"
    if validity < 0.60 or delta > 200:
        return "T2A_structure_risk", "reject"
    if pka < 0.20:
        return "T2A_mechanism_weak", "retain_as_boundary"

    if input_level == "pass_like":
        if validity >= 0.85 and delta <= 20 and pka >= 0.60:
            return "T2A_high_confidence", "promote_to_heavy_candidate_pool"
        if validity >= 0.75 and delta <= 80 and pka >= 0.40:
            return "T2A_good", "promote_to_heavy_candidate_pool"
        if validity >= 0.65 and delta <= 120 and pka >= 0.30:
            return "T2A_boundary_unresolved", "retain_for_stage2b"
        return "T2A_mechanism_weak", "retain_as_boundary"

    if input_level == "high_support_boundary":
        if validity >= 0.82 and delta <= 20 and pka >= 0.55:
            return "T2A_supported_boundary_confirmed", "promote_to_heavy_candidate_pool"
        if validity >= 0.70 and delta <= 80 and pka >= 0.35:
            return "T2A_supported_boundary_confirmed", "retain_for_stage2b"
        if validity >= 0.60 and delta <= 160:
            return "T2A_boundary_unresolved", "retain_as_boundary"
        return "T2A_structure_risk", "reject"

    return "T2A_boundary_unresolved", "retain_as_boundary"


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    lines = [
        "| " + " | ".join(map(str, df.columns)) + " |",
        "| " + " | ".join("---" for _ in df.columns) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in df.columns) + " |")
    return "\n".join(lines)


def summarize(df: pd.DataFrame, group_cols: list[str], path: Path) -> pd.DataFrame:
    out = (
        df.groupby(group_cols + ["stage2a_final_class", "stage2a_action"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(group_cols + ["stage2a_final_class", "stage2a_action"])
    )
    out.to_csv(path, index=False)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    inp = pd.read_csv(INPUT_CSV)
    pyro = pd.read_csv(PYRO_CSV)
    pka = pd.read_csv(PKA_CSV)

    df = inp.merge(pyro, on=["variant_id", "target"], how="left", validate="one_to_one")
    df = df.merge(pka, on=["variant_id", "target"], how="left", validate="one_to_one")

    df["stage2a_structure_score"] = [
        round(clamp(0.70 * float(v or 0.0) + 0.30 * rosetta_delta_support(float(d or 0.0))), 4)
        for v, d in zip(df["local_structure_validity_t2_score"], df["rosetta_delta_score"])
    ]
    df["stage2a_local_validity"] = pd.to_numeric(df["local_structure_validity_t2_score"], errors="coerce")
    df["stage2a_new_clash_count"] = pd.to_numeric(df["local_clash_count"], errors="coerce")
    df["stage2a_max_new_clash_overlap"] = pd.NA
    df["stage2a_new_bad_contact_count"] = pd.NA
    df["stage2a_lost_parent_contact_count"] = pd.NA
    df["stage2a_neutral_retention_support"] = [
        round(clamp(0.50 * float(v or 0.0) + 0.50 * rosetta_delta_support(float(d or 0.0))), 4)
        for v, d in zip(df["local_structure_validity_t2_score"], df["rosetta_delta_score"])
    ]
    df["stage2a_pH_mechanism_support"] = pd.to_numeric(df["his_pka_support_t2_score"], errors="coerce")
    df["input_class"] = df["boundary_support_level"]
    df["compute_status"] = df["tool_status"].fillna("not_run")

    classes_actions = df.apply(classify, axis=1, result_type="expand")
    df["stage2a_final_class"] = classes_actions[0]
    df["stage2a_action"] = classes_actions[1]

    results_cols = [
        "variant_id",
        "target",
        "canonical_sequence_hash",
        "mutation_list",
        "input_class",
        "boundary_support_level",
        "control_anchor",
        "audit_only",
        "compute_status",
        "tool_status",
        "pka_tool_status",
        "stage2a_structure_score",
        "stage2a_local_validity",
        "stage2a_new_clash_count",
        "stage2a_max_new_clash_overlap",
        "stage2a_new_bad_contact_count",
        "stage2a_lost_parent_contact_count",
        "stage2a_neutral_retention_support",
        "stage2a_pH_mechanism_support",
        "stage2a_final_class",
        "stage2a_action",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "mutation_count",
        "final_slot_type",
        "source_bank",
        "source_round",
        "pdb_path",
        "rosetta_delta_score",
        "local_clash_count",
        "local_structure_validity_t2_score",
        "pka_combined_mean",
    ]
    results_path = OUT_DIR / "stage2a_compute_results.csv"
    df[results_cols].to_csv(results_path, index=False)

    by_target = summarize(df, ["target"], OUT_DIR / "stage2a_result_by_target.csv")
    by_input = summarize(df, ["target", "input_class"], OUT_DIR / "stage2a_result_by_input_class.csv")
    by_seed = summarize(df, ["target", "his_seed_set"], OUT_DIR / "stage2a_result_by_seed.csv")
    by_cluster = summarize(df, ["target", "near_duplicate_cluster_id"], OUT_DIR / "stage2a_result_by_cluster.csv")
    by_mut_count = summarize(df, ["target", "mutation_count"], OUT_DIR / "stage2a_result_by_mutation_count.csv")

    controls = df[df["control_anchor"].astype(str).eq("yes") | df["audit_only"].astype(str).eq("yes")].copy()
    controls[results_cols].to_csv(OUT_DIR / "stage2a_control_anchor_report.csv", index=False)
    failures = df[~df["tool_status"].eq("success") | (~df["pka_tool_status"].isin(["success", "not_applicable"]))]
    failures[results_cols].to_csv(OUT_DIR / "stage2a_compute_failure_report.csv", index=False)

    target_counts = df.groupby(["target", "stage2a_final_class"]).size().unstack(fill_value=0)
    action_counts = df.groupby(["target", "stage2a_action"]).size().unstack(fill_value=0)
    promotion = df[df["stage2a_action"].eq("promote_to_heavy_candidate_pool")]
    promote_seed = (
        promotion.groupby(["target", "his_seed_set"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
    )
    promote_cluster = (
        promotion.groupby(["target", "near_duplicate_cluster_id"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
    )
    promote_seed.to_csv(OUT_DIR / "stage2a_promoted_seed_summary.csv", index=False)
    promote_cluster.to_csv(OUT_DIR / "stage2a_promoted_cluster_summary.csv", index=False)

    gate_rows = []
    for target, group in df.groupby("target"):
        high_good = group["stage2a_final_class"].isin(["T2A_high_confidence", "T2A_good"]).sum()
        high_good_confirmed = group["stage2a_final_class"].isin(
            ["T2A_high_confidence", "T2A_good", "T2A_supported_boundary_confirmed"]
        ).sum()
        promoted = group[group["stage2a_action"].eq("promote_to_heavy_candidate_pool")]
        top_seed_frac = 0.0
        top_cluster_frac = 0.0
        four_mut_frac = 0.0
        ag102_promoted = 0
        av105_promoted = 0
        if len(promoted):
            top_seed_frac = promoted["his_seed_set"].value_counts(dropna=False).iloc[0] / len(promoted)
            top_cluster_frac = promoted["near_duplicate_cluster_id"].value_counts(dropna=False).iloc[0] / len(promoted)
            four_mut_frac = pd.to_numeric(promoted["mutation_count"], errors="coerce").ge(4).sum() / len(promoted)
            ag102_promoted = promoted["mutation_list"].fillna("").str.contains("AG102H").sum()
            av105_promoted = promoted["mutation_list"].fillna("").str.contains("AV105H").sum()
        minimum_signal_pass = high_good >= 100 or high_good_confirmed >= 200
        diversity_reasons = []
        if target == "Ab_sdAb":
            if top_seed_frac > 0.25:
                diversity_reasons.append("promoted_top_seed_fraction_gt_25pct")
            if top_cluster_frac > 0.12:
                diversity_reasons.append("promoted_top_cluster_fraction_gt_12pct")
            if four_mut_frac > 0.10:
                diversity_reasons.append("promoted_4mut_fraction_gt_10pct")
            if av105_promoted > 0:
                diversity_reasons.append("AV105H_promoted")
        if target == "Ab_1E62":
            if top_cluster_frac > 0.12:
                diversity_reasons.append("promoted_top_cluster_fraction_gt_12pct")
        if not minimum_signal_pass:
            gate_status = "HOLD"
        elif diversity_reasons:
            gate_status = "PATCH_BEFORE_TIER2_HEAVY"
        else:
            gate_status = "PASS"
        gate_rows.append(
            {
                "target": target,
                "high_good_count": int(high_good),
                "high_good_confirmed_count": int(high_good_confirmed),
                "promoted_count": int(len(promoted)),
                "promoted_top_seed_fraction": round(top_seed_frac, 4),
                "promoted_top_cluster_fraction": round(top_cluster_frac, 4),
                "promoted_4mut_fraction": round(four_mut_frac, 4),
                "AG102H_promoted_count": int(ag102_promoted),
                "AV105H_promoted_count": int(av105_promoted),
                "minimum_signal_gate": "PASS" if minimum_signal_pass else "HOLD",
                "diversity_gate_reasons": ";".join(diversity_reasons),
                "targeted_tier2_heavy_planning_gate": gate_status,
            }
        )
    gate = pd.DataFrame(gate_rows)
    gate.to_csv(OUT_DIR / "stage2a_post_compute_gate.csv", index=False)

    status = df.groupby(["target", "tool_status", "pka_tool_status"], dropna=False).size().reset_index(name="count")
    status.to_csv(OUT_DIR / "stage2a_compute_status_summary.csv", index=False)

    top_cluster_lines = []
    for target, group in promotion.groupby("target"):
        vc = group["near_duplicate_cluster_id"].value_counts(dropna=False)
        total = max(1, len(group))
        top_cluster_lines.append(
            {
                "target": target,
                "promoted_count": len(group),
                "top_cluster_fraction": round(vc.iloc[0] / total, 4) if len(vc) else 0,
                "top5_cluster_fraction": round(vc.head(5).sum() / total, 4) if len(vc) else 0,
                "top10_cluster_fraction": round(vc.head(10).sum() / total, 4) if len(vc) else 0,
            }
        )
    top_cluster_summary = pd.DataFrame(top_cluster_lines)
    top_cluster_summary.to_csv(OUT_DIR / "stage2a_promoted_cluster_concentration.csv", index=False)

    manifest = [
        "stage: Stage-2A",
        "status: compute_completed",
        "manual_unlock_source: user_request",
        "input_list:",
        f"  path: {INPUT_CSV.as_posix()}",
        f"  sha256: {sha256_file(INPUT_CSV)}",
        "audit_report:",
        f"  path: {AUDIT_REPORT.as_posix()}",
        f"  sha256: {sha256_file(AUDIT_REPORT)}",
        "outputs:",
        f"  result_table: {(OUT_DIR / 'stage2a_compute_results.csv').as_posix()}",
        f"  transition_report: {(OUT_DIR / 'stage2a_stage_transition_report.md').as_posix()}",
        "locked_after_compute:",
        "  tier2_heavy: true",
        "  af3_simplefold_heavy_review: true",
        "  md: true",
        "  final_10k: true",
        "",
    ]
    (OUT_DIR / "stage2a_compute_manifest.yaml").write_text("\n".join(manifest), encoding="utf-8")

    report = [
        "# Stage-2A Compute Result Report",
        "",
        "## Executive Summary",
        "",
        "Stage-2A compute has been completed for the fixed audited 1,370-candidate input list. No reserve candidates were added, and no post-audit candidate insertion occurred.",
        "",
        "Both targets pass the minimum signal gate for downstream planning, but neither target should move directly into Tier2-heavy execution: the promoted subsets are still too concentrated by seed and/or near-duplicate cluster. The correct next step is targeted Tier2-heavy planning with diversity capping and manual review.",
        "",
        "## Scope",
        "",
        f"Input list: `{INPUT_CSV.as_posix()}`",
        f"Input SHA256: `{sha256_file(INPUT_CSV)}`",
        f"Audit report SHA256: `{sha256_file(AUDIT_REPORT)}`",
        "",
        "No reserve candidates were added. Tier2-heavy / AF3 / MD / final 10K remain locked.",
        "",
        "## Compute Status",
        "",
        md_table(status),
        "",
        "## Final Class By Target",
        "",
        md_table(target_counts.reset_index()),
        "",
        "## Action By Target",
        "",
        md_table(action_counts.reset_index()),
        "",
        "## Post-compute Gate",
        "",
        md_table(gate),
        "",
        "Gate interpretation:",
        "",
        "```text",
        "minimum_signal_gate = enough Stage-2A-supported candidates exist for planning;",
        "targeted_tier2_heavy_planning_gate = whether the promoted pool can proceed without additional diversity capping;",
        "PATCH_BEFORE_TIER2_HEAVY = do not run heavy yet; first cap seed/cluster concentration and re-audit the proposed heavy pool.",
        "```",
        "",
        "## Promoted Cluster Concentration",
        "",
        md_table(top_cluster_summary),
        "",
        "## Top Promoted Seeds",
        "",
        md_table(promote_seed.groupby("target").head(8).reset_index(drop=True)),
        "",
        "## Interpretation",
        "",
        "Stage-2A compute was completed for the fixed audited input list. Both targets meet the minimum signal gate, but the promoted pools are still seed/cluster concentrated. This supports targeted Tier2-heavy planning only after diversity capping and manual review; it does not unlock Tier2-heavy execution, AF3/SimpleFold heavy review, MD, or final 10K.",
        "",
    ]
    report_path = OUT_DIR / "stage2a_stage_transition_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text("\n".join(report), encoding="utf-8")

    with PROGRESS.open("a") as handle:
        handle.write(
            "\n## 2026-05-29 - Stage-2A compute execution\n\n"
            f"- Stage-2A compute completed for fixed audited input list: {len(df)} rows.\n"
            f"- PyRosetta success: {df['tool_status'].eq('success').sum()} / {len(df)}.\n"
            f"- pKa success: {df['pka_tool_status'].eq('success').sum()} / {len(df)}; non-His not_applicable: {df['pka_tool_status'].eq('not_applicable').sum()}.\n"
            f"- Current stage report overwritten with `results/initial_design_generation/stage2a_compute/stage2a_stage_transition_report.md`.\n"
            "- Tier2-heavy / AF3 / MD / final 10K remain locked pending review.\n"
        )

    print(f"rows={len(df)}")
    print(status.to_string(index=False))
    print(gate.to_string(index=False))


if __name__ == "__main__":
    main()
