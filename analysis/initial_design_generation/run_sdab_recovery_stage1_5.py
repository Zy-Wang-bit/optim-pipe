#!/usr/bin/env python3
"""Recalibrate an sdAb recovery-round Stage-1 result set and update the bank."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.initial_design_generation.build_manual_spot_check_structure_features import (
    load_pka_detail,
    parent_rows,
    structure_feature_row,
)
from analysis.initial_design_generation.run_stage1_5_stage2a_recalibration import add_refined_labels
from analysis.initial_design_generation.run_sdab_recovery_loop import markdown_table, supported_boundary_status


ROOT = Path(__file__).resolve().parents[2]
ROUND = ROOT / "results/initial_design_generation/sdab_recovery_loop/round_01"
BASE_STAGE15 = ROOT / "results/initial_design_generation/stage1_5_stage2a/stage1_refined_structure_risk.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round-dir", default=str(ROUND))
    parser.add_argument("--base-stage15", default=str(BASE_STAGE15))
    parser.add_argument("--candidate-csv", default=str(ROUND / "sdab_stage1_minibatch_list.csv"))
    parser.add_argument("--pyrosetta-results", default=str(ROUND / "stage1/tier2b_full_pyrosetta_results.csv"))
    parser.add_argument("--pka-summary", default=str(ROUND / "stage1/tier2b_full_pka_summary.csv"))
    parser.add_argument("--pka-detail", default=str(ROUND / "stage1/tier2b_full_pka_detail.csv"))
    parser.add_argument("--current-bank", default=str(ROUND / "sdab_candidate_bank.csv"))
    return parser.parse_args()


def csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def num(series: pd.Series | Any) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def build_stage1_table(args: argparse.Namespace) -> pd.DataFrame:
    candidates = pd.read_csv(args.candidate_csv, low_memory=False)
    pyro = pd.read_csv(args.pyrosetta_results, low_memory=False)
    pka = pd.read_csv(args.pka_summary, low_memory=False)
    base = pd.read_csv(args.base_stage15, low_memory=False)
    parent = base[(base["target"].eq("Ab_sdAb")) & (base["control_type"].eq("parent_wt_anchor"))].copy()
    if len(parent) != 1:
        raise RuntimeError(f"Expected one Ab_sdAb parent anchor row, got {len(parent)}")

    pyro = pyro.rename(
        columns={
            "tool_status": "pyrosetta_tool_status",
            "feature_missing_reason": "pyrosetta_feature_missing_reason",
            "feature_confidence": "pyrosetta_feature_confidence",
        }
    )
    df = candidates.merge(pyro, on=["variant_id", "target"], how="left", validate="one_to_one")
    df = df.merge(pka, on=["variant_id", "target"], how="left", validate="one_to_one")
    df["control_type"] = ""
    df["tier2_class"] = df["cheap_gate_status"].map(
        {
            "stage1_minibatch_eligible": "T2_good_candidate",
            "tier1_recovery_candidate": "T2_release_possible_but_neutral_risky",
        }
    ).fillna("T2_release_possible_but_neutral_risky")
    df["tier1_review_class"] = df["cheap_gate_status"].map(
        {
            "stage1_minibatch_eligible": "A_tier1_priority_candidate",
            "tier1_recovery_candidate": "B_rescue_enriched_candidate",
        }
    ).fillna("B_rescue_enriched_candidate")
    df["rosetta_delta"] = df["rosetta_delta_score"]
    df["local_validity_score"] = df["local_structure_validity_t2_score"]
    df["neutral_retention_t2_score"] = df["neutral_retention_score"]
    df["acidic_release_mechanism_t2_score"] = df["acidic_release_support_score"]
    df["global_weakening_risk_t2_score"] = df["global_weakening_risk_score"]
    df["mpnn_compatibility_t2_score"] = df["primary_generation_route"].astype(str).str.contains("MPNN", regex=False).map({True: 0.7, False: 0.45})
    df["foldx_neutral_interface_t2_score"] = df.get("foldx_proxy_feature", 0.0)
    df["glycan_coverage_risk_t2_score"] = 0.0
    df["selection_class_from_tier1_or_stage1"] = df["tier1_review_class"] + "|" + df["tier2_class"]

    parent = parent.copy()
    parent["tier1_review_class"] = parent.get("tier1_review_class", "")
    combined = pd.concat([parent, df], ignore_index=True, sort=False)
    return combined


def recalibrate(combined: pd.DataFrame, pka_detail_path: Path) -> pd.DataFrame:
    parents = parent_rows(combined)
    parent = parents.get("Ab_sdAb")
    if parent is None:
        raise RuntimeError("Missing Ab_sdAb parent anchor for structure feature extraction")
    pka_detail = load_pka_detail(pka_detail_path)
    parent_cache: dict[str, dict[str, object]] = {}
    feature_rows = []
    variants = combined[~combined["control_type"].eq("parent_wt_anchor")].copy()
    for _, row in variants.iterrows():
        feature_rows.append(structure_feature_row(row, row, parent, parent_cache, pka_detail))
    features = pd.DataFrame(feature_rows)
    merged = variants.merge(
        features.drop(columns=[c for c in ["target", "mutation_list"] if c in features.columns]),
        on="variant_id",
        how="left",
        suffixes=("", "_feature"),
    )
    for col in features.columns:
        if col in {"variant_id", "target", "mutation_list"}:
            continue
        feature_col = f"{col}_feature"
        if feature_col in merged.columns:
            merged[col] = merged[feature_col]
            merged = merged.drop(columns=[feature_col])
    refined = add_refined_labels(merged)
    refined["supported_boundary_status"] = refined.apply(supported_boundary_status, axis=1)
    refined["bank_eligibility"] = "reject"
    refined.loc[refined["refined_structure_risk_class"].eq("T2_pass_like_structure"), "bank_eligibility"] = "pass_like"
    refined.loc[
        refined["supported_boundary_status"].astype(str).str.startswith("supported_boundary"),
        "bank_eligibility",
    ] = "supported_boundary"
    refined.loc[refined["refined_structure_risk_class"].eq("T2_severe_structure_risk"), "bank_eligibility"] = "reject_severe"
    return refined


def update_bank(current_bank: pd.DataFrame, refined: pd.DataFrame, round_label: str) -> pd.DataFrame:
    current = current_bank.copy()
    current["bank_source"] = current.get("bank_source", "pre_recovery_stage2a")
    current["canonical_bank_key"] = current.get("sequence", current["variant_id"]).fillna(current["variant_id"]).astype(str)
    additions = refined[refined["bank_eligibility"].isin(["pass_like", "supported_boundary"])].copy()
    additions["bank_source"] = f"{round_label}_stage1_5"
    additions["canonical_bank_key"] = additions.get("sequence", additions["variant_id"]).fillna(additions["variant_id"]).astype(str)
    bank = pd.concat([current, additions], ignore_index=True, sort=False)
    bank = bank.sort_values(
        ["bank_eligibility", "refined_structure_risk_class", "neutral_retention_score"],
        ascending=[True, True, False],
    )
    bank = bank.drop_duplicates("canonical_bank_key", keep="first").reset_index(drop=True)
    return bank


def apply_bank_cluster_cap(bank: pd.DataFrame, max_cluster_fraction: float = 0.12) -> tuple[pd.DataFrame, pd.DataFrame]:
    if bank.empty or "near_duplicate_cluster_id" not in bank.columns:
        return bank, bank.head(0).copy()
    kept = bank.copy()
    dropped_rows = []
    while len(kept) > 1000:
        counts = kept.groupby("near_duplicate_cluster_id", dropna=False).size()
        if counts.empty:
            break
        top_cluster = counts.idxmax()
        top_fraction = int(counts.max()) / max(1, len(kept))
        if top_fraction <= max_cluster_fraction:
            break
        cluster_rows = kept[kept["near_duplicate_cluster_id"].eq(top_cluster)].copy()
        cluster_rows["_is_pass_like"] = cluster_rows["refined_structure_risk_class"].eq("T2_pass_like_structure").astype(int)
        cluster_rows["_neutral"] = num(cluster_rows.get("neutral_retention_score", pd.Series(0, index=cluster_rows.index))).fillna(0)
        cluster_rows = cluster_rows.sort_values(
            ["_is_pass_like", "_neutral", "variant_id"],
            ascending=[True, True, True],
        )
        drop_idx = cluster_rows.index[0]
        row = kept.loc[[drop_idx]].copy()
        row["bank_drop_reason"] = f"near_duplicate_cluster_cap>{max_cluster_fraction:.2f}"
        dropped_rows.append(row)
        kept = kept.drop(index=drop_idx).reset_index(drop=True)
    dropped = pd.concat(dropped_rows, ignore_index=True, sort=False) if dropped_rows else bank.head(0).copy()
    return kept.reset_index(drop=True), dropped


def gate_metrics(bank: pd.DataFrame, refined: pd.DataFrame) -> dict[str, Any]:
    total = len(bank)
    pass_like = int(bank["refined_structure_risk_class"].eq("T2_pass_like_structure").sum()) if total else 0
    boundary = int(bank["bank_eligibility"].eq("supported_boundary").sum()) if total else 0
    severe = int(bank["refined_structure_risk_class"].eq("T2_severe_structure_risk").sum()) if total else 0
    unsupported = int(bank["supported_boundary_status"].eq("unsupported_boundary").sum()) if "supported_boundary_status" in bank else 0
    top_seed = int(bank.groupby("his_seed_set", dropna=False).size().max()) if total else 0
    top_cluster = int(bank.groupby("near_duplicate_cluster_id", dropna=False).size().max()) if total and "near_duplicate_cluster_id" in bank else 0
    four_mut = int(num(bank["mutation_count"]).ge(4).sum()) if total else 0
    round_total = len(refined)
    round_severe = int(refined["refined_structure_risk_class"].eq("T2_severe_structure_risk").sum())
    round_pass_like = int(refined["refined_structure_risk_class"].eq("T2_pass_like_structure").sum())
    return {
        "bank_unique_count": total,
        "bank_pass_like_count": pass_like,
        "bank_pass_like_fraction": pass_like / max(1, total),
        "bank_supported_boundary_count": boundary,
        "bank_boundary_fraction": boundary / max(1, total),
        "bank_severe_count": severe,
        "bank_unsupported_boundary_count": unsupported,
        "bank_top_seed_fraction": top_seed / max(1, total),
        "bank_top_near_duplicate_cluster_fraction": top_cluster / max(1, total),
        "bank_4mut_fraction": four_mut / max(1, total),
        "round_total": round_total,
        "round_pass_like_count": round_pass_like,
        "round_pass_like_fraction": round_pass_like / max(1, round_total),
        "round_severe_count": round_severe,
        "round_severe_fraction": round_severe / max(1, round_total),
        "round_supported_boundary_count": int(refined["bank_eligibility"].eq("supported_boundary").sum()),
    }


def gate_table(metrics: dict[str, Any]) -> pd.DataFrame:
    minimal = {
        "bank_unique_count>=600": metrics["bank_unique_count"] >= 600,
        "bank_pass_like_fraction>=0.15": metrics["bank_pass_like_fraction"] >= 0.15,
        "bank_severe_count==0": metrics["bank_severe_count"] == 0,
        "bank_unsupported_boundary_count==0": metrics["bank_unsupported_boundary_count"] == 0,
        "bank_top_seed_fraction<=0.30": metrics["bank_top_seed_fraction"] <= 0.30,
        "bank_top_near_duplicate_cluster_fraction<=0.15": metrics["bank_top_near_duplicate_cluster_fraction"] <= 0.15,
        "bank_4mut_fraction<=0.15": metrics["bank_4mut_fraction"] <= 0.15,
    }
    full = {
        "bank_unique_count>=1000": metrics["bank_unique_count"] >= 1000,
        "bank_pass_like_count>=200_or_fraction>=0.20": metrics["bank_pass_like_count"] >= 200 or metrics["bank_pass_like_fraction"] >= 0.20,
        "bank_severe_count==0": metrics["bank_severe_count"] == 0,
        "bank_unsupported_boundary_count==0": metrics["bank_unsupported_boundary_count"] == 0,
        "bank_boundary_fraction<=0.70": metrics["bank_boundary_fraction"] <= 0.70,
        "bank_top_seed_fraction<=0.25": metrics["bank_top_seed_fraction"] <= 0.25,
        "bank_top_near_duplicate_cluster_fraction<=0.12": metrics["bank_top_near_duplicate_cluster_fraction"] <= 0.12,
        "bank_4mut_fraction<=0.10": metrics["bank_4mut_fraction"] <= 0.10,
    }
    rows = []
    for scope, checks in [("minimal_sync_gate", minimal), ("full_stage2a_gate", full)]:
        for check, ok in checks.items():
            rows.append({"scope": scope, "check": check, "status": "PASS" if ok else "FAIL"})
    return pd.DataFrame(rows)


def write_reports(round_dir: Path, refined: pd.DataFrame, bank: pd.DataFrame, metrics: dict[str, Any], gates: pd.DataFrame, round_label: str) -> None:
    metric_table = pd.DataFrame([{"metric": k, "value": v} for k, v in metrics.items()])
    round_classes = refined.groupby(["refined_structure_risk_class", "bank_eligibility"], dropna=False).size().reset_index(name="count")
    bank_classes = bank.groupby(["refined_structure_risk_class", "bank_eligibility"], dropna=False).size().reset_index(name="count")
    minimal_pass = gates[gates["scope"].eq("minimal_sync_gate")]["status"].eq("PASS").all()
    full_pass = gates[gates["scope"].eq("full_stage2a_gate")]["status"].eq("PASS").all()
    lines = [
        f"# sdAb Recovery {round_label} Iteration Report",
        "",
        "## Decision",
        "",
        f"Minimal sync gate: `{'PASS' if minimal_pass else 'NOT_PASSED'}`",
        f"Full Stage-2A gate: `{'PASS' if full_pass else 'NOT_PASSED'}`",
        "",
        "This report does not unlock Stage-2A compute automatically.",
        "",
        f"## {round_label} Refined Classes",
        "",
        markdown_table(round_classes),
        "",
        "## Cumulative Bank Classes",
        "",
        markdown_table(bank_classes),
        "",
        "## Gate Metrics",
        "",
        markdown_table(metric_table),
        "",
        "## Gate Checks",
        "",
        markdown_table(gates),
    ]
    text = "\n".join(lines) + "\n"
    (round_dir / "sdab_recovery_iteration_report.md").write_text(text, encoding="utf-8")
    (round_dir / "sdab_stage2a_gate_report.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    round_dir = Path(args.round_dir)
    round_label = round_dir.name
    combined = build_stage1_table(args)
    csv(combined, round_dir / f"stage1/sdab_{round_label}_stage1_combined_input.csv")
    refined = recalibrate(combined, Path(args.pka_detail))
    csv(refined, round_dir / "sdab_stage1_5_refined_results.csv")
    current_bank = pd.read_csv(args.current_bank, low_memory=False)
    bank = update_bank(current_bank, refined, round_label)
    bank, dropped = apply_bank_cluster_cap(bank)
    csv(dropped, round_dir / "sdab_candidate_bank_diversity_dropped.csv")
    csv(bank, round_dir / "sdab_stage2a_candidate_bank.csv")
    # Preserve the plan's expected cumulative-bank filename as an updated alias.
    csv(bank, round_dir / "sdab_candidate_bank.csv")
    metrics = gate_metrics(bank, refined)
    gates = gate_table(metrics)
    csv(pd.DataFrame([{"metric": k, "value": v} for k, v in metrics.items()]), round_dir / "sdab_stage2a_gate_metrics.csv")
    csv(gates, round_dir / "sdab_stage2a_gate_checks.csv")
    write_reports(round_dir, refined, bank, metrics, gates, round_label)
    print(
        {
            "status": "complete",
            "round": round_label,
            "round_rows": len(refined),
            "bank_rows": len(bank),
            "minimal_sync_gate": "PASS" if gates[gates["scope"].eq("minimal_sync_gate")]["status"].eq("PASS").all() else "NOT_PASSED",
            "full_stage2a_gate": "PASS" if gates[gates["scope"].eq("full_stage2a_gate")]["status"].eq("PASS").all() else "NOT_PASSED",
        }
    )


if __name__ == "__main__":
    main()
