#!/usr/bin/env python3
"""Build Tier2 Stage-1 diagnostic tables and calibration notes.

This is a table-only Stage 1.5 audit. It does not launch new structural tools.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_INPUT = Path("results/initial_design_generation/tier2_staged/full_stage1/tier2b_stage1_results.csv")
DEFAULT_OUT = Path("results/initial_design_generation/tier2_stage1_diagnostics")

STRONG_GOOD = {"T2_strong_candidate", "T2_good_candidate"}

PRIORITY_SEEDS = {
    "Ab_1E62": {
        "LQ35H;LY38H",
        "LK24H;LY38H",
        "LK24H;LQ35H",
        "LY31H;LQ35H",
    },
    "Ab_sdAb": {
        "AD110H;AY111H",
    },
}

SECONDARY_SEEDS = {
    "Ab_1E62": set(),
    "Ab_sdAb": {
        "AQ100H;AD110H",
        "AQ100H;AY111H",
    },
}

REVIEW_ROUTE_PRIORITY = {
    "Ab_1E62": {
        "His_plus_rescue",
        "His_rule",
        "ProteinMPNN_seeded_rescue",
        "wetlab_informed_expansion",
        "structure_or_interface_guided",
    },
    "Ab_sdAb": {
        "His_plus_rescue",
        "ProteinMPNN_seeded_rescue",
        "wetlab_informed_expansion",
        "structure_or_interface_guided",
    },
}


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(pd.NA, index=df.index, dtype="Float64")
    return pd.to_numeric(df[col], errors="coerce")


def bool_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].astype(str).str.lower().isin({"true", "1", "yes", "y"})


def safe_col(df: pd.DataFrame, col: str, default: str = "missing") -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index)
    return df[col].fillna(default).astype(str).replace({"": default, "nan": default})


def add_rate_columns(grouped: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    out = grouped.copy()
    total = out.groupby(group_cols, dropna=False)["count"].transform("sum")
    out["fraction_within_group"] = (out["count"] / total).round(6)
    strong_good = out["tier2_class"].isin(STRONG_GOOD)
    sg = out[strong_good].groupby(group_cols, dropna=False)["count"].sum().rename("strong_good_count")
    risk = out[out["tier2_class"].eq("T2_structure_risk")].groupby(group_cols, dropna=False)["count"].sum().rename("structure_risk_count")
    totals = out.groupby(group_cols, dropna=False)["count"].sum().rename("total_count")
    joined = pd.concat([totals, sg, risk], axis=1).fillna(0).reset_index()
    joined["strong_good_rate"] = (joined["strong_good_count"] / joined["total_count"]).round(6)
    joined["structure_risk_rate"] = (joined["structure_risk_count"] / joined["total_count"]).round(6)
    return out.merge(joined[group_cols + ["total_count", "strong_good_rate", "structure_risk_rate"]], on=group_cols, how="left")


def class_table(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    grouped = df.groupby(group_cols + ["tier2_class"], dropna=False).size().reset_index(name="count")
    return add_rate_columns(grouped, group_cols).sort_values(group_cols + ["count"], ascending=[True] * len(group_cols) + [False])


def summarize_features(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    features = [
        "neutral_retention_t2_score",
        "acidic_release_mechanism_t2_score",
        "global_weakening_risk_t2_score",
        "local_structure_validity_t2_score",
        "local_clash_count",
        "rosetta_delta_score",
        "his_pka_support_t2_score",
        "foldx_ddG_vs_parent_pH7_4",
        "foldx_neutral_interface_t2_score",
        "mutation_count",
        "His_count",
        "rescue_count",
    ]
    rows: list[dict[str, object]] = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = dict(zip(group_cols, keys))
        base["count"] = len(sub)
        for col in features:
            vals = numeric(sub, col).dropna()
            if vals.empty:
                base[f"{col}_median"] = pd.NA
                base[f"{col}_p25"] = pd.NA
                base[f"{col}_p75"] = pd.NA
            else:
                base[f"{col}_median"] = round(float(vals.median()), 6)
                base[f"{col}_p25"] = round(float(vals.quantile(0.25)), 6)
                base[f"{col}_p75"] = round(float(vals.quantile(0.75)), 6)
        rows.append(base)
    return pd.DataFrame(rows)


def add_structure_risk_drivers(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    local = numeric(out, "local_structure_validity_t2_score")
    clash = numeric(out, "local_clash_count").fillna(0)
    delta = numeric(out, "rosetta_delta_score").fillna(0)
    pyro_bad = ~safe_col(out, "pyrosetta_tool_status", "missing").eq("success")
    low_local = local.lt(0.35)
    high_clash = clash.ge(10)
    high_delta = delta.gt(600)

    labels: list[str] = []
    for a, b, c, d in zip(pyro_bad, low_local.fillna(False), high_clash, high_delta):
        parts = []
        if a:
            parts.append("pyrosetta_not_success")
        if b:
            parts.append("local_validity_lt_0.35")
        if c:
            parts.append("local_clash_ge_10")
        if d:
            parts.append("rosetta_delta_gt_600")
        labels.append(";".join(parts) if parts else "no_structure_risk_trigger")
    out["structure_risk_driver"] = labels

    relaxed_low_local = local.lt(0.30)
    relaxed_high_clash = clash.ge(12)
    relaxed_high_delta = delta.gt(800)
    relaxed_trigger = pyro_bad | relaxed_low_local.fillna(False) | relaxed_high_clash | relaxed_high_delta
    out["structure_risk_under_relaxed_thresholds"] = relaxed_trigger
    out["structure_risk_near_threshold_only"] = out["tier2_class"].eq("T2_structure_risk") & ~relaxed_trigger
    return out


def seed_policy(target: str, seed: str) -> str:
    if seed in PRIORITY_SEEDS.get(target, set()):
        return "priority_seed"
    if seed in SECONDARY_SEEDS.get(target, set()):
        return "secondary_seed"
    if target == "Ab_sdAb" and "AE108H" in seed:
        return "deprioritize_ae108h_seed"
    if seed in {"", "missing", "none", "nan"}:
        return "no_his_seed"
    return "unprioritized_seed"


def add_calibration_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Add decision-oriented Stage-1.5 labels without changing Tier2 classes."""

    out = df.copy()
    local = numeric(out, "local_structure_validity_t2_score")
    clash = numeric(out, "local_clash_count").fillna(0)
    delta = numeric(out, "rosetta_delta_score").fillna(0)
    pyro_bad = ~safe_col(out, "pyrosetta_tool_status", "missing").eq("success")
    low_local = local.lt(0.35).fillna(False)
    high_clash = clash.ge(10)
    high_delta = delta.gt(600)

    severe = (
        pyro_bad
        | local.lt(0.30).fillna(False)
        | clash.ge(12)
        | delta.gt(800)
        | (low_local & high_clash & high_delta)
    )
    local_clash_boundary = low_local & high_clash & ~severe
    delta_only_boundary = high_delta & ~low_local & ~high_clash & delta.le(800)
    risk = out["tier2_class"].eq("T2_structure_risk")

    out["stage1_5_structure_risk_band"] = "not_structure_risk"
    out.loc[risk & severe, "stage1_5_structure_risk_band"] = "severe_structure_risk"
    out.loc[risk & local_clash_boundary, "stage1_5_structure_risk_band"] = "local_clash_boundary_review"
    out.loc[risk & delta_only_boundary, "stage1_5_structure_risk_band"] = "rosetta_delta_only_boundary_review"
    out.loc[
        risk
        & out["stage1_5_structure_risk_band"].eq("not_structure_risk"),
        "stage1_5_structure_risk_band",
    ] = "other_boundary_review"

    out["stage1_5_seed_policy"] = [
        seed_policy(target, seed)
        for target, seed in zip(safe_col(out, "target"), safe_col(out, "his_seed_set"))
    ]
    out["stage1_5_route_supported"] = [
        route in REVIEW_ROUTE_PRIORITY.get(target, set())
        for target, route in zip(safe_col(out, "target"), safe_col(out, "primary_generation_route"))
    ]
    out["stage1_5_mutation_order_policy"] = "standard_review"
    mut_count = numeric(out, "mutation_count")
    out.loc[out["target"].eq("Ab_1E62") & mut_count.le(3), "stage1_5_mutation_order_policy"] = "preferred_order"
    out.loc[
        out["target"].eq("Ab_1E62")
        & mut_count.eq(4)
        & out["stage1_5_seed_policy"].eq("priority_seed"),
        "stage1_5_mutation_order_policy",
    ] = "allowed_only_with_priority_seed"
    out.loc[out["target"].eq("Ab_1E62") & mut_count.gt(4), "stage1_5_mutation_order_policy"] = "hold_high_order"
    out.loc[out["target"].eq("Ab_sdAb") & mut_count.le(3), "stage1_5_mutation_order_policy"] = "preferred_order"
    out.loc[out["target"].eq("Ab_sdAb") & mut_count.eq(4), "stage1_5_mutation_order_policy"] = "caution_4mut"
    out.loc[out["target"].eq("Ab_sdAb") & mut_count.gt(4), "stage1_5_mutation_order_policy"] = "hold_high_order"

    foldx_success = safe_col(out, "foldx_tool_status").eq("success")
    out["stage1_5_foldx_neutral_retained_proxy"] = (
        bool_series(out, "tier2_foldx_subset_selected")
        & foldx_success
        & numeric(out, "foldx_ddG_vs_parent_pH7_4").le(2.0)
    )

    neutral = numeric(out, "neutral_retention_t2_score").fillna(0)
    release = numeric(out, "acidic_release_mechanism_t2_score").fillna(0)
    weakening = numeric(out, "global_weakening_risk_t2_score").fillna(1)
    local_score = numeric(out, "local_structure_validity_t2_score").fillna(0)
    support = (
        neutral
        + release
        + local_score
        - weakening
        + out["stage1_5_route_supported"].astype(int) * 0.25
        + out["stage1_5_foldx_neutral_retained_proxy"].astype(int) * 0.50
    )
    support += out["stage1_5_seed_policy"].map(
        {
            "priority_seed": 0.75,
            "secondary_seed": 0.35,
            "unprioritized_seed": 0.0,
            "no_his_seed": -0.25,
            "deprioritize_ae108h_seed": -0.50,
        }
    ).fillna(0)
    out["stage1_5_review_support_score"] = support.round(6)

    reviewable_boundary = risk & out["stage1_5_structure_risk_band"].isin(
        {
            "rosetta_delta_only_boundary_review",
            "local_clash_boundary_review",
            "other_boundary_review",
        }
    )
    supported_seed = out["stage1_5_seed_policy"].isin({"priority_seed", "secondary_seed"})
    good_nonrisk = out["tier2_class"].isin(STRONG_GOOD)

    out["stage1_5_recommended_action"] = "hold_for_now"
    out.loc[out["tier2_class"].eq("T2_control_anchor"), "stage1_5_recommended_action"] = "keep_control_anchor"
    out.loc[good_nonrisk, "stage1_5_recommended_action"] = "keep_for_targeted_pool_review"
    out.loc[
        good_nonrisk & supported_seed & out["stage1_5_route_supported"],
        "stage1_5_recommended_action",
    ] = "priority_for_targeted_pool_review"
    out.loc[
        risk & out["stage1_5_structure_risk_band"].eq("severe_structure_risk"),
        "stage1_5_recommended_action",
    ] = "hold_severe_structure_risk"
    out.loc[
        reviewable_boundary & supported_seed & out["stage1_5_route_supported"],
        "stage1_5_recommended_action",
    ] = "manual_boundary_review"
    out.loc[
        reviewable_boundary
        & out["target"].eq("Ab_sdAb")
        & out["stage1_5_foldx_neutral_retained_proxy"]
        & out["stage1_5_seed_policy"].ne("deprioritize_ae108h_seed"),
        "stage1_5_recommended_action",
    ] = "manual_boundary_review"
    out.loc[
        out["stage1_5_seed_policy"].eq("deprioritize_ae108h_seed")
        & out["tier2_class"].ne("T2_strong_candidate"),
        "stage1_5_recommended_action",
    ] = "hold_ae108h_broad_expansion"
    out.loc[
        risk & out["stage1_5_structure_risk_band"].eq("severe_structure_risk"),
        "stage1_5_recommended_action",
    ] = "hold_severe_structure_risk"
    out.loc[out["tier2_class"].eq("T2_control_anchor"), "stage1_5_recommended_action"] = "keep_control_anchor"
    return out


def calibration_summary(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["target", "stage1_5_structure_risk_band", "stage1_5_recommended_action"]
    out = df.groupby(cols, dropna=False).size().reset_index(name="count")
    target_totals = df.groupby("target", dropna=False).size().rename("target_total")
    out = out.merge(target_totals, on="target", how="left")
    out["fraction_of_target"] = (out["count"] / out["target_total"]).round(6)
    return out.sort_values(["target", "stage1_5_structure_risk_band", "count"], ascending=[True, True, False])


def boundary_review_candidates(df: pd.DataFrame, per_target: int = 120) -> pd.DataFrame:
    sub = df[df["stage1_5_recommended_action"].eq("manual_boundary_review")].copy()
    if sub.empty:
        return sub
    sub = sub.sort_values(
        [
            "target",
            "stage1_5_review_support_score",
            "stage1_5_foldx_neutral_retained_proxy",
            "local_structure_validity_t2_score",
            "variant_id",
        ],
        ascending=[True, False, False, False, True],
    )
    picked = sub.groupby(["target", "his_seed_set", "primary_generation_route"], dropna=False).head(8)
    picked = picked.groupby("target", dropna=False).head(per_target)
    cols = [
        "target",
        "variant_id",
        "tier2_class",
        "stage1_5_structure_risk_band",
        "stage1_5_recommended_action",
        "stage1_5_seed_policy",
        "stage1_5_route_supported",
        "stage1_5_mutation_order_policy",
        "stage1_5_review_support_score",
        "primary_generation_route",
        "his_seed_set",
        "mutation_count",
        "His_count",
        "rescue_signature",
        "neutral_retention_t2_score",
        "acidic_release_mechanism_t2_score",
        "global_weakening_risk_t2_score",
        "local_structure_validity_t2_score",
        "local_clash_count",
        "rosetta_delta_score",
        "structure_risk_driver",
        "stage1_5_foldx_neutral_retained_proxy",
        "foldx_tool_status",
        "foldx_ddG_vs_parent_pH7_4",
        "pdb_path",
        "mutation_list",
    ]
    return picked[[c for c in cols if c in picked.columns]]


def targeted_pool_review_candidates(df: pd.DataFrame, per_target: int = 200) -> pd.DataFrame:
    sub = df[df["stage1_5_recommended_action"].isin(
        {"priority_for_targeted_pool_review", "keep_for_targeted_pool_review"}
    )].copy()
    if sub.empty:
        return sub
    sub = sub.sort_values(
        [
            "target",
            "stage1_5_recommended_action",
            "stage1_5_review_support_score",
            "tier2_class",
            "variant_id",
        ],
        ascending=[True, True, False, True, True],
    )
    picked = sub.groupby(["target", "his_seed_set", "primary_generation_route"], dropna=False).head(10)
    picked = picked.groupby("target", dropna=False).head(per_target)
    cols = [
        "target",
        "variant_id",
        "tier2_class",
        "stage1_5_recommended_action",
        "stage1_5_seed_policy",
        "stage1_5_mutation_order_policy",
        "stage1_5_review_support_score",
        "primary_generation_route",
        "his_seed_set",
        "mutation_count",
        "His_count",
        "rescue_signature",
        "neutral_retention_t2_score",
        "acidic_release_mechanism_t2_score",
        "global_weakening_risk_t2_score",
        "local_structure_validity_t2_score",
        "local_clash_count",
        "rosetta_delta_score",
        "foldx_tool_status",
        "foldx_ddG_vs_parent_pH7_4",
        "mutation_list",
    ]
    return picked[[c for c in cols if c in picked.columns]]


def foldx_concordance(df: pd.DataFrame) -> pd.DataFrame:
    sub = df[bool_series(df, "tier2_foldx_subset_selected")].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["foldx_success"] = safe_col(sub, "foldx_tool_status").eq("success")
    sub["foldx_neutral_retained_proxy"] = numeric(sub, "foldx_ddG_vs_parent_pH7_4").le(2.0)
    sub["is_strong_good"] = sub["tier2_class"].isin(STRONG_GOOD)
    rows = []
    for keys, group in sub.groupby(["target", "tier2_class"], dropna=False):
        target, tier2_class = keys
        ddg = numeric(group, "foldx_ddG_vs_parent_pH7_4")
        rows.append(
            {
                "target": target,
                "tier2_class": tier2_class,
                "foldx_selected_count": len(group),
                "foldx_success_count": int(group["foldx_success"].sum()),
                "foldx_success_rate": round(float(group["foldx_success"].mean()), 6),
                "foldx_ddG_vs_parent_pH7_4_median": round(float(ddg.median()), 6) if ddg.notna().any() else pd.NA,
                "foldx_ddG_vs_parent_pH7_4_p75": round(float(ddg.quantile(0.75)), 6) if ddg.notna().any() else pd.NA,
                "foldx_neutral_retained_proxy_rate": round(float(group["foldx_neutral_retained_proxy"].mean()), 6),
            }
        )
    return pd.DataFrame(rows)


def route_seed_risk_report(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["target", "primary_generation_route", "his_seed_set"]
    rows = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        target, route, seed = keys
        total = len(sub)
        strong_good = int(sub["tier2_class"].isin(STRONG_GOOD).sum())
        structure_risk = int(sub["tier2_class"].eq("T2_structure_risk").sum())
        rows.append(
            {
                "target": target,
                "primary_generation_route": route,
                "his_seed_set": seed,
                "total_count": total,
                "strong_good_count": strong_good,
                "strong_good_rate": round(strong_good / total, 6) if total else 0,
                "structure_risk_count": structure_risk,
                "structure_risk_rate": round(structure_risk / total, 6) if total else 0,
                "median_local_structure_validity": round(float(numeric(sub, "local_structure_validity_t2_score").median()), 6),
                "median_rosetta_delta": round(float(numeric(sub, "rosetta_delta_score").median()), 6),
                "median_neutral_retention_t2": round(float(numeric(sub, "neutral_retention_t2_score").median()), 6),
                "median_acidic_release_t2": round(float(numeric(sub, "acidic_release_mechanism_t2_score").median()), 6),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(["target", "strong_good_rate", "total_count"], ascending=[True, False, False])


def structure_risk_driver_summary(df: pd.DataFrame) -> pd.DataFrame:
    risk = df[df["tier2_class"].eq("T2_structure_risk")].copy()
    if risk.empty:
        return pd.DataFrame(
            columns=[
                "target",
                "structure_risk_driver",
                "count",
                "fraction_of_structure_risk",
                "fraction_of_stage1",
            ]
        )
    stage_totals = df.groupby("target", dropna=False).size().rename("stage1_total")
    risk_totals = risk.groupby("target", dropna=False).size().rename("structure_risk_total")
    out = risk.groupby(["target", "structure_risk_driver"], dropna=False).size().reset_index(name="count")
    out = out.merge(stage_totals, on="target", how="left")
    out = out.merge(risk_totals, on="target", how="left")
    out["fraction_of_structure_risk"] = (out["count"] / out["structure_risk_total"]).round(6)
    out["fraction_of_stage1"] = (out["count"] / out["stage1_total"]).round(6)
    return out.sort_values(["target", "count"], ascending=[True, False])


def manual_spot_check_panel(df: pd.DataFrame) -> pd.DataFrame:
    specs = {
        "Ab_1E62": {
            "T2_strong_candidate": 20,
            "T2_good_candidate": 20,
            "T2_structure_risk": 20,
            "T2_release_possible_but_neutral_risky": 10,
        },
        "Ab_sdAb": {
            "T2_strong_candidate": 20,
            "T2_good_candidate": 20,
            "T2_structure_risk": 30,
            "T2_release_possible_but_neutral_risky": 10,
            "T2_control_anchor": 10,
        },
    }
    rank_cols = [
        "target",
        "tier2_class",
        "tier2_stage1_bucket",
        "tier1_review_class",
        "primary_generation_route",
        "his_seed_set",
        "mutation_count",
        "His_count",
        "rescue_signature",
        "near_duplicate_cluster_id",
        "neutral_retention_t2_score",
        "acidic_release_mechanism_t2_score",
        "global_weakening_risk_t2_score",
        "local_structure_validity_t2_score",
        "local_clash_count",
        "rosetta_delta_score",
        "his_pka_support_t2_score",
        "foldx_tool_status",
        "foldx_ddG_vs_parent_pH7_4",
        "structure_risk_driver",
        "pdb_path",
    ]
    rows = []
    sort_df = df.copy()
    sort_df["manual_rank_score"] = (
        numeric(sort_df, "neutral_retention_t2_score").fillna(0)
        + numeric(sort_df, "acidic_release_mechanism_t2_score").fillna(0)
        + numeric(sort_df, "local_structure_validity_t2_score").fillna(0)
        - numeric(sort_df, "global_weakening_risk_t2_score").fillna(1)
    )
    for target, cls_spec in specs.items():
        for cls, n in cls_spec.items():
            sub = sort_df[sort_df["target"].eq(target) & sort_df["tier2_class"].eq(cls)].copy()
            if sub.empty:
                continue
            sub = sub.sort_values(
                ["manual_rank_score", "near_duplicate_cluster_id", "variant_id"],
                ascending=[False, True, True],
            )
            picked = sub.groupby("near_duplicate_cluster_id", dropna=False).head(2).head(n).copy()
            picked["spot_check_reason"] = f"{target}:{cls}"
            rows.append(picked[["variant_id", "spot_check_reason"] + [c for c in rank_cols if c in picked.columns]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)


def render_top_table(df: pd.DataFrame, cols: Iterable[str], n: int = 12) -> str:
    cols = [c for c in cols if c in df.columns]
    if not cols or df.empty:
        return "_No rows._"
    view = df[cols].head(n).copy()
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, low_memory=False)
    df = add_structure_risk_drivers(df)
    df = add_calibration_labels(df)

    outputs: dict[str, pd.DataFrame] = {}
    outputs["tier2_stage1_calibrated_results.csv"] = df
    outputs["tier2_class_by_target.csv"] = class_table(df, ["target"])
    outputs["tier2_class_by_tier1_class.csv"] = class_table(df, ["target", "tier1_review_class"])
    outputs["tier2_class_by_generation_route.csv"] = class_table(df, ["target", "primary_generation_route"])
    outputs["tier2_class_by_selection_bucket.csv"] = class_table(df, ["target", "tier2_stage1_bucket"])
    outputs["tier2_class_by_his_count.csv"] = class_table(df, ["target", "His_count"])
    outputs["tier2_class_by_his_seed_set.csv"] = class_table(df, ["target", "his_seed_set"])
    outputs["tier2_class_by_mutation_order.csv"] = class_table(df, ["target", "mutation_count"])
    outputs["tier2_class_by_rescue_objective.csv"] = class_table(df, ["target", "rescue_signature"])
    outputs["tier2_class_by_near_duplicate_cluster.csv"] = class_table(df, ["target", "near_duplicate_cluster_id"])
    outputs["tier2_feature_summary_by_class.csv"] = summarize_features(df, ["target", "tier2_class"])
    outputs["foldx_subset_concordance_by_class.csv"] = foldx_concordance(df)
    outputs["tier2_class_transition_matrix.csv"] = class_table(df, ["target", "tier1_review_class"])
    outputs["tier2_route_seed_risk_report.csv"] = route_seed_risk_report(df)
    outputs["tier2_structure_risk_driver_summary.csv"] = structure_risk_driver_summary(df)
    outputs["tier2_stage1_5_calibration_summary.csv"] = calibration_summary(df)
    outputs["tier2_boundary_review_candidates.csv"] = boundary_review_candidates(df)
    outputs["tier2_targeted_pool_review_candidates.csv"] = targeted_pool_review_candidates(df)
    outputs["manual_spot_check_panel.csv"] = manual_spot_check_panel(df)

    for name, table in outputs.items():
        write_csv(table, out_dir / name)

    target_summary = outputs["tier2_class_by_target.csv"]
    route_summary = (
        outputs["tier2_class_by_generation_route.csv"]
        .drop_duplicates(["target", "primary_generation_route"])
        .sort_values(["target", "strong_good_rate"], ascending=[True, False])
    )
    seed_report = outputs["tier2_route_seed_risk_report.csv"]
    risk_driver = outputs["tier2_structure_risk_driver_summary.csv"]
    calibration = outputs["tier2_stage1_5_calibration_summary.csv"]
    boundary_review = outputs["tier2_boundary_review_candidates.csv"]
    targeted_review = outputs["tier2_targeted_pool_review_candidates.csv"]
    near_threshold = df[df["structure_risk_near_threshold_only"]].groupby("target").size().reset_index(name="near_threshold_structure_risk_count")
    total_risk = df[df["tier2_class"].eq("T2_structure_risk")].groupby("target").size().reset_index(name="structure_risk_count")
    threshold = total_risk.merge(near_threshold, on="target", how="left").fillna(0)
    threshold["near_threshold_fraction_of_structure_risk"] = (
        threshold["near_threshold_structure_risk_count"] / threshold["structure_risk_count"]
    ).round(6)
    write_csv(threshold, out_dir / "tier2_threshold_sensitivity_summary.csv")

    report = [
        "# Tier2 Stage-1 Diagnostic Summary",
        "",
        "Status: `stage1_5_diagnostics_complete`",
        "",
        "This audit uses the existing 7,000-row Full Stage-1 output only. It does not run Tier2-C expansion, Tier2-heavy, AF3, SimpleFold, MD, or final selection.",
        "",
        "## Class Summary",
        render_top_table(target_summary, ["target", "tier2_class", "count", "fraction_within_group", "strong_good_rate", "structure_risk_rate"], 20),
        "",
        "## Route-Level Signal",
        "Rows below are one row per target / route, sorted by strong+good rate.",
        "",
        render_top_table(route_summary, ["target", "primary_generation_route", "total_count", "strong_good_rate", "structure_risk_rate"], 20),
        "",
        "## Structure-Risk Drivers",
        render_top_table(risk_driver, ["target", "structure_risk_driver", "count", "fraction_of_structure_risk", "fraction_of_stage1"], 30),
        "",
        "## Threshold Sensitivity",
        "Near-threshold rows are currently `T2_structure_risk` but would not trigger under a coarse relaxed check: local validity <0.30, clashes >=12, or Rosetta delta >800.",
        "",
        render_top_table(threshold, ["target", "structure_risk_count", "near_threshold_structure_risk_count", "near_threshold_fraction_of_structure_risk"], 10),
        "",
        "## Stage-1.5 Calibration Labels",
        "These labels do not change the original Tier2 class. They split `T2_structure_risk` into severe risk versus reviewable boundary cases and tag rows that should remain held.",
        "",
        render_top_table(calibration, ["target", "stage1_5_structure_risk_band", "stage1_5_recommended_action", "count", "fraction_of_target"], 40),
        "",
        "### Boundary Review Queue",
        f"Rows marked `manual_boundary_review`: {len(boundary_review)}. This is a capped review queue, not an expansion list.",
        "",
        "### Targeted Pool Review Queue",
        f"Rows marked for targeted pool review in the capped queue: {len(targeted_review)}. Tier2-C expansion remains HOLD until these labels are manually accepted.",
        "",
        "## Top Route/Seed Combinations By Strong+Good Rate",
        "Minimum `total_count >= 10`; shown separately by target to avoid one target hiding the other.",
        "",
        "### Ab_1E62",
        render_top_table(seed_report[seed_report["target"].eq("Ab_1E62") & seed_report["total_count"].ge(10)], ["target", "primary_generation_route", "his_seed_set", "total_count", "strong_good_rate", "structure_risk_rate", "median_local_structure_validity", "median_rosetta_delta"], 20),
        "",
        "### Ab_sdAb",
        render_top_table(seed_report[seed_report["target"].eq("Ab_sdAb") & seed_report["total_count"].ge(10)], ["target", "primary_generation_route", "his_seed_set", "total_count", "strong_good_rate", "structure_risk_rate", "median_local_structure_validity", "median_rosetta_delta"], 20),
        "",
        "## Deliverables",
        "",
        "- `tier2_class_by_target.csv`",
        "- `tier2_class_by_tier1_class.csv`",
        "- `tier2_class_by_generation_route.csv`",
        "- `tier2_class_by_selection_bucket.csv`",
        "- `tier2_class_by_his_count.csv`",
        "- `tier2_class_by_his_seed_set.csv`",
        "- `tier2_class_by_mutation_order.csv`",
        "- `tier2_class_by_rescue_objective.csv`",
        "- `tier2_class_by_near_duplicate_cluster.csv`",
        "- `tier2_feature_summary_by_class.csv`",
        "- `foldx_subset_concordance_by_class.csv`",
        "- `tier2_class_transition_matrix.csv`",
        "- `tier2_route_seed_risk_report.csv`",
        "- `tier2_structure_risk_driver_summary.csv`",
        "- `tier2_stage1_calibrated_results.csv`",
        "- `tier2_stage1_5_calibration_summary.csv`",
        "- `tier2_threshold_sensitivity_summary.csv`",
        "- `tier2_boundary_review_candidates.csv`",
        "- `tier2_targeted_pool_review_candidates.csv`",
        "- `manual_spot_check_panel.csv`",
    ]
    (out_dir / "tier2_stage1_diagnostic_summary.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    notes = [
        "# Tier2 Threshold Calibration Notes",
        "",
        "No thresholds were changed in this run.",
        "",
        "Current `T2_structure_risk` decomposition follows the implemented Full Stage-1 collector:",
        "",
        "```text",
        "pyrosetta_tool_status != success",
        "OR local_structure_validity_t2_score < 0.35",
        "OR local_clash_count >= 10",
        "OR rosetta_delta_score > 600",
        "```",
        "",
        "Immediate interpretation rules:",
        "",
        "- If structure-risk rows are dominated by severe Rosetta delta / clashes, treat them as real Tier2 risk and do not expand that source broadly.",
        "- If a large fraction is near-threshold only, inspect the spot-check panel before relaxing any threshold.",
        "- FoldX subset is a supporting check only; it was not run for all 7,000 rows.",
        "- sdAb should retain stricter expansion criteria unless route/seed diagnostics show a clear high-yield region.",
        "- `stage1_5_*` labels are decision-support labels only. They do not unlock Tier2-C, Tier2-heavy, or final selection.",
        "- `severe_structure_risk` should stay held. Boundary labels are reviewable only when seed, route, neutral-retention, and diversity support are coherent.",
        "- sdAb `AE108H` broad expansion remains held unless later evidence overrides this rule.",
        "",
        "Suggested decision boundary after reviewing this report:",
        "",
        "```text",
        "Tier2-C expansion: HOLD until route/seed risk is interpreted.",
        "Tier2-heavy pilot compute: HOLD until manual spot-check and threshold interpretation are complete.",
        "Final 10K selection: NO-GO.",
        "```",
    ]
    (out_dir / "tier2_threshold_calibration_notes.md").write_text("\n".join(notes) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
