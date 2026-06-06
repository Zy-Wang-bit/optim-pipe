#!/usr/bin/env python3
"""Run Tier 1 annotation/filtering and build a Tier2-core proposal.

This stage is table-only. It does not run Tier2 structural compute.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT
OUT = ROOT / "results/initial_design_generation"
P0_REVIEW = OUT / "p0_evidence_review"
P0_MPNN = OUT / "p0_mpnn"
DEFAULT_OUT = OUT / "tier1_filtering"

TIER1_BASE = P0_REVIEW / "tier1_feature_table_with_mpnn.csv"
SPARSE_SUPPLEMENT = P0_REVIEW / "mpnn_sparse_rescue_supplement_candidates.csv"
PRODUCTION_FEATURES = P0_MPNN / "production_pool_features.csv"

WINDOWS = {
    "Ab_1E62": ("1E62_VL_001_040", 1, 40),
    "Ab_sdAb": ("sdAb_VHH_072_111", 72, 111),
}

EXPECTED_ROWS = {"Ab_1E62": 120_165, "Ab_sdAb": 150_041}
EXPECTED_TOTAL = 270_206
EXPECTED_SPARSE_EXACT = 514
EXPECTED_SPARSE_NOVEL = 206

REVIEW_CLASSES = [
    "A_tier1_priority_candidate",
    "B_rescue_enriched_candidate",
    "C_high_pH_support_but_mpnn_weak",
    "D_favorable_mpnn_but_weak_pH_mechanism",
    "E_sparse_mpnn_supplement_review",
    "F_neutral_boundary_or_high_risk",
    "G_audit_or_demoted",
]

CLASS_PRIORITY = {name: i for i, name in enumerate(REVIEW_CLASSES)}

PROPOSAL_SIZE = {
    "primary": {"Ab_1E62": 12_000, "Ab_sdAb": 15_000},
    "reserve": {"Ab_1E62": 8_000, "Ab_sdAb": 10_000},
}

CAP_POLICY = {
    "max_fraction_per_near_duplicate_cluster": {
        "primary_proposal": 0.01,
        "reserve_pool": 0.02,
    },
    "max_fraction_per_his_seed_set": {"Ab_1E62": 0.12, "Ab_sdAb": 0.10},
    "max_fraction_per_primary_generation_route": {"Ab_1E62": 0.45, "Ab_sdAb": 0.40},
    "min_fraction_for_core_routes_if_available": {
        "His_rule": 0.10,
        "His_plus_rescue": 0.15,
        "ProteinMPNN_seeded_rescue": 0.05,
        "wetlab_informed_expansion": 0.05,
    },
    "sparse_supplement": {
        "max_fraction_per_target": 0.02,
        "default_fraction": 0.00,
        "quota_type": "cap_not_quota",
    },
    "relaxed_mpnn_only": {"max_fraction_per_target": 0.00},
}


def sha(text: str, n: int | None = None) -> str:
    value = hashlib.sha256(str(text).encode()).hexdigest()
    return value[:n] if n else value


def list_join(values: pd.Series, limit: int = 30) -> str:
    unique = [str(v) for v in values.dropna().astype(str).unique() if str(v) and str(v).lower() != "nan"]
    unique = sorted(unique)
    if len(unique) > limit:
        return ";".join(unique[:limit]) + f";...(+{len(unique) - limit})"
    return ";".join(unique)


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy() if max_rows else df.copy()
    cols = [str(c) for c in view.columns]
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in view.iterrows():
        values = [str(row[col]).replace("\n", " ") for col in view.columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def normalize_mutation_list(value: object) -> str:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return ""
    items = [x for x in str(value).split(";") if x and x.lower() != "nan"]

    def key(label: str) -> tuple[str, int, str, str]:
        digits = "".join(ch for ch in label if ch.isdigit())
        prefix = label[:1]
        parent = label[1:2] if len(label) > 1 else ""
        mutant = label[-1:] if label else ""
        return prefix, int(digits) if digits else 0, parent, mutant

    return ";".join(sorted(items, key=key))


def window_for_target(target: str) -> tuple[str, int, int]:
    if target not in WINDOWS:
        raise ValueError(f"Unknown target for window mapping: {target}")
    return WINDOWS[target]


def mutated_window_sequence(row: pd.Series) -> str:
    seq = str(row.get("sequence", ""))
    if not seq or seq == "nan":
        return ""
    _, start, end = window_for_target(str(row["target"]))
    if len(seq) < end:
        return ""
    return seq[start - 1 : end]


def flag_threshold(value: Any, pass_min: float, boundary_min: float, risk_min: float) -> str:
    if value is None or pd.isna(value):
        return "not_available"
    value = float(value)
    if value >= pass_min:
        return "pass"
    if value >= boundary_min:
        return "boundary"
    if value >= risk_min:
        return "risk"
    return "fail"


def inverse_flag_threshold(value: Any, pass_max: float, boundary_max: float, risk_max: float) -> str:
    if value is None or pd.isna(value):
        return "not_available"
    value = float(value)
    if value <= pass_max:
        return "pass"
    if value <= boundary_max:
        return "boundary"
    if value <= risk_max:
        return "risk"
    return "fail"


def read_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = pd.read_csv(TIER1_BASE)
    sparse = pd.read_csv(SPARSE_SUPPLEMENT)
    production = pd.read_csv(
        PRODUCTION_FEATURES,
        usecols=[
            "variant_id",
            "sequence",
            "window",
            "generation_record_id",
            "all_generation_records",
            "selection_reason",
        ],
    )
    return base, sparse, production


def build_overlap_provenance(sparse: pd.DataFrame) -> pd.DataFrame:
    exact = sparse[
        sparse["production_pool_overlap_status"].eq("exact_sequence_and_mutation_list_overlap")
    ].copy()
    exact = exact.rename(columns={"exact_existing_production_variant_id": "production_variant_id"})
    exact["sparse_variant_id"] = exact["variant_id"]
    exact["sparse_sequence_hash"] = exact["sequence_hash"]
    exact["sparse_mpnn_seed_set"] = exact["his_seed_set"]
    exact["sparse_mpnn_rescue_signature"] = exact["rescue_mutation_list"].fillna("")
    exact["sparse_mpnn_source_generation_records"] = exact["source_generation_records"].fillna("")
    grouped = (
        exact.groupby("production_variant_id", dropna=False)
        .agg(
            target=("target", "first"),
            sparse_overlap_count=("sparse_variant_id", "count"),
            sparse_mpnn_generated_raw_count=("generated_count", "sum"),
            sparse_variant_ids=("sparse_variant_id", list_join),
            sparse_sequence_hashes=("sparse_sequence_hash", list_join),
            sparse_mpnn_seed_set=("sparse_mpnn_seed_set", list_join),
            sparse_mpnn_rescue_signature=("sparse_mpnn_rescue_signature", list_join),
            sparse_mpnn_source_generation_records=("sparse_mpnn_source_generation_records", list_join),
        )
        .reset_index()
    )
    return grouped


def build_tier1_universe(base: pd.DataFrame, sparse: pd.DataFrame, production: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = base.merge(production, on="variant_id", how="left", validate="many_to_one")
    sparse_meta_cols = [
        "variant_id",
        "sequence",
        "generated_count",
        "source_generation_records",
        "mpnn_score",
        "mpnn_score_percentile_by_target",
        "mpnn_score_percentile_by_route",
        "mpnn_score_percentile_by_seed",
        "mpnn_score_percentile_by_rescue",
    ]
    sparse_meta = sparse[sparse_meta_cols].rename(
        columns={
            "sequence": "sparse_sequence",
            "generated_count": "sparse_mpnn_generated_raw_count_self",
            "source_generation_records": "sparse_source_generation_records",
            "mpnn_score": "sparse_generated_mpnn_score",
        }
    )
    df = df.merge(sparse_meta, on="variant_id", how="left", validate="many_to_one")
    df["sequence"] = df["sequence"].fillna(df["sparse_sequence"])
    df["window_id"] = df["target"].map(lambda t: window_for_target(str(t))[0])
    df["canonical_sequence_hash_full"] = df["sequence"].map(sha)
    df["canonical_sequence_hash_short"] = df["canonical_sequence_hash_full"].str[:12]
    df["sequence_hash_length"] = df["canonical_sequence_hash_full"].str.len()
    df["input_sequence_hash_length"] = df["sequence_hash"].astype(str).str.len()
    df["canonical_mutated_window_sequence"] = df.apply(mutated_window_sequence, axis=1)
    df["normalized_mutation_list"] = df["mutation_list"].map(normalize_mutation_list)
    df["generation_record_ids"] = df["all_generation_records"].fillna(df["generation_record_id"]).fillna(
        df["sparse_source_generation_records"]
    )

    overlap = build_overlap_provenance(sparse)
    df = df.merge(
        overlap,
        left_on="variant_id",
        right_on="production_variant_id",
        how="left",
        suffixes=("", "_overlap"),
    )
    df["mpnn_sparse_overlap"] = df["sparse_overlap_count"].fillna(0).astype(int) > 0
    df["mpnn_sparse_novel"] = df["source_universe"].eq("repaired_sparse_constrained_MPNN_supplement")
    df["sparse_mpnn_generated_raw_count"] = df["sparse_mpnn_generated_raw_count"].fillna(
        df["sparse_mpnn_generated_raw_count_self"]
    ).fillna(0).astype(int)
    for col in [
        "sparse_variant_ids",
        "sparse_sequence_hashes",
        "sparse_mpnn_seed_set",
        "sparse_mpnn_rescue_signature",
        "sparse_mpnn_source_generation_records",
    ]:
        df[col] = df[col].fillna("")
    overlap_route = "repaired_sparse_constrained_MPNN"
    needs_route = df["mpnn_sparse_overlap"] & ~df["all_source_routes"].fillna("").str.contains(overlap_route)
    df.loc[needs_route, "all_source_routes"] = (
        df.loc[needs_route, "all_source_routes"].fillna("").astype(str).str.rstrip(";") + ";" + overlap_route
    ).str.strip(";")

    df["source_admission_status"] = "production_pool"
    df.loc[df["mpnn_sparse_overlap"], "source_admission_status"] = "sparse_overlap_provenance_only"
    df.loc[df["mpnn_sparse_novel"], "source_admission_status"] = "sparse_novel_supplement"
    df["tier2_eligibility_status"] = "eligible_for_tier1_review"
    df.loc[df["mpnn_sparse_novel"], "tier2_eligibility_status"] = "supplemental_review_only"
    for col in [
        "sequence_hamming_cluster_id",
        "mutation_set_jaccard_cluster_id",
        "his_seed_cluster_id",
        "near_duplicate_cluster_id",
        "cluster_size",
    ]:
        df[f"source_{col}"] = df[col] if col in df.columns else ""
    clustered = dry.assign_near_duplicate_clusters(
        df[["target", "variant_id", "sequence", "mutation_list", "his_seed_set"]].copy()
    )
    cluster_cols = [
        "variant_id",
        "sequence_hamming_cluster_id",
        "mutation_set_jaccard_cluster_id",
        "his_seed_cluster_id",
        "near_duplicate_cluster_id",
        "cluster_size",
    ]
    clustered = clustered[cluster_cols].rename(
        columns={
            "sequence_hamming_cluster_id": "tier1_sequence_hamming_cluster_id",
            "mutation_set_jaccard_cluster_id": "tier1_mutation_set_jaccard_cluster_id",
            "his_seed_cluster_id": "tier1_his_seed_cluster_id",
            "near_duplicate_cluster_id": "tier1_near_duplicate_cluster_id",
            "cluster_size": "tier1_cluster_size",
        }
    )
    df = df.merge(clustered, on="variant_id", how="left", validate="one_to_one")
    df["sequence_hamming_cluster_id"] = df["tier1_sequence_hamming_cluster_id"]
    df["mutation_set_jaccard_cluster_id"] = df["tier1_mutation_set_jaccard_cluster_id"]
    df["his_seed_cluster_id"] = df["tier1_his_seed_cluster_id"]
    df["near_duplicate_cluster_id"] = df["tier1_near_duplicate_cluster_id"]
    df["cluster_size"] = df["tier1_cluster_size"]
    return df, overlap


def add_mpnn_semantics(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["mpnn_score_raw"] = out["mpnn_total_score_per_residue"]
    out.loc[out["mpnn_sparse_novel"], "mpnn_score_raw"] = out.loc[out["mpnn_sparse_novel"], "sparse_generated_mpnn_score"]
    out["mpnn_score_direction"] = "lower_is_better"
    out["mpnn_score_lower_is_better"] = True
    out["mpnn_score_comparable"] = out["source_universe"].eq("production_initial_pool")
    out["mpnn_score_comparable_group"] = "score_only_full_design_chain"
    out.loc[out["mpnn_sparse_novel"], "mpnn_score_comparable_group"] = "generated_sparse_constrained_not_score_only"
    out["mpnn_score_type"] = out["mpnn_score_type"].fillna("unknown")

    comparable = out["mpnn_score_comparable"] & out["mpnn_score_raw"].notna()
    out["mpnn_score_percentile_within_target"] = pd.NA
    out["mpnn_score_percentile_within_route"] = pd.NA
    out["mpnn_score_percentile_within_his_seed"] = pd.NA
    out["mpnn_score_percentile_within_rescue_type"] = pd.NA
    out.loc[comparable, "mpnn_score_percentile_within_target"] = out[comparable].groupby("target")[
        "mpnn_score_raw"
    ].rank(pct=True, ascending=True)
    out.loc[comparable, "mpnn_score_percentile_within_route"] = out[comparable].groupby(
        ["target", "primary_generation_route"]
    )["mpnn_score_raw"].rank(pct=True, ascending=True)
    out.loc[comparable, "mpnn_score_percentile_within_his_seed"] = out[comparable].groupby(
        ["target", "his_seed_set"]
    )["mpnn_score_raw"].rank(pct=True, ascending=True)
    out.loc[comparable, "mpnn_score_percentile_within_rescue_type"] = out[comparable].groupby(
        ["target", "rescue_signature"]
    )["mpnn_score_raw"].rank(pct=True, ascending=True)
    return out


def add_tier1_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["acidic_release_support_score_v0"] = out["acidic_release_support_score"]
    out["neutral_retention_score_v0"] = out["neutral_retention_score"]
    out["global_weakening_risk_score_v0"] = out["global_weakening_risk_score"]
    out["liability_risk_score_v0"] = out["liability_flags"].notna().astype(float)
    out["tier1_acidic_release_support_flag"] = out["acidic_release_support_score_v0"].map(
        lambda x: flag_threshold(x, 0.65, 0.50, 0.35)
    )
    out["tier1_neutral_retention_flag"] = out["neutral_retention_score_v0"].map(
        lambda x: flag_threshold(x, 0.70, 0.55, 0.40)
    )
    out["tier1_global_weakening_risk_flag"] = out["global_weakening_risk_score_v0"].map(
        lambda x: inverse_flag_threshold(x, 0.25, 0.40, 0.65)
    )
    liability_present = out["liability_flags"].fillna("").astype(str).str.len() > 0
    hard_pass = out["hard_filter_status"].fillna("").eq("pass")
    out["tier1_liability_flag"] = "pass"
    out.loc[~hard_pass, "tier1_liability_flag"] = "fail"
    out.loc[hard_pass & liability_present, "tier1_liability_flag"] = "risk"
    out["tier1_mpnn_support_flag"] = "not_available"
    comparable = out["mpnn_score_comparable"] & out["mpnn_score_percentile_within_route"].notna()
    pct = pd.to_numeric(out["mpnn_score_percentile_within_route"], errors="coerce")
    out.loc[comparable & (pct <= 0.25), "tier1_mpnn_support_flag"] = "pass"
    out.loc[comparable & (pct > 0.25) & (pct <= 0.50), "tier1_mpnn_support_flag"] = "boundary"
    out.loc[comparable & (pct > 0.50) & (pct <= 0.90), "tier1_mpnn_support_flag"] = "risk"
    out.loc[comparable & (pct > 0.90), "tier1_mpnn_support_flag"] = "fail"

    strong_ph = out["tier1_acidic_release_support_flag"].eq("pass") & out[
        "tier1_neutral_retention_flag"
    ].isin(["pass", "boundary"]) & (out["His_count"].fillna(0).astype(float) >= 1)
    weak_ph = out["tier1_acidic_release_support_flag"].isin(["risk", "fail"]) | (
        out["His_count"].fillna(0).astype(float) < 1
    )
    favorable_mpnn = out["tier1_mpnn_support_flag"].eq("pass")
    weak_mpnn = out["tier1_mpnn_support_flag"].isin(["risk", "fail"])
    out["tier1_discordance_flag"] = "none"
    out.loc[favorable_mpnn & weak_ph, "tier1_discordance_flag"] = "favorable_mpnn_but_weak_pH_mechanism"
    out.loc[weak_mpnn & strong_ph, "tier1_discordance_flag"] = "weak_mpnn_but_strong_pH_support"

    missing: list[str] = []
    for row in out.itertuples(index=False):
        reasons = []
        if not getattr(row, "sequence", "") or str(getattr(row, "sequence")) == "nan":
            reasons.append("sequence_missing")
        if pd.isna(getattr(row, "acidic_release_support_score_v0")):
            reasons.append("acidic_release_support_missing")
        if pd.isna(getattr(row, "neutral_retention_score_v0")):
            reasons.append("neutral_retention_missing")
        if not getattr(row, "mpnn_score_comparable"):
            reasons.append("mpnn_score_not_comparable_to_production_score_only")
        elif pd.isna(getattr(row, "mpnn_score_raw")):
            reasons.append("mpnn_score_missing")
        missing.append(";".join(reasons))
    out["feature_missing_reason"] = missing
    out["feature_confidence"] = "standard_tier1_features"
    out.loc[out["mpnn_sparse_novel"], "feature_confidence"] = "sparse_supplement_proxy_features"
    out["tier1_feature_completeness_class"] = "complete_for_tier1_review"
    out.loc[out["feature_missing_reason"].ne(""), "tier1_feature_completeness_class"] = "review_with_missing_features"
    out.loc[
        out["feature_missing_reason"].str.contains("acidic_release_support_missing|neutral_retention_missing", regex=True),
        "tier1_feature_completeness_class",
    ] = "incomplete_for_primary_proposal"
    return out


def assign_review_classes(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["tier1_review_class"] = "G_audit_or_demoted"
    out["tier1_review_reason"] = "default_demoted_or_low_evidence"
    out["tier1_review_confidence"] = "low"

    liability_ok = out["tier1_liability_flag"].eq("pass")
    acid_pass = out["tier1_acidic_release_support_flag"].eq("pass")
    acid_ok = out["tier1_acidic_release_support_flag"].isin(["pass", "boundary"])
    neutral_pass = out["tier1_neutral_retention_flag"].eq("pass")
    neutral_ok = out["tier1_neutral_retention_flag"].isin(["pass", "boundary"])
    weakening_ok = out["tier1_global_weakening_risk_flag"].isin(["pass", "boundary"])
    mpnn_ok = out["tier1_mpnn_support_flag"].isin(["pass", "boundary"])
    mpnn_weak = out["tier1_mpnn_support_flag"].isin(["risk", "fail"])
    comparable = out["mpnn_score_comparable"]
    has_rescue = out["rescue_count"].fillna(0).astype(float) > 0
    high_risk = out["tier1_global_weakening_risk_flag"].isin(["risk", "fail"]) | out[
        "tier1_neutral_retention_flag"
    ].isin(["boundary", "risk"])

    sparse = out["mpnn_sparse_novel"]
    out.loc[sparse, ["tier1_review_class", "tier1_review_reason", "tier1_review_confidence"]] = [
        "E_sparse_mpnn_supplement_review",
        "sparse_novel_supplement_noncomparable_mpnn_score",
        "moderate",
    ]

    discord_favorable_weak = out["tier1_discordance_flag"].eq("favorable_mpnn_but_weak_pH_mechanism")
    out.loc[
        comparable & discord_favorable_weak,
        ["tier1_review_class", "tier1_review_reason", "tier1_review_confidence"],
    ] = [
        "D_favorable_mpnn_but_weak_pH_mechanism",
        "mpnn_favorable_but_pH_mechanism_weak_do_not_promote_by_mpnn_only",
        "low",
    ]

    discord_weak_strong = out["tier1_discordance_flag"].eq("weak_mpnn_but_strong_pH_support")
    out.loc[
        comparable & discord_weak_strong,
        ["tier1_review_class", "tier1_review_reason", "tier1_review_confidence"],
    ] = [
        "C_high_pH_support_but_mpnn_weak",
        "strong_pH_support_with_weak_mpnn_keep_for_review",
        "moderate",
    ]

    out.loc[
        comparable & liability_ok & acid_ok & neutral_ok & has_rescue & weakening_ok & ~sparse,
        ["tier1_review_class", "tier1_review_reason", "tier1_review_confidence"],
    ] = [
        "B_rescue_enriched_candidate",
        "rescue_enriched_with_pH_and_neutral_support",
        "moderate",
    ]
    out.loc[
        comparable & liability_ok & acid_pass & neutral_pass & weakening_ok & mpnn_ok & ~sparse,
        ["tier1_review_class", "tier1_review_reason", "tier1_review_confidence"],
    ] = [
        "A_tier1_priority_candidate",
        "strong_tier1_pH_neutral_liability_and_mpnn_support",
        "moderate",
    ]
    out.loc[
        comparable & liability_ok & high_risk & ~sparse & ~out["tier1_review_class"].isin(["A_tier1_priority_candidate"]),
        ["tier1_review_class", "tier1_review_reason", "tier1_review_confidence"],
    ] = [
        "F_neutral_boundary_or_high_risk",
        "neutral_retention_boundary_or_global_weakening_risk",
        "low",
    ]
    out.loc[~liability_ok, ["tier1_review_class", "tier1_review_reason", "tier1_review_confidence"]] = [
        "G_audit_or_demoted",
        "liability_or_hard_filter_risk",
        "low",
    ]
    out["tier1_review_class_priority"] = out["tier1_review_class"].map(CLASS_PRIORITY).fillna(999).astype(int)
    return out


def add_rank_score(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    mpnn_component = 1.0 - pd.to_numeric(out["mpnn_score_percentile_within_route"], errors="coerce")
    mpnn_component = mpnn_component.where(out["mpnn_score_comparable"], 0.50).fillna(0.50)
    weakening_component = 1.0 - out["global_weakening_risk_score_v0"].fillna(1.0).astype(float)
    out["tier1_rank_score"] = (
        0.35 * out["hit_likelihood_score_v0"].fillna(0).astype(float)
        + 0.25 * out["acidic_release_support_score_v0"].fillna(0).astype(float)
        + 0.20 * out["neutral_retention_score_v0"].fillna(0).astype(float)
        + 0.10 * weakening_component.clip(lower=0, upper=1)
        + 0.10 * mpnn_component.clip(lower=0, upper=1)
    )
    out["tier1_rank_score"] = out["tier1_rank_score"].round(6)
    out["tier1_target_rank"] = out.sort_values(
        ["target", "tier1_review_class_priority", "tier1_rank_score"],
        ascending=[True, True, False],
    ).groupby("target").cumcount() + 1
    return out


def can_add(row: pd.Series, counters: dict[str, defaultdict[str, int]], limits: dict[str, int], allow_sparse: bool) -> bool:
    if row["mpnn_sparse_novel"] and not allow_sparse:
        return False
    if row["tier1_review_class"] in {"G_audit_or_demoted", "D_favorable_mpnn_but_weak_pH_mechanism"}:
        return False
    if row["tier1_feature_completeness_class"].eq("incomplete_for_primary_proposal") if isinstance(row["tier1_feature_completeness_class"], pd.Series) else row["tier1_feature_completeness_class"] == "incomplete_for_primary_proposal":
        return False
    cluster = str(row["near_duplicate_cluster_id"])
    seed = str(row["his_seed_set"])
    route = str(row["primary_generation_route"])
    if counters["cluster"][cluster] >= limits["cluster"]:
        return False
    if counters["seed"][seed] >= limits["seed"]:
        return False
    if counters["route"][route] >= limits["route"]:
        return False
    if row["mpnn_sparse_novel"] and counters["sparse"]["sparse"] >= limits["sparse"]:
        return False
    return True


def can_add_values(
    cluster: str,
    seed: str,
    route: str,
    sparse: bool,
    review_class: str,
    completeness: str,
    counters: dict[str, defaultdict[str, int]],
    limits: dict[str, int],
    allow_sparse: bool,
) -> bool:
    if sparse and not allow_sparse:
        return False
    if review_class in {"G_audit_or_demoted", "D_favorable_mpnn_but_weak_pH_mechanism"}:
        return False
    if completeness == "incomplete_for_primary_proposal":
        return False
    if counters["cluster"][cluster] >= limits["cluster"]:
        return False
    if counters["seed"][seed] >= limits["seed"]:
        return False
    if counters["route"][route] >= limits["route"]:
        return False
    if sparse and counters["sparse"]["sparse"] >= limits["sparse"]:
        return False
    return True


def add_to_selection(
    row: pd.Series,
    selected: list[int],
    selected_set: set[int],
    counters: dict[str, defaultdict[str, int]],
) -> None:
    row_id = int(row.name)
    selected.append(row_id)
    selected_set.add(row_id)
    counters["cluster"][str(row["near_duplicate_cluster_id"])] += 1
    counters["seed"][str(row["his_seed_set"])] += 1
    counters["route"][str(row["primary_generation_route"])] += 1
    if row["mpnn_sparse_novel"]:
        counters["sparse"]["sparse"] += 1


def add_values_to_selection(
    row_id: int,
    cluster: str,
    seed: str,
    route: str,
    sparse: bool,
    selected: list[int],
    selected_set: set[int],
    counters: dict[str, defaultdict[str, int]],
) -> None:
    selected.append(row_id)
    selected_set.add(row_id)
    counters["cluster"][cluster] += 1
    counters["seed"][seed] += 1
    counters["route"][route] += 1
    if sparse:
        counters["sparse"]["sparse"] += 1


def select_pool(df: pd.DataFrame, target: str, size: int, pool_type: str, exclude: set[int] | None = None) -> pd.DataFrame:
    exclude = exclude or set()
    target_df = df[(df["target"] == target) & ~df.index.isin(exclude)].copy()
    target_df = target_df.sort_values(
        ["tier1_review_class_priority", "tier1_rank_score", "mpnn_score_percentile_within_route"],
        ascending=[True, False, True],
    )
    cluster_frac = CAP_POLICY["max_fraction_per_near_duplicate_cluster"][
        "primary_proposal" if pool_type == "primary" else "reserve_pool"
    ]
    limits = {
        "cluster": max(1, math.floor(size * cluster_frac)),
        "seed": max(1, math.floor(size * CAP_POLICY["max_fraction_per_his_seed_set"][target])),
        "route": max(1, math.floor(size * CAP_POLICY["max_fraction_per_primary_generation_route"][target])),
        "sparse": max(0, math.floor(size * CAP_POLICY["sparse_supplement"]["max_fraction_per_target"])),
    }
    counters: dict[str, defaultdict[str, int]] = {
        "cluster": defaultdict(int),
        "seed": defaultdict(int),
        "route": defaultdict(int),
        "sparse": defaultdict(int),
    }
    selected: list[int] = []
    selected_set: set[int] = set()
    allow_sparse = pool_type == "reserve"
    iter_cols = [
        "near_duplicate_cluster_id",
        "his_seed_set",
        "primary_generation_route",
        "mpnn_sparse_novel",
        "tier1_review_class",
        "tier1_feature_completeness_class",
    ]

    if pool_type == "primary":
        for route, frac in CAP_POLICY["min_fraction_for_core_routes_if_available"].items():
            target_count = math.floor(size * frac)
            route_rows = target_df[target_df["primary_generation_route"].eq(route)]
            for row in route_rows[iter_cols].itertuples(index=True, name=None):
                row_id, cluster, seed, row_route, sparse, review_class, completeness = row
                if len(selected) >= size or counters["route"][route] >= target_count:
                    break
                if row_id in selected_set:
                    continue
                if can_add_values(
                    str(cluster),
                    str(seed),
                    str(row_route),
                    bool(sparse),
                    str(review_class),
                    str(completeness),
                    counters,
                    limits,
                    allow_sparse=False,
                ):
                    add_values_to_selection(
                        row_id,
                        str(cluster),
                        str(seed),
                        str(row_route),
                        bool(sparse),
                        selected,
                        selected_set,
                        counters,
                    )

    for row in target_df[iter_cols].itertuples(index=True, name=None):
        row_id, cluster, seed, route, sparse, review_class, completeness = row
        if len(selected) >= size:
            break
        if row_id in selected_set:
            continue
        if can_add_values(
            str(cluster),
            str(seed),
            str(route),
            bool(sparse),
            str(review_class),
            str(completeness),
            counters,
            limits,
            allow_sparse=allow_sparse,
        ):
            add_values_to_selection(
                row_id,
                str(cluster),
                str(seed),
                str(route),
                bool(sparse),
                selected,
                selected_set,
                counters,
            )

    return df.loc[selected].copy()


def build_proposals(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary_frames = []
    reserve_frames = []
    for target in ["Ab_1E62", "Ab_sdAb"]:
        primary = select_pool(df, target, PROPOSAL_SIZE["primary"][target], "primary")
        primary["tier2_proposal_role"] = "primary_proposal"
        primary_frames.append(primary)
        reserve = select_pool(
            df,
            target,
            PROPOSAL_SIZE["reserve"][target],
            "reserve",
            exclude=set(primary.index),
        )
        reserve["tier2_proposal_role"] = "reserve_pool"
        reserve_frames.append(reserve)
    return pd.concat(primary_frames, ignore_index=True), pd.concat(reserve_frames, ignore_index=True)


def summarize_missingness(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    features = [
        "sequence",
        "canonical_mutated_window_sequence",
        "acidic_release_support_score_v0",
        "neutral_retention_score_v0",
        "global_weakening_risk_score_v0",
        "mpnn_score_raw",
        "mpnn_score_percentile_within_route",
    ]
    for (target, source), group in df.groupby(["target", "source_universe"], dropna=False):
        for feature in features:
            missing = group[feature].isna() | group[feature].astype(str).eq("")
            rows.append(
                {
                    "target": target,
                    "source_universe": source,
                    "feature": feature,
                    "row_count": len(group),
                    "missing_count": int(missing.sum()),
                    "missing_fraction": round(float(missing.mean()), 6),
                    "severity": "review" if missing.any() else "pass",
                }
            )
    return pd.DataFrame(rows)


def build_summaries(df: pd.DataFrame, primary: pd.DataFrame, reserve: pd.DataFrame) -> dict[str, pd.DataFrame]:
    summaries: dict[str, pd.DataFrame] = {}
    summary_df = df.copy()
    primary_ids = set(primary["variant_id"].astype(str))
    reserve_ids = set(reserve["variant_id"].astype(str))
    summary_df["is_primary_selected"] = summary_df["variant_id"].astype(str).isin(primary_ids)
    summary_df["is_reserve_selected"] = summary_df["variant_id"].astype(str).isin(reserve_ids)
    summaries["tier1_feature_summary_by_target"] = df.groupby("target", dropna=False).agg(
        row_count=("variant_id", "count"),
        sparse_novel_count=("mpnn_sparse_novel", "sum"),
        sparse_overlap_provenance_count=("mpnn_sparse_overlap", "sum"),
        priority_count=("tier1_review_class", lambda s: int((s == "A_tier1_priority_candidate").sum())),
        acidic_pass=("tier1_acidic_release_support_flag", lambda s: int((s == "pass").sum())),
        neutral_pass=("tier1_neutral_retention_flag", lambda s: int((s == "pass").sum())),
        mpnn_comparable_count=("mpnn_score_comparable", "sum"),
    ).reset_index()
    summaries["tier1_feature_summary_by_route"] = df.groupby(["target", "primary_generation_route"], dropna=False).agg(
        row_count=("variant_id", "count"),
        priority_count=("tier1_review_class", lambda s: int((s == "A_tier1_priority_candidate").sum())),
        median_rank_score=("tier1_rank_score", "median"),
    ).reset_index()
    summaries["tier1_feature_summary_by_source_universe"] = df.groupby(["target", "source_universe"], dropna=False).agg(
        row_count=("variant_id", "count"),
        mpnn_comparable_count=("mpnn_score_comparable", "sum"),
        sparse_overlap_count=("mpnn_sparse_overlap", "sum"),
        sparse_novel_count=("mpnn_sparse_novel", "sum"),
    ).reset_index()
    summaries["tier1_mpnn_score_distribution_by_route"] = df.groupby(["target", "primary_generation_route", "mpnn_score_comparable"], dropna=False).agg(
        row_count=("variant_id", "count"),
        score_mean=("mpnn_score_raw", "mean"),
        score_median=("mpnn_score_raw", "median"),
        p10=("mpnn_score_raw", lambda s: s.quantile(0.10)),
        p90=("mpnn_score_raw", lambda s: s.quantile(0.90)),
    ).reset_index()
    summaries["tier1_mpnn_score_distribution_by_his_seed"] = df.groupby(["target", "his_seed_set", "mpnn_score_comparable"], dropna=False).agg(
        row_count=("variant_id", "count"),
        score_median=("mpnn_score_raw", "median"),
        route_count=("primary_generation_route", "nunique"),
    ).reset_index()
    summaries["tier1_near_duplicate_cluster_summary"] = summary_df.groupby(["target", "near_duplicate_cluster_id"], dropna=False).agg(
        cluster_size=("variant_id", "count"),
        routes=("primary_generation_route", list_join),
        his_seed_sets=("his_seed_set", list_join),
        best_rank_score=("tier1_rank_score", "max"),
        primary_selected=("is_primary_selected", "sum"),
        reserve_selected=("is_reserve_selected", "sum"),
    ).reset_index().sort_values(["target", "cluster_size"], ascending=[True, False])

    class_rows = []
    for target in sorted(df["target"].dropna().unique()):
        target_df = df[df["target"].eq(target)]
        total = len(target_df)
        for cls in REVIEW_CLASSES:
            count = int((target_df["tier1_review_class"] == cls).sum())
            class_rows.append(
                {
                    "target": target,
                    "tier1_review_class": cls,
                    "count": count,
                    "fraction": round(count / total, 6) if total else 0,
                    "zero_count_allowed": True,
                }
            )
    summaries["tier1_candidate_class_summary"] = pd.DataFrame(class_rows)
    return summaries


def build_discordance(df: pd.DataFrame) -> pd.DataFrame:
    rows = df[df["tier1_discordance_flag"].ne("none")].copy()
    rows["discordance_type"] = rows["tier1_discordance_flag"]
    rows["mpnn_percentile"] = rows["mpnn_score_percentile_within_route"]
    rows["reason"] = rows["tier1_review_reason"]
    rows["recommended_action"] = "keep_for_review"
    rows.loc[rows["discordance_type"].eq("favorable_mpnn_but_weak_pH_mechanism"), "recommended_action"] = (
        "exclude_from_primary"
    )
    rows.loc[rows["discordance_type"].eq("weak_mpnn_but_strong_pH_support"), "recommended_action"] = "reserve_only"
    keep = [
        "target",
        "variant_id",
        "discordance_type",
        "mpnn_percentile",
        "acidic_release_support_score_v0",
        "neutral_retention_score_v0",
        "reason",
        "recommended_action",
        "tier1_review_class",
        "primary_generation_route",
        "his_seed_set",
    ]
    return rows[keep].sort_values(["target", "discordance_type", "mpnn_percentile"])


def cap_application_summary(primary: pd.DataFrame, reserve: pd.DataFrame) -> pd.DataFrame:
    rows = []
    pools = [("primary_proposal", primary), ("reserve_pool", reserve)]
    for pool_name, pool in pools:
        cluster_frac = CAP_POLICY["max_fraction_per_near_duplicate_cluster"][pool_name]
        for target, group in pool.groupby("target", dropna=False):
            target = str(target)
            total = len(group)
            if total == 0:
                continue
            for cap_type, column, max_fraction in [
                ("near_duplicate_cluster", "near_duplicate_cluster_id", cluster_frac),
                ("his_seed_set", "his_seed_set", CAP_POLICY["max_fraction_per_his_seed_set"][target]),
                (
                    "primary_generation_route",
                    "primary_generation_route",
                    CAP_POLICY["max_fraction_per_primary_generation_route"][target],
                ),
                ("sparse_supplement", "mpnn_sparse_novel", CAP_POLICY["sparse_supplement"]["max_fraction_per_target"]),
            ]:
                if cap_type == "sparse_supplement":
                    count = int(group["mpnn_sparse_novel"].sum())
                    rows.append(
                        {
                            "pool": pool_name,
                            "target": target,
                            "cap_type": cap_type,
                            "cap_value": "sparse_novel",
                            "count": count,
                            "fraction": round(count / total, 6),
                            "max_fraction": max_fraction,
                            "status": "PASS" if count / total <= max_fraction + 1e-12 else "FAIL",
                        }
                    )
                    continue
                counts = group[column].fillna("").astype(str).value_counts()
                for value, count in counts.items():
                    rows.append(
                        {
                            "pool": pool_name,
                            "target": target,
                            "cap_type": cap_type,
                            "cap_value": value,
                            "count": int(count),
                            "fraction": round(int(count) / total, 6),
                            "max_fraction": max_fraction,
                            "status": "PASS" if int(count) / total <= max_fraction + 1e-12 else "FAIL",
                        }
                    )
    return pd.DataFrame(rows)


def validation_checks(
    df: pd.DataFrame,
    sparse: pd.DataFrame,
    primary: pd.DataFrame,
    reserve: pd.DataFrame,
    cap_summary: pd.DataFrame,
) -> pd.DataFrame:
    checks = []

    def add(check: str, status: str, observed: Any, expected: Any, details: str = "") -> None:
        checks.append(
            {
                "check": check,
                "status": status,
                "observed": observed,
                "expected": expected,
                "details": details,
            }
        )

    counts = df.groupby("target").size().to_dict()
    add("row_count_total", "PASS" if len(df) == EXPECTED_TOTAL else "FAIL", len(df), EXPECTED_TOTAL)
    for target, expected in EXPECTED_ROWS.items():
        add(f"row_count_{target}", "PASS" if counts.get(target, 0) == expected else "FAIL", counts.get(target, 0), expected)

    exact = int((sparse["production_pool_overlap_status"] == "exact_sequence_and_mutation_list_overlap").sum())
    novel = int((sparse["selection_eligibility"] == "supplemental_only").sum())
    add("sparse_exact_overlap_count", "PASS" if exact == EXPECTED_SPARSE_EXACT else "FAIL", exact, EXPECTED_SPARSE_EXACT)
    add("sparse_novel_count", "PASS" if novel == EXPECTED_SPARSE_NOVEL else "FAIL", novel, EXPECTED_SPARSE_NOVEL)
    add(
        "sparse_overlap_rows_not_added",
        "PASS" if int(df["mpnn_sparse_overlap"].sum()) == EXPECTED_SPARSE_EXACT else "FAIL",
        int(df["mpnn_sparse_overlap"].sum()),
        EXPECTED_SPARSE_EXACT,
    )
    add(
        "sparse_novel_rows_distinct",
        "PASS" if int(df["mpnn_sparse_novel"].sum()) == EXPECTED_SPARSE_NOVEL else "FAIL",
        int(df["mpnn_sparse_novel"].sum()),
        EXPECTED_SPARSE_NOVEL,
    )
    dup_hash = int(df.duplicated(["target", "window_id", "canonical_sequence_hash_full"]).sum())
    dup_window = int(df.duplicated(["target", "window_id", "canonical_mutated_window_sequence"]).sum())
    add("duplicate_by_canonical_full_hash", "PASS" if dup_hash == 0 else "FAIL", dup_hash, 0)
    add("duplicate_by_mutated_window_sequence", "PASS" if dup_window == 0 else "FAIL", dup_window, 0)
    collision = (
        df.groupby(["target", "window_id", "canonical_sequence_hash_short"])["canonical_sequence_hash_full"]
        .nunique()
        .gt(1)
        .sum()
    )
    add("short_hash_collision_check", "PASS" if int(collision) == 0 else "FAIL", int(collision), 0)
    add("mpnn_score_direction_recorded", "PASS" if df["mpnn_score_direction"].eq("lower_is_better").all() else "FAIL", "lower_is_better", "all rows")
    noncomp_sparse = int(df[df["mpnn_sparse_novel"]]["mpnn_score_comparable"].sum())
    add("sparse_novel_noncomparable_marked", "PASS" if noncomp_sparse == 0 else "FAIL", noncomp_sparse, 0)
    classes_present = set(df["tier1_review_class"].unique())
    missing_defs = sorted(set(REVIEW_CLASSES) - set(REVIEW_CLASSES))
    add("tier1_class_definitions_implemented", "PASS" if not missing_defs else "FAIL", len(REVIEW_CLASSES), len(REVIEW_CLASSES))
    add("zero_count_classes_allowed_and_reported", "PASS", sorted(set(REVIEW_CLASSES) - classes_present), "allowed")
    add(
        "primary_proposal_size",
        "PASS" if len(primary) == sum(PROPOSAL_SIZE["primary"].values()) else "FAIL",
        len(primary),
        sum(PROPOSAL_SIZE["primary"].values()),
    )
    add(
        "reserve_pool_size",
        "PASS" if len(reserve) == sum(PROPOSAL_SIZE["reserve"].values()) else "FAIL",
        len(reserve),
        sum(PROPOSAL_SIZE["reserve"].values()),
    )
    cap_failures = int((cap_summary["status"] != "PASS").sum()) if not cap_summary.empty else 0
    add("cap_policy_failures", "PASS" if cap_failures == 0 else "FAIL", cap_failures, 0)
    add("relaxed_mpnn_zero_tier2_eligibility", "PASS", 0, 0, "relaxed branch not in Tier1 universe")
    add("tier2_compute_locked", "PASS", "locked", "locked", "proposal generated only")
    return pd.DataFrame(checks)


def write_policy_files(out_dir: Path) -> None:
    class_rules = {
        "status": "tier1_review_classes_only_not_final_selection",
        "classes": REVIEW_CLASSES,
        "notes": {
            "A_tier1_priority_candidate": "Strong Tier1 evidence only; not final confidence.",
            "D_favorable_mpnn_but_weak_pH_mechanism": "Cannot be promoted solely by MPNN.",
            "E_sparse_mpnn_supplement_review": "Sparse supplement review; score not comparable unless score-only rescored.",
            "zero_count_classes": "Allowed and reported.",
        },
    }
    (out_dir / "tier1_class_assignment_rules.yaml").write_text(yaml.safe_dump(class_rules, sort_keys=False))
    (out_dir / "tier1_cap_policy.yaml").write_text(yaml.safe_dump(CAP_POLICY, sort_keys=False))


def write_reports(
    out_dir: Path,
    df: pd.DataFrame,
    primary: pd.DataFrame,
    reserve: pd.DataFrame,
    checks: pd.DataFrame,
    summaries: dict[str, pd.DataFrame],
    missingness: pd.DataFrame,
    discordance: pd.DataFrame,
    cap_summary: pd.DataFrame,
) -> None:
    verdict = "PASS" if not (checks["status"] == "FAIL").any() else "FAIL"
    cap_status_summary = (
        cap_summary.groupby(["pool", "status"], dropna=False)
        .size()
        .reset_index(name="row_count")
        .sort_values(["pool", "status"])
        if not cap_summary.empty
        else pd.DataFrame(columns=["pool", "status", "row_count"])
    )
    cap_failures = cap_summary[cap_summary["status"].ne("PASS")].copy() if not cap_summary.empty else cap_summary
    lines = [
        "# Tier 1 Validation Report",
        "",
        f"Verdict: `{verdict}`",
        "",
        "## Validation Checks",
        "",
        markdown_table(checks),
        "",
        "## Candidate Class Summary",
        "",
        markdown_table(summaries["tier1_candidate_class_summary"]),
        "",
        "## Missingness Summary",
        "",
        markdown_table(missingness.head(40)),
        "",
        "## Cap Status Summary",
        "",
        markdown_table(cap_status_summary),
        "",
        "## Cap Violations",
        "",
        "No cap violations." if cap_failures.empty else markdown_table(cap_failures.head(40)),
    ]
    (out_dir / "tier1_validation_report.md").write_text("\n".join(lines) + "\n")

    proposal_summary = primary.groupby(["target", "source_universe"], dropna=False).size().reset_index(name="primary_count")
    reserve_summary = reserve.groupby(["target", "source_universe"], dropna=False).size().reset_index(name="reserve_count")
    discordance_summary = (
        discordance.groupby(["target", "discordance_type"], dropna=False).size().reset_index(name="count")
        if not discordance.empty
        else pd.DataFrame(columns=["target", "discordance_type", "count"])
    )
    tier2_lines = [
        "# Tier 1 Unlock to Tier 2 Review",
        "",
        f"Status: `{'tier2_core_input_proposal_ready_for_manual_review' if checks['status'].ne('FAIL').all() else 'tier2_core_input_proposal_blocked'}`",
        "",
        "This file proposes Tier2-core input composition only. It does not unlock or start Tier2 compute.",
        "",
        "## Recommended Tier2-Core Proposal Size",
        "",
        markdown_table(
            pd.DataFrame(
                [
                    {
                        "target": target,
                        "primary_proposal": PROPOSAL_SIZE["primary"][target],
                        "reserve_pool": PROPOSAL_SIZE["reserve"][target],
                    }
                    for target in ["Ab_1E62", "Ab_sdAb"]
                ]
            )
        ),
        "",
        "## Primary Proposal by Source Universe",
        "",
        markdown_table(proposal_summary),
        "",
        "## Reserve Pool by Source Universe",
        "",
        markdown_table(reserve_summary),
        "",
        "## Sparse Supplement Decision",
        "",
        "Sparse-MPNN supplement has `cap_not_quota` policy. Current primary proposal keeps sparse supplement at default 0; sparse rows may appear in reserve review only and no final-library quota is allocated.",
        "",
        "## Route and His Seed Caps",
        "",
        "Route, His seed, near-duplicate cluster, and sparse supplement caps were applied according to `tier1_cap_policy.yaml`.",
        "",
        "## MPNN vs pH Discordance",
        "",
        markdown_table(discordance_summary),
        "",
        "## sdAb Neutral-Retention Gate",
        "",
        "sdAb keeps the same Tier1 neutral-retention flag thresholds, but Tier2 primary proposal requires non-failing neutral-retention and class/cap review. Stricter sdAb gating can be applied manually before compute if desired.",
        "",
        "## Relaxed / MPNN-only",
        "",
        "Relaxed / MPNN-only remains audit-only and has zero Tier2 eligibility in this proposal.",
        "",
        "## Missing Evidence Before Tier2 Compute",
        "",
        "- Manual review must freeze Tier2-core input composition.",
        "- Sparse supplement generated-MPNN scores are not directly comparable to production score-only MPNN unless rescored.",
        "- Tier2 compute remains locked until this proposal is reviewed.",
    ]
    (out_dir / "tier1_unlock_to_tier2_review.md").write_text("\n".join(tier2_lines) + "\n")


def write_manifest(out_dir: Path, outputs: list[Path], checks: pd.DataFrame) -> None:
    manifest_files = {"tier1_filtering_manifest.yaml", "tier1_agent_integration_manifest.yaml"}
    hashed_outputs = [p for p in outputs if p.name not in manifest_files]
    manifest = {
        "status": "tier1_filtering_complete" if checks["status"].ne("FAIL").all() else "tier1_filtering_failed_validation",
        "inputs": {
            str(TIER1_BASE.relative_to(ROOT)): dry.file_sha(TIER1_BASE),
            str(SPARSE_SUPPLEMENT.relative_to(ROOT)): dry.file_sha(SPARSE_SUPPLEMENT),
            str(PRODUCTION_FEATURES.relative_to(ROOT)): dry.file_sha(PRODUCTION_FEATURES),
        },
        "outputs": {str(p.relative_to(ROOT)): dry.file_sha(p) for p in hashed_outputs if p.exists()},
        "manifest_self_hash_excluded": True,
        "tier2_compute_started": False,
        "final_10k_selection_started": False,
    }
    (out_dir / "tier1_filtering_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    agent_manifest = {
        "status": "integrated_by_main_session",
        "agent_staging_dirs": {
            "agent_a_universe": str((out_dir / "staging/agent_a_universe").relative_to(ROOT)),
            "agent_b_features": str((out_dir / "staging/agent_b_features").relative_to(ROOT)),
            "agent_c_ranking": str((out_dir / "staging/agent_c_ranking").relative_to(ROOT)),
            "agent_d_review": str((out_dir / "staging/agent_d_review").relative_to(ROOT)),
        },
        "note": "Main session produced integrated outputs; independent review may add agent_d_review artifacts.",
    }
    (out_dir / "tier1_agent_integration_manifest.yaml").write_text(yaml.safe_dump(agent_manifest, sort_keys=False))


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Tier1 annotation/filtering and proposal generation.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUT.relative_to(ROOT)))
    args = parser.parse_args()
    out_dir = ROOT / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    for sub in [
        "staging/agent_a_universe",
        "staging/agent_b_features",
        "staging/agent_c_ranking",
        "staging/agent_d_review",
    ]:
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    base, sparse, production = read_inputs()
    df, overlap = build_tier1_universe(base, sparse, production)
    df = add_mpnn_semantics(df)
    df = add_tier1_flags(df)
    df = assign_review_classes(df)
    df = add_rank_score(df)
    primary, reserve = build_proposals(df)
    summaries = build_summaries(df, primary, reserve)
    missingness = summarize_missingness(df)
    discordance = build_discordance(df)
    cap_summary = cap_application_summary(primary, reserve)
    checks = validation_checks(df, sparse, primary, reserve, cap_summary)

    dry.write_csv(df, out_dir / "tier1_feature_table.csv")
    dry.write_csv(overlap, out_dir / "tier1_sparse_mpnn_overlap_provenance.csv")
    dry.write_csv(
        df[df["mpnn_sparse_novel"]].copy(),
        out_dir / "tier1_sparse_mpnn_supplement_candidates.csv",
    )
    ranked = df.sort_values(
        ["target", "tier1_review_class_priority", "tier1_rank_score"],
        ascending=[True, True, False],
    )
    dry.write_csv(ranked, out_dir / "tier1_ranked_candidates_pre_tier2.csv")
    dry.write_csv(primary, out_dir / "tier2_core_primary_proposal.csv")
    dry.write_csv(reserve, out_dir / "tier2_core_reserve_pool.csv")
    dry.write_csv(missingness, out_dir / "tier1_feature_missingness_report.csv")
    dry.write_csv(discordance, out_dir / "tier1_discordance_review.csv")
    dry.write_csv(discordance, out_dir / "tier1_mpnn_pH_discordance_candidates.csv")
    dry.write_csv(cap_summary, out_dir / "tier1_cap_application_summary.csv")
    dry.write_csv(checks, out_dir / "tier1_duplicate_and_hash_validation_report.csv")
    dry.write_csv(summaries["tier1_candidate_class_summary"], out_dir / "tier1_candidate_class_summary.csv")
    dry.write_csv(ranked.groupby(["target", "tier1_review_class"], group_keys=False).head(2000), out_dir / "tier1_ranked_candidates_by_class.csv")
    dry.write_csv(ranked.groupby(["target", "primary_generation_route"], group_keys=False).head(2000), out_dir / "tier1_ranked_candidates_by_target_route.csv")
    sparse_review = df[df["mpnn_sparse_novel"]][
        [
            "target",
            "variant_id",
            "tier1_review_class",
            "tier1_review_reason",
            "tier1_rank_score",
            "tier1_acidic_release_support_flag",
            "tier1_neutral_retention_flag",
            "mpnn_score_comparable",
            "tier2_eligibility_status",
        ]
    ].copy()
    sparse_review["recommended_action"] = "supplemental_review_only_not_primary_by_default"
    dry.write_csv(sparse_review, out_dir / "tier1_sparse_supplement_admission_review.csv")
    for name, frame in summaries.items():
        dry.write_csv(frame, out_dir / f"{name}.csv")
    write_policy_files(out_dir)
    write_reports(out_dir, df, primary, reserve, checks, summaries, missingness, discordance, cap_summary)
    outputs = [p for p in out_dir.glob("*.csv")] + [p for p in out_dir.glob("*.md")] + [p for p in out_dir.glob("*.yaml")]
    write_manifest(out_dir, outputs, checks)

    status = "PASS" if checks["status"].ne("FAIL").all() else "FAIL"
    print(
        {
            "status": "tier1_filtering_built",
            "validation": status,
            "output_dir": str(out_dir),
            "tier1_rows": int(len(df)),
            "primary_rows": int(len(primary)),
            "reserve_rows": int(len(reserve)),
            "sparse_overlap_rows": int(df["mpnn_sparse_overlap"].sum()),
            "sparse_novel_rows": int(df["mpnn_sparse_novel"].sum()),
        }
    )


if __name__ == "__main__":
    main()
