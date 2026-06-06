#!/usr/bin/env python3
"""Build P0 evidence review and sparse-MPNN supplement artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT
OUT = ROOT / "results/initial_design_generation"
P0 = OUT / "p0_mpnn"
SPARSE = OUT / "p0_mpnn_constrained_sparse_full"
RUNNER = OUT / "p0_mpnn_runner_inputs_constrained_sparse_full"
PRODUCTION_FEATURES = P0 / "production_pool_features.csv"
EVIDENCE_OUT = OUT / "p0_evidence_review"


def sha(text: str, n: int = 12) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:n]


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


def parse_header(header: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in header.split(","):
        if "=" in part:
            key, value = part.strip().split("=", 1)
            out[key.strip()] = value.strip()
    return out


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def split_mutations(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [item for item in str(value).split(";") if item and item.lower() != "nan"]


def is_his_mutation(label: str) -> bool:
    return bool(re.match(r"[A-Z][A-Z]\d+H$", label))


def rescue_mutation_list(mutation_list: object) -> str:
    return ";".join(m for m in split_mutations(mutation_list) if not is_his_mutation(m))


def load_sparse_pool(p0_dir: Path) -> pd.DataFrame:
    frames = []
    for target in ["1E62", "sdAb"]:
        path = p0_dir / f"constrained_mpnn_rescue_pool_{target}.csv"
        if path.exists():
            frames.append(pd.read_csv(path))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def sparse_score_records(runner_dir: Path) -> pd.DataFrame:
    rows: list[dict] = []
    for fasta in sorted((runner_dir / "mpnn_out").glob("**/seqs/*.fa")):
        records = parse_fasta(fasta)
        if not records:
            continue
        job_name = fasta.stem
        parent_meta = parse_header(records[0][0])
        parent_score = float(parent_meta.get("score", "nan"))
        parent_global = float(parent_meta.get("global_score", "nan"))
        for header, sequence in records[1:]:
            if not header.startswith("T="):
                continue
            meta = parse_header(header)
            sample = str(meta.get("sample", ""))
            score = float(meta.get("score", "nan"))
            global_score = float(meta.get("global_score", "nan"))
            rows.append(
                {
                    "generation_record_id": job_name,
                    "mpnn_sample_id": sample,
                    "sequence_hash": sha(sequence, 12),
                    "mpnn_score": score,
                    "mpnn_global_score": global_score,
                    "mpnn_native_score": parent_score,
                    "mpnn_native_global_score": parent_global,
                    "mpnn_parent_delta": score - parent_score,
                    "mpnn_seq_recovery": float(meta.get("seq_recovery", "nan")),
                    "mpnn_temperature_from_fasta": meta.get("T", ""),
                    "mpnn_model_version": parent_meta.get("model_name", ""),
                    "mpnn_git_hash": parent_meta.get("git_hash", ""),
                    "mpnn_random_seed": parent_meta.get("seed", ""),
                }
            )
    return pd.DataFrame(rows)


def list_join(values: pd.Series, limit: int = 20) -> str:
    unique = [str(v) for v in values.dropna().astype(str).unique() if str(v) and str(v).lower() != "nan"]
    unique = sorted(unique)
    if len(unique) > limit:
        return ";".join(unique[:limit]) + f";...(+{len(unique) - limit})"
    return ";".join(unique)


def build_sparse_supplement(production: pd.DataFrame, sparse_pool: pd.DataFrame, score: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pool = sparse_pool.copy()
    pool["mpnn_sample_id"] = pool["mpnn_sample_id"].astype(str)
    score["mpnn_sample_id"] = score["mpnn_sample_id"].astype(str)
    pool = pool.merge(
        score,
        on=["generation_record_id", "mpnn_sample_id", "sequence_hash"],
        how="left",
        validate="many_to_one",
    )
    pool["rescue_mutation_list"] = pool["mutation_list"].map(rescue_mutation_list)
    pool["rescue_subset_id"] = pool["generation_record_id"].str.extract(r"(subset\d+)$", expand=False).fillna("")
    grouped = pool.groupby(["target", "sequence_hash", "mutation_list"], dropna=False)
    supp = grouped.agg(
        sequence=("sequence", "first"),
        generated_count=("variant_id", "count"),
        hard_filter_status=("hard_filter_status", lambda s: "pass" if (s == "pass").all() else "mixed_or_fail"),
        his_seed_set=("his_seed_set", "first"),
        rescue_mutation_list=("rescue_mutation_list", "first"),
        rescue_subset_id=("rescue_subset_id", list_join),
        mutation_count=("mutation_count", "first"),
        His_count=("His_count", "first"),
        rescue_count=("rescue_count", "first"),
        mpnn_score=("mpnn_score", "min"),
        mpnn_score_mean=("mpnn_score", "mean"),
        mpnn_global_score=("mpnn_global_score", "min"),
        mpnn_parent_delta=("mpnn_parent_delta", "min"),
        mpnn_parent_delta_mean=("mpnn_parent_delta", "mean"),
        mpnn_seq_recovery_mean=("mpnn_seq_recovery", "mean"),
        mpnn_temperature=("mpnn_temperature", list_join),
        mpnn_model_version=("mpnn_model_version_y", lambda s: list_join(s) if "mpnn_model_version_y" in pool else ""),
        source_generation_records=("generation_record_id", list_join),
    ).reset_index()
    supp["source_branch"] = "repaired_sparse_constrained_MPNN"
    supp["new_canonical_nxs_t_status"] = "pass"
    supp["rescue_count_status"] = "pass"
    supp["production_sequence_hash_10"] = supp["sequence_hash"].astype(str).str[:10]
    supp["variant_id"] = supp.apply(
        lambda r: f"mpnn_sparse_{str(r.target).replace('Ab_', '')}_{r.sequence_hash}", axis=1
    )
    supp["mpnn_score_percentile_by_target"] = supp.groupby("target")["mpnn_score"].rank(pct=True, ascending=True)
    supp["mpnn_score_percentile_by_route"] = supp.groupby(["target", "source_branch"])["mpnn_score"].rank(
        pct=True, ascending=True
    )
    supp["mpnn_score_percentile_by_seed"] = supp.groupby(["target", "his_seed_set"])["mpnn_score"].rank(
        pct=True, ascending=True
    )
    supp["mpnn_score_percentile_by_rescue"] = supp.groupby(["target", "rescue_mutation_list"])["mpnn_score"].rank(
        pct=True, ascending=True
    )

    prod = production.copy()
    prod["production_sequence_hash_10"] = prod["sequence_hash"].astype(str)
    exact = prod.groupby(["target", "production_sequence_hash_10", "mutation_list"], dropna=False).agg(
        exact_existing_production_variant_id=("variant_id", list_join),
        exact_existing_routes=("all_source_routes", list_join),
        exact_existing_count=("variant_id", "count"),
    ).reset_index()
    seq = prod.groupby(["target", "production_sequence_hash_10"], dropna=False).agg(
        sequence_existing_production_variant_id=("variant_id", list_join),
        sequence_existing_count=("variant_id", "count"),
    ).reset_index()
    mut = prod.groupby(["target", "mutation_list"], dropna=False).agg(
        mutation_existing_production_variant_id=("variant_id", list_join),
        mutation_existing_routes=("all_source_routes", list_join),
        mutation_existing_count=("variant_id", "count"),
    ).reset_index()
    supp = supp.merge(exact, on=["target", "production_sequence_hash_10", "mutation_list"], how="left")
    supp = supp.merge(seq, on=["target", "production_sequence_hash_10"], how="left")
    supp = supp.merge(mut, on=["target", "mutation_list"], how="left")
    supp["existing_production_variant_id"] = supp["exact_existing_production_variant_id"].fillna(
        supp["sequence_existing_production_variant_id"]
    ).fillna(supp["mutation_existing_production_variant_id"])
    supp["production_pool_overlap_status"] = "novel_sequence_and_mutation_list"
    supp.loc[supp["mutation_existing_count"].notna(), "production_pool_overlap_status"] = "mutation_list_overlap_only"
    supp.loc[supp["sequence_existing_count"].notna(), "production_pool_overlap_status"] = "sequence_hash_overlap_only"
    supp.loc[supp["exact_existing_count"].notna(), "production_pool_overlap_status"] = (
        "exact_sequence_and_mutation_list_overlap"
    )
    supp["all_source_routes"] = supp["mutation_existing_routes"].fillna("")
    supp["all_source_routes"] = supp["all_source_routes"].where(
        supp["all_source_routes"].eq(""),
        supp["all_source_routes"] + ";repaired_sparse_constrained_MPNN",
    )
    supp.loc[supp["all_source_routes"].eq(""), "all_source_routes"] = "repaired_sparse_constrained_MPNN"
    supp["selection_eligibility"] = "supplemental_only"
    supp.loc[
        supp["production_pool_overlap_status"].eq("exact_sequence_and_mutation_list_overlap"),
        "selection_eligibility",
    ] = "provenance_only_existing_production"

    clustered = dry.assign_near_duplicate_clusters(
        supp[
            [
                "target",
                "variant_id",
                "sequence",
                "mutation_list",
                "his_seed_set",
            ]
        ].copy()
    )
    supp = supp.merge(
        clustered[
            [
                "variant_id",
                "sequence_hamming_cluster_id",
                "mutation_set_jaccard_cluster_id",
                "his_seed_cluster_id",
                "near_duplicate_cluster_id",
                "cluster_size",
            ]
        ],
        on="variant_id",
        how="left",
    )
    overlap = supp.groupby(["target", "production_pool_overlap_status", "selection_eligibility"], dropna=False).agg(
        candidate_count=("variant_id", "count"),
        generated_raw_count=("generated_count", "sum"),
    ).reset_index()
    near_dup = supp.groupby(["target", "near_duplicate_cluster_id"], dropna=False).agg(
        cluster_size=("variant_id", "count"),
        generated_raw_count=("generated_count", "sum"),
        his_seed_sets=("his_seed_set", list_join),
        rescue_mutations=("rescue_mutation_list", list_join),
        overlap_statuses=("production_pool_overlap_status", list_join),
        eligibility=("selection_eligibility", list_join),
    ).reset_index()
    target_counts = supp.groupby("target").size().to_dict()
    near_dup["fraction_within_target"] = near_dup.apply(
        lambda r: r["cluster_size"] / target_counts.get(r["target"], 1), axis=1
    )
    near_dup = near_dup.sort_values(["target", "cluster_size"], ascending=[True, False])
    return supp, overlap, near_dup


def heuristic_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    non_his = out["rescue_count"].fillna(0).astype(float)
    mut_count = out["mutation_count"].fillna(0).astype(float)
    his_count = out["His_count"].fillna(0).astype(float)
    foldx_proxy = (1.0 - 0.12 * mut_count - 0.05 * non_his).clip(lower=0.0)
    neutral = (1.0 - 0.10 * non_his - 0.07 * his_count).clip(lower=0.0)
    acidic = (0.20 + 0.25 * his_count).clip(upper=1.0)
    out["neutral_retention_score"] = neutral.round(4)
    out["acidic_release_support_score"] = acidic.round(4)
    out["global_weakening_risk_score"] = (0.08 * mut_count + 0.06 * non_his).clip(upper=1.0).round(4)
    out["foldx_proxy_feature"] = foldx_proxy.round(4)
    out["hit_likelihood_score_v0"] = (0.35 * neutral + 0.35 * acidic + 0.30 * foldx_proxy).round(4)
    out["score_confidence_flag"] = "sparse_mpnn_supplement_proxy_only"
    return out


def build_tier1_feature_table(production: pd.DataFrame, supplement: pd.DataFrame) -> pd.DataFrame:
    prod_cols = [
        "target",
        "variant_id",
        "sequence_hash",
        "mutation_list",
        "mutation_count",
        "His_count",
        "his_seed_set",
        "primary_generation_route",
        "all_source_routes",
        "rescue_mutation_list",
        "rescue_signature",
        "rescue_count",
        "hard_filter_status",
        "liability_flags",
        "near_duplicate_cluster_id",
        "cluster_size",
        "hit_likelihood_score_v0",
        "neutral_retention_score",
        "acidic_release_support_score",
        "global_weakening_risk_score",
        "display_or_expression_risk_score",
        "glycan_or_epitope_risk_score",
        "foldx_proxy_feature",
        "interface_hotspot_proxy",
        "buildability_light_status",
        "mpnn_score_status",
        "mpnn_total_score_per_residue",
        "mpnn_parent_delta",
        "mpnn_route_percentile",
        "mpnn_his_seed_percentile",
        "mpnn_rescue_type_percentile",
        "mpnn_score_confidence",
    ]
    prod = production[prod_cols].copy()
    prod["source_universe"] = "production_initial_pool"
    prod["source_branch"] = "rule_constrained_generation"
    prod["production_pool_overlap_status"] = "current_production_pool"
    prod["selection_eligibility"] = "tier1_candidate"
    prod["mpnn_score_type"] = "score_only_full_design_chain"
    prod["p0_evidence_class"] = "production_rule_candidate_with_mpnn_score"

    novel = supplement[
        supplement["selection_eligibility"].eq("supplemental_only")
        & supplement["production_pool_overlap_status"].ne("exact_sequence_and_mutation_list_overlap")
    ].copy()
    novel = heuristic_scores(novel)
    supp = pd.DataFrame(
        {
            "target": novel["target"],
            "variant_id": novel["variant_id"],
            "sequence_hash": novel["sequence_hash"],
            "mutation_list": novel["mutation_list"],
            "mutation_count": novel["mutation_count"],
            "His_count": novel["His_count"],
            "his_seed_set": novel["his_seed_set"],
            "primary_generation_route": "repaired_sparse_constrained_MPNN",
            "all_source_routes": novel["all_source_routes"],
            "rescue_mutation_list": novel["rescue_mutation_list"],
            "rescue_signature": novel["rescue_mutation_list"],
            "rescue_count": novel["rescue_count"],
            "hard_filter_status": novel["hard_filter_status"],
            "liability_flags": "",
            "near_duplicate_cluster_id": novel["near_duplicate_cluster_id"],
            "cluster_size": novel["cluster_size"],
            "hit_likelihood_score_v0": novel["hit_likelihood_score_v0"],
            "neutral_retention_score": novel["neutral_retention_score"],
            "acidic_release_support_score": novel["acidic_release_support_score"],
            "global_weakening_risk_score": novel["global_weakening_risk_score"],
            "display_or_expression_risk_score": "",
            "glycan_or_epitope_risk_score": "",
            "foldx_proxy_feature": novel["foldx_proxy_feature"],
            "interface_hotspot_proxy": "sparse_mpnn_rescue_supplement",
            "buildability_light_status": "pending_for_supplement",
            "mpnn_score_status": "generated_by_mpnn_with_score",
            "mpnn_total_score_per_residue": novel["mpnn_score"],
            "mpnn_parent_delta": novel["mpnn_parent_delta"],
            "mpnn_route_percentile": novel["mpnn_score_percentile_by_route"],
            "mpnn_his_seed_percentile": novel["mpnn_score_percentile_by_seed"],
            "mpnn_rescue_type_percentile": novel["mpnn_score_percentile_by_rescue"],
            "mpnn_score_confidence": "generated_sparse_constrained_score",
            "source_universe": "repaired_sparse_constrained_MPNN_supplement",
            "source_branch": "repaired_sparse_constrained_MPNN",
            "production_pool_overlap_status": novel["production_pool_overlap_status"],
            "selection_eligibility": novel["selection_eligibility"],
            "mpnn_score_type": "generated_candidate_score",
            "p0_evidence_class": "deduplicated_sparse_mpnn_supplement_candidate",
        }
    )
    combined = pd.concat([prod, supp], ignore_index=True)
    high_ph = (
        (combined["acidic_release_support_score"].fillna(0) >= 0.65)
        & (combined["neutral_retention_score"].fillna(0) >= 0.70)
        & (combined["His_count"].fillna(0) >= 1)
    )
    favorable = combined["mpnn_route_percentile"].fillna(1.0) <= 0.25
    unfavorable = combined["mpnn_route_percentile"].fillna(0.0) >= 0.90
    weak_ph = (combined["acidic_release_support_score"].fillna(0) < 0.50) | (
        combined["His_count"].fillna(0) < 1
    )
    combined["p0_evidence_flag"] = "standard_review"
    combined.loc[high_ph & favorable, "p0_evidence_flag"] = "high_pH_support_and_favorable_mpnn"
    combined.loc[high_ph & unfavorable, "p0_evidence_flag"] = "high_pH_support_but_unfavorable_mpnn"
    combined.loc[weak_ph & favorable, "p0_evidence_flag"] = "favorable_mpnn_but_weak_pH_mechanism"
    combined["tier1_review_status"] = "eligible_for_p0_review"
    combined.loc[
        combined["source_universe"].eq("repaired_sparse_constrained_MPNN_supplement"),
        "tier1_review_status",
    ] = "eligible_supplement_after_dedup_review"
    return combined


def build_percentile_table(tier1: pd.DataFrame) -> pd.DataFrame:
    rows = tier1[
        [
            "target",
            "variant_id",
            "source_universe",
            "primary_generation_route",
            "his_seed_set",
            "rescue_signature",
            "mpnn_score_type",
            "mpnn_total_score_per_residue",
        ]
    ].copy()
    rows["mpnn_percentile_by_target"] = rows.groupby(["mpnn_score_type", "target"])[
        "mpnn_total_score_per_residue"
    ].rank(pct=True, ascending=True)
    rows["mpnn_percentile_by_target_route"] = rows.groupby(
        ["mpnn_score_type", "target", "primary_generation_route"]
    )["mpnn_total_score_per_residue"].rank(pct=True, ascending=True)
    rows["mpnn_percentile_by_target_seed"] = rows.groupby(["mpnn_score_type", "target", "his_seed_set"])[
        "mpnn_total_score_per_residue"
    ].rank(pct=True, ascending=True)
    rows["mpnn_percentile_by_target_rescue"] = rows.groupby(["mpnn_score_type", "target", "rescue_signature"])[
        "mpnn_total_score_per_residue"
    ].rank(pct=True, ascending=True)
    return rows


def write_sparse_report(out_dir: Path, supplement: pd.DataFrame, overlap: pd.DataFrame, near_dup: pd.DataFrame) -> None:
    summary = supplement.groupby("target", dropna=False).agg(
        dedup_candidate_count=("variant_id", "count"),
        generated_raw_count=("generated_count", "sum"),
        novel_candidate_count=(
            "production_pool_overlap_status",
            lambda s: int((s == "novel_sequence_and_mutation_list").sum()),
        ),
        supplemental_eligible_count=("selection_eligibility", lambda s: int((s == "supplemental_only").sum())),
        exact_overlap_count=(
            "production_pool_overlap_status",
            lambda s: int((s == "exact_sequence_and_mutation_list_overlap").sum()),
        ),
        largest_cluster_size=("cluster_size", "max"),
    ).reset_index()
    lines = [
        "# Sparse-MPNN Supplement Dedup Report",
        "",
        "Status: `deduplicated_sparse_supplement_ready_for_p0_review`",
        "",
        "The repaired sparse constrained-MPNN raw pool is not merged directly. This report deduplicates it and checks overlap against the current production pool.",
        "",
        "## Target Summary",
        "",
        markdown_table(summary),
        "",
        "## Production Overlap",
        "",
        markdown_table(overlap),
        "",
        "## Top Near-Duplicate Clusters",
        "",
        markdown_table(near_dup, max_rows=30),
        "",
        "## Decision Boundary",
        "",
        "- `supplemental_only` rows may enter P0 / Tier 1 prioritization universe.",
        "- `provenance_only_existing_production` rows should not create duplicate candidate rows.",
        "- No sparse-MPNN final-library quota is reserved.",
    ]
    (out_dir / "mpnn_sparse_supplement_dedup_report.md").write_text("\n".join(lines) + "\n")


def write_evidence_report(out_dir: Path, tier1: pd.DataFrame, supplement: pd.DataFrame, overlap: pd.DataFrame) -> None:
    tier_summary = tier1.groupby(["target", "source_universe"], dropna=False).agg(
        row_count=("variant_id", "count"),
        high_pH_support=("p0_evidence_flag", lambda s: int(s.astype(str).str.contains("high_pH_support").sum())),
        favorable_mpnn=("mpnn_route_percentile", lambda s: int((s.fillna(1) <= 0.25).sum())),
    ).reset_index()
    supp_summary = supplement.groupby("target", dropna=False).agg(
        dedup_sparse_count=("variant_id", "count"),
        supplemental_eligible=("selection_eligibility", lambda s: int((s == "supplemental_only").sum())),
        exact_overlap=("production_pool_overlap_status", lambda s: int((s == "exact_sequence_and_mutation_list_overlap").sum())),
    ).reset_index()
    lines = [
        "# P0 Evidence Review Report",
        "",
        "Status: `p0_evidence_review_complete_tier1_annotation_ready`",
        "",
        "## Decisions Applied",
        "",
        "- Score-only ProteinMPNN features are accepted as P0 evidence features.",
        "- Original broad constrained-MPNN generated candidates remain rejected.",
        "- Repaired sparse constrained-MPNN is admitted only as a deduplicated supplemental rescue source.",
        "- Relaxed / MPNN-only remains audit-only with zero preallocated final quota.",
        "- Tier 2 and final 10K remain locked.",
        "",
        "## Tier 1 Feature Table Summary",
        "",
        markdown_table(tier_summary),
        "",
        "## Sparse Supplement Summary",
        "",
        markdown_table(supp_summary),
        "",
        "## Sparse Supplement Overlap",
        "",
        markdown_table(overlap),
        "",
        "## Required Next Step Before Tier 2",
        "",
        "Manual review must approve Tier 2-core input composition, including per-target/route/cluster caps and whether any sparse-MPNN supplement rows are admitted.",
    ]
    (out_dir / "p0_evidence_review_report.md").write_text("\n".join(lines) + "\n")


def write_shortlist_policy(out_dir: Path) -> None:
    lines = [
        "# Tier 1 Shortlist Policy",
        "",
        "Status: `policy_defined_tier2_still_locked`",
        "",
        "## Candidate Classes",
        "",
        "1. High-confidence rule-generated candidates: high pH-support, acceptable neutral-retention, no hard liabilities, and non-poor route/seed-specific MPNN score.",
        "2. Rescue-enriched candidates: strong neutral-retention rescue logic with acceptable pH-support and controlled rescue count.",
        "3. Sparse-MPNN supplement candidates: deduplicated, production-novel, hard-filter passing rows from repaired sparse constrained-MPNN.",
        "4. Audit-only candidates: relaxed / MPNN-only rows; not eligible for Tier 2 unless manually unlocked later.",
        "",
        "## ProteinMPNN Use",
        "",
        "- Use target/route/seed/rescue-specific percentiles.",
        "- Do not use one global cutoff.",
        "- Do not reject His-trigger candidates solely for weaker MPNN score.",
        "- Do not promote candidates with no pH-support solely because MPNN score is favorable.",
        "",
        "## Sparse Supplement Caps",
        "",
        "- P0 / Tier 1 universe: all deduplicated novel sparse-MPNN candidates may be annotated.",
        "- Tier 2-core recommendation cap: 1E62 <= 100-250; sdAb <= 50-150.",
        "- Final 10K: no preallocated quota; default 0-2% only if later evidence supports it.",
        "",
        "## Tier 2 Unlock Requirements",
        "",
        "Tier 2-core remains locked until manual review approves:",
        "",
        "- P0 evidence review report;",
        "- Tier 1 feature table;",
        "- sparse supplement dedup/overlap report;",
        "- per-target/route/cluster caps;",
        "- missing-feature policy;",
        "- demotion rule for high-MPNN-score but weak pH mechanism candidates.",
    ]
    (out_dir / "tier1_shortlist_policy.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build P0 evidence review artifacts.")
    parser.add_argument("--output-dir", default=str(EVIDENCE_OUT.relative_to(ROOT)))
    args = parser.parse_args()
    out_dir = ROOT / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    production = pd.read_csv(PRODUCTION_FEATURES)
    sparse_pool = load_sparse_pool(SPARSE)
    sparse_scores = sparse_score_records(RUNNER)
    supplement, overlap, near_dup = build_sparse_supplement(production, sparse_pool, sparse_scores)
    tier1 = build_tier1_feature_table(production, supplement)
    percentile = build_percentile_table(tier1)

    dry.write_csv(supplement, out_dir / "mpnn_sparse_rescue_supplement_candidates.csv")
    dry.write_csv(overlap, out_dir / "mpnn_sparse_vs_production_overlap.csv")
    dry.write_csv(near_dup, out_dir / "mpnn_sparse_near_duplicate_summary.csv")
    dry.write_csv(tier1, out_dir / "tier1_feature_table_with_mpnn.csv")
    dry.write_csv(percentile, out_dir / "mpnn_score_percentile_by_target_route_seed_rescue.csv")

    low_score_high = tier1[
        (tier1["mpnn_route_percentile"].fillna(1) <= 0.10)
        & (tier1["acidic_release_support_score"].fillna(0) >= 0.65)
        & (tier1["neutral_retention_score"].fillna(0) >= 0.70)
        & (tier1["His_count"].fillna(0) >= 1)
    ].copy()
    high_score_weak = tier1[
        (tier1["mpnn_route_percentile"].fillna(0) >= 0.90)
        & ((tier1["acidic_release_support_score"].fillna(0) < 0.50) | (tier1["His_count"].fillna(0) < 1))
    ].copy()
    dry.write_csv(low_score_high.head(2000), out_dir / "mpnn_low_score_high_pH_support_candidates.csv")
    dry.write_csv(high_score_weak.head(2000), out_dir / "mpnn_high_score_no_pH_mechanism_candidates.csv")

    write_sparse_report(out_dir, supplement, overlap, near_dup)
    write_evidence_report(out_dir, tier1, supplement, overlap)
    write_shortlist_policy(out_dir)

    print(
        {
            "status": "p0_evidence_review_built",
            "output_dir": str(out_dir),
            "production_rows": int(len(production)),
            "sparse_raw_rows": int(len(sparse_pool)),
            "sparse_dedup_rows": int(len(supplement)),
            "tier1_feature_rows": int(len(tier1)),
            "sparse_supplemental_rows": int((supplement["selection_eligibility"] == "supplemental_only").sum()),
        }
    )


if __name__ == "__main__":
    main()
