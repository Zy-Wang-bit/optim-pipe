#!/usr/bin/env python3
"""Build Backfill v2 expanded Tier2 candidate pools.

This stage constructs larger evidence-traceable candidate mother pools from
the frozen Tier2 candidate bank templates. It performs table scans and
constrained local expansion only. It does not run structure prediction,
Rosetta/FoldX/MD, glycan modeling, or final library selection.
"""

from __future__ import annotations

import hashlib
import itertools
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.initial_design_generation import build_expanded_tier2_candidate_pool as v1

TASK_DIR = ROOT / ".tasks/active/initial-design-generation"
BANK_DIR = ROOT / "results/initial_design_generation/tier2_candidate_bank"
OUT_DIR = ROOT / "results/initial_design_generation/expanded_tier2_candidate_pool_v2"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"

BANK = BANK_DIR / "tier2_candidate_bank.csv"
TEMPLATES = BANK_DIR / "tier2_candidate_bank_expansion_templates.csv"

SOURCE_FILES = [
    ("stage2a_candidate_list", ROOT / "results/initial_design_generation/stage1_5_stage2a/stage2a_candidate_list.csv", 0),
    (
        "sdab_recovery_passlike_supplement",
        ROOT / "results/initial_design_generation/sdab_recovery_loop/round_02_passlike_supplement/sdab_candidate_bank.csv",
        1,
    ),
    ("tier2_candidate_snapshot", ROOT / "results/initial_design_generation/tier2_staged/tier2_candidate_snapshot.csv", 2),
    ("tier2_core_reserve_pool", ROOT / "results/initial_design_generation/tier1_filtering/tier2_core_reserve_pool.csv", 3),
    (
        "tier1_ranked_production_pool",
        ROOT / "results/initial_design_generation/tier1_filtering/tier1_ranked_candidates_pre_tier2.csv",
        4,
    ),
]
PRODUCTION_POOL = ROOT / "results/initial_design_generation/production_initial_pool/production_initial_pool_candidates_all.csv"

READ_COLS = [
    "target",
    "variant_id",
    "sequence",
    "sequence_hash",
    "canonical_sequence_hash_full",
    "canonical_recovery_sequence_hash",
    "canonical_sequence_hash_short",
    "canonical_mutated_window_sequence",
    "mutation_list",
    "normalized_mutation_list",
    "mutation_count",
    "His_count",
    "his_seed_set",
    "rescue_mutation_list",
    "rescue_signature",
    "rescue_count",
    "primary_generation_route",
    "all_source_routes",
    "near_duplicate_cluster_id",
    "tier1_near_duplicate_cluster_id",
    "source_near_duplicate_cluster_id",
    "hard_filter_status",
    "forbidden_pair_status",
    "exact_duplicate_status",
    "buildability_light_status",
    "selection_eligibility",
    "tier2_eligibility_status",
    "tier1_review_class",
    "bank_eligibility",
    "stage2a_list_action",
    "t2_class_current",
    "tier2_class",
    "mpnn_total_score_per_residue",
    "mpnn_score_percentile_within_target",
    "mpnn_score_percentile_by_target",
    "neutral_retention_score",
    "acidic_release_support_score",
    "global_weakening_risk_score",
    "display_or_expression_risk_score",
    "glycan_or_epitope_risk_score",
    "hit_likelihood_score_v0",
    "tier1_rank_score",
]

AA_ALPHABET = tuple("ADEFGIKLMNPQRSTVWY")  # no Cys, no His

WINDOW_RULES = {
    "Ab_1E62": {
        "chain": "L",
        "window_start": 1,
        "window_end": 40,
        "hard_protect_positions": {23},
        "excluded_positions": {23},
        "local_expansion_allowed_statuses": {
            "active_primary_template",
            "active_primary_backfill_template",
        },
    },
    "Ab_sdAb": {
        "chain": "A",
        "window_start": 72,
        "window_end": 111,
        "hard_protect_positions": {96},
        "excluded_positions": {96, 102, 105},
        "local_expansion_allowed_statuses": {
            "secondary_template_complex_weak",
            "secondary_template_complex_unchecked",
        },
    },
}

TARGET_CONFIG = {
    "Ab_1E62": {
        "target_count": 20000,
        "seed_cap": 5000,
        "cluster_cap": 1000,
        "four_mut_cap": 12000,
        "seed_only_cap": 5000,
        "neutral_boundary_cap": 10000,
        "high_or_supported_min": 8000,
        "ay111h_containing_cap": None,
        "template_status_quota": {
            "active_primary_template": 10000,
            "active_primary_backfill_template": 12000,
            "limited_boundary_template": 1000,
        },
        "template_cap": {
            "active_primary_template": 1200,
            "active_primary_backfill_template": 1000,
            "limited_boundary_template": 200,
        },
    },
    "Ab_sdAb": {
        "target_count": 20000,
        "seed_cap": 6000,
        "cluster_cap": 1000,
        "four_mut_cap": 3000,
        "seed_only_cap": 6000,
        "neutral_boundary_cap": None,
        "high_or_supported_min": None,
        "ay111h_containing_cap": 7000,
        "template_status_quota": {
            "secondary_template_complex_weak": 9000,
            "secondary_template_complex_unchecked": 7500,
            "low_priority_secondary_template": 3500,
            "boundary_representative_only": 500,
        },
        "template_cap": {
            "secondary_template_complex_weak": 650,
            "secondary_template_complex_unchecked": 650,
            "low_priority_secondary_template": 500,
            "boundary_representative_only": 50,
        },
    },
}

MATCH_TIER_RANK = {
    "template_seed": 0,
    "A_strong": 1,
    "B_medium": 2,
    "local_expansion_A": 3,
    "local_expansion_B": 4,
    "C_seed_only": 5,
}

SOURCE_CLASS_SCORE = {
    "A_tier1_priority_candidate": 0.30,
    "B_rescue_enriched_candidate": 0.24,
    "E_sparse_mpnn_supplement_review": 0.10,
    "D_favorable_mpnn_but_weak_pH_mechanism": 0.02,
    "F_neutral_boundary_or_high_risk": -0.08,
}


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def short_hash(text: str, n: int = 10) -> str:
    return sha256_text(text)[:n]


def md_table(df: pd.DataFrame) -> str:
    return v1.md_table(df)


def as_bool(value: object) -> bool:
    return v1.as_bool(value)


def sha256_file(path: Path) -> str:
    return v1.sha256_file(path)


def parse_mutation(mut: object) -> tuple[str, str, int, str] | None:
    match = re.match(r"^([A-Z])([A-Z]?)(\d+)([A-Z])$", str(mut))
    if not match:
        return None
    chain, wt, pos, aa = match.groups()
    return chain, wt, int(pos), aa


def mutation_position(mut: object) -> str:
    parsed = parse_mutation(mut)
    if not parsed:
        return str(mut)
    chain, _wt, pos, _aa = parsed
    return f"{chain}{pos}"


def split_mutations(mutation_list: object) -> tuple[list[str], list[str]]:
    if not isinstance(mutation_list, str) or not mutation_list or mutation_list.lower() == "nan":
        return [], []
    muts = [m.strip() for m in mutation_list.split(";") if m.strip() and m.strip().lower() != "nan"]
    his = [m for m in muts if m.endswith("H")]
    non_his = [m for m in muts if not m.endswith("H")]
    return his, non_his


def mutation_count_class(n: int | float | object) -> str:
    try:
        value = int(float(n))
    except (TypeError, ValueError):
        return "unknown"
    if value <= 2:
        return "low_order"
    if value == 3:
        return "mid_order"
    return "four_mut"


def canonical_key(row: pd.Series | dict[str, Any]) -> str:
    target = str(row.get("target"))
    for col in ["canonical_sequence_hash_full", "canonical_recovery_sequence_hash", "sequence_hash"]:
        val = row.get(col)
        if isinstance(val, str) and val and val.lower() != "nan":
            return f"{target}|{val}"
    seq = row.get("sequence")
    if isinstance(seq, str) and seq and seq.lower() != "nan":
        return f"{target}|{sha256_text(seq)}"
    variant_id = row.get("variant_id")
    if isinstance(variant_id, str) and variant_id:
        return f"{target}|variant|{variant_id}"
    return f"{target}|mutation|{row.get('mutation_list', '')}"


def selected_seed_contains_ay111h(row: pd.Series | dict[str, Any]) -> bool:
    text = f"{row.get('his_seed_set', '')};{row.get('mutation_list', '')}"
    return "AY111H" in text


def selected_seed_contains_sdab_risk(row: pd.Series | dict[str, Any]) -> bool:
    if str(row.get("target")) != "Ab_sdAb":
        return False
    text = f"{row.get('his_seed_set', '')};{row.get('mutation_list', '')}"
    return "AG102H" in text or "AV105H" in text


def source_inventory() -> pd.DataFrame:
    rows = []
    for source_pool, path, _rank in SOURCE_FILES:
        if not path.exists():
            continue
        rows.append(
            {
                "source_pool": source_pool,
                "path": str(path.relative_to(ROOT)),
                "rows": sum(1 for _ in path.open()) - 1,
                "sha256": sha256_file(path),
                "used_for_selection": True,
            }
        )
    if PRODUCTION_POOL.exists():
        rows.append(
            {
                "source_pool": "production_initial_pool_270k_reference",
                "path": str(PRODUCTION_POOL.relative_to(ROOT)),
                "rows": sum(1 for _ in PRODUCTION_POOL.open()) - 1,
                "sha256": sha256_file(PRODUCTION_POOL),
                "used_for_selection": False,
            }
        )
    return pd.DataFrame(rows)


def read_source(path: Path, source_pool: str, source_rank: int) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    use = [c for c in READ_COLS if c in cols]
    df = pd.read_csv(path, usecols=use, low_memory=False)
    for col in READ_COLS:
        if col not in df.columns:
            df[col] = pd.NA
    cluster_cols = ["near_duplicate_cluster_id", "tier1_near_duplicate_cluster_id", "source_near_duplicate_cluster_id"]
    df["near_duplicate_cluster_id"] = df[cluster_cols].bfill(axis=1).iloc[:, 0]
    df["source_pool"] = source_pool
    df["source_rank"] = source_rank
    return df


def load_candidate_sources() -> pd.DataFrame:
    frames = [read_source(path, source_pool, rank) for source_pool, path, rank in SOURCE_FILES]
    frames = [frame for frame in frames if not frame.empty]
    df = pd.concat(frames, ignore_index=True)
    numeric_cols = [
        "mutation_count",
        "His_count",
        "mpnn_total_score_per_residue",
        "mpnn_score_percentile_within_target",
        "mpnn_score_percentile_by_target",
        "neutral_retention_score",
        "acidic_release_support_score",
        "global_weakening_risk_score",
        "display_or_expression_risk_score",
        "glycan_or_epitope_risk_score",
        "hit_likelihood_score_v0",
        "tier1_rank_score",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["mpnn_score_percentile_within_target"] = df["mpnn_score_percentile_within_target"].fillna(
        df["mpnn_score_percentile_by_target"]
    )
    df["canonical_unique_key"] = df.apply(canonical_key, axis=1)
    source_lists = (
        df.groupby("canonical_unique_key")["source_pool"]
        .agg(lambda s: ";".join(sorted(set(map(str, s)))))
        .reset_index(name="source_pool_list")
    )
    occurrence = df.groupby("canonical_unique_key").size().reset_index(name="source_occurrence_count")
    df["_source_quality_sort"] = -df["tier1_rank_score"].fillna(-999)
    df = df.sort_values(["source_rank", "_source_quality_sort", "variant_id"], ascending=[True, True, True])
    df = df.drop_duplicates("canonical_unique_key", keep="first")
    df = df.merge(source_lists, on="canonical_unique_key", how="left")
    df = df.merge(occurrence, on="canonical_unique_key", how="left")
    return df


def canonical_nxs_sites(seq: str) -> set[int]:
    sites = set()
    for idx in range(len(seq) - 2):
        if seq[idx] == "N" and seq[idx + 1] != "P" and seq[idx + 2] in {"S", "T"}:
            sites.add(idx + 1)
    return sites


def parent_sequences_from_sources(sources: pd.DataFrame) -> dict[str, str]:
    parents: dict[str, str] = {}
    for target, group in sources.dropna(subset=["sequence"]).groupby("target", dropna=False):
        for _, row in group.iterrows():
            seq = str(row["sequence"])
            muts = split_mutations(row.get("mutation_list"))[0] + split_mutations(row.get("mutation_list"))[1]
            if not seq or seq.lower() == "nan" or not muts:
                continue
            chars = list(seq)
            ok = True
            for mut in muts:
                parsed = parse_mutation(mut)
                if not parsed:
                    ok = False
                    break
                _chain, wt, pos, aa = parsed
                if not (1 <= pos <= len(chars)):
                    ok = False
                    break
                if chars[pos - 1] == aa:
                    chars[pos - 1] = wt
                elif chars[pos - 1] == wt:
                    pass
                else:
                    # Some imported rows may have inconsistent coordinates; try another.
                    ok = False
                    break
            if ok:
                parents[str(target)] = "".join(chars)
                break
    return parents


def apply_mutations(parent: str, mutations: list[str]) -> str | None:
    chars = list(parent)
    for mut in mutations:
        parsed = parse_mutation(mut)
        if not parsed:
            return None
        _chain, wt, pos, aa = parsed
        if not (1 <= pos <= len(chars)):
            return None
        if chars[pos - 1] != wt:
            return None
        chars[pos - 1] = aa
    return "".join(chars)


def source_rejection_counts(df: pd.DataFrame, reason: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    return (
        df.groupby(["target", "source_pool"], dropna=False)
        .size()
        .reset_index(name="count")
        .assign(template_status="prefilter", reject_reason=reason)
    )


def filter_candidates(df: pd.DataFrame, templates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    allowed_seeds = (
        templates[templates["allowed_backfill"].map(as_bool)]
        .groupby("target")["seed_pattern"]
        .apply(lambda s: set(map(str, s)))
        .to_dict()
    )
    mask = pd.Series(True, index=df.index)
    reasons: list[pd.DataFrame] = []
    allowed_pair = {f"{target}||{seed}" for target, seeds in allowed_seeds.items() for seed in seeds}
    pair = df["target"].astype(str) + "||" + df["his_seed_set"].astype(str)

    failed = ~pair.isin(allowed_pair)
    reasons.append(source_rejection_counts(df[mask & failed], "no_template_seed_match"))
    mask &= ~failed

    failed = df["hard_filter_status"].notna() & ~df["hard_filter_status"].fillna("pass").eq("pass")
    reasons.append(source_rejection_counts(df[mask & failed], "hard_filter_not_pass"))
    mask &= ~failed

    failed = df["forbidden_pair_status"].notna() & ~df["forbidden_pair_status"].fillna("pass").eq("pass")
    reasons.append(source_rejection_counts(df[mask & failed], "forbidden_pair_not_pass"))
    mask &= ~failed

    failed = df["buildability_light_status"].notna() & ~df["buildability_light_status"].fillna("pass").eq("pass")
    reasons.append(source_rejection_counts(df[mask & failed], "buildability_not_pass"))
    mask &= ~failed

    failed = df["mutation_count"].fillna(99).gt(4)
    reasons.append(source_rejection_counts(df[mask & failed], "mutation_count_gt_4"))
    mask &= ~failed

    failed = df["tier1_review_class"].fillna("").eq("G_audit_or_demoted")
    reasons.append(source_rejection_counts(df[mask & failed], "tier1_audit_or_demoted"))
    mask &= ~failed

    risk_text = df["his_seed_set"].fillna("").astype(str) + ";" + df["mutation_list"].fillna("").astype(str)
    failed = df["target"].eq("Ab_sdAb") & (
        risk_text.str.contains("AG102H", regex=False) | risk_text.str.contains("AV105H", regex=False)
    )
    reasons.append(source_rejection_counts(df[mask & failed], "sdAb_AG102H_or_AV105H_main_excluded"))
    mask &= ~failed

    out = df[mask].copy()
    rejection = pd.concat([item for item in reasons if not item.empty], ignore_index=True) if reasons else pd.DataFrame()
    if not rejection.empty:
        rejection = (
            rejection.groupby(["target", "source_pool", "template_status", "reject_reason"], dropna=False)["count"]
            .sum()
            .reset_index()
            .sort_values(["target", "source_pool", "count"], ascending=[True, True, False])
        )
    return out, rejection


def template_match_for_candidate(row: pd.Series, tpl: pd.Series) -> dict[str, Any]:
    _his, cand_non_his = split_mutations(row.get("mutation_list"))
    _tpl_his, tpl_non_his = split_mutations(tpl.get("mutation_pattern"))
    cand_non_his_set = set(cand_non_his)
    cand_non_his_pos = {mutation_position(m) for m in cand_non_his}
    tpl_non_his_pos = {mutation_position(m) for m in tpl_non_his}
    rescue_pattern = str(tpl.get("rescue_pattern", "none"))
    rescue_pos = mutation_position(rescue_pattern) if rescue_pattern and rescue_pattern.lower() != "nan" else "none"

    exact_rescue = rescue_pattern in cand_non_his_set
    same_rescue_position = rescue_pos != "none" and rescue_pos in cand_non_his_pos
    overlap_count = len(cand_non_his_pos & tpl_non_his_pos)
    tpl_count_class = mutation_count_class(len(split_mutations(tpl.get("mutation_pattern"))[0]) + len(tpl_non_his))
    cand_count_class = mutation_count_class(row.get("mutation_count"))
    route_text = f"{row.get('primary_generation_route', '')};{row.get('all_source_routes', '')}"
    route_compatible = (
        "ProteinMPNN_seeded_rescue" in route_text
        or "His_plus_rescue" in route_text
        or "wetlab_informed" in route_text
        or "structure_or_interface" in route_text
        or pd.isna(row.get("primary_generation_route"))
    )

    if exact_rescue and (tpl_count_class == cand_count_class or overlap_count >= 1) and route_compatible:
        match_tier = "A_strong"
        match_level = "exact_rescue"
        match_score = 4
    elif exact_rescue or same_rescue_position or overlap_count >= 1:
        match_tier = "B_medium"
        match_level = "same_rescue_position" if not exact_rescue else "exact_rescue"
        match_score = 3 if exact_rescue else 2
    else:
        match_tier = "C_seed_only"
        match_level = "seed_only"
        match_score = 0
    return {
        "template_id": tpl["template_id"],
        "template_status": tpl["template_status"],
        "recommended_expansion_role": tpl["recommended_expansion_role"],
        "complex_diagnostic_required": as_bool(tpl["complex_diagnostic_required"]),
        "template_priority": int(tpl["priority_level"]),
        "template_match_tier": match_tier,
        "template_match_level": match_level,
        "template_match_score": match_score,
        "same_position_overlap_count": overlap_count,
        "evidence_inheritance_reason": (
            f"{match_tier};exact_rescue={exact_rescue};same_rescue_position={same_rescue_position};"
            f"position_overlap={overlap_count};route_compatible={route_compatible}"
        ),
    }


def template_options(row: pd.Series, by_target_seed: dict[tuple[str, str], pd.DataFrame]) -> list[dict[str, Any]]:
    key = (str(row["target"]), str(row.get("his_seed_set")))
    if key not in by_target_seed:
        return []
    options = [template_match_for_candidate(row, tpl) for _, tpl in by_target_seed[key].iterrows()]
    return sorted(
        options,
        key=lambda x: (
            MATCH_TIER_RANK[x["template_match_tier"]],
            -int(x["template_match_score"]),
            int(x["template_priority"]),
            str(x["template_id"]),
        ),
    )


def quality_tier(row: pd.Series | dict[str, Any]) -> str:
    cls = str(row.get("tier1_review_class", ""))
    bank = str(row.get("bank_eligibility", ""))
    t2 = str(row.get("t2_class_current", row.get("tier2_class", "")))
    if bank == "pass_like" or cls == "A_tier1_priority_candidate" or t2 in {"T2_strong_candidate", "T2_pass_like_structure"}:
        return "high_support"
    if cls in {"B_rescue_enriched_candidate", "E_sparse_mpnn_supplement_review"} or t2 == "T2_good_candidate":
        return "supported"
    if cls == "D_favorable_mpnn_but_weak_pH_mechanism":
        return "mpnn_favorable_pH_weak"
    if cls == "F_neutral_boundary_or_high_risk":
        return "neutral_boundary_or_high_risk"
    if cls == "local_expansion_constrained":
        return "supported"
    return "proxy_support"


def candidate_score(row: pd.Series, options: list[dict[str, Any]]) -> float:
    best = options[0] if options else {}
    score = 0.0
    score += SOURCE_CLASS_SCORE.get(str(row.get("tier1_review_class", "")), 0.0)
    score += max(0, 6 - MATCH_TIER_RANK.get(str(best.get("template_match_tier")), 6)) * 0.16
    if row.get("bank_eligibility") == "pass_like":
        score += 0.30
    elif row.get("bank_eligibility") == "supported_boundary":
        score += 0.08
    t2 = str(row.get("t2_class_current", row.get("tier2_class", "")))
    if t2 in {"T2_strong_candidate", "T2_pass_like_structure"}:
        score += 0.25
    elif t2 == "T2_good_candidate":
        score += 0.15
    elif t2 == "T2_release_possible_but_neutral_risky":
        score -= 0.10
    if row.get("stage2a_list_action") == "stage2a_include_priority":
        score += 0.12
    for col, weight, direction in [
        ("neutral_retention_score", 0.14, 1),
        ("acidic_release_support_score", 0.18, 1),
        ("global_weakening_risk_score", 0.12, -1),
        ("display_or_expression_risk_score", 0.06, -1),
        ("glycan_or_epitope_risk_score", 0.05, -1),
        ("tier1_rank_score", 0.12, 1),
    ]:
        val = row.get(col)
        if pd.notna(val):
            score += float(val) * weight * direction
    pct = row.get("mpnn_score_percentile_within_target")
    if pd.notna(pct):
        score += max(0.0, 1.0 - float(pct)) * 0.18
    mut_count = row.get("mutation_count")
    if pd.notna(mut_count):
        if int(mut_count) <= 2:
            score += 0.05
        elif int(mut_count) == 4:
            score -= 0.04
    score += max(0, 5 - int(row.get("source_rank", 5))) * 0.02
    return score


def build_backfill_candidates(candidates: pd.DataFrame, templates: pd.DataFrame) -> pd.DataFrame:
    backfill_templates = templates[templates["allowed_backfill"].map(as_bool)].copy()
    by_target_seed = {
        key: group.copy()
        for key, group in backfill_templates.groupby(["target", "seed_pattern"], dropna=False)
    }
    rows: list[dict[str, Any]] = []
    for _, row in candidates.iterrows():
        opts = template_options(row, by_target_seed)
        if not opts:
            continue
        item = row.to_dict()
        item["_template_options"] = opts
        item["candidate_quality_tier"] = quality_tier(row)
        item["expanded_pool_score"] = candidate_score(row, opts)
        item["source_variant_id"] = row["variant_id"]
        rows.append(item)
    return pd.DataFrame(rows)


def local_expansion_position_order(target: str, seed_positions: set[int], base_positions: set[int]) -> list[int]:
    rule = WINDOW_RULES[target]
    positions = []
    for pos in range(rule["window_start"], rule["window_end"] + 1):
        if pos in rule["excluded_positions"] or pos in seed_positions or pos in base_positions:
            continue
        positions.append(pos)
    # Deterministic center-out order around known windows.
    center = 30 if target == "Ab_1E62" else 100
    return sorted(positions, key=lambda p: (abs(p - center), p))


def mutation_list_string(muts: list[str]) -> str:
    return ";".join(sorted(muts, key=lambda m: (parse_mutation(m)[2] if parse_mutation(m) else 999, m)))


def generate_local_expansion_candidates(templates: pd.DataFrame, sources: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, str, str]:
    parents = parent_sequences_from_sources(sources)
    parent_nxs = {target: canonical_nxs_sites(seq) for target, seq in parents.items()}
    rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    config_lines = [
        "stage: backfill_v2_local_expansion",
        "mode: constrained_template_local_expansion",
        "no_heavy_compute: true",
        "alphabet: ADEFGIKLMNPQRSTVWY",
        "hard_rules:",
        "  new_cys: forbidden",
        "  new_canonical_nxs_t: forbidden",
        "  mutation_count_max: 4",
        "  sdab_AG102H_AV105H_main: forbidden",
        "  sdab_complex_diagnostic_required: true",
    ]
    plan_lines = [
        "# Local Expansion Plan",
        "",
        "Local expansion is used only to fill Backfill v2 target counts after existing source-pool backfill is insufficient.",
        "It keeps the His seed fixed, stays inside the selected 40 aa window, blocks hard-protect positions, and does not run heavy compute.",
        "",
    ]
    seen_hashes: set[str] = set()

    for _, tpl in templates.sort_values(["target", "priority_level", "template_id"]).iterrows():
        target = str(tpl["target"])
        status = str(tpl["template_status"])
        if target not in WINDOW_RULES or status not in WINDOW_RULES[target]["local_expansion_allowed_statuses"]:
            continue
        if target not in parents:
            audit_rows.append(
                {
                    "target": target,
                    "template_id": tpl["template_id"],
                    "template_status": status,
                    "generated": 0,
                    "reason": "missing_parent_sequence",
                }
            )
            continue
        parent = parents[target]
        _his, base_non_his = split_mutations(tpl["mutation_pattern"])
        seed_his, _seed_non = split_mutations(tpl["seed_pattern"])
        base_muts = seed_his + base_non_his
        seed_positions = {parse_mutation(m)[2] for m in seed_his if parse_mutation(m)}
        base_positions = {parse_mutation(m)[2] for m in base_non_his if parse_mutation(m)}
        if not seed_his:
            continue
        max_per_template = 650 if target == "Ab_sdAb" else 1200
        if status in {"active_primary_template", "secondary_template_complex_weak", "secondary_template_complex_unchecked"}:
            max_per_template = 750 if target == "Ab_sdAb" else 1200
        positions = local_expansion_position_order(target, seed_positions, base_positions)
        generated = 0

        def add_candidate(muts: list[str], tier: str, reason: str) -> None:
            nonlocal generated
            if generated >= max_per_template:
                return
            mut_str = mutation_list_string(muts)
            if target == "Ab_sdAb" and ("AG102H" in mut_str or "AV105H" in mut_str):
                return
            if any(m.endswith("C") for m in muts):
                return
            seq = apply_mutations(parent, mut_str.split(";"))
            if not seq:
                return
            if canonical_nxs_sites(seq) - parent_nxs[target]:
                return
            seq_hash_full = sha256_text(seq)
            key = f"{target}|{seq_hash_full}"
            if key in seen_hashes:
                return
            seen_hashes.add(key)
            muts_parsed = [parse_mutation(m) for m in mut_str.split(";")]
            non_his_positions = [p[2] for p in muts_parsed if p and p[3] != "H"]
            cluster_sig = ";".join([str(p) for p in sorted(non_his_positions)]) or "seed_only"
            rows.append(
                {
                    "target": target,
                    "variant_id": f"{target.replace('Ab_', '')}_local_v2_{short_hash(target + seq_hash_full)}",
                    "source_variant_id": pd.NA,
                    "sequence": seq,
                    "sequence_hash": seq_hash_full[:10],
                    "canonical_sequence_hash_full": seq_hash_full,
                    "canonical_recovery_sequence_hash": pd.NA,
                    "canonical_unique_key": key,
                    "mutation_list": mut_str,
                    "normalized_mutation_list": mut_str,
                    "mutation_count": len(mut_str.split(";")),
                    "His_count": sum(1 for m in mut_str.split(";") if m.endswith("H")),
                    "his_seed_set": tpl["seed_pattern"],
                    "near_duplicate_cluster_id": f"{target}_lex_{short_hash(str(tpl['seed_pattern']) + cluster_sig, 12)}",
                    "assigned_template_id": tpl["template_id"],
                    "assigned_template_status": status,
                    "assigned_expansion_role": tpl["recommended_expansion_role"],
                    "template_match_tier": tier,
                    "template_match_level": "local_expansion",
                    "expanded_pool_role": "constrained_local_expansion",
                    "candidate_quality_tier": "supported",
                    "expanded_pool_score": (
                        (1.95 if tier == "local_expansion_A" else 1.65)
                        + (0.45 if len(mut_str.split(";")) <= 3 else -0.25)
                    ),
                    "source_pool": "constrained_local_expansion_v2",
                    "source_pool_list": "constrained_local_expansion_v2",
                    "source_occurrence_count": 1,
                    "hard_filter_status": "pass",
                    "forbidden_pair_status": "pass",
                    "buildability_light_status": "pass",
                    "tier1_review_class": "local_expansion_constrained",
                    "bank_eligibility": pd.NA,
                    "stage2a_list_action": pd.NA,
                    "t2_class_current": pd.NA,
                    "tier2_class": pd.NA,
                    "selection_eligibility": "local_expansion_candidate",
                    "tier2_eligibility_status": pd.NA,
                    "mpnn_total_score_per_residue": pd.NA,
                    "mpnn_score_percentile_within_target": pd.NA,
                    "neutral_retention_score": pd.NA,
                    "acidic_release_support_score": pd.NA,
                    "global_weakening_risk_score": pd.NA,
                    "display_or_expression_risk_score": pd.NA,
                    "glycan_or_epitope_risk_score": pd.NA,
                    "hit_likelihood_score_v0": pd.NA,
                    "tier1_rank_score": pd.NA,
                    "complex_diagnostic_required": bool(tpl["complex_diagnostic_required"]) or target == "Ab_sdAb",
                    "same_position_overlap_count": len(base_positions),
                    "evidence_inheritance_reason": reason,
                }
            )
            generated += 1

        # 1E62 needs enough <=3 mutation candidates to keep the 4-mut fraction
        # below the v2 cap. This is a bounded single-rescue-position expansion:
        # fixed His seed plus one non-His local rescue-style mutation.
        if target == "Ab_1E62":
            low_order_generated = 0
            for pos in positions[:24]:
                for aa in AA_ALPHABET:
                    wt = parent[pos - 1]
                    if aa == wt:
                        continue
                    add_candidate(
                        seed_his + [f"{WINDOW_RULES[target]['chain']}{wt}{pos}{aa}"],
                        "local_expansion_B",
                        "local_expansion_B;single_non_his_rescue_position_low_order;keeps_His_seed_fixed",
                    )
                    low_order_generated += 1
                    if generated >= max_per_template or low_order_generated >= 300:
                        break
                if generated >= max_per_template or low_order_generated >= 300:
                    break

        # A-tier local expansion: keep exact template rescue mutations and add
        # one or two bounded non-His local mutations.
        for extra_count in [1, 2]:
            if len(base_muts) + extra_count > 4:
                continue
            for pos_combo in itertools.combinations(positions[:18], extra_count):
                for aa_combo in itertools.product(AA_ALPHABET, repeat=extra_count):
                    extra = []
                    skip = False
                    for pos, aa in zip(pos_combo, aa_combo):
                        wt = parent[pos - 1]
                        if aa == wt:
                            skip = True
                            break
                        extra.append(f"{WINDOW_RULES[target]['chain']}{wt}{pos}{aa}")
                    if skip:
                        continue
                    add_candidate(
                        base_muts + extra,
                        "local_expansion_A",
                        "local_expansion_A;keeps_exact_template_rescue;adds_bounded_non_his_context",
                    )
                    if generated >= max_per_template:
                        break
                if generated >= max_per_template:
                    break
            if generated >= max_per_template:
                break

        # B-tier local expansion: same rescue position, alternate rescue AA,
        # with one additional bounded non-His mutation if capacity remains.
        if generated < max_per_template and base_non_his:
            rescue = base_non_his[0]
            parsed = parse_mutation(rescue)
            if parsed:
                chain, wt, pos, original_aa = parsed
                for alt in AA_ALPHABET:
                    if alt in {wt, original_aa}:
                        continue
                    replacement = f"{chain}{wt}{pos}{alt}"
                    if len(seed_his) + 2 <= 4:
                        for extra_pos in positions[:18]:
                            if extra_pos == pos:
                                continue
                            for aa in AA_ALPHABET:
                                wt2 = parent[extra_pos - 1]
                                if aa == wt2:
                                    continue
                                add_candidate(
                                    seed_his + [replacement, f"{WINDOW_RULES[target]['chain']}{wt2}{extra_pos}{aa}"],
                                    "local_expansion_B",
                                    "local_expansion_B;same_rescue_position_alternate_rescue;adds_bounded_non_his_context",
                                )
                                if generated >= max_per_template:
                                    break
                            if generated >= max_per_template:
                                break
                    if generated >= max_per_template:
                        break

        audit_rows.append(
            {
                "target": target,
                "template_id": tpl["template_id"],
                "template_status": status,
                "generated": generated,
                "reason": "generated" if generated else "no_valid_local_combinations",
            }
        )

    # Controlled 1E62 low-frequency His-seed expansion. The strict template
    # seed universe has too little <=3-mut capacity to meet the 20K target
    # while respecting the 4-mut cap. This branch is bounded to known His
    # control/anchor positions inside VL:1-40 and remains template-bank
    # secondary/backfill evidence, not a final-selection claim.
    target = "Ab_1E62"
    if target in parents:
        parent = parents[target]
        known_his_positions = [24, 26, 28, 31, 33, 34, 35, 38]
        template_seed_set = set(templates[templates["target"].eq(target)]["seed_pattern"].astype(str))
        lowfreq_seed_pairs = []
        for p1, p2 in itertools.combinations(known_his_positions, 2):
            seed = mutation_list_string(
                [
                    f"L{parent[p1 - 1]}{p1}H",
                    f"L{parent[p2 - 1]}{p2}H",
                ]
            )
            if seed not in template_seed_set:
                lowfreq_seed_pairs.append((seed, {p1, p2}))
        generated_total = 0
        for seed, seed_positions in lowfreq_seed_pairs:
            template_id = f"tpl_1E62_lowfreq_seed_{short_hash(seed, 10)}"
            generated = 0
            positions = local_expansion_position_order(target, seed_positions, set())
            for pos in positions[:24]:
                for aa in AA_ALPHABET:
                    wt = parent[pos - 1]
                    if aa == wt:
                        continue
                    muts = seed.split(";") + [f"L{wt}{pos}{aa}"]
                    mut_str = mutation_list_string(muts)
                    seq = apply_mutations(parent, mut_str.split(";"))
                    if not seq or canonical_nxs_sites(seq) - parent_nxs[target]:
                        continue
                    seq_hash_full = sha256_text(seq)
                    key = f"{target}|{seq_hash_full}"
                    if key in seen_hashes:
                        continue
                    seen_hashes.add(key)
                    cluster_sig = ";".join(str(parse_mutation(m)[2]) for m in mut_str.split(";") if parse_mutation(m))
                    rows.append(
                        {
                            "target": target,
                            "variant_id": f"1E62_lowfreq_local_v2_{short_hash(target + seq_hash_full)}",
                            "source_variant_id": pd.NA,
                            "sequence": seq,
                            "sequence_hash": seq_hash_full[:10],
                            "canonical_sequence_hash_full": seq_hash_full,
                            "canonical_recovery_sequence_hash": pd.NA,
                            "canonical_unique_key": key,
                            "mutation_list": mut_str,
                            "normalized_mutation_list": mut_str,
                            "mutation_count": len(mut_str.split(";")),
                            "His_count": sum(1 for m in mut_str.split(";") if m.endswith("H")),
                            "his_seed_set": seed,
                            "near_duplicate_cluster_id": f"Ab_1E62_lfh_{short_hash(seed + cluster_sig, 12)}",
                            "assigned_template_id": template_id,
                            "assigned_template_status": "active_primary_backfill_template",
                            "assigned_expansion_role": "low_frequency_seed_diversity_template",
                            "template_match_tier": "local_expansion_B",
                            "template_match_level": "low_frequency_seed_local_expansion",
                            "expanded_pool_role": "constrained_local_expansion",
                            "candidate_quality_tier": "supported",
                            "expanded_pool_score": 2.05,
                            "source_pool": "constrained_low_frequency_seed_local_expansion_v2",
                            "source_pool_list": "constrained_low_frequency_seed_local_expansion_v2",
                            "source_occurrence_count": 1,
                            "hard_filter_status": "pass",
                            "forbidden_pair_status": "pass",
                            "buildability_light_status": "pass",
                            "tier1_review_class": "local_expansion_constrained",
                            "bank_eligibility": pd.NA,
                            "stage2a_list_action": pd.NA,
                            "t2_class_current": pd.NA,
                            "tier2_class": pd.NA,
                            "selection_eligibility": "local_expansion_candidate",
                            "tier2_eligibility_status": pd.NA,
                            "mpnn_total_score_per_residue": pd.NA,
                            "mpnn_score_percentile_within_target": pd.NA,
                            "neutral_retention_score": pd.NA,
                            "acidic_release_support_score": pd.NA,
                            "global_weakening_risk_score": pd.NA,
                            "display_or_expression_risk_score": pd.NA,
                            "glycan_or_epitope_risk_score": pd.NA,
                            "hit_likelihood_score_v0": pd.NA,
                            "tier1_rank_score": pd.NA,
                            "complex_diagnostic_required": False,
                            "same_position_overlap_count": 0,
                            "evidence_inheritance_reason": "local_expansion_B;low_frequency_His_seed_control_anchor_derived;single_non_his_rescue_position_low_order",
                        }
                    )
                    generated += 1
                    generated_total += 1
                    if generated >= 240:
                        break
                if generated >= 240:
                    break
            audit_rows.append(
                {
                    "target": target,
                    "template_id": template_id,
                    "template_status": "active_primary_backfill_template",
                    "generated": generated,
                    "reason": "low_frequency_His_seed_control_anchor_derived",
                }
            )
            if generated_total >= 6000:
                break

    plan_lines += [
        "## Local Expansion Scope",
        "",
        "- 1E62: enabled for active primary and active primary backfill templates.",
        "- 1E62: additionally enables a bounded low-frequency His-seed branch from known VL:1-40 control/anchor His positions when strict template seeds cannot supply enough <=3-mut capacity.",
        "- sdAb: enabled only for complex-weak / complex-unchecked secondary templates and every row remains `complex_diagnostic_required`.",
        "- Boundary templates are not used as local expansion drivers.",
        "",
    ]
    return pd.DataFrame(rows), pd.DataFrame(audit_rows), "\n".join(plan_lines) + "\n", "\n".join(config_lines) + "\n"


def bank_seed_rows(bank: pd.DataFrame, templates: pd.DataFrame, parents: dict[str, str]) -> pd.DataFrame:
    tpl = templates[["template_id", "variant_id"]]
    out = bank.merge(tpl, on="variant_id", how="left", validate="one_to_one")
    seqs = []
    keys = []
    sequence_hashes = []
    for _, row in out.iterrows():
        seq = pd.NA
        seq_hash = pd.NA
        key = f"{row['target']}|bank|{row['variant_id']}"
        parent = parents.get(str(row["target"]))
        if parent:
            muts = [m for m in str(row["mutation_list"]).split(";") if m and m.lower() != "nan"]
            seq_try = apply_mutations(parent, muts)
            if seq_try:
                seq = seq_try
                seq_hash = sha256_text(seq_try)
                key = f"{row['target']}|{seq_hash}"
        seqs.append(seq)
        sequence_hashes.append(seq_hash)
        keys.append(key)
    out["sequence"] = seqs
    out["canonical_sequence_hash_full"] = sequence_hashes
    out["sequence_hash"] = [s[:10] if isinstance(s, str) else pd.NA for s in sequence_hashes]
    out["canonical_unique_key"] = keys
    out["source_pool"] = "tier2_candidate_bank"
    out["source_pool_list"] = "tier2_candidate_bank"
    out["source_occurrence_count"] = 1
    out["expanded_pool_role"] = "bank_template_seed"
    out["template_match_tier"] = "template_seed"
    out["template_match_level"] = "template_seed"
    out["candidate_quality_tier"] = "bank_reviewed"
    out["expanded_pool_score"] = 10.0 - out["bank_priority"].fillna(9) * 0.01
    out["source_variant_id"] = out["variant_id"]
    out["same_position_overlap_count"] = pd.NA
    out["evidence_inheritance_reason"] = "frozen_tier2_candidate_bank_template_seed"
    return out


def prepare_selection_candidates(existing_backfill: pd.DataFrame, local_expansion: pd.DataFrame) -> pd.DataFrame:
    if existing_backfill.empty:
        existing = pd.DataFrame()
    else:
        existing = existing_backfill.copy()
        for col in ["assigned_template_id", "assigned_template_status", "assigned_expansion_role"]:
            existing[col] = pd.NA
        existing["template_match_tier"] = existing["_template_options"].map(lambda opts: opts[0]["template_match_tier"])
        existing["template_match_level"] = existing["_template_options"].map(lambda opts: opts[0]["template_match_level"])
        existing["same_position_overlap_count"] = existing["_template_options"].map(lambda opts: opts[0]["same_position_overlap_count"])
        existing["evidence_inheritance_reason"] = existing["_template_options"].map(lambda opts: opts[0]["evidence_inheritance_reason"])
        existing["expanded_pool_role"] = "evidence_guided_backfill"
    if local_expansion.empty:
        local = pd.DataFrame()
    else:
        local = local_expansion.copy()
        local["_template_options"] = local.apply(
            lambda row: [
                {
                    "template_id": row["assigned_template_id"],
                    "template_status": row["assigned_template_status"],
                    "recommended_expansion_role": row["assigned_expansion_role"],
                    "complex_diagnostic_required": bool(row["complex_diagnostic_required"]),
                    "template_priority": 99,
                    "template_match_tier": row["template_match_tier"],
                    "template_match_level": row["template_match_level"],
                    "template_match_score": 3 if row["template_match_tier"] == "local_expansion_A" else 2,
                    "same_position_overlap_count": row["same_position_overlap_count"],
                    "evidence_inheritance_reason": row["evidence_inheritance_reason"],
                }
            ],
            axis=1,
        )
    if existing.empty:
        return local
    if local.empty:
        return existing
    all_cols = sorted(set(existing.columns) | set(local.columns))
    for df in [existing, local]:
        for col in all_cols:
            if col not in df.columns:
                df[col] = pd.NA
    return pd.concat([existing[all_cols], local[all_cols]], ignore_index=True)


def select_pool(bank_rows: pd.DataFrame, candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected: list[dict[str, Any]] = []
    selected_keys: set[str] = set()
    selected_variant_keys: set[tuple[str, str]] = set()
    counters = {
        "seed": defaultdict(Counter),
        "cluster": defaultdict(Counter),
        "template_status": defaultdict(Counter),
        "template": defaultdict(Counter),
        "four_mut": Counter(),
        "ay111h": Counter(),
        "seed_only": Counter(),
        "neutral": Counter(),
    }
    target_counts: Counter[str] = Counter()
    reject_counts: Counter[tuple[str, str, str, str]] = Counter()

    def can_add(item: dict[str, Any], option: dict[str, Any]) -> tuple[bool, str]:
        target = str(item["target"])
        cfg = TARGET_CONFIG[target]
        seed = str(item.get("his_seed_set", ""))
        cluster = str(item.get("near_duplicate_cluster_id", "missing_cluster"))
        status = str(option["template_status"])
        template_id = str(option["template_id"])
        match_tier = str(option.get("template_match_tier", item.get("template_match_tier", "C_seed_only")))
        mut_count = int(float(item.get("mutation_count", 99))) if pd.notna(item.get("mutation_count")) else 99
        quality = str(item.get("candidate_quality_tier", ""))

        if target_counts[target] >= cfg["target_count"]:
            return False, "target_count_full"
        if selected_seed_contains_sdab_risk(item):
            return False, "sdAb_AG102H_or_AV105H_main_excluded"
        if counters["seed"][target][seed] >= cfg["seed_cap"]:
            return False, "seed_cap_full"
        if counters["cluster"][target][cluster] >= cfg["cluster_cap"]:
            return False, "cluster_cap_full"
        if status not in cfg["template_status_quota"]:
            return False, "template_status_not_allowed"
        if counters["template_status"][target][status] >= cfg["template_status_quota"][status]:
            return False, "template_status_quota_full"
        if counters["template"][target][template_id] >= cfg["template_cap"].get(status, 100):
            return False, "template_cap_full"
        if mut_count == 4 and counters["four_mut"][target] >= cfg["four_mut_cap"]:
            return False, "four_mut_cap_full"
        if match_tier == "C_seed_only" and counters["seed_only"][target] >= cfg["seed_only_cap"]:
            return False, "seed_only_cap_full"
        if cfg["neutral_boundary_cap"] is not None and quality == "neutral_boundary_or_high_risk":
            if counters["neutral"][target] >= cfg["neutral_boundary_cap"]:
                return False, "neutral_boundary_cap_full"
        if cfg["ay111h_containing_cap"] is not None and selected_seed_contains_ay111h(item):
            if counters["ay111h"][target] >= cfg["ay111h_containing_cap"]:
                return False, "AY111H_containing_cap_full"
        return True, "selected"

    def add(item: dict[str, Any], option: dict[str, Any], role: str) -> None:
        target = str(item["target"])
        seed = str(item.get("his_seed_set", ""))
        cluster = str(item.get("near_duplicate_cluster_id", "missing_cluster"))
        mut_count = int(float(item.get("mutation_count", 99))) if pd.notna(item.get("mutation_count")) else 99
        match_tier = str(option.get("template_match_tier", item.get("template_match_tier", "C_seed_only")))
        quality = str(item.get("candidate_quality_tier", ""))
        item = dict(item)
        item["assigned_template_id"] = option["template_id"]
        item["assigned_template_status"] = option["template_status"]
        item["assigned_expansion_role"] = option["recommended_expansion_role"]
        item["template_match_tier"] = match_tier
        item["template_match_level"] = option.get("template_match_level", item.get("template_match_level"))
        item["same_position_overlap_count"] = option.get("same_position_overlap_count", item.get("same_position_overlap_count"))
        item["evidence_inheritance_reason"] = option.get(
            "evidence_inheritance_reason", item.get("evidence_inheritance_reason")
        )
        item["complex_diagnostic_required"] = option.get(
            "complex_diagnostic_required", item.get("complex_diagnostic_required", False)
        )
        item["expanded_pool_role"] = role if role != "candidate_default" else item.get("expanded_pool_role", "evidence_guided_backfill")
        selected.append(item)
        target_counts[target] += 1
        selected_keys.add(str(item["canonical_unique_key"]))
        selected_variant_keys.add((target, str(item.get("variant_id"))))
        counters["seed"][target][seed] += 1
        counters["cluster"][target][cluster] += 1
        counters["template_status"][target][str(option["template_status"])] += 1
        counters["template"][target][str(option["template_id"])] += 1
        if mut_count == 4:
            counters["four_mut"][target] += 1
        if selected_seed_contains_ay111h(item):
            counters["ay111h"][target] += 1
        if match_tier == "C_seed_only":
            counters["seed_only"][target] += 1
        if quality == "neutral_boundary_or_high_risk":
            counters["neutral"][target] += 1

    for _, row in bank_rows.sort_values(["target", "bank_priority", "variant_id"]).iterrows():
        option = {
            "template_id": row["template_id"],
            "template_status": row["template_status"],
            "recommended_expansion_role": row["recommended_expansion_role"],
            "complex_diagnostic_required": as_bool(row.get("complex_diagnostic_required", False)),
            "template_match_tier": "template_seed",
            "template_match_level": "template_seed",
            "same_position_overlap_count": pd.NA,
            "evidence_inheritance_reason": "frozen_tier2_candidate_bank_template_seed",
        }
        ok, reason = can_add(row.to_dict(), option)
        if ok:
            add(row.to_dict(), option, "bank_template_seed")
        else:
            reject_counts[(row["target"], "tier2_candidate_bank", str(row["template_status"]), reason)] += 1

    candidates = candidates.copy()
    candidates["match_tier_order"] = candidates["template_match_tier"].map(MATCH_TIER_RANK).fillna(9).astype(int)
    candidates["_four_mut_sort"] = pd.to_numeric(candidates["mutation_count"], errors="coerce").fillna(99).eq(4).astype(int)
    candidates = candidates.sort_values(
        ["target", "_four_mut_sort", "match_tier_order", "expanded_pool_score", "source_rank", "variant_id"],
        ascending=[True, True, True, False, True, True],
    )
    for _, row in candidates.iterrows():
        target = str(row["target"])
        if target_counts[target] >= TARGET_CONFIG[target]["target_count"]:
            reject_counts[(target, str(row["source_pool"]), "any", "target_count_full")] += 1
            continue
        if str(row["canonical_unique_key"]) in selected_keys or (target, str(row.get("variant_id"))) in selected_variant_keys:
            reject_counts[(target, str(row["source_pool"]), "any", "duplicate_selected")] += 1
            continue
        option_reasons = []
        for option in row["_template_options"]:
            ok, reason = can_add(row.to_dict(), option)
            if ok:
                add(row.to_dict(), option, "candidate_default")
                break
            option_reasons.append((option, reason))
        else:
            if option_reasons:
                option, reason = option_reasons[0]
                reject_counts[(target, str(row["source_pool"]), str(option["template_status"]), reason)] += 1
            else:
                reject_counts[(target, str(row["source_pool"]), "none", "no_template_option")] += 1

    selected_df = pd.DataFrame(selected)
    selected_df = selected_df.drop(columns=["_four_mut_sort"], errors="ignore")
    rejection = pd.DataFrame(
        [
            {
                "target": target,
                "source_pool": source_pool,
                "template_status": template_status,
                "reject_reason": reason,
                "count": count,
            }
            for (target, source_pool, template_status, reason), count in reject_counts.items()
        ]
    )
    if not rejection.empty:
        rejection = rejection.sort_values(["target", "count"], ascending=[True, False])
    return selected_df, rejection


def rule_text() -> str:
    return """stage: backfill_v2_20k_candidate_pool
no_heavy_compute: true
target_counts:
  Ab_1E62: 20000
  Ab_sdAb: 20000
match_tiers:
  A_strong:
    description: same His seed plus exact template rescue or equivalent position-set evidence
  B_medium:
    description: same His seed plus partial rescue or same-position evidence
  C_seed_only:
    description: same His seed only; capped diversity support
caps:
  Ab_1E62:
    seed_only_fraction_max: 0.25
    top_seed_fraction_max: 0.25
    top_cluster_fraction_max: 0.05
    four_mut_fraction_max: 0.60
    neutral_boundary_or_high_risk_fraction_max: 0.50
  Ab_sdAb:
    seed_only_fraction_max: 0.30
    top_seed_fraction_max: 0.30
    top_cluster_fraction_max: 0.05
    four_mut_fraction_max: 0.15
    AY111H_containing_fraction_max: 0.35
    AG102H_AV105H_main_candidates: 0
locks:
  final_10k_or_15k: true
  broad_tier2_heavy: true
  new_af3_simplefold_pyrosetta_foldx: true
  md: true
  glycan_low_risk_claim: true
  sdab_structure_confirmed_upgrade: true
local_expansion:
  Ab_1E62:
    low_frequency_His_seed_control_anchor_branch: enabled_when_strict_template_seed_low_order_capacity_is_insufficient
    allowed_His_positions: [24, 26, 28, 31, 33, 34, 35, 38]
    mutation_count: 3
    status: active_primary_backfill_template
    final_selection_claim: false
"""


def write_outputs(
    selected: pd.DataFrame,
    existing_backfill: pd.DataFrame,
    local_expansion: pd.DataFrame,
    local_audit: pd.DataFrame,
    rejection: pd.DataFrame,
    source_inv: pd.DataFrame,
    local_plan: str,
    local_config: str,
) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "backfill_v2_template_match_rules.yaml").write_text(rule_text(), encoding="utf-8")
    (OUT_DIR / "local_expansion_plan.md").write_text(local_plan, encoding="utf-8")
    (OUT_DIR / "local_expansion_config.yaml").write_text(local_config, encoding="utf-8")
    local_expansion.drop(columns=["_template_options"], errors="ignore").to_csv(
        OUT_DIR / "local_expansion_candidates.csv", index=False
    )
    local_audit.to_csv(OUT_DIR / "local_expansion_audit_report.csv", index=False)

    keep = [
        "target",
        "variant_id",
        "source_variant_id",
        "sequence",
        "sequence_hash",
        "canonical_sequence_hash_full",
        "canonical_recovery_sequence_hash",
        "canonical_unique_key",
        "mutation_list",
        "normalized_mutation_list",
        "mutation_count",
        "His_count",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "assigned_template_id",
        "assigned_template_status",
        "assigned_expansion_role",
        "template_match_tier",
        "template_match_level",
        "same_position_overlap_count",
        "evidence_inheritance_reason",
        "expanded_pool_role",
        "candidate_quality_tier",
        "expanded_pool_score",
        "source_pool",
        "source_pool_list",
        "source_occurrence_count",
        "hard_filter_status",
        "forbidden_pair_status",
        "buildability_light_status",
        "tier1_review_class",
        "bank_eligibility",
        "stage2a_list_action",
        "t2_class_current",
        "tier2_class",
        "selection_eligibility",
        "tier2_eligibility_status",
        "mpnn_total_score_per_residue",
        "mpnn_score_percentile_within_target",
        "neutral_retention_score",
        "acidic_release_support_score",
        "global_weakening_risk_score",
        "display_or_expression_risk_score",
        "glycan_or_epitope_risk_score",
        "hit_likelihood_score_v0",
        "tier1_rank_score",
        "complex_diagnostic_required",
    ]
    for col in keep:
        if col not in selected.columns:
            selected[col] = pd.NA
        if col not in existing_backfill.columns:
            existing_backfill[col] = pd.NA
    selected[keep].to_csv(OUT_DIR / "expanded_tier2_candidate_pool_v2.csv", index=False)
    selected[selected["expanded_pool_role"].isin(["evidence_guided_backfill", "constrained_local_expansion"])][keep].to_csv(
        OUT_DIR / "evidence_guided_backfill_candidates_v2.csv", index=False
    )
    existing_backfill.drop(columns=["_template_options"], errors="ignore").to_csv(
        OUT_DIR / "existing_source_backfill_availability_v2.csv", index=False
    )
    rejection.to_csv(OUT_DIR / "backfill_v2_rejection_report.csv", index=False)
    source_inv.to_csv(OUT_DIR / "expanded_pool_v2_source_inventory.csv", index=False)

    by_template = (
        selected.groupby(
            ["target", "assigned_template_id", "assigned_template_status", "assigned_expansion_role"],
            dropna=False,
        )
        .agg(
            count=("variant_id", "size"),
            bank_seed_count=("expanded_pool_role", lambda s: int((s == "bank_template_seed").sum())),
            backfill_count=("expanded_pool_role", lambda s: int((s == "evidence_guided_backfill").sum())),
            local_expansion_count=("expanded_pool_role", lambda s: int((s == "constrained_local_expansion").sum())),
            A_count=("template_match_tier", lambda s: int(s.isin(["A_strong", "local_expansion_A"]).sum())),
            B_count=("template_match_tier", lambda s: int(s.isin(["B_medium", "local_expansion_B"]).sum())),
            C_count=("template_match_tier", lambda s: int((s == "C_seed_only").sum())),
        )
        .reset_index()
        .sort_values(["target", "count"], ascending=[True, False])
    )
    by_template.to_csv(OUT_DIR / "backfill_v2_by_template_summary.csv", index=False)

    by_seed = (
        selected.groupby(["target", "his_seed_set"], dropna=False)
        .agg(
            count=("variant_id", "size"),
            four_mut_count=("mutation_count", lambda s: int((pd.to_numeric(s, errors="coerce") == 4).sum())),
            cluster_count=("near_duplicate_cluster_id", "nunique"),
            seed_only_count=("template_match_tier", lambda s: int((s == "C_seed_only").sum())),
        )
        .reset_index()
        .sort_values(["target", "count"], ascending=[True, False])
    )
    by_seed.to_csv(OUT_DIR / "expanded_pool_v2_by_seed.csv", index=False)

    by_cluster = (
        selected.groupby(["target", "near_duplicate_cluster_id"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
    )
    by_cluster.to_csv(OUT_DIR / "expanded_pool_v2_by_cluster.csv", index=False)

    by_tier = (
        selected.groupby(["target", "template_match_tier", "expanded_pool_role"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "template_match_tier", "count"], ascending=[True, True, False])
    )
    by_tier.to_csv(OUT_DIR / "expanded_pool_v2_by_match_tier.csv", index=False)

    by_source = (
        selected.groupby(["target", "source_pool", "template_match_tier"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "source_pool", "count"], ascending=[True, True, False])
    )
    by_source.to_csv(OUT_DIR / "expanded_pool_v2_by_source.csv", index=False)
    (OUT_DIR / "expanded_pool_v2_by_source.md").write_text(md_table(by_source) + "\n", encoding="utf-8")

    by_target = audit_rows(selected)
    by_target.to_csv(OUT_DIR / "expanded_pool_v2_audit_checks.csv", index=False)
    (OUT_DIR / "expanded_pool_v2_by_target.md").write_text(md_table(by_target) + "\n", encoding="utf-8")
    (OUT_DIR / "expanded_pool_v2_by_seed.md").write_text(md_table(by_seed) + "\n", encoding="utf-8")
    (OUT_DIR / "expanded_pool_v2_by_cluster.md").write_text(md_table(by_cluster.groupby("target", group_keys=False).head(30)) + "\n", encoding="utf-8")
    (OUT_DIR / "expanded_pool_v2_by_match_tier.md").write_text(md_table(by_tier) + "\n", encoding="utf-8")
    (OUT_DIR / "local_expansion_audit_report.md").write_text(md_table(local_audit) + "\n", encoding="utf-8")
    write_report(selected, rejection, source_inv, by_target, by_template, by_seed, by_cluster, by_tier, by_source, local_audit)


def audit_rows(selected: pd.DataFrame) -> pd.DataFrame:
    rows = []
    duplicate_count = int(selected["canonical_unique_key"].duplicated().sum()) if "canonical_unique_key" in selected else 0
    for target, cfg in TARGET_CONFIG.items():
        sub = selected[selected["target"].eq(target)].copy()
        mut = pd.to_numeric(sub["mutation_count"], errors="coerce")
        seed_counts = sub["his_seed_set"].value_counts() if len(sub) else pd.Series(dtype=int)
        cluster_counts = sub["near_duplicate_cluster_id"].value_counts() if len(sub) else pd.Series(dtype=int)
        top_seed_count = int(seed_counts.max()) if len(seed_counts) else 0
        top_cluster_count = int(cluster_counts.max()) if len(cluster_counts) else 0
        seed_only_count = int(sub["template_match_tier"].eq("C_seed_only").sum()) if len(sub) else 0
        neutral_count = int(sub["candidate_quality_tier"].eq("neutral_boundary_or_high_risk").sum()) if len(sub) else 0
        high_supported_count = int(sub["candidate_quality_tier"].isin(["high_support", "supported", "bank_reviewed"]).sum()) if len(sub) else 0
        ay111h_count = int(sub.apply(selected_seed_contains_ay111h, axis=1).sum()) if len(sub) else 0
        sdab_risk_count = int(sub.apply(selected_seed_contains_sdab_risk, axis=1).sum()) if len(sub) else 0
        complex_diag = int(sub["complex_diagnostic_required"].fillna(False).astype(bool).sum()) if len(sub) else 0

        verdict = "PASS"
        reasons = []
        if len(sub) < cfg["target_count"]:
            verdict = "PATCH"
            reasons.append("below_20k_target")
        if duplicate_count:
            verdict = "FAIL"
            reasons.append("canonical_duplicate_count_nonzero")
        if top_seed_count > cfg["seed_cap"]:
            verdict = "FAIL"
            reasons.append("top_seed_above_cap")
        if top_cluster_count > cfg["cluster_cap"]:
            verdict = "FAIL"
            reasons.append("top_cluster_above_cap")
        if int((mut == 4).sum()) > cfg["four_mut_cap"]:
            verdict = "FAIL"
            reasons.append("four_mut_above_cap")
        if seed_only_count > cfg["seed_only_cap"]:
            verdict = "FAIL"
            reasons.append("seed_only_above_cap")
        if cfg["neutral_boundary_cap"] is not None and neutral_count > cfg["neutral_boundary_cap"]:
            verdict = "FAIL"
            reasons.append("neutral_boundary_or_high_risk_above_cap")
        if cfg["high_or_supported_min"] is not None and high_supported_count < cfg["high_or_supported_min"]:
            verdict = "FAIL"
            reasons.append("high_support_plus_supported_below_min")
        if cfg["ay111h_containing_cap"] is not None and ay111h_count > cfg["ay111h_containing_cap"]:
            verdict = "FAIL"
            reasons.append("AY111H_containing_above_cap")
        if target == "Ab_sdAb":
            if sdab_risk_count:
                verdict = "FAIL"
                reasons.append("AG102H_or_AV105H_main_rows_present")
            if complex_diag != len(sub):
                verdict = "FAIL"
                reasons.append("sdAb_complex_diagnostic_required_not_retained")
        rows.append(
            {
                "target": target,
                "expanded_count": len(sub),
                "target_count": cfg["target_count"],
                "bank_seed_count": int(sub["expanded_pool_role"].eq("bank_template_seed").sum()) if len(sub) else 0,
                "existing_backfill_count": int(sub["expanded_pool_role"].eq("evidence_guided_backfill").sum()) if len(sub) else 0,
                "local_expansion_count": int(sub["expanded_pool_role"].eq("constrained_local_expansion").sum()) if len(sub) else 0,
                "A_or_local_A_count": int(sub["template_match_tier"].isin(["A_strong", "local_expansion_A"]).sum()) if len(sub) else 0,
                "B_or_local_B_count": int(sub["template_match_tier"].isin(["B_medium", "local_expansion_B"]).sum()) if len(sub) else 0,
                "seed_only_count": seed_only_count,
                "seed_only_cap": cfg["seed_only_cap"],
                "seed_only_fraction": seed_only_count / len(sub) if len(sub) else 0.0,
                "top_seed_count": top_seed_count,
                "seed_cap": cfg["seed_cap"],
                "top_seed_fraction": top_seed_count / len(sub) if len(sub) else 0.0,
                "top_cluster_count": top_cluster_count,
                "cluster_cap": cfg["cluster_cap"],
                "top_cluster_fraction": top_cluster_count / len(sub) if len(sub) else 0.0,
                "four_mut_count": int((mut == 4).sum()),
                "four_mut_cap": cfg["four_mut_cap"],
                "four_mut_fraction": float((mut == 4).mean()) if len(sub) else 0.0,
                "neutral_boundary_or_high_risk_count": neutral_count,
                "neutral_boundary_or_high_risk_cap": cfg["neutral_boundary_cap"],
                "high_support_plus_supported_count": high_supported_count,
                "high_support_plus_supported_min": cfg["high_or_supported_min"],
                "AY111H_containing_count": ay111h_count if target == "Ab_sdAb" else pd.NA,
                "AY111H_containing_cap": cfg["ay111h_containing_cap"] if target == "Ab_sdAb" else pd.NA,
                "AG102H_AV105H_main_count": sdab_risk_count if target == "Ab_sdAb" else pd.NA,
                "complex_diagnostic_required_count": complex_diag,
                "canonical_duplicate_count_global": duplicate_count,
                "verdict": verdict,
                "reasons": ";".join(reasons) if reasons else "all_v2_gates_pass",
            }
        )
    return pd.DataFrame(rows)


def write_report(
    selected: pd.DataFrame,
    rejection: pd.DataFrame,
    source_inv: pd.DataFrame,
    by_target: pd.DataFrame,
    by_template: pd.DataFrame,
    by_seed: pd.DataFrame,
    by_cluster: pd.DataFrame,
    by_tier: pd.DataFrame,
    by_source: pd.DataFrame,
    local_audit: pd.DataFrame,
) -> None:
    overall = "EXPANDED_TIER2_CANDIDATE_POOL_V2_BUILT__FINAL_SELECTION_LOCKED"
    if (by_target["verdict"] == "FAIL").any():
        overall = "EXPANDED_TIER2_CANDIDATE_POOL_V2_BUILT_WITH_FAILING_GATES__PATCH_REQUIRED"
    elif (by_target["verdict"] == "PATCH").any():
        overall = "EXPANDED_TIER2_CANDIDATE_POOL_V2_BUILT_WITH_PATCH__FINAL_SELECTION_LOCKED"

    report = [
        "# Expanded Tier2 Candidate Pool v2",
        "",
        "## Executive Summary",
        "",
        f"Verdict: `{overall}`.",
        "",
        "This stage used frozen Tier2 candidate bank templates to construct 20K-per-target expanded candidate mother pools. It used existing source-pool backfill first and constrained local expansion only where existing A/B evidence inheritance was insufficient.",
        "",
        "No AF3, SimpleFold, PyRosetta, FoldX, MD, explicit glycan modeling, or final library selection was run.",
        "",
        "The output is a candidate mother pool for later final 10K-15K planning, not a synthesis-ready final library.",
        "",
        "## Source Inventory",
        "",
        md_table(source_inv),
        "",
        "## Audit Checks",
        "",
        md_table(by_target),
        "",
        "## Match Tier Distribution",
        "",
        md_table(by_tier),
        "",
        "## Source Distribution",
        "",
        md_table(by_source.groupby("target", group_keys=False).head(20)),
        "",
        "## Seed Summary",
        "",
        md_table(by_seed),
        "",
        "## Top Cluster Summary",
        "",
        md_table(by_cluster.groupby("target", group_keys=False).head(20)),
        "",
        "## Template Summary",
        "",
        md_table(by_template.groupby("target", group_keys=False).head(30)),
        "",
        "## Local Expansion Audit",
        "",
        md_table(local_audit.groupby("target", group_keys=False).head(30)),
        "",
        "## Rejection Summary",
        "",
        md_table(
            rejection.groupby(["target", "reject_reason"], dropna=False)["count"]
            .sum()
            .reset_index()
            .sort_values(["target", "count"], ascending=[True, False])
            .head(50)
        ),
        "",
        "## Interpretation",
        "",
        "- 1E62 remains the primary expanded branch.",
        "- sdAb remains the secondary / exploratory branch; all sdAb rows retain `complex_diagnostic_required` and are not structure-confirmed.",
        "- A/B/local-expansion rows are explicitly separated from C seed-only rows; seed-only support remains capped.",
        "- Local expansion is constrained by the fixed 40 aa window, fixed His seed, hard-protect positions, no new Cys, no new canonical N-X-S/T, and mutation-count limits.",
        "- Glycan remains unchecked; this report must not be cited as low glycan risk.",
        "",
        "## Next Gate",
        "",
        "Allowed now:",
        "",
        "- Review the v2 expanded pool family, match-tier, seed, cluster, source, and local-expansion balance.",
        "- Use this v2 pool as the candidate universe for later final candidate-pool construction after explicit review.",
        "- Plan parent complex baseline diagnostics, especially before upgrading sdAb confidence.",
        "",
        "Still locked:",
        "",
        "- Final 10K / 15K library selection.",
        "- Broad Tier2-heavy or broad AF3/SimpleFold compute.",
        "- MD or constant-pH MD.",
        "- Treating sdAb as structure-confirmed.",
        "- Reporting glycan low risk.",
        "",
        "## Output Files",
        "",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool_v2/expanded_tier2_candidate_pool_v2.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool_v2/evidence_guided_backfill_candidates_v2.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool_v2/local_expansion_candidates.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool_v2/backfill_v2_by_template_summary.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool_v2/backfill_v2_rejection_report.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool_v2/expanded_pool_v2_audit_report.md`",
        "",
    ]
    text = "\n".join(report)
    (OUT_DIR / "expanded_pool_v2_audit_report.md").write_text(text, encoding="utf-8")
    (OUT_DIR / "expanded_pool_v2_next_gate_report.md").write_text("\n".join(report[: report.index("## Output Files")]) + "\n", encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(text, encoding="utf-8")

    manifest = [
        "stage: expanded_tier2_candidate_pool_v2",
        f"status: {overall}",
        f"bank: {BANK.relative_to(ROOT)}",
        f"bank_sha256: {sha256_file(BANK)}",
        f"templates: {TEMPLATES.relative_to(ROOT)}",
        f"templates_sha256: {sha256_file(TEMPLATES)}",
        f"expanded_rows: {len(selected)}",
        "targets:",
    ]
    for _, row in by_target.iterrows():
        manifest += [
            f"  {row['target']}:",
            f"    expanded_count: {row['expanded_count']}",
            f"    verdict: {row['verdict']}",
            f"    reasons: {row['reasons']}",
        ]
    manifest += [
        "locks:",
        "  final_10k_or_15k: true",
        "  broad_tier2_heavy: true",
        "  new_af3_simplefold_pyrosetta_foldx: true",
        "  md: true",
        "  sdab_structure_confirmed_upgrade: true",
        "  glycan_low_risk_claim: true",
    ]
    (OUT_DIR / "expanded_pool_v2_manifest.yaml").write_text("\n".join(manifest) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    templates = pd.read_csv(TEMPLATES)
    bank = pd.read_csv(BANK)
    source_inv = source_inventory()
    sources = load_candidate_sources()
    parents = parent_sequences_from_sources(sources)
    filtered, prefilter_rejections = filter_candidates(sources, templates)
    existing_backfill = build_backfill_candidates(filtered, templates)
    local_expansion, local_audit, local_plan, local_config = generate_local_expansion_candidates(templates, sources)
    candidates = prepare_selection_candidates(existing_backfill, local_expansion)
    bank_rows = bank_seed_rows(bank, templates, parents)
    selected, cap_rejections = select_pool(bank_rows, candidates)
    selected = selected.drop(columns=["_template_options"], errors="ignore")

    if prefilter_rejections.empty:
        rejection = cap_rejections
    elif cap_rejections.empty:
        rejection = prefilter_rejections
    else:
        rejection = pd.concat([prefilter_rejections, cap_rejections], ignore_index=True)
        rejection = (
            rejection.groupby(["target", "source_pool", "template_status", "reject_reason"], dropna=False)["count"]
            .sum()
            .reset_index()
            .sort_values(["target", "count"], ascending=[True, False])
        )

    write_outputs(selected, existing_backfill, local_expansion, local_audit, rejection, source_inv, local_plan, local_config)
    print(f"wrote {OUT_DIR / 'expanded_tier2_candidate_pool_v2.csv'} rows={len(selected)}")
    print(
        selected.groupby(["target", "expanded_pool_role"], dropna=False)
        .size()
        .reset_index(name="count")
        .to_string(index=False)
    )
    print(f"wrote {OUT_DIR / 'expanded_pool_v2_audit_report.md'}")


if __name__ == "__main__":
    main()
