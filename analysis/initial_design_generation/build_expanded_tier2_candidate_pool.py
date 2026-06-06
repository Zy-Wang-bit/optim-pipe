#!/usr/bin/env python3
"""Build an expanded Tier2 candidate pool by template-guided backfill.

This stage returns to existing production / reserve / recovery pools and finds
variants similar to the evidence-backed Tier2 candidate bank templates. It does
not generate new variants and does not run AF3, SimpleFold, PyRosetta, MD, or
final library selection.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"
BANK_DIR = ROOT / "results/initial_design_generation/tier2_candidate_bank"
OUT_DIR = ROOT / "results/initial_design_generation/expanded_tier2_candidate_pool"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"

BANK = BANK_DIR / "tier2_candidate_bank.csv"
TEMPLATES = BANK_DIR / "tier2_candidate_bank_expansion_templates.csv"

SOURCE_FILES = [
    (
        "stage2a_candidate_list",
        ROOT / "results/initial_design_generation/stage1_5_stage2a/stage2a_candidate_list.csv",
        0,
    ),
    (
        "sdab_recovery_passlike_supplement",
        ROOT / "results/initial_design_generation/sdab_recovery_loop/round_02_passlike_supplement/sdab_candidate_bank.csv",
        1,
    ),
    (
        "tier2_candidate_snapshot",
        ROOT / "results/initial_design_generation/tier2_staged/tier2_candidate_snapshot.csv",
        2,
    ),
    (
        "tier2_core_reserve_pool",
        ROOT / "results/initial_design_generation/tier1_filtering/tier2_core_reserve_pool.csv",
        3,
    ),
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
    "mutation_list",
    "normalized_mutation_list",
    "mutation_count",
    "His_count",
    "his_seed_set",
    "near_duplicate_cluster_id",
    "tier1_near_duplicate_cluster_id",
    "source_near_duplicate_cluster_id",
    "hard_filter_status",
    "forbidden_pair_status",
    "exact_duplicate_status",
    "selection_eligibility",
    "tier2_eligibility_status",
    "tier1_review_class",
    "bank_eligibility",
    "stage2a_list_action",
    "t2_class_current",
    "tier2_class",
    "mpnn_total_score_per_residue",
    "mpnn_score_percentile_within_target",
    "neutral_retention_score",
    "acidic_release_support_score",
    "global_weakening_risk_score",
    "display_or_expression_risk_score",
    "glycan_or_epitope_risk_score",
    "hit_likelihood_score_v0",
    "tier1_rank_score",
]

TARGET_CONFIG = {
    "Ab_1E62": {
        "target_count": 12000,
        "minimum_count": 10000,
        "seed_cap": 2700,
        "cluster_cap": 500,
        "four_mut_cap": 9600,
        "ay111h_containing_cap": None,
        "template_status_quota": {
            "active_primary_template": 4800,
            "active_primary_backfill_template": 6600,
            "limited_boundary_template": 600,
        },
        "template_cap": {
            "active_primary_template": 400,
            "active_primary_backfill_template": 350,
            "limited_boundary_template": 100,
        },
    },
    "Ab_sdAb": {
        "target_count": 8000,
        "minimum_count": 5000,
        "seed_cap": 2000,
        "cluster_cap": 450,
        "four_mut_cap": 800,
        "ay111h_containing_cap": 2400,
        "template_status_quota": {
            "secondary_template_complex_weak": 3200,
            "secondary_template_complex_unchecked": 2800,
            "low_priority_secondary_template": 2000,
            "boundary_representative_only": 4,
        },
        "template_cap": {
            "secondary_template_complex_weak": 220,
            "secondary_template_complex_unchecked": 220,
            "low_priority_secondary_template": 180,
            "boundary_representative_only": 4,
        },
    },
}

SOURCE_CLASS_SCORE = {
    "A_tier1_priority_candidate": 0.30,
    "B_rescue_enriched_candidate": 0.24,
    "E_sparse_mpnn_supplement_review": 0.10,
    "D_favorable_mpnn_but_weak_pH_mechanism": 0.02,
    "F_neutral_boundary_or_high_risk": -0.05,
}


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


def mutation_position(mut: str) -> str:
    match = re.match(r"^([A-Z])([A-Z]?)(\d+)([A-Z])$", str(mut))
    if not match:
        return str(mut)
    return f"{match.group(1)}{match.group(3)}"


def split_mutations(mutation_list: object) -> tuple[list[str], list[str]]:
    if not isinstance(mutation_list, str) or not mutation_list:
        return [], []
    muts = [m.strip() for m in mutation_list.split(";") if m.strip() and m.strip().lower() != "nan"]
    his = [m for m in muts if m.endswith("H")]
    non_his = [m for m in muts if not m.endswith("H")]
    return his, non_his


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def canonical_key(row: pd.Series) -> str:
    for col in ["canonical_sequence_hash_full", "canonical_recovery_sequence_hash", "sequence_hash"]:
        val = row.get(col)
        if isinstance(val, str) and val and val.lower() != "nan":
            return f"{row['target']}|{val}"
    variant_id = row.get("variant_id")
    if isinstance(variant_id, str) and variant_id:
        return f"{row['target']}|variant|{variant_id}"
    return f"{row['target']}|mutation|{row.get('mutation_list', '')}"


def read_source(path: Path, source_pool: str, source_rank: int) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    use = [c for c in READ_COLS if c in cols]
    df = pd.read_csv(path, usecols=use)
    for col in READ_COLS:
        if col not in df.columns:
            df[col] = pd.NA
    df["source_pool"] = source_pool
    df["source_rank"] = source_rank
    return df


def source_inventory() -> pd.DataFrame:
    rows = []
    for source_pool, path, _rank in SOURCE_FILES:
        if path.exists():
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


def load_candidate_sources() -> pd.DataFrame:
    frames = [read_source(path, source_pool, rank) for source_pool, path, rank in SOURCE_FILES]
    frames = [f for f in frames if not f.empty]
    df = pd.concat(frames, ignore_index=True)
    for col in [
        "mutation_count",
        "His_count",
        "mpnn_total_score_per_residue",
        "mpnn_score_percentile_within_target",
        "neutral_retention_score",
        "acidic_release_support_score",
        "global_weakening_risk_score",
        "display_or_expression_risk_score",
        "glycan_or_epitope_risk_score",
        "hit_likelihood_score_v0",
        "tier1_rank_score",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

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


def filter_candidates(df: pd.DataFrame, templates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    allowed_seeds = (
        templates[templates["allowed_backfill"].map(as_bool)]
        .groupby("target")["seed_pattern"]
        .apply(lambda s: set(map(str, s)))
        .to_dict()
    )
    reason_frames: list[pd.DataFrame] = []
    mask = pd.Series(True, index=df.index)

    allowed_pair = {
        f"{target}||{seed}"
        for target, seed_set in allowed_seeds.items()
        for seed in seed_set
    }
    pair = df["target"].astype(str) + "||" + df["his_seed_set"].astype(str)
    no_seed = ~pair.isin(allowed_pair)
    reason_frames.append(reject_counts_for(df[mask & no_seed], "no_template_seed_match"))
    mask &= ~no_seed

    hard_fail = df["hard_filter_status"].notna() & ~df["hard_filter_status"].fillna("pass").eq("pass")
    reason_frames.append(reject_counts_for(df[mask & hard_fail], "hard_filter_not_pass"))
    mask &= ~hard_fail

    forbidden_fail = df["forbidden_pair_status"].notna() & ~df["forbidden_pair_status"].fillna("pass").eq("pass")
    reason_frames.append(reject_counts_for(df[mask & forbidden_fail], "forbidden_pair_not_pass"))
    mask &= ~forbidden_fail

    mut_fail = df["mutation_count"].fillna(99).gt(4)
    reason_frames.append(reject_counts_for(df[mask & mut_fail], "mutation_count_gt_4"))
    mask &= ~mut_fail

    audit_demoted = df["tier1_review_class"].fillna("").eq("G_audit_or_demoted")
    reason_frames.append(reject_counts_for(df[mask & audit_demoted], "tier1_audit_or_demoted"))
    mask &= ~audit_demoted

    risk_text = df["his_seed_set"].fillna("").astype(str) + ";" + df["mutation_list"].fillna("").astype(str)
    risk = df["target"].eq("Ab_sdAb") & (
        risk_text.str.contains("AG102H", regex=False)
        | risk_text.str.contains("AV105H", regex=False)
    )
    reason_frames.append(reject_counts_for(df[mask & risk], "sdAb_AG102H_or_AV105H_main_excluded"))
    mask &= ~risk

    filtered = df[mask].copy()
    reject_summary = pd.concat([f for f in reason_frames if not f.empty], ignore_index=True) if reason_frames else pd.DataFrame()
    if not reject_summary.empty:
        reject_summary = (
            reject_summary.groupby(["target", "source_pool", "template_status", "reject_reason"], dropna=False)["count"]
            .sum()
            .reset_index(name="count")
            .sort_values(["target", "source_pool", "count"], ascending=[True, True, False])
        )
    return filtered, reject_summary


def reject_counts_for(df: pd.DataFrame, reason: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    return (
        df.groupby(["target", "source_pool"], dropna=False)
        .size()
        .reset_index(name="count")
        .assign(template_status="prefilter", reject_reason=reason)
    )


def template_options(row: pd.Series, by_target_seed: dict[tuple[str, str], pd.DataFrame]) -> list[dict[str, object]]:
    key = (row["target"], str(row.get("his_seed_set")))
    if key not in by_target_seed:
        return []
    _his, non_his = split_mutations(row.get("mutation_list"))
    non_his_set = set(non_his)
    non_his_pos = {mutation_position(m) for m in non_his}
    options = []
    for _, tpl in by_target_seed[key].iterrows():
        rescue_pattern = str(tpl.get("rescue_pattern", "none"))
        if rescue_pattern == "none" or rescue_pattern.lower() == "nan":
            match_level = "seed_only"
            match_score = 0
        elif rescue_pattern in non_his_set:
            match_level = "exact_rescue"
            match_score = 3
        elif mutation_position(rescue_pattern) in non_his_pos:
            match_level = "same_rescue_position"
            match_score = 2
        else:
            match_level = "seed_only"
            match_score = 0
        options.append(
            {
                "template_id": tpl["template_id"],
                "template_status": tpl["template_status"],
                "recommended_expansion_role": tpl["recommended_expansion_role"],
                "complex_diagnostic_required": as_bool(tpl["complex_diagnostic_required"]),
                "template_priority": int(tpl["priority_level"]),
                "match_level": match_level,
                "match_score": match_score,
            }
        )
    return sorted(options, key=lambda x: (-int(x["match_score"]), int(x["template_priority"]), str(x["template_id"])))


def quality_tier(row: pd.Series) -> str:
    cls = str(row.get("tier1_review_class", ""))
    bank = str(row.get("bank_eligibility", ""))
    t2 = str(row.get("t2_class_current", row.get("tier2_class", "")))
    if bank == "pass_like" or cls == "A_tier1_priority_candidate" or t2 == "T2_strong_candidate":
        return "high_support"
    if cls in {"B_rescue_enriched_candidate", "E_sparse_mpnn_supplement_review"} or t2 == "T2_good_candidate":
        return "supported"
    if cls == "D_favorable_mpnn_but_weak_pH_mechanism":
        return "mpnn_favorable_pH_weak"
    if cls == "F_neutral_boundary_or_high_risk":
        return "neutral_boundary_or_high_risk"
    return "proxy_support"


def candidate_score(row: pd.Series, options: list[dict[str, object]]) -> float:
    best_match = max([o["match_score"] for o in options], default=0)
    score = 0.0
    score += SOURCE_CLASS_SCORE.get(str(row.get("tier1_review_class", "")), 0.0)
    if row.get("bank_eligibility") == "pass_like":
        score += 0.30
    elif row.get("bank_eligibility") == "supported_boundary":
        score += 0.08
    t2 = str(row.get("t2_class_current", row.get("tier2_class", "")))
    if t2 == "T2_strong_candidate":
        score += 0.25
    elif t2 == "T2_good_candidate":
        score += 0.15
    elif t2 == "T2_release_possible_but_neutral_risky":
        score -= 0.10
    if row.get("stage2a_list_action") == "stage2a_include_priority":
        score += 0.12
    score += float(row.get("neutral_retention_score", 0) or 0) * 0.15
    score += float(row.get("acidic_release_support_score", 0) or 0) * 0.18
    score -= float(row.get("global_weakening_risk_score", 0) or 0) * 0.12
    pct = row.get("mpnn_score_percentile_within_target")
    if pd.notna(pct):
        score += max(0.0, 1.0 - float(pct)) * 0.18
    rank = row.get("tier1_rank_score")
    if pd.notna(rank):
        score += float(rank) * 0.12
    mut_count = row.get("mutation_count")
    if pd.notna(mut_count):
        if int(mut_count) <= 2:
            score += 0.05
        elif int(mut_count) == 4:
            score -= 0.05
    score += best_match * 0.05
    score += max(0, 5 - int(row.get("source_rank", 5))) * 0.02
    return score


def selected_seed_contains_ay111h(row: pd.Series) -> bool:
    text = f"{row.get('his_seed_set', '')};{row.get('mutation_list', '')}"
    return "AY111H" in text


def bank_seed_rows(bank: pd.DataFrame, templates: pd.DataFrame) -> pd.DataFrame:
    tpl = templates[["template_id", "variant_id"]]
    out = bank.merge(tpl, on="variant_id", how="left", validate="one_to_one")
    out["canonical_unique_key"] = out.apply(lambda r: f"{r['target']}|bank|{r['variant_id']}", axis=1)
    out["source_pool"] = "tier2_candidate_bank"
    out["source_pool_list"] = "tier2_candidate_bank"
    out["source_occurrence_count"] = 1
    out["expanded_pool_role"] = "bank_template_seed"
    out["template_match_level"] = "template_seed"
    out["candidate_quality_tier"] = "bank_reviewed"
    out["expanded_pool_score"] = 10.0 - out["bank_priority"].fillna(9) * 0.01
    out["source_variant_id"] = out["variant_id"]
    return out


def build_backfill_candidates(candidates: pd.DataFrame, templates: pd.DataFrame) -> pd.DataFrame:
    backfill_templates = templates[templates["allowed_backfill"].map(as_bool)].copy()
    by_target_seed = {
        key: group.copy()
        for key, group in backfill_templates.groupby(["target", "seed_pattern"], dropna=False)
    }
    rows = []
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


def select_pool(bank_rows: pd.DataFrame, candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected: list[dict[str, object]] = []
    selected_keys: set[str] = set()
    selected_variant_keys: set[tuple[str, str]] = set()
    counters = {
        "seed": defaultdict(Counter),
        "cluster": defaultdict(Counter),
        "template_status": defaultdict(Counter),
        "template": defaultdict(Counter),
        "four_mut": Counter(),
        "ay111h": Counter(),
    }
    target_counts: Counter[str] = Counter()
    reject_counts: Counter[tuple[str, str, str, str]] = Counter()

    def can_add(item: dict[str, object], option: dict[str, object], is_bank: bool = False) -> tuple[bool, str]:
        target = str(item["target"])
        cfg = TARGET_CONFIG[target]
        seed = str(item.get("his_seed_set", ""))
        cluster = str(item.get("near_duplicate_cluster_id", "missing_cluster"))
        status = str(option["template_status"])
        template_id = str(option["template_id"])
        mut_count = int(float(item.get("mutation_count", 99))) if pd.notna(item.get("mutation_count")) else 99

        if target_counts[target] >= cfg["target_count"]:
            return False, "target_count_full"
        if counters["seed"][target][seed] >= cfg["seed_cap"]:
            return False, "seed_cap_full"
        if counters["cluster"][target][cluster] >= cfg["cluster_cap"]:
            return False, "cluster_cap_full"
        if status in cfg["template_status_quota"] and counters["template_status"][target][status] >= cfg["template_status_quota"][status]:
            return False, "template_status_quota_full"
        if status not in cfg["template_status_quota"]:
            return False, "template_status_not_allowed"
        if counters["template"][target][template_id] >= cfg["template_cap"].get(status, 100):
            return False, "template_cap_full"
        if mut_count == 4 and counters["four_mut"][target] >= cfg["four_mut_cap"]:
            return False, "four_mut_cap_full"
        if cfg["ay111h_containing_cap"] is not None and selected_seed_contains_ay111h(pd.Series(item)):
            if counters["ay111h"][target] >= cfg["ay111h_containing_cap"]:
                return False, "AY111H_containing_cap_full"
        return True, "selected"

    def add(item: dict[str, object], option: dict[str, object], role: str) -> None:
        target = str(item["target"])
        seed = str(item.get("his_seed_set", ""))
        cluster = str(item.get("near_duplicate_cluster_id", "missing_cluster"))
        mut_count = int(float(item.get("mutation_count", 99))) if pd.notna(item.get("mutation_count")) else 99
        item = dict(item)
        item["assigned_template_id"] = option["template_id"]
        item["assigned_template_status"] = option["template_status"]
        item["assigned_expansion_role"] = option["recommended_expansion_role"]
        item["template_match_level"] = option.get("match_level", item.get("template_match_level", "seed_only"))
        item["complex_diagnostic_required"] = option.get("complex_diagnostic_required", item.get("complex_diagnostic_required", False))
        item["expanded_pool_role"] = role
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
        if selected_seed_contains_ay111h(pd.Series(item)):
            counters["ay111h"][target] += 1

    # Seed reviewed bank rows first.
    for _, row in bank_rows.sort_values(["target", "bank_priority", "variant_id"]).iterrows():
        option = {
            "template_id": row["template_id"],
            "template_status": row["template_status"],
            "recommended_expansion_role": row["recommended_expansion_role"],
            "complex_diagnostic_required": as_bool(row.get("complex_diagnostic_required", False)),
            "match_level": "template_seed",
        }
        ok, reason = can_add(row.to_dict(), option, is_bank=True)
        if ok:
            add(row.to_dict(), option, "bank_template_seed")
        else:
            reject_counts[(row["target"], "tier2_candidate_bank", str(row["template_status"]), reason)] += 1

    candidates = candidates.sort_values(
        ["target", "expanded_pool_score", "source_rank", "variant_id"],
        ascending=[True, False, True, True],
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
                add(row.to_dict(), option, "evidence_guided_backfill")
                break
            option_reasons.append((option, reason))
        else:
            if option_reasons:
                option, reason = option_reasons[0]
                reject_counts[(target, str(row["source_pool"]), str(option["template_status"]), reason)] += 1
            else:
                reject_counts[(target, str(row["source_pool"]), "none", "no_template_option")] += 1

    selected_df = pd.DataFrame(selected)
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
    ).sort_values(["target", "count"], ascending=[True, False])
    return selected_df, rejection


def write_outputs(selected: pd.DataFrame, rejection: pd.DataFrame, source_inv: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
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
        "template_match_level",
        "expanded_pool_role",
        "candidate_quality_tier",
        "expanded_pool_score",
        "source_pool",
        "source_pool_list",
        "source_occurrence_count",
        "hard_filter_status",
        "forbidden_pair_status",
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
    selected[keep].to_csv(OUT_DIR / "expanded_tier2_candidate_pool.csv", index=False)
    selected[selected["expanded_pool_role"].eq("evidence_guided_backfill")][keep].to_csv(
        OUT_DIR / "evidence_guided_backfill_candidates.csv", index=False
    )
    rejection.to_csv(OUT_DIR / "backfill_rejection_report.csv", index=False)
    source_inv.to_csv(OUT_DIR / "expanded_pool_source_inventory.csv", index=False)

    by_template = (
        selected.groupby(["target", "assigned_template_id", "assigned_template_status", "assigned_expansion_role"], dropna=False)
        .agg(
            count=("variant_id", "size"),
            bank_seed_count=("expanded_pool_role", lambda s: int((s == "bank_template_seed").sum())),
            backfill_count=("expanded_pool_role", lambda s: int((s == "evidence_guided_backfill").sum())),
            exact_rescue_count=("template_match_level", lambda s: int((s == "exact_rescue").sum())),
            same_position_count=("template_match_level", lambda s: int((s == "same_rescue_position").sum())),
            seed_only_count=("template_match_level", lambda s: int((s == "seed_only").sum())),
        )
        .reset_index()
        .sort_values(["target", "count"], ascending=[True, False])
    )
    by_template.to_csv(OUT_DIR / "backfill_by_template_summary.csv", index=False)

    by_seed = (
        selected.groupby(["target", "his_seed_set"], dropna=False)
        .agg(
            count=("variant_id", "size"),
            four_mut_count=("mutation_count", lambda s: int((pd.to_numeric(s, errors="coerce") == 4).sum())),
            cluster_count=("near_duplicate_cluster_id", "nunique"),
        )
        .reset_index()
        .sort_values(["target", "count"], ascending=[True, False])
    )
    by_seed.to_csv(OUT_DIR / "expanded_pool_by_seed.csv", index=False)

    by_cluster = (
        selected.groupby(["target", "near_duplicate_cluster_id"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
    )
    by_cluster.to_csv(OUT_DIR / "expanded_pool_by_cluster.csv", index=False)

    by_role = (
        selected.groupby(["target", "assigned_template_status", "expanded_pool_role", "candidate_quality_tier"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "assigned_template_status", "count"], ascending=[True, True, False])
    )
    by_role.to_csv(OUT_DIR / "expanded_pool_by_role.csv", index=False)

    write_reports(selected, rejection, source_inv, by_template, by_seed, by_cluster, by_role)


def audit_rows(selected: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target, cfg in TARGET_CONFIG.items():
        sub = selected[selected["target"].eq(target)].copy()
        mut = pd.to_numeric(sub["mutation_count"], errors="coerce")
        seed_counts = sub["his_seed_set"].value_counts() if len(sub) else pd.Series(dtype=int)
        cluster_counts = sub["near_duplicate_cluster_id"].value_counts() if len(sub) else pd.Series(dtype=int)
        top_seed_count = int(seed_counts.max()) if len(seed_counts) else 0
        top_cluster_count = int(cluster_counts.max()) if len(cluster_counts) else 0
        top_seed_frac = top_seed_count / len(sub) if len(sub) else 0.0
        top_cluster_frac = top_cluster_count / len(sub) if len(sub) else 0.0
        ay111h_count = int(sub.apply(selected_seed_contains_ay111h, axis=1).sum()) if len(sub) else 0
        verdict = "PASS"
        reasons = []
        if len(sub) < cfg["minimum_count"]:
            verdict = "FAIL"
            reasons.append("below_minimum_count")
        elif len(sub) < cfg["target_count"]:
            reasons.append("below_preferred_target_count_within_acceptable_range")
        if top_seed_count > cfg["seed_cap"]:
            verdict = "PATCH" if verdict == "PASS" else verdict
            reasons.append("top_seed_count_above_absolute_cap")
        if top_cluster_count > cfg["cluster_cap"]:
            verdict = "PATCH" if verdict == "PASS" else verdict
            reasons.append("top_cluster_count_above_absolute_cap")
        if int((mut == 4).sum()) > cfg["four_mut_cap"]:
            verdict = "PATCH" if verdict == "PASS" else verdict
            reasons.append("four_mut_above_cap")
        if cfg["ay111h_containing_cap"] is not None and ay111h_count > cfg["ay111h_containing_cap"]:
            verdict = "PATCH" if verdict == "PASS" else verdict
            reasons.append("AY111H_containing_above_cap")
        rows.append(
            {
                "target": target,
                "expanded_count": len(sub),
                "target_count": cfg["target_count"],
                "minimum_count": cfg["minimum_count"],
                "bank_seed_count": int(sub["expanded_pool_role"].eq("bank_template_seed").sum()) if len(sub) else 0,
                "backfill_count": int(sub["expanded_pool_role"].eq("evidence_guided_backfill").sum()) if len(sub) else 0,
                "top_seed_count": top_seed_count,
                "seed_cap": cfg["seed_cap"],
                "top_seed_fraction": float(top_seed_frac),
                "top_cluster_count": top_cluster_count,
                "cluster_cap": cfg["cluster_cap"],
                "top_cluster_fraction": float(top_cluster_frac),
                "four_mut_count": int((mut == 4).sum()),
                "four_mut_cap": cfg["four_mut_cap"],
                "four_mut_fraction": float((mut == 4).mean()) if len(sub) else 0.0,
                "AY111H_containing_count": ay111h_count if target == "Ab_sdAb" else pd.NA,
                "AY111H_containing_cap": cfg["ay111h_containing_cap"] if target == "Ab_sdAb" else pd.NA,
                "AY111H_containing_fraction": float(ay111h_count / len(sub)) if target == "Ab_sdAb" and len(sub) else pd.NA,
                "complex_diagnostic_required_count": int(sub["complex_diagnostic_required"].fillna(False).astype(bool).sum()) if len(sub) else 0,
                "verdict": verdict,
                "reasons": ";".join(reasons) if reasons else "all_caps_pass",
            }
        )
    return pd.DataFrame(rows)


def write_reports(
    selected: pd.DataFrame,
    rejection: pd.DataFrame,
    source_inv: pd.DataFrame,
    by_template: pd.DataFrame,
    by_seed: pd.DataFrame,
    by_cluster: pd.DataFrame,
    by_role: pd.DataFrame,
) -> None:
    audit = audit_rows(selected)
    audit.to_csv(OUT_DIR / "expanded_pool_audit_checks.csv", index=False)

    report = [
        "# Expanded Tier2 Candidate Pool",
        "",
        "## Executive Summary",
        "",
        "Verdict: `EXPANDED_TIER2_CANDIDATE_POOL_BUILT__FINAL_10K_LOCKED`.",
        "",
        "This stage used the frozen Tier2 candidate bank templates to backfill similar candidates from existing production / Tier1 / reserve / Stage-2A / sdAb recovery pools. It did not generate new variants and did not run new AF3, SimpleFold, PyRosetta, MD, glycan modeling, or final library selection.",
        "",
        "The expanded pool is a larger candidate pool for later final candidate construction, not a final 10K library.",
        "",
        "## Source Inventory",
        "",
        md_table(source_inv),
        "",
        "## Audit Checks",
        "",
        md_table(audit),
        "",
        "## Role / Quality Distribution",
        "",
        md_table(by_role),
        "",
        "## Seed Summary",
        "",
        md_table(by_seed),
        "",
        "## Top Cluster Summary",
        "",
        md_table(by_cluster.groupby("target", group_keys=False).head(20)),
        "",
        "## Template Backfill Summary",
        "",
        md_table(by_template.groupby("target", group_keys=False).head(30)),
        "",
        "## Rejection Summary",
        "",
        md_table(rejection.groupby(["target", "reject_reason"], dropna=False)["count"].sum().reset_index().sort_values(["target", "count"], ascending=[True, False]).head(40)),
        "",
        "## Interpretation",
        "",
        "- 1E62 remains the primary expanded branch. Its pool is built from active primary templates, active primary backfill templates, and a limited boundary-template quota.",
        "- sdAb remains a secondary branch. Backfill is allowed, but all sdAb rows retain complex-diagnostic requirements and are not structure-confirmed.",
        "- Boundary templates are capped and retained for diversity / traceability, not as broad expansion drivers.",
        "- Glycan remains unchecked because this stage used existing feature tables and did not run explicit glycan review.",
        "",
        "## Next Gate",
        "",
        "Allowed now:",
        "",
        "- Audit expanded pool family balance, seed balance, cluster balance, and source composition.",
        "- Plan parent complex baseline diagnostic, especially for sdAb.",
        "- Use the expanded pool as the input universe for later final candidate pool construction after explicit review.",
        "",
        "Still locked:",
        "",
        "- Final 10K library selection.",
        "- Broad Tier2-heavy or broad AF3/SimpleFold compute.",
        "- MD or constant-pH MD.",
        "- Treating sdAb as structure-confirmed.",
        "- Reporting glycan low risk.",
        "",
        "## Output Files",
        "",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool/expanded_tier2_candidate_pool.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool/evidence_guided_backfill_candidates.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool/backfill_by_template_summary.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool/backfill_rejection_report.csv`",
        "- `results/initial_design_generation/expanded_tier2_candidate_pool/expanded_pool_audit_report.md`",
        "",
    ]
    (OUT_DIR / "expanded_pool_audit_report.md").write_text("\n".join(report), encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text("\n".join(report), encoding="utf-8")

    manifest = [
        "stage: expanded_tier2_candidate_pool",
        "status: complete",
        f"bank: {BANK.relative_to(ROOT)}",
        f"bank_sha256: {sha256_file(BANK)}",
        f"templates: {TEMPLATES.relative_to(ROOT)}",
        f"templates_sha256: {sha256_file(TEMPLATES)}",
        f"expanded_rows: {len(selected)}",
        f"backfill_rows: {int(selected['expanded_pool_role'].eq('evidence_guided_backfill').sum())}",
        "locks:",
        "  final_10k: true",
        "  broad_tier2_heavy: true",
        "  md: true",
        "  sdab_structure_confirmed_upgrade: true",
        "  glycan_low_risk_claim: true",
    ]
    (OUT_DIR / "expanded_tier2_candidate_pool_manifest.yaml").write_text("\n".join(manifest) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    templates = pd.read_csv(TEMPLATES)
    bank = pd.read_csv(BANK)
    source_inv = source_inventory()
    sources = load_candidate_sources()
    filtered, prefilter_rejections = filter_candidates(sources, templates)
    backfill_candidates = build_backfill_candidates(filtered, templates)
    bank_rows = bank_seed_rows(bank, templates)
    selected, cap_rejections = select_pool(bank_rows, backfill_candidates)
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
    write_outputs(selected, rejection, source_inv)
    print(f"wrote {OUT_DIR / 'expanded_tier2_candidate_pool.csv'} rows={len(selected)}")
    print(f"wrote {OUT_DIR / 'evidence_guided_backfill_candidates.csv'} rows={int(selected['expanded_pool_role'].eq('evidence_guided_backfill').sum())}")
    print(f"wrote {OUT_DIR / 'expanded_pool_audit_report.md'}")


if __name__ == "__main__":
    main()
