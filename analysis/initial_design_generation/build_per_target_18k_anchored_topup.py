#!/usr/bin/env python3
"""Build per-target 18K candidate drafts by anchored top-up from frozen 20K pools.

This step preserves the existing per-target 15K drafts as anchors and adds only
main candidates from the same target's expanded Tier2 v2 mother pool. It does
not generate new variants or run structural computation.
"""

import csv
import os
import sys
from collections import Counter


SCRIPT_DIR = os.path.dirname(__file__)
sys.path.insert(0, SCRIPT_DIR)

import build_per_target_15k_candidate_pools as base  # noqa: E402


ROOT = base.ROOT
EXPANDED_V2 = base.EXPANDED_V2
OUT_DIR = os.path.join(
    ROOT,
    "results/initial_design_generation/per_target_18k_candidate_pools",
)
TASK_REPORT = os.path.join(
    ROOT,
    ".tasks/active/initial-design-generation/current_stage_report.md",
)

ANCHOR_15K = {
    "Ab_1E62": os.path.join(
        ROOT,
        "results/initial_design_generation/per_target_15k_candidate_pools/"
        "final_candidate_pool_1E62_15k_draft.csv",
    ),
    "Ab_sdAb": os.path.join(
        ROOT,
        "results/initial_design_generation/per_target_15k_candidate_pools/"
        "final_candidate_pool_sdAb_15k_draft.csv",
    ),
}

TARGETS = ["Ab_1E62", "Ab_sdAb"]
TARGET_SIZE = 18000
CONTROL_ANCHOR_COUNT = 20

TARGET_CONFIG = {
    "Ab_1E62": {
        "short_name": "1E62",
        "top_seed_cap": int(TARGET_SIZE * 0.28),
        "top_cluster_cap": int(TARGET_SIZE * 0.06),
        "seed_only_cap": int(TARGET_SIZE * 0.20),
        "neutral_boundary_cap": int(TARGET_SIZE * 0.25),
        "local_expansion_cap": int(TARGET_SIZE * 0.65),
        "four_mut_cap": int(TARGET_SIZE * 0.65),
        "ay111h_cap": None,
        "exclude_ag102h_av105h_main": False,
    },
    "Ab_sdAb": {
        "short_name": "sdAb",
        "top_seed_cap": int(TARGET_SIZE * 0.30),
        "top_cluster_cap": int(TARGET_SIZE * 0.06),
        "seed_only_cap": int(TARGET_SIZE * 0.20),
        "neutral_boundary_cap": int(TARGET_SIZE * 0.10),
        "local_expansion_cap": int(TARGET_SIZE * 0.40),
        "four_mut_cap": int(TARGET_SIZE * 0.15),
        "ay111h_cap": int(TARGET_SIZE * 0.35),
        "exclude_ag102h_av105h_main": True,
    },
}


def row_key(row):
    for field in ("canonical_unique_key", "canonical_sequence_hash_full", "variant_id"):
        value = str(row.get(field, "") or "").strip()
        if value:
            if field == "canonical_sequence_hash_full":
                return "%s|%s" % (row.get("target", ""), value)
            return value
    return ""


def read_rows(path):
    fields, rows = base.read_csv(path)
    return fields, rows


def write_csv(path, fields, rows):
    base.write_csv(path, fields, rows)


def is_topup(row):
    return str(row.get("topup_delta_flag", "")).lower() == "true"


def decorate_anchor(row, target):
    out = dict(row)
    out["anchored_15k_status"] = "anchor_from_frozen_15k"
    out["topup_delta_flag"] = "False"
    out["per_target_18k_selection_reason"] = "preserved_from_frozen_15k_anchor"
    out["not_final_order_list"] = "True"
    out["manual_review_required"] = "True"
    if target == "Ab_sdAb":
        out["complex_diagnostic_required"] = "True"
    return out


def decorate_topup(row, target):
    out = dict(row)
    out["candidate_type"] = "main_candidate"
    out["evidence_class"] = row.get("assigned_template_status", "") or row.get("assigned_expansion_role", "")
    out["quality_tier"] = row.get("candidate_quality_tier", "")
    out["local_expansion_flag"] = "True" if base.is_local(row) else "False"
    out["seed_only_flag"] = "True" if base.is_seed_only(row) else "False"
    out["glycan_status"] = "glycan_unchecked_no_explicit_review"
    out["buildability_status"] = "pending_audit"
    out["anchored_15k_status"] = "new_18k_topup"
    out["topup_delta_flag"] = "True"
    out["selection_reason"] = "selected_from_expanded_v2_per_target_18k_anchored_topup"
    out["per_target_18k_selection_reason"] = "cap_aware_topup_from_unselected_20k_mother_pool"
    out["not_final_order_list"] = "True"
    out["manual_review_required"] = "True"
    if target == "Ab_sdAb":
        out["complex_diagnostic_required"] = "True"
    return out


def make_counts_from_rows(rows):
    counts = base.make_counts()
    for row in rows:
        base.add_count(row, counts)
    return counts


def cap_violations(rows, cfg):
    violations = []
    seed = Counter(r.get("his_seed_set", "") or "NA" for r in rows)
    cluster = Counter(r.get("near_duplicate_cluster_id", "") or "NA" for r in rows)
    if seed and seed.most_common(1)[0][1] > cfg["top_seed_cap"]:
        violations.append("top_seed_cap_exceeded")
    if cluster and cluster.most_common(1)[0][1] > cfg["top_cluster_cap"]:
        violations.append("top_cluster_cap_exceeded")
    checks = [
        ("seed_only_cap_exceeded", sum(1 for r in rows if base.is_seed_only(r)), cfg["seed_only_cap"]),
        ("neutral_boundary_cap_exceeded", sum(1 for r in rows if base.is_neutral_boundary(r)), cfg["neutral_boundary_cap"]),
        ("local_expansion_cap_exceeded", sum(1 for r in rows if base.is_local(r)), cfg["local_expansion_cap"]),
        ("four_mut_cap_exceeded", sum(1 for r in rows if base.is_four_mut(r)), cfg["four_mut_cap"]),
    ]
    if cfg.get("ay111h_cap") is not None:
        checks.append(
            ("ay111h_cap_exceeded", sum(1 for r in rows if base.contains_ay111h(r)), cfg["ay111h_cap"])
        )
    for name, value, cap in checks:
        if value > cap:
            violations.append(name)
    return violations


def selection_score(row, target):
    score = base.selection_score(row, target)
    quality = row.get("candidate_quality_tier", "")
    tier = row.get("template_match_tier", "")
    source = row.get("source_pool", "")
    if quality in ("high_support", "supported"):
        score += 120
    if tier in ("A_strong", "B_medium", "local_expansion_A"):
        score += 80
    if source in ("stage2a_candidate_list", "sdab_recovery_passlike_supplement"):
        score += 40
    if base.is_neutral_boundary(row):
        score -= 120
    if base.is_seed_only(row):
        score -= 80
    return score


def select_target(target, source_rows, anchor_rows, parents):
    cfg = TARGET_CONFIG[target]
    anchors = [decorate_anchor(r, target) for r in anchor_rows if r.get("target") == target]
    selected_keys = set()
    selected_ids = set()
    selected_hashes = set()
    for row in anchors:
        key = row_key(row)
        if key:
            selected_keys.add(key)
        if row.get("variant_id"):
            selected_ids.add(row.get("variant_id"))
        if row.get("canonical_sequence_hash_full"):
            selected_hashes.add(row.get("canonical_sequence_hash_full"))

    rejection = Counter()
    eligible = []
    for row in source_rows:
        if row.get("target") != target:
            continue
        key = row_key(row)
        if (
            (key and key in selected_keys)
            or (row.get("variant_id") and row.get("variant_id") in selected_ids)
            or (
                row.get("canonical_sequence_hash_full")
                and row.get("canonical_sequence_hash_full") in selected_hashes
            )
        ):
            rejection["already_in_15k_anchor_or_duplicate"] += 1
            continue
        if cfg.get("exclude_ag102h_av105h_main") and base.contains_ag102h_av105h(row):
            rejection["ag102h_av105h_main_excluded"] += 1
            continue
        failures = base.row_hard_failures(row, parents)
        if failures:
            rejection["topup_hard_filter_rejected"] += 1
            continue
        eligible.append(row)

    eligible.sort(key=lambda r: (-selection_score(r, target), r.get("variant_id", "")))
    counts = make_counts_from_rows(anchors)
    topups = []
    for row in eligible:
        if len(anchors) + len(topups) >= TARGET_SIZE:
            break
        ok, reason = base.can_add(row, counts, cfg)
        if not ok:
            rejection[reason] += 1
            continue
        decorated = decorate_topup(row, target)
        key = row_key(decorated)
        if (
            (key and key in selected_keys)
            or (decorated.get("variant_id") and decorated.get("variant_id") in selected_ids)
            or (
                decorated.get("canonical_sequence_hash_full")
                and decorated.get("canonical_sequence_hash_full") in selected_hashes
            )
        ):
            rejection["topup_duplicate_rejected"] += 1
            continue
        topups.append(decorated)
        base.add_count(decorated, counts)
        if key:
            selected_keys.add(key)
        if decorated.get("variant_id"):
            selected_ids.add(decorated.get("variant_id"))
        if decorated.get("canonical_sequence_hash_full"):
            selected_hashes.add(decorated.get("canonical_sequence_hash_full"))

    final_rows = anchors + topups
    final_rows.sort(
        key=lambda r: (
            0 if str(r.get("topup_delta_flag", "")).lower() != "true" else 1,
            0 if r.get("candidate_type") == "main_candidate" else 1,
            -selection_score(r, target),
            r.get("variant_id", ""),
        )
    )
    short = cfg["short_name"]
    for idx, row in enumerate(final_rows, 1):
        row["per_target_18k_row_id"] = "%s_18K_%05d" % (short, idx)
        row["draft_status"] = "draft_not_final_order_list"
    return final_rows, topups, rejection, eligible


def audit_target(target, rows, parents):
    total = len(rows)
    main = [r for r in rows if r.get("candidate_type") == "main_candidate"]
    controls = [r for r in rows if r.get("candidate_type") == "control_anchor"]
    topups = [r for r in rows if is_topup(r)]
    failures = []
    for row in rows:
        fail = base.row_hard_failures(row, parents)
        if fail:
            failures.append({
                "target": target,
                "variant_id": row.get("variant_id", ""),
                "candidate_type": row.get("candidate_type", ""),
                "failure": ";".join(fail),
            })

    counters = {}
    for field in ("variant_id", "canonical_unique_key", "canonical_sequence_hash_full"):
        c = Counter(r.get(field, "") for r in rows if r.get(field, ""))
        counters[field] = sum(1 for v in c.values() if v > 1)

    seed = Counter(r.get("his_seed_set", "") or "NA" for r in rows)
    cluster = Counter(r.get("near_duplicate_cluster_id", "") or "NA" for r in rows)
    summary = {
        "target": target,
        "selected_count": total,
        "main_count": len(main),
        "control_anchor_count": len(controls),
        "topup_count": len(topups),
        "canonical_duplicate_count": counters["canonical_unique_key"],
        "sequence_duplicate_count": counters["canonical_sequence_hash_full"],
        "variant_id_duplicate_count": counters["variant_id"],
        "hard_failure_count": len(failures),
        "top_seed": "",
        "top_seed_count": 0,
        "top_seed_fraction": "0",
        "top_cluster": "",
        "top_cluster_count": 0,
        "top_cluster_fraction": "0",
        "four_mut_count": sum(1 for r in rows if base.is_four_mut(r)),
        "four_mut_fraction": "%.6f" % (sum(1 for r in rows if base.is_four_mut(r)) / total if total else 0),
        "local_expansion_count": sum(1 for r in rows if base.is_local(r)),
        "local_expansion_fraction": "%.6f" % (sum(1 for r in rows if base.is_local(r)) / total if total else 0),
        "seed_only_count": sum(1 for r in rows if base.is_seed_only(r)),
        "seed_only_fraction": "%.6f" % (sum(1 for r in rows if base.is_seed_only(r)) / total if total else 0),
        "neutral_boundary_or_high_risk_count": sum(1 for r in rows if base.is_neutral_boundary(r)),
        "neutral_boundary_or_high_risk_fraction": "%.6f" % (sum(1 for r in rows if base.is_neutral_boundary(r)) / total if total else 0),
        "glycan_low_risk_claim_count": sum(1 for r in rows if "low" in str(r.get("glycan_status", "")).lower()),
        "sdAb_complex_diagnostic_required_count": "",
        "AG102H_AV105H_main_count": "",
        "AY111H_containing_count": "",
        "AY111H_containing_fraction": "",
    }
    if seed:
        summary["top_seed"], summary["top_seed_count"] = seed.most_common(1)[0]
        summary["top_seed_fraction"] = "%.6f" % (summary["top_seed_count"] / total)
    if cluster:
        summary["top_cluster"], summary["top_cluster_count"] = cluster.most_common(1)[0]
        summary["top_cluster_fraction"] = "%.6f" % (summary["top_cluster_count"] / total)
    if target == "Ab_sdAb":
        diag = sum(1 for r in rows if base.truthy(r.get("complex_diagnostic_required", "")))
        summary["sdAb_complex_diagnostic_required_count"] = diag
        summary["AG102H_AV105H_main_count"] = sum(1 for r in main if base.contains_ag102h_av105h(r))
        ay = sum(1 for r in rows if base.contains_ay111h(r))
        summary["AY111H_containing_count"] = ay
        summary["AY111H_containing_fraction"] = "%.6f" % (ay / total if total else 0)

    verdict = "PASS"
    reasons = []
    if total != TARGET_SIZE:
        verdict = "PATCH" if total >= 17000 else "FAIL"
        reasons.append("selected_count_not_18000")
    if len(main) != TARGET_SIZE - CONTROL_ANCHOR_COUNT:
        verdict = "PATCH" if verdict == "PASS" else verdict
        reasons.append("main_count_not_17980")
    if len(controls) != CONTROL_ANCHOR_COUNT:
        verdict = "PATCH" if verdict == "PASS" else verdict
        reasons.append("control_anchor_count_not_20")
    if counters["canonical_unique_key"] or counters["canonical_sequence_hash_full"] or counters["variant_id"]:
        verdict = "FAIL"
        reasons.append("duplicate_detected")
    if failures:
        verdict = "FAIL"
        reasons.append("hard_failure_detected")
    cap_fail = cap_violations(rows, TARGET_CONFIG[target])
    if cap_fail:
        verdict = "PATCH" if verdict == "PASS" else verdict
        reasons.extend(cap_fail)
    if target == "Ab_sdAb":
        if summary["sdAb_complex_diagnostic_required_count"] != total:
            verdict = "FAIL"
            reasons.append("sdAb_complex_diagnostic_not_100pct")
        if summary["AG102H_AV105H_main_count"]:
            verdict = "FAIL"
            reasons.append("AG102H_AV105H_main_present")
    if summary["glycan_low_risk_claim_count"]:
        verdict = "FAIL"
        reasons.append("glycan_low_risk_claim_present")
    summary["verdict"] = verdict
    summary["reasons"] = ";".join(reasons) if reasons else "all_per_target_18k_gates_pass"
    return summary, failures


def topup_delta_summary(target, topups, rejection):
    row = {
        "target": target,
        "topup_count": len(topups),
        "topup_from_unselected_mother_pool_count": len(topups),
        "topup_duplicate_rejected_count": rejection.get("topup_duplicate_rejected", 0),
        "topup_cap_rejected_count": sum(
            v for k, v in rejection.items()
            if k.endswith("_cap_full") or k in ("seed_cap_full", "cluster_cap_full")
        ),
        "topup_hard_filter_rejected_count": rejection.get("topup_hard_filter_rejected", 0),
        "topup_seed_cap_rejected_count": rejection.get("seed_cap_full", 0),
        "topup_cluster_cap_rejected_count": rejection.get("cluster_cap_full", 0),
        "topup_local_expansion_cap_rejected_count": rejection.get("local_expansion_cap_full", 0),
        "topup_seed_only_cap_rejected_count": rejection.get("seed_only_cap_full", 0),
        "topup_neutral_boundary_or_high_risk_count": sum(1 for r in topups if base.is_neutral_boundary(r)),
        "topup_4mut_count": sum(1 for r in topups if base.is_four_mut(r)),
        "topup_source_pool_distribution": ";".join(
            "%s:%s" % kv for kv in Counter(r.get("source_pool", "") or "NA" for r in topups).most_common()
        ),
        "topup_evidence_class_distribution": ";".join(
            "%s:%s" % kv for kv in Counter(r.get("evidence_class", "") or "NA" for r in topups).most_common()
        ),
        "topup_quality_tier_distribution": ";".join(
            "%s:%s" % kv for kv in Counter(r.get("quality_tier", "") or "NA" for r in topups).most_common()
        ),
        "topup_template_match_tier_distribution": ";".join(
            "%s:%s" % kv for kv in Counter(r.get("template_match_tier", "") or "NA" for r in topups).most_common()
        ),
    }
    return row


def md_table(rows, columns):
    return base.md_table(rows, columns)


def write_target_audit(path, target, rows, topups, summary, failures, rejection):
    with open(path, "w") as h:
        h.write("# %s per-target 18K anchored top-up audit\n\n" % TARGET_CONFIG[target]["short_name"])
        h.write("Verdict: `%s`.\n\n" % summary["verdict"])
        h.write("This preserves the frozen 15K draft and adds main candidates from the same target's frozen 20K mother pool.\n")
        h.write("No new variants, backfill, AF3, SimpleFold, PyRosetta, FoldX, MD, glycan modeling, or wet-lab order generation were run.\n\n")
        cols = [
            "target", "selected_count", "main_count", "control_anchor_count", "topup_count",
            "canonical_duplicate_count", "sequence_duplicate_count", "variant_id_duplicate_count",
            "hard_failure_count", "top_seed", "top_seed_count", "top_seed_fraction",
            "top_cluster", "top_cluster_count", "top_cluster_fraction", "four_mut_count",
            "four_mut_fraction", "local_expansion_count", "local_expansion_fraction",
            "seed_only_count", "seed_only_fraction", "neutral_boundary_or_high_risk_count",
            "neutral_boundary_or_high_risk_fraction", "glycan_low_risk_claim_count",
            "sdAb_complex_diagnostic_required_count", "AG102H_AV105H_main_count",
            "AY111H_containing_count", "AY111H_containing_fraction", "verdict", "reasons",
        ]
        h.write("## Final 18K Summary\n\n")
        h.write(md_table([summary], cols))
        h.write("\n\n## Top-up Delta Summary\n\n")
        delta = topup_delta_summary(target, topups, rejection)
        delta_cols = [
            "target", "topup_count", "topup_from_unselected_mother_pool_count",
            "topup_duplicate_rejected_count", "topup_cap_rejected_count",
            "topup_hard_filter_rejected_count", "topup_seed_cap_rejected_count",
            "topup_cluster_cap_rejected_count", "topup_local_expansion_cap_rejected_count",
            "topup_seed_only_cap_rejected_count", "topup_neutral_boundary_or_high_risk_count",
            "topup_4mut_count", "topup_source_pool_distribution",
            "topup_evidence_class_distribution", "topup_quality_tier_distribution",
            "topup_template_match_tier_distribution",
        ]
        h.write(md_table([delta], delta_cols))
        h.write("\n\n## Final Evidence Class Distribution\n\n")
        h.write(md_table(base.by_count(rows, ["candidate_type", "evidence_class", "quality_tier", "template_match_tier"])[:50], ["candidate_type", "evidence_class", "quality_tier", "template_match_tier", "count"]))
        h.write("\n\n## Top-up Source Distribution\n\n")
        h.write(md_table(base.by_count(topups, ["source_pool", "quality_tier", "template_match_tier"])[:50], ["source_pool", "quality_tier", "template_match_tier", "count"]))
        h.write("\n\n## Final Seed Distribution\n\n")
        h.write(md_table(base.top_rows(rows, "his_seed_set", 25), ["his_seed_set", "count", "fraction"]))
        h.write("\n\n## Final Cluster Distribution\n\n")
        h.write(md_table(base.top_rows(rows, "near_duplicate_cluster_id", 25), ["near_duplicate_cluster_id", "count", "fraction"]))
        h.write("\n\n## Rejection / Cap Pressure Summary\n\n")
        rej = [{"reject_reason": k, "count": v} for k, v in rejection.most_common()]
        h.write(md_table(rej, ["reject_reason", "count"]) if rej else "No rejected rows under target construction.\n")
        h.write("\n\n## Locked Items\n\n")
        h.write("- These 18K files are pre-synthesis candidate pools, not final order lists.\n")
        h.write("- Glycan remains unchecked; no low-risk glycan claim is made.\n")
        if target == "Ab_sdAb":
            h.write("- sdAb remains secondary / complex-diagnostic-required for every row.\n")
            h.write("- AG102H / AV105H are excluded from main candidates.\n")


def write_comparison(path, summaries, delta_rows):
    overall = "PASS_PER_TARGET_18K_ANCHORED_TOPUP_BUILT__FINAL_ORDER_LOCKED"
    if any(s["verdict"] != "PASS" for s in summaries):
        overall = "PATCH_OR_FAIL_REVIEW_REQUIRED__FINAL_ORDER_LOCKED"
    with open(path, "w") as h:
        h.write("# Per-target 18K anchored top-up candidate draft comparison\n\n")
        h.write("Verdict: `%s`.\n\n" % overall)
        h.write("Each antibody now has its own 18K pre-synthesis candidate draft. The existing 15K drafts were preserved as anchors, and the added rows come only from the same target's frozen 20K expanded Tier2 mother pool.\n\n")
        cols = [
            "target", "selected_count", "main_count", "control_anchor_count", "topup_count",
            "top_seed", "top_seed_fraction", "top_cluster_fraction", "four_mut_fraction",
            "local_expansion_fraction", "seed_only_fraction", "neutral_boundary_or_high_risk_fraction",
            "sdAb_complex_diagnostic_required_count", "AG102H_AV105H_main_count",
            "glycan_low_risk_claim_count", "verdict", "reasons",
        ]
        h.write("## Final 18K Summary\n\n")
        h.write(md_table(summaries, cols))
        h.write("\n\n## Top-up Delta Summary\n\n")
        delta_cols = [
            "target", "topup_count", "topup_duplicate_rejected_count", "topup_cap_rejected_count",
            "topup_hard_filter_rejected_count", "topup_seed_cap_rejected_count",
            "topup_cluster_cap_rejected_count", "topup_local_expansion_cap_rejected_count",
            "topup_seed_only_cap_rejected_count", "topup_neutral_boundary_or_high_risk_count",
            "topup_4mut_count", "topup_source_pool_distribution",
            "topup_quality_tier_distribution", "topup_template_match_tier_distribution",
        ]
        h.write(md_table(delta_rows, delta_cols))
        h.write("\n\n## Interpretation\n\n")
        h.write("- 18K expansion used anchored top-up from the frozen 20K mother pools; no new backfill or generation was introduced.\n")
        h.write("- 1E62 top-up is expected to increase local-expansion and 4-mut exposure because the remaining admissible 1E62 space is dominated by constrained local expansion.\n")
        h.write("- sdAb remains secondary / exploratory and complex-diagnostic-required; the larger pool does not upgrade it to structure-confirmed.\n")
        h.write("- Glycan remains unchecked and must not be reported as low risk.\n")
        h.write("- Final synthesis / wet-lab order generation remains locked pending vendor-level buildability and order-format checks.\n")
        h.write("\n## Output Files\n\n")
        h.write("- `final_candidate_pool_1E62_18k_draft.csv`\n")
        h.write("- `final_candidate_pool_sdAb_18k_draft.csv`\n")
        h.write("- `final_candidate_pool_1E62_18k_audit.md`\n")
        h.write("- `final_candidate_pool_sdAb_18k_audit.md`\n")
        h.write("- `per_target_18k_topup_delta_summary.csv`\n")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    source_fields, source_rows = read_rows(EXPANDED_V2)
    parents = base.derive_parents(source_rows)
    anchor_fields = []
    anchor_rows_by_target = {}
    for target, path in ANCHOR_15K.items():
        fields, rows = read_rows(path)
        anchor_fields.extend(fields)
        anchor_rows_by_target[target] = rows

    required = [
        "variant_id", "target", "sequence", "mutation_list", "his_seed_set",
        "near_duplicate_cluster_id", "candidate_type", "evidence_class",
        "quality_tier", "template_match_tier", "source_pool",
        "local_expansion_flag", "seed_only_flag", "complex_diagnostic_required",
        "glycan_status", "buildability_status", "selection_reason",
        "anchored_15k_status", "topup_delta_flag", "per_target_18k_selection_reason",
        "per_target_18k_row_id", "draft_status", "not_final_order_list",
        "manual_review_required",
    ]
    output_fields = []
    for field in required + anchor_fields + source_fields:
        if field not in output_fields:
            output_fields.append(field)

    summaries = []
    delta_rows = []
    all_failures = []
    for target in TARGETS:
        rows, topups, rejection, _eligible = select_target(
            target,
            source_rows,
            anchor_rows_by_target[target],
            parents,
        )
        summary, failures = audit_target(target, rows, parents)
        summaries.append(summary)
        all_failures.extend(failures)
        delta_rows.append(topup_delta_summary(target, topups, rejection))
        short = TARGET_CONFIG[target]["short_name"]
        write_csv(
            os.path.join(OUT_DIR, "final_candidate_pool_%s_18k_draft.csv" % short),
            output_fields,
            rows,
        )
        write_csv(
            os.path.join(OUT_DIR, "final_candidate_pool_%s_18k_topup_delta_rows.csv" % short),
            output_fields,
            topups,
        )
        write_target_audit(
            os.path.join(OUT_DIR, "final_candidate_pool_%s_18k_audit.md" % short),
            target,
            rows,
            topups,
            summary,
            failures,
            rejection,
        )

    write_csv(
        os.path.join(OUT_DIR, "per_target_18k_audit_summary.csv"),
        list(summaries[0].keys()),
        summaries,
    )
    write_csv(
        os.path.join(OUT_DIR, "per_target_18k_topup_delta_summary.csv"),
        list(delta_rows[0].keys()),
        delta_rows,
    )
    write_csv(
        os.path.join(OUT_DIR, "per_target_18k_hard_failure_rows.csv"),
        ["target", "variant_id", "candidate_type", "failure"],
        all_failures,
    )
    manifest = [
        {"input": "expanded_v2", "path": os.path.relpath(EXPANDED_V2, ROOT), "rows": len(source_rows), "sha256": base.sha256_file(EXPANDED_V2)},
    ]
    for target, path in ANCHOR_15K.items():
        manifest.append({
            "input": "%s_15k_anchor" % target,
            "path": os.path.relpath(path, ROOT),
            "rows": len(anchor_rows_by_target[target]),
            "sha256": base.sha256_file(path),
        })
    write_csv(
        os.path.join(OUT_DIR, "per_target_18k_input_manifest.csv"),
        ["input", "path", "rows", "sha256"],
        manifest,
    )
    comparison = os.path.join(OUT_DIR, "final_candidate_pool_18k_by_target_comparison.md")
    write_comparison(comparison, summaries, delta_rows)
    with open(comparison) as src, open(TASK_REPORT, "w") as dst:
        dst.write(src.read())
    if any(s["verdict"] != "PASS" for s in summaries):
        print("PATCH_OR_FAIL_REVIEW_REQUIRED__FINAL_ORDER_LOCKED")
    else:
        print("PASS_PER_TARGET_18K_ANCHORED_TOPUP_BUILT__FINAL_ORDER_LOCKED")
    print(os.path.relpath(OUT_DIR, ROOT))


if __name__ == "__main__":
    main()
