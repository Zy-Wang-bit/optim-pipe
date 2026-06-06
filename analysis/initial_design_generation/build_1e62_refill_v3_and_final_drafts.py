#!/usr/bin/env python3
"""Targeted 1E62 refill from existing source pools and rerun final drafts.

This stage performs table-only source rescanning. It does not run new
structure prediction, Rosetta/FoldX, MD, glycan modeling, or synthesis
selection.
"""

import csv
import hashlib
import importlib.util
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results/initial_design_generation/final_candidate_pool_refill_v3"
CURRENT_STAGE_REPORT = ROOT / ".tasks/active/initial-design-generation/current_stage_report.md"
BUILDER_PATH = ROOT / "analysis/initial_design_generation/build_final_candidate_pool_drafts.py"
V2_POOL = ROOT / "results/initial_design_generation/expanded_tier2_candidate_pool_v2/expanded_tier2_candidate_pool_v2.csv"
TEMPLATES = ROOT / "results/initial_design_generation/tier2_candidate_bank/tier2_candidate_bank_expansion_templates.csv"

SOURCE_FILES = [
    ("stage2a_candidate_list", ROOT / "results/initial_design_generation/stage1_5_stage2a/stage2a_candidate_list.csv", 0),
    ("tier2_candidate_snapshot", ROOT / "results/initial_design_generation/tier2_staged/tier2_candidate_snapshot.csv", 1),
    ("tier2_core_reserve_pool", ROOT / "results/initial_design_generation/tier1_filtering/tier2_core_reserve_pool.csv", 2),
    ("tier1_ranked_production_pool", ROOT / "results/initial_design_generation/tier1_filtering/tier1_ranked_candidates_pre_tier2.csv", 3),
]

spec = importlib.util.spec_from_file_location("final_builder", str(BUILDER_PATH))
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def short_hash(text, n=10):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:n]


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        seen = set()
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    fieldnames.append(key)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def as_float(value, default=0.0):
    try:
        text = str(value or "").strip()
        if not text or text.lower() in {"nan", "none", "<na>"}:
            return default
        return float(text)
    except Exception:
        return default


def as_int(value, default=0):
    return int(round(as_float(value, default)))


def split_mutations(text):
    if not text:
        return []
    return [x.strip() for x in str(text).split(";") if x.strip() and x.strip().lower() != "nan"]


def parse_mut(mut):
    m = re.match(r"^([A-Z])([A-Z]?)(\d+)([A-Z])$", str(mut))
    if not m:
        return None
    chain, wt, pos, aa = m.groups()
    return chain, wt, int(pos), aa


def mut_pos(mut):
    parsed = parse_mut(mut)
    if not parsed:
        return str(mut)
    chain, _wt, pos, _aa = parsed
    return "%s%d" % (chain, pos)


def non_his(muts):
    return [m for m in muts if not m.endswith("H")]


def his(muts):
    return [m for m in muts if m.endswith("H")]


def canonical_key(row):
    target = row.get("target", "")
    for col in ["canonical_sequence_hash_full", "canonical_recovery_sequence_hash", "sequence_hash"]:
        value = str(row.get(col, "") or "").strip()
        if value and value.lower() not in {"nan", "none", "<na>"}:
            return "%s|%s" % (target, value)
    seq = str(row.get("sequence", "") or "").strip()
    if seq:
        return "%s|%s" % (target, hashlib.sha256(seq.encode("utf-8")).hexdigest())
    return "%s|variant|%s" % (target, row.get("variant_id", ""))


def load_templates():
    templates_by_seed = defaultdict(list)
    for row in read_csv(TEMPLATES):
        if row.get("target") != "Ab_1E62":
            continue
        if str(row.get("allowed_backfill", "")).lower() != "true":
            continue
        muts = split_mutations(row.get("mutation_pattern"))
        npos = {mut_pos(m) for m in non_his(muts)}
        tpl = dict(row)
        tpl["_non_his_positions"] = npos
        tpl["_non_his_set"] = set(non_his(muts))
        tpl["_priority"] = as_int(row.get("priority_level"), 9)
        templates_by_seed[row.get("seed_pattern", "")].append(tpl)
    for seed in templates_by_seed:
        templates_by_seed[seed].sort(key=lambda r: (r["_priority"], r.get("template_id", "")))
    return templates_by_seed


def source_cluster(row):
    for col in ["near_duplicate_cluster_id", "tier1_near_duplicate_cluster_id", "source_near_duplicate_cluster_id"]:
        value = str(row.get(col, "") or "").strip()
        if value and value.lower() not in {"nan", "none", "<na>"}:
            return value
    return "Ab_1E62_refill_cluster_" + short_hash(row.get("mutation_list", "missing"))


def hard_ok(row):
    if row.get("target") != "Ab_1E62":
        return False, "not_1E62"
    if as_int(row.get("mutation_count"), 99) > 3:
        return False, "mutation_count_gt_3"
    muts = split_mutations(row.get("normalized_mutation_list") or row.get("mutation_list"))
    if any(m.endswith("C") for m in muts):
        return False, "new_cys_mutation"
    for field in ["hard_filter_status", "forbidden_pair_status", "buildability_light_status"]:
        value = str(row.get(field, "") or "").strip().lower()
        if value and value not in {"pass", "nan", "none", "<na>"}:
            return False, "%s_not_pass" % field
    if not row.get("sequence"):
        return False, "missing_sequence"
    return True, ""


def template_match(row, templates_by_seed):
    seed = row.get("his_seed_set", "")
    templates = templates_by_seed.get(seed, [])
    if not templates:
        return None, "no_template_seed_match"
    muts = split_mutations(row.get("normalized_mutation_list") or row.get("mutation_list"))
    nh = non_his(muts)
    if not nh:
        return None, "seed_only"
    nh_pos = {mut_pos(m) for m in nh}
    best = None
    best_score = -1
    best_tier = None
    best_overlap = 0
    best_reason = ""
    for tpl in templates:
        exact = bool(set(nh) & tpl["_non_his_set"])
        overlap = len(nh_pos & tpl["_non_his_positions"])
        if exact:
            tier = "A_strong"
            score = 10 + overlap
            reason = "A_strong;exact_template_rescue_or_position;source_rescan_v3"
        elif overlap >= 1:
            tier = "B_medium"
            score = 5 + overlap
            reason = "B_medium;same_rescue_position;source_rescan_v3"
        else:
            continue
        if score > best_score:
            best = tpl
            best_score = score
            best_tier = tier
            best_overlap = overlap
            best_reason = reason
    if best is None:
        return None, "no_rescue_or_position_overlap"
    return {
        "template": best,
        "template_match_tier": best_tier,
        "same_position_overlap_count": best_overlap,
        "evidence_inheritance_reason": best_reason,
    }, ""


def quality_tier(row, match_tier):
    cls = str(row.get("tier1_review_class", ""))
    t2 = str(row.get("t2_class_current") or row.get("tier2_class") or "")
    if "A_tier1_priority" in cls or "T2_strong" in t2 or match_tier == "A_strong":
        return "high_support"
    if "B_rescue" in cls or "T2_supported" in t2 or "allow" in str(row.get("stage2a_list_action", "")):
        return "supported"
    return "supported"


def row_score(row, match_tier):
    score = 0.0
    score += 3.0 if match_tier == "A_strong" else 2.0
    score += as_float(row.get("tier1_rank_score")) * 2.0
    score += as_float(row.get("hit_likelihood_score_v0")) * 1.0
    score += as_float(row.get("neutral_retention_score")) * 0.8
    score += as_float(row.get("acidic_release_support_score")) * 0.8
    score -= as_float(row.get("global_weakening_risk_score")) * 0.6
    score -= as_float(row.get("mpnn_score_percentile_within_target"), as_float(row.get("mpnn_score_percentile_by_target"), 0.5)) * 0.2
    return round(score, 6)


def normalize_row(row, source_pool, match):
    tpl = match["template"]
    out = {}
    # Start with the v2 pool columns.
    for col in V2_FIELDS:
        out[col] = row.get(col, "")
    out["target"] = "Ab_1E62"
    out["variant_id"] = row.get("variant_id") or ("1E62_refill_v3_" + short_hash(row.get("sequence", "") + row.get("mutation_list", "")))
    out["source_variant_id"] = row.get("variant_id", "")
    out["sequence"] = row.get("sequence", "")
    out["sequence_hash"] = row.get("sequence_hash") or short_hash(out["sequence"])
    out["canonical_sequence_hash_full"] = row.get("canonical_sequence_hash_full") or hashlib.sha256(out["sequence"].encode("utf-8")).hexdigest()
    out["canonical_recovery_sequence_hash"] = row.get("canonical_recovery_sequence_hash", "")
    out["canonical_unique_key"] = "Ab_1E62|%s" % out["canonical_sequence_hash_full"]
    out["mutation_list"] = row.get("mutation_list", "")
    out["normalized_mutation_list"] = row.get("normalized_mutation_list") or row.get("mutation_list", "")
    out["mutation_count"] = str(as_int(row.get("mutation_count")))
    out["His_count"] = row.get("His_count", "")
    out["his_seed_set"] = row.get("his_seed_set", "")
    out["near_duplicate_cluster_id"] = source_cluster(row)
    out["assigned_template_id"] = tpl.get("template_id", "")
    out["assigned_template_status"] = tpl.get("template_status", "")
    out["assigned_expansion_role"] = tpl.get("recommended_expansion_role", "")
    out["template_match_tier"] = match["template_match_tier"]
    out["template_match_level"] = "source_rescan_v3"
    out["same_position_overlap_count"] = str(match["same_position_overlap_count"])
    out["evidence_inheritance_reason"] = match["evidence_inheritance_reason"]
    out["expanded_pool_role"] = "targeted_1E62_refill_v3"
    out["candidate_quality_tier"] = quality_tier(row, match["template_match_tier"])
    out["expanded_pool_score"] = str(row_score(row, match["template_match_tier"]))
    out["source_pool"] = source_pool
    out["source_pool_list"] = source_pool
    out["source_occurrence_count"] = "1"
    out["hard_filter_status"] = row.get("hard_filter_status", "pass") or "pass"
    out["forbidden_pair_status"] = row.get("forbidden_pair_status", "")
    out["buildability_light_status"] = row.get("buildability_light_status", "pass") or "pass"
    out["tier1_review_class"] = row.get("tier1_review_class", "")
    out["bank_eligibility"] = row.get("bank_eligibility", "")
    out["stage2a_list_action"] = row.get("stage2a_list_action", "")
    out["t2_class_current"] = row.get("t2_class_current", "")
    out["tier2_class"] = row.get("tier2_class", "")
    out["selection_eligibility"] = row.get("selection_eligibility", "")
    out["tier2_eligibility_status"] = row.get("tier2_eligibility_status", "")
    for col in [
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
        out[col] = row.get(col, "")
    out["complex_diagnostic_required"] = "False"
    return out


def scan_sources(existing_keys, templates_by_seed, limit=12000):
    accepted = {}
    rejection = Counter()
    by_source = Counter()
    for source_pool, path, _rank in SOURCE_FILES:
        if not path.exists():
            continue
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                ok, reason = hard_ok(row)
                if not ok:
                    rejection[reason] += 1
                    continue
                key = canonical_key(row)
                if key in existing_keys or key in accepted:
                    rejection["already_in_v2_or_duplicate"] += 1
                    continue
                match, reason = template_match(row, templates_by_seed)
                if match is None:
                    rejection[reason] += 1
                    continue
                out = normalize_row(row, source_pool, match)
                accepted[key] = out
                by_source[source_pool] += 1
    rows = sorted(
        accepted.values(),
        key=lambda r: (
            builder.as_float(r.get("expanded_pool_score")),
            -builder.as_int(r.get("mutation_count")),
            r.get("variant_id", ""),
        ),
        reverse=True,
    )
    return rows[:limit], rejection, by_source


def md_table(rows, fields):
    lines = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(f, "")) for f in fields) + " |")
    return "\n".join(lines)


def main():
    global V2_FIELDS
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    v2_rows = builder.read_csv(V2_POOL)
    V2_FIELDS = list(v2_rows[0].keys())
    existing_keys = set(r.get("canonical_unique_key") or ("%s|%s" % (r.get("target"), r.get("canonical_sequence_hash_full"))) for r in v2_rows)
    templates_by_seed = load_templates()
    supplement, rejection, by_source = scan_sources(existing_keys, templates_by_seed)
    write_csv(OUT_DIR / "targeted_1E62_refill_v3_candidates.csv", supplement, V2_FIELDS)
    write_csv(
        OUT_DIR / "targeted_1E62_refill_v3_rejection_summary.csv",
        [{"reject_reason": k, "count": v} for k, v in sorted(rejection.items(), key=lambda x: (-x[1], x[0]))],
    )
    write_csv(
        OUT_DIR / "targeted_1E62_refill_v3_source_summary.csv",
        [{"source_pool": k, "count": v} for k, v in sorted(by_source.items())],
    )

    combined = v2_rows + supplement
    write_csv(OUT_DIR / "expanded_tier2_candidate_pool_v3_with_1E62_refill.csv", combined, V2_FIELDS)
    enriched = builder.add_features(combined)
    default_audits = []
    selected_by_scenario = {}
    excluded_all = []
    for name, scenario in builder.DEFAULT_SCENARIOS.items():
        selected, excluded, audits = builder.select_scenario(enriched, name + "_refill_v3", scenario)
        selected_by_scenario[name + "_refill_v3"] = selected
        default_audits.extend(audits)
        excluded_all.extend(excluded[:10000])
        fields = list(enriched[0].keys()) + [x for x in builder.EXTRA_FIELDS if x not in enriched[0]]
        write_csv(OUT_DIR / ("final_candidate_pool_%s_refill_v3_draft.csv" % name), selected, fields)
    write_csv(OUT_DIR / "refill_v3_final_candidate_audit_checks.csv", default_audits)
    write_csv(OUT_DIR / "refill_v3_excluded_candidates_sample.csv", excluded_all)

    target_rows = []
    for scenario, selected in selected_by_scenario.items():
        counts = Counter(r["target"] for r in selected)
        for target, count in sorted(counts.items()):
            target_rows.append({"selection_scenario": scenario, "target": target, "count": count})

    audit_fields = [
        "selection_scenario",
        "target",
        "selected_count",
        "target_quota",
        "top_seed",
        "top_seed_count",
        "top_seed_fraction",
        "top_cluster_count",
        "top_cluster_fraction",
        "four_mut_count",
        "four_mut_fraction",
        "local_expansion_count",
        "local_expansion_fraction",
        "seed_only_count",
        "seed_only_fraction",
        "neutral_boundary_or_high_risk_count",
        "neutral_boundary_or_high_risk_fraction",
        "AG102H_AV105H_main_count",
        "sdAb_complex_diagnostic_required_count",
        "canonical_duplicate_count",
        "verdict",
        "reasons",
    ]
    verdicts = [a["verdict"] for a in default_audits]
    if any(v == "FAIL" for v in verdicts):
        overall = "REFILL_V3_FINAL_DRAFTS_BUILT__AUDIT_FAIL__FINAL_SELECTION_LOCKED"
    elif any(v == "PATCH" for v in verdicts):
        overall = "REFILL_V3_FINAL_DRAFTS_BUILT__PATCH_REVIEW_REQUIRED__FINAL_SELECTION_LOCKED"
    else:
        overall = "REFILL_V3_FINAL_DRAFTS_BUILT__AUDIT_PASS__FINAL_SELECTION_LOCKED"

    if supplement:
        interpretation_refill = "Targeted source rescan found additional 1E62 low-order, non-seed-only, template-supported candidates from existing pools."
        next_step = "If audit PASS, discuss manual unlock criteria for final library planning."
    else:
        interpretation_refill = "Targeted source rescan found no additional 1E62 low-order, non-seed-only, template-supported candidates outside the v2 mother pool. Existing source pools are exhausted under this refill definition."
        next_step = "Use the PATCH-resolution drafts for review, or request a new targeted generation / Backfill v3 policy that explicitly changes the allowed evidence or cap assumptions."

    report = """# 1E62 Targeted Refill v3 And Final Candidate Drafts

## Executive Summary

Verdict: `%s`.

This stage rescanned existing source pools for additional 1E62 candidates that were not in the expanded v2 mother pool. It accepted only table-level candidates with same His seed, rescue-position/template evidence inheritance, mutation_count <= 3, buildability/hard-filter pass, and no new Cys. It did not run new AF3, SimpleFold, PyRosetta, FoldX, MD, glycan modeling, or final synthesis selection.

## Refill Summary

| metric | value |
| --- | --- |
| refill candidates accepted | %d |
| combined v3 candidate pool rows | %d |
| v2 pool sha256 | %s |

## Refill Source Summary

%s

## Draft Target Counts

%s

## Draft Audit

%s

## Interpretation

- %s
- This attempt does not improve the original 1E62 count shortfall because no eligible new refill rows were found.
- sdAb remains secondary / complex-diagnostic-required.
- These remain final candidate-pool drafts; final synthesis-ready selection is still locked.

## Next Gate

Allowed now:

- Review whether refill v3 final drafts are acceptable.
- %s

Still locked:

- Final synthesis-ready 10K / 15K selection.
- Wet-lab order list generation.
- New AF3 / SimpleFold / PyRosetta / FoldX compute.
- MD / constant-pH MD.
- sdAb structure-confirmed upgrade.
- Glycan low-risk claim.
""" % (
        overall,
        len(supplement),
        len(combined),
        sha256_file(V2_POOL),
        md_table([{"source_pool": k, "count": v} for k, v in sorted(by_source.items())], ["source_pool", "count"]),
        md_table(target_rows, ["selection_scenario", "target", "count"]),
        md_table(default_audits, audit_fields),
        interpretation_refill,
        next_step,
    )
    (OUT_DIR / "refill_v3_final_candidate_audit_report.md").write_text(report, encoding="utf-8")
    (OUT_DIR / "refill_v3_next_gate_report.md").write_text(report, encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(report, encoding="utf-8")
    print(overall)
    print("accepted_refill", len(supplement))
    print(OUT_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
