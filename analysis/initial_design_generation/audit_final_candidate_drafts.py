#!/usr/bin/env python3
"""Audit final candidate-pool drafts without unlocking final synthesis.

Inputs are the current refill-v3 10K/15K draft tables. This script only
summarizes and audits existing tables; it does not generate variants or run
structure/MD/glycan compute.
"""

import csv
import hashlib
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results/initial_design_generation/final_candidate_draft_audit"
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"

DRAFTS = {
    "10k": ROOT / "results/initial_design_generation/final_candidate_pool_refill_v3/final_candidate_pool_10k_refill_v3_draft.csv",
    "15k": ROOT / "results/initial_design_generation/final_candidate_pool_refill_v3/final_candidate_pool_15k_refill_v3_draft.csv",
}
CONTROL_PANEL = ROOT / "results/initial_design_generation/final_candidate_pool_planning/final_candidate_pool_control_anchor_panel.csv"
REFILL_REPORT = ROOT / "results/initial_design_generation/final_candidate_pool_refill_v3/refill_v3_final_candidate_audit_report.md"
PATCH_RESOLUTION_REPORT = ROOT / "results/initial_design_generation/final_candidate_pool_patch_resolution/patch_resolution_audit_report.md"


NOMINAL = {
    "10k": {"total": 10000, "targets": {"Ab_1E62": 6000, "Ab_sdAb": 4000}},
    "15k": {"total": 15000, "targets": {"Ab_1E62": 8500, "Ab_sdAb": 6500}},
}


def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


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
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def md_table(rows, fields):
    lines = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(field, "")) for field in fields) + " |")
    return "\n".join(lines)


def as_bool(value):
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def as_int(value, default=0):
    try:
        text = str(value or "").strip()
        if not text or text.lower() in {"nan", "none", "<na>"}:
            return default
        return int(round(float(text)))
    except Exception:
        return default


def fraction(count, total):
    if not total:
        return 0.0
    return round(count / total, 6)


def norm(value, fallback):
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none", "<na>"}:
        return fallback
    return text


def top_counter(counter):
    if not counter:
        return "", 0
    return counter.most_common(1)[0]


def contains_ag102h_av105h(row):
    text = (row.get("his_seed_set", "") or "") + ";" + (row.get("mutation_list", "") or "")
    return "AG102H" in text or "AV105H" in text


def contains_ay111h(row):
    text = (row.get("his_seed_set", "") or "") + ";" + (row.get("mutation_list", "") or "")
    return "AY111H" in text


def is_local(row):
    return str(row.get("is_local_expansion", "")).lower() == "true" or row.get("source_pool", "").startswith("constrained")


def is_seed_only(row):
    return str(row.get("is_seed_only", "")).lower() == "true" or row.get("template_match_tier") == "C_seed_only"


def is_four_mut(row):
    return as_int(row.get("mutation_count")) >= 4


def is_neutral_boundary(row):
    return str(row.get("is_neutral_boundary_or_high_risk", "")).lower() == "true" or row.get("candidate_quality_tier") == "neutral_boundary_or_high_risk"


def hard_status_bad(row):
    for field in ["hard_filter_status", "forbidden_pair_status", "buildability_light_status"]:
        value = str(row.get(field, "") or "").strip().lower()
        if value and value not in {"pass", "nan", "none", "<na>"}:
            return True
    return False


def summarize_group(rows, scenario, target):
    subset = [r for r in rows if r.get("target") == target]
    total = len(subset)
    seed = Counter(norm(r.get("his_seed_set"), "missing_seed") for r in subset)
    cluster = Counter(norm(r.get("near_duplicate_cluster_id"), "missing_cluster") for r in subset)
    top_seed, top_seed_count = top_counter(seed)
    top_cluster, top_cluster_count = top_counter(cluster)
    four = sum(is_four_mut(r) for r in subset)
    local = sum(is_local(r) for r in subset)
    seed_only = sum(is_seed_only(r) for r in subset)
    neutral = sum(is_neutral_boundary(r) for r in subset)
    agav = sum(contains_ag102h_av105h(r) for r in subset)
    ay111h = sum(contains_ay111h(r) for r in subset)
    sdab_diag = sum(as_bool(r.get("complex_diagnostic_required")) for r in subset)
    build_bad = sum(hard_status_bad(r) for r in subset)
    nominal = NOMINAL[scenario]["targets"].get(target, 0)
    return {
        "scenario": scenario,
        "target": target,
        "selected_count": total,
        "nominal_quota": nominal,
        "shortfall": max(0, nominal - total),
        "top_seed": top_seed,
        "top_seed_count": top_seed_count,
        "top_seed_fraction": fraction(top_seed_count, total),
        "top_cluster": top_cluster,
        "top_cluster_count": top_cluster_count,
        "top_cluster_fraction": fraction(top_cluster_count, total),
        "four_mut_count": four,
        "four_mut_fraction": fraction(four, total),
        "local_expansion_count": local,
        "local_expansion_fraction": fraction(local, total),
        "seed_only_count": seed_only,
        "seed_only_fraction": fraction(seed_only, total),
        "neutral_boundary_or_high_risk_count": neutral,
        "neutral_boundary_or_high_risk_fraction": fraction(neutral, total),
        "AY111H_containing_count": ay111h if target == "Ab_sdAb" else "",
        "AY111H_containing_fraction": fraction(ay111h, total) if target == "Ab_sdAb" else "",
        "AG102H_AV105H_main_count": agav if target == "Ab_sdAb" else "",
        "sdAb_complex_diagnostic_required_count": sdab_diag if target == "Ab_sdAb" else "",
        "buildability_or_hard_status_bad_count": build_bad,
    }


def counter_rows(rows, scenario, fields, top_n=None):
    c = Counter(tuple(norm(row.get(field), "missing") for field in fields) for row in rows)
    items = c.most_common(top_n) if top_n else sorted(c.items())
    out = []
    for key, count in items:
        rec = {"scenario": scenario}
        for field, value in zip(fields, key):
            rec[field] = value
        rec["count"] = count
        out.append(rec)
    return out


def audit_scenario(scenario, rows):
    keys = [r.get("canonical_unique_key") or (r.get("target", "") + "|" + r.get("canonical_sequence_hash_full", "")) for r in rows]
    dupes = len(keys) - len(set(keys))
    target_counts = Counter(r.get("target") for r in rows)
    total = len(rows)
    nominal_total = NOMINAL[scenario]["total"]
    shortfall = max(0, nominal_total - total)
    agav = sum(contains_ag102h_av105h(r) for r in rows if r.get("target") == "Ab_sdAb")
    sdab_rows = [r for r in rows if r.get("target") == "Ab_sdAb"]
    sdab_diag = sum(as_bool(r.get("complex_diagnostic_required")) for r in sdab_rows)
    build_bad = sum(hard_status_bad(r) for r in rows)
    glycan_low_claim = 0
    for r in rows:
        for field in ["glycan_status", "glycan_risk_status", "glycan_claim"]:
            if str(r.get(field, "")).strip().lower() == "low":
                glycan_low_claim += 1

    verdict = "REVIEWABLE_PATCH"
    reasons = []
    if dupes:
        verdict = "FAIL"
        reasons.append("canonical_duplicate_present")
    if agav:
        verdict = "FAIL"
        reasons.append("sdAb_AG102H_or_AV105H_main_present")
    if sdab_diag != len(sdab_rows):
        verdict = "FAIL"
        reasons.append("sdAb_complex_diagnostic_required_not_retained")
    if build_bad:
        verdict = "FAIL"
        reasons.append("hard_or_buildability_status_bad")
    if glycan_low_claim:
        verdict = "FAIL"
        reasons.append("glycan_low_risk_claim_present")
    if shortfall and verdict != "FAIL":
        reasons.append("below_nominal_size")
    if any(target_counts.get(t, 0) < q for t, q in NOMINAL[scenario]["targets"].items()) and verdict != "FAIL":
        reasons.append("target_quota_shortfall")
    if not reasons:
        verdict = "PASS_REVIEWABLE"
        reasons.append("all_audit_checks_pass")
    return {
        "scenario": scenario,
        "selected_count": total,
        "nominal_total": nominal_total,
        "shortfall": shortfall,
        "Ab_1E62_count": target_counts.get("Ab_1E62", 0),
        "Ab_sdAb_count": target_counts.get("Ab_sdAb", 0),
        "canonical_duplicate_count": dupes,
        "sdAb_AG102H_AV105H_main_count": agav,
        "sdAb_complex_diagnostic_required_count": sdab_diag,
        "sdAb_selected_count": len(sdab_rows),
        "buildability_or_hard_status_bad_count": build_bad,
        "glycan_low_risk_claim_count": glycan_low_claim,
        "verdict": verdict,
        "reasons": ";".join(reasons),
    }


def write_single_audit(scenario, rows, group_summary, scenario_audit, control_summary):
    fields = [
        "scenario",
        "target",
        "selected_count",
        "nominal_quota",
        "shortfall",
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
        "AY111H_containing_count",
        "AY111H_containing_fraction",
        "AG102H_AV105H_main_count",
        "sdAb_complex_diagnostic_required_count",
        "buildability_or_hard_status_bad_count",
    ]
    source = counter_rows(rows, scenario, ["target", "source_pool", "template_match_tier"])
    quality = counter_rows(rows, scenario, ["target", "candidate_quality_tier", "template_match_tier"])
    mut = counter_rows(rows, scenario, ["target", "mutation_count"])
    seed_top = counter_rows(rows, scenario, ["target", "his_seed_set"], top_n=20)
    cluster_top = counter_rows(rows, scenario, ["target", "near_duplicate_cluster_id"], top_n=20)

    write_csv(OUT_DIR / ("final_candidate_draft_%s_by_source.csv" % scenario), source)
    write_csv(OUT_DIR / ("final_candidate_draft_%s_by_quality.csv" % scenario), quality)
    write_csv(OUT_DIR / ("final_candidate_draft_%s_by_mutation_count.csv" % scenario), mut)
    write_csv(OUT_DIR / ("final_candidate_draft_%s_top_seed.csv" % scenario), seed_top)
    write_csv(OUT_DIR / ("final_candidate_draft_%s_top_cluster.csv" % scenario), cluster_top)

    text = """# Final Candidate Draft %s Audit

## Verdict

`%s`

Reasons: `%s`

## Scenario Summary

%s

## Target Audit

%s

## Source / Match Tier Summary

%s

## Quality / Match Tier Summary

%s

## Mutation Count Summary

%s

## Top Seed Summary

%s

## Top Cluster Summary

%s

## Control / Anchor Panel

%s

## Interpretation

- This is a reviewable draft, not a synthesis-ready final library.
- sdAb remains secondary and complex-diagnostic-required.
- Glycan remains unchecked; no low-risk glycan claim is made.
- New heavy compute and MD remain locked.
""" % (
        scenario,
        scenario_audit["verdict"],
        scenario_audit["reasons"],
        md_table([scenario_audit], ["scenario", "selected_count", "nominal_total", "shortfall", "Ab_1E62_count", "Ab_sdAb_count", "canonical_duplicate_count", "verdict", "reasons"]),
        md_table(group_summary, fields),
        md_table(source[:25], ["scenario", "target", "source_pool", "template_match_tier", "count"]),
        md_table(quality, ["scenario", "target", "candidate_quality_tier", "template_match_tier", "count"]),
        md_table(mut, ["scenario", "target", "mutation_count", "count"]),
        md_table(seed_top, ["scenario", "target", "his_seed_set", "count"]),
        md_table(cluster_top, ["scenario", "target", "near_duplicate_cluster_id", "count"]),
        md_table(control_summary, ["target", "control_anchor_count", "top_control_seed", "top_control_seed_count", "risk_control_AG102H_AV105H_count"]),
    )
    (OUT_DIR / ("final_candidate_draft_%s_audit.md" % scenario)).write_text(text, encoding="utf-8")


def summarize_controls():
    if not CONTROL_PANEL.exists():
        return []
    rows = read_csv(CONTROL_PANEL)
    out = []
    for target in sorted(set(r.get("target") for r in rows)):
        subset = [r for r in rows if r.get("target") == target]
        seeds = Counter(norm(r.get("his_seed_set"), "missing_seed") for r in subset)
        top_seed, top_count = top_counter(seeds)
        risk = sum(contains_ag102h_av105h(r) for r in subset)
        out.append({
            "target": target,
            "control_anchor_count": len(subset),
            "top_control_seed": top_seed,
            "top_control_seed_count": top_count,
            "risk_control_AG102H_AV105H_count": risk,
        })
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_audits = []
    all_group = []
    control_summary = summarize_controls()
    source_inputs = []

    for scenario, path in DRAFTS.items():
        rows = read_csv(path)
        scenario_audit = audit_scenario(scenario, rows)
        group_summary = [summarize_group(rows, scenario, target) for target in ["Ab_1E62", "Ab_sdAb"]]
        all_audits.append(scenario_audit)
        all_group.extend(group_summary)
        write_single_audit(scenario, rows, group_summary, scenario_audit, control_summary)
        source_inputs.append({"scenario": scenario, "path": str(path.relative_to(ROOT)), "rows": len(rows), "sha256": sha256_file(path)})

    if CONTROL_PANEL.exists():
        source_inputs.append({"scenario": "controls", "path": str(CONTROL_PANEL.relative_to(ROOT)), "rows": len(read_csv(CONTROL_PANEL)), "sha256": sha256_file(CONTROL_PANEL)})
    if REFILL_REPORT.exists():
        source_inputs.append({"scenario": "refill_v3_report", "path": str(REFILL_REPORT.relative_to(ROOT)), "rows": "", "sha256": sha256_file(REFILL_REPORT)})
    if PATCH_RESOLUTION_REPORT.exists():
        source_inputs.append({"scenario": "patch_resolution_report", "path": str(PATCH_RESOLUTION_REPORT.relative_to(ROOT)), "rows": "", "sha256": sha256_file(PATCH_RESOLUTION_REPORT)})

    write_csv(OUT_DIR / "final_candidate_draft_audit_summary.csv", all_audits)
    write_csv(OUT_DIR / "final_candidate_draft_target_summary.csv", all_group)
    write_csv(OUT_DIR / "final_candidate_draft_control_anchor_summary.csv", control_summary)
    write_csv(OUT_DIR / "final_candidate_draft_input_manifest.csv", source_inputs)

    overall = "FINAL_DRAFT_AUDIT_COMPLETE__15K_PRIMARY_REVIEW_RECOMMENDED__FINAL_SELECTION_LOCKED"
    comparison = """# Final Candidate Draft Comparison

## Executive Summary

Verdict: `%s`.

The 10K and 15K refill-v3 drafts were audited as reviewable PATCH drafts. This stage did not generate new candidates, did not run AF3, SimpleFold, PyRosetta, FoldX, MD, or glycan modeling, and did not unlock final synthesis selection.

Recommendation: use the 15K draft as the primary review object because it provides 14,220 candidates within the intended 10K-15K practical capacity range. Do not repeat the same refill scan. If exact 15K is required, choose a new policy: limited 1E62 rule relaxation or new targeted 1E62 generation.

## Inputs

%s

## Scenario Comparison

%s

## Target-Level Audit

%s

## Control / Anchor Panel

%s

## Decision Options

1. Accept the 15K reviewable PATCH draft as the main final candidate draft basis.
2. If exact 10K / 15K size is mandatory, define a new 1E62 policy: limited rule relaxation or new targeted generation.
3. Do not continue same-rule source rescan; it has already returned zero eligible refill candidates.

## Still Locked

- Final synthesis-ready library selection.
- Wet-lab order list generation.
- New AF3 / SimpleFold / PyRosetta / FoldX compute.
- MD / constant-pH MD.
- sdAb structure-confirmed upgrade.
- Glycan low-risk claim.
""" % (
        overall,
        md_table(source_inputs, ["scenario", "path", "rows", "sha256"]),
        md_table(all_audits, ["scenario", "selected_count", "nominal_total", "shortfall", "Ab_1E62_count", "Ab_sdAb_count", "canonical_duplicate_count", "sdAb_AG102H_AV105H_main_count", "sdAb_complex_diagnostic_required_count", "sdAb_selected_count", "glycan_low_risk_claim_count", "verdict", "reasons"]),
        md_table(all_group, ["scenario", "target", "selected_count", "nominal_quota", "shortfall", "top_seed", "top_seed_fraction", "top_cluster_fraction", "four_mut_fraction", "local_expansion_fraction", "seed_only_fraction", "neutral_boundary_or_high_risk_fraction", "AY111H_containing_fraction", "AG102H_AV105H_main_count", "sdAb_complex_diagnostic_required_count"]),
        md_table(control_summary, ["target", "control_anchor_count", "top_control_seed", "top_control_seed_count", "risk_control_AG102H_AV105H_count"]),
    )
    (OUT_DIR / "final_candidate_draft_comparison.md").write_text(comparison, encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(comparison, encoding="utf-8")
    print(overall)
    print(OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
