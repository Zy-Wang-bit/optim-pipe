#!/usr/bin/env python3
"""Build per-target 15K candidate-pool drafts from expanded v2 pools.

New corrected objective:

- Ab_1E62 gets its own 15K pre-synthesis candidate pool.
- Ab_sdAb gets its own 15K pre-synthesis candidate pool.
- Controls / anchors are counted inside each target's 15K pool.

No new variants are generated and no heavy computation is run.
"""

import csv
import hashlib
import os
import re
from collections import Counter, defaultdict


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
EXPANDED_V2 = os.path.join(
    ROOT,
    "results/initial_design_generation/expanded_tier2_candidate_pool_v2/"
    "expanded_tier2_candidate_pool_v2.csv",
)
CONTROL_PANEL = os.path.join(
    ROOT,
    "results/initial_design_generation/final_candidate_pool_planning/"
    "final_candidate_pool_control_anchor_panel.csv",
)
OUT_DIR = os.path.join(
    ROOT,
    "results/initial_design_generation/per_target_15k_candidate_pools",
)
TASK_REPORT = os.path.join(
    ROOT,
    ".tasks/active/initial-design-generation/current_stage_report.md",
)

LOOKUP_SOURCES = [
    EXPANDED_V2,
    os.path.join(
        ROOT,
        "results/initial_design_generation/final_candidate_pool_refill_v3/"
        "final_candidate_pool_15k_refill_v3_draft.csv",
    ),
    os.path.join(
        ROOT,
        "results/initial_design_generation/tier2_staged/tier2_candidate_snapshot.csv",
    ),
    os.path.join(
        ROOT,
        "results/initial_design_generation/stage1_5_stage2a/stage2a_candidate_list.csv",
    ),
    os.path.join(
        ROOT,
        "results/initial_design_generation/tier1_filtering/"
        "tier1_ranked_candidates_pre_tier2.csv",
    ),
    os.path.join(
        ROOT,
        "results/initial_design_generation/production_initial_pool/"
        "production_initial_pool_candidates_all.csv",
    ),
]

TARGETS = ["Ab_1E62", "Ab_sdAb"]
TARGET_SIZE = 15000
ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")

TARGET_CONFIG = {
    "Ab_1E62": {
        "short_name": "1E62",
        "top_seed_cap": int(15000 * 0.28),
        "top_cluster_cap": int(15000 * 0.06),
        "seed_only_cap": int(15000 * 0.20),
        "neutral_boundary_cap": int(15000 * 0.25),
        # 1E62 needs a small policy relaxation to reach 15K from the
        # existing 20K mother pool. The reviewed policy allows this before
        # any new targeted generation. Keep it below the 65% hard-like range
        # used for 4-mut tolerance.
        "local_expansion_cap": 9500,
        "four_mut_cap": int(15000 * 0.65),
        "ay111h_cap": None,
        "exclude_ag102h_av105h_main": False,
    },
    "Ab_sdAb": {
        "short_name": "sdAb",
        "top_seed_cap": int(15000 * 0.30),
        "top_cluster_cap": int(15000 * 0.06),
        "seed_only_cap": int(15000 * 0.20),
        "neutral_boundary_cap": int(15000 * 0.10),
        "local_expansion_cap": int(15000 * 0.40),
        "four_mut_cap": int(15000 * 0.15),
        "ay111h_cap": int(15000 * 0.35),
        "exclude_ag102h_av105h_main": True,
    },
}


def read_csv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], [dict(row) for row in reader]


def write_csv(path, fields, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def truthy(value):
    return str(value).strip().lower() in ("1", "true", "yes", "y")


def as_int(value, default=0):
    try:
        if value is None or str(value).strip() == "":
            return default
        return int(float(value))
    except Exception:
        return default


def as_float(value, default=0.0):
    try:
        if value is None or str(value).strip() == "":
            return default
        return float(value)
    except Exception:
        return default


def parse_mutations(mutation_list):
    muts = []
    for token in str(mutation_list or "").split(";"):
        token = token.strip()
        if not token:
            continue
        match = re.match(r"^([A-Z])([A-Z])(\d+)([A-Z])$", token)
        if not match:
            muts.append((token, None, None, None, None, "parse_error"))
            continue
        chain, old, pos, new = match.groups()
        muts.append((token, chain, old, int(pos), new, "ok"))
    return muts


def motif_positions(seq):
    out = set()
    for i in range(max(0, len(seq) - 2)):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in ("S", "T"):
            out.add(i + 1)
    return out


def reverse_to_parent(row):
    seq = row.get("sequence", "")
    if not seq:
        return None
    parent = list(seq)
    for token, chain, old, pos, new, status in parse_mutations(row.get("mutation_list", "")):
        if status != "ok" or pos < 1 or pos > len(parent) or parent[pos - 1] != new:
            return None
        parent[pos - 1] = old
    return "".join(parent)


def derive_parents(rows):
    parents = {}
    for target in TARGETS:
        for row in rows:
            if row.get("target") != target:
                continue
            p = reverse_to_parent(row)
            if p:
                parents[target] = p
                break
    return parents


def apply_mutations(parent_seq, mutation_list):
    if not parent_seq:
        return "", "missing_parent_sequence"
    seq = list(parent_seq)
    problems = []
    for token, chain, old, pos, new, status in parse_mutations(mutation_list):
        if status != "ok":
            problems.append("%s:%s" % (token, status))
            continue
        if pos < 1 or pos > len(seq):
            problems.append("%s:position_out_of_range" % token)
            continue
        if seq[pos - 1] != old:
            problems.append("%s:parent_mismatch_%s_%s" % (token, old, seq[pos - 1]))
            continue
        seq[pos - 1] = new
    status = "reconstructed_from_parent_mutation_list"
    if problems:
        status = "mutation_reconstruction_warning:" + "|".join(problems)
    return "".join(seq), status


def build_sequence_lookup(paths):
    lookup = {}
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            if "variant_id" not in fields or "sequence" not in fields:
                continue
            for row in reader:
                vid = row.get("variant_id", "")
                seq = row.get("sequence", "")
                if vid and seq and vid not in lookup:
                    lookup[vid] = (seq, os.path.relpath(path, ROOT))
    return lookup


def is_local(row):
    return row.get("template_match_tier", "").startswith("local_expansion") or "local_expansion" in row.get("source_pool", "")


def is_seed_only(row):
    return row.get("template_match_tier", "") == "C_seed_only"


def is_four_mut(row):
    return as_int(row.get("mutation_count", ""), 0) >= 4


def is_neutral_boundary(row):
    return row.get("candidate_quality_tier", "") == "neutral_boundary_or_high_risk"


def contains_ay111h(row):
    return "AY111H" in row.get("mutation_list", "")


def contains_ag102h_av105h(row):
    muts = row.get("mutation_list", "")
    return "AG102H" in muts or "AV105H" in muts


def bad_status(value):
    text = str(value or "").strip().lower()
    if text == "":
        return False
    return text not in ("pass", "passed", "ok", "true", "eligible", "accepted")


def row_hard_failures(row, parents):
    failures = []
    seq = row.get("sequence", "")
    target = row.get("target", "")
    parent = parents.get(target, "")
    if not seq:
        failures.append("missing_sequence")
        return failures
    illegal = sorted(set(seq) - ALLOWED_AA)
    if illegal:
        failures.append("illegal_aa:" + "".join(illegal))
    if parent and len(seq) != len(parent):
        failures.append("sequence_length_mismatch")
    if parent and seq.count("C") > parent.count("C"):
        failures.append("new_cys_by_sequence_count")
    for token, chain, old, pos, new, status in parse_mutations(row.get("mutation_list", "")):
        if status != "ok":
            failures.append("mutation_parse_error:" + token)
        elif new == "C":
            failures.append("new_cys_by_mutation:" + token)
    if parent:
        new_motifs = motif_positions(seq) - motif_positions(parent)
        if new_motifs:
            failures.append("new_canonical_NXS_T_motif:" + ",".join(map(str, sorted(new_motifs))))
    if target == "Ab_sdAb":
        muts = set(row.get("mutation_list", "").split(";"))
        if "AV105H" in muts and "AD110H" in muts:
            failures.append("forbidden_pair_AV105H_AD110H")
    for field in ("hard_filter_status", "buildability_light_status", "forbidden_pair_status"):
        if bad_status(row.get(field, "")):
            failures.append("bad_%s:%s" % (field, row.get(field, "")))
    class_text = " ".join(
        str(row.get(field, ""))
        for field in ("candidate_quality_tier", "tier1_review_class", "tier2_class", "t2_class_current")
    ).lower()
    if any(term in class_text for term in ("severe", "unsupported", "reject", "structure_risk", "mechanism_weak")):
        failures.append("severe_or_unsupported_class")
    return failures


def normalize_control(row, lookup, parents):
    out = dict(row)
    if row.get("variant_id") in lookup:
        seq, source = lookup[row.get("variant_id")]
        status = "sequence_lookup_by_variant_id"
    else:
        seq, status = apply_mutations(parents.get(row.get("target", ""), ""), row.get("mutation_list", ""))
        source = "parent_sequence_plus_mutation_list"
    full_hash = sha256_text(seq) if seq else ""
    out["sequence"] = seq
    out["sequence_hash"] = full_hash[:10]
    out["canonical_sequence_hash_full"] = full_hash
    out["canonical_unique_key"] = "%s|%s" % (row.get("target", ""), full_hash)
    out["candidate_type"] = "control_anchor"
    out["evidence_class"] = "control_anchor"
    out["quality_tier"] = "control_anchor"
    out["source_pool"] = "control_anchor_panel"
    out["local_expansion_flag"] = "False"
    out["seed_only_flag"] = "False"
    out["glycan_status"] = "glycan_unchecked_no_explicit_review"
    out["buildability_status"] = "pending_audit"
    out["selection_reason"] = "control_anchor_counted_inside_target_15k"
    out["sequence_source"] = source
    out["sequence_reconstruction_status"] = status
    out["complex_diagnostic_required"] = "True" if row.get("target") == "Ab_sdAb" else ""
    out["not_final_order_list"] = "True"
    out["manual_review_required"] = "True"
    return out


def selection_score(row, target):
    score = as_float(row.get("expanded_pool_score", ""), 0.0)
    quality = row.get("candidate_quality_tier", "")
    tier = row.get("template_match_tier", "")
    status = row.get("assigned_template_status", "")
    source = row.get("source_pool", "")

    score += {
        "bank_reviewed": 1200,
        "high_support": 950,
        "supported": 820,
        "mpnn_favorable_pH_weak": 650,
        "neutral_boundary_or_high_risk": 260,
    }.get(quality, 0)
    score += {
        "template_seed": 1100,
        "A_strong": 980,
        "local_expansion_A": 900,
        "B_medium": 800,
        "local_expansion_B": 680,
        "C_seed_only": 260,
    }.get(tier, 0)
    score += {
        "active_primary_template": 180,
        "active_primary_backfill_template": 90,
        "limited_boundary_template": -80,
        "secondary_template_complex_weak": 120,
        "secondary_template_complex_unchecked": 90,
        "low_priority_secondary_template": -120,
        "boundary_representative_only": -180,
    }.get(status, 0)
    if source == "tier2_candidate_bank":
        score += 250
    elif source in ("stage2a_candidate_list", "sdab_recovery_passlike_supplement"):
        score += 100
    elif source == "tier2_candidate_snapshot":
        score += 55
    elif source == "tier1_ranked_production_pool":
        score += 30
    if is_seed_only(row):
        score -= 160
    if is_neutral_boundary(row):
        score -= 180
    if is_local(row):
        score -= 50
    if is_four_mut(row):
        score -= 30 if target == "Ab_1E62" else 80
    score -= max(0, as_int(row.get("mutation_count", ""), 0) - 4) * 40
    return score


def can_add(row, counts, cfg):
    seed = row.get("his_seed_set", "") or "NA"
    cluster = row.get("near_duplicate_cluster_id", "") or "NA"
    if counts["seed"][seed] >= cfg["top_seed_cap"]:
        return False, "seed_cap_full"
    if counts["cluster"][cluster] >= cfg["top_cluster_cap"]:
        return False, "cluster_cap_full"
    if is_seed_only(row) and counts["seed_only"] >= cfg["seed_only_cap"]:
        return False, "seed_only_cap_full"
    if is_neutral_boundary(row) and counts["neutral"] >= cfg["neutral_boundary_cap"]:
        return False, "neutral_boundary_cap_full"
    if is_local(row) and counts["local"] >= cfg["local_expansion_cap"]:
        return False, "local_expansion_cap_full"
    if is_four_mut(row) and counts["four_mut"] >= cfg["four_mut_cap"]:
        return False, "four_mut_cap_full"
    if cfg.get("ay111h_cap") is not None and contains_ay111h(row) and counts["ay111h"] >= cfg["ay111h_cap"]:
        return False, "ay111h_cap_full"
    return True, "selected"


def add_count(row, counts):
    counts["seed"][row.get("his_seed_set", "") or "NA"] += 1
    counts["cluster"][row.get("near_duplicate_cluster_id", "") or "NA"] += 1
    if is_seed_only(row):
        counts["seed_only"] += 1
    if is_neutral_boundary(row):
        counts["neutral"] += 1
    if is_local(row):
        counts["local"] += 1
    if is_four_mut(row):
        counts["four_mut"] += 1
    if contains_ay111h(row):
        counts["ay111h"] += 1


def make_counts():
    return {
        "seed": Counter(),
        "cluster": Counter(),
        "seed_only": 0,
        "neutral": 0,
        "local": 0,
        "four_mut": 0,
        "ay111h": 0,
    }


def decorate_main(row, target):
    out = dict(row)
    out["candidate_type"] = "main_candidate"
    out["evidence_class"] = row.get("assigned_template_status", "") or row.get("assigned_expansion_role", "")
    out["quality_tier"] = row.get("candidate_quality_tier", "")
    out["local_expansion_flag"] = "True" if is_local(row) else "False"
    out["seed_only_flag"] = "True" if is_seed_only(row) else "False"
    out["glycan_status"] = "glycan_unchecked_no_explicit_review"
    out["buildability_status"] = "pending_audit"
    out["selection_reason"] = "selected_from_expanded_v2_per_target_15k"
    out["not_final_order_list"] = "True"
    out["manual_review_required"] = "True"
    if target == "Ab_sdAb":
        out["complex_diagnostic_required"] = "True"
    return out


def select_target(target, source_rows, control_rows, lookup, parents):
    cfg = TARGET_CONFIG[target]
    target_controls = [normalize_control(r, lookup, parents) for r in control_rows if r.get("target") == target]
    control_hashes = set(r.get("canonical_sequence_hash_full", "") for r in target_controls if r.get("canonical_sequence_hash_full"))
    control_keys = set(r.get("canonical_unique_key", "") for r in target_controls if r.get("canonical_unique_key"))
    control_ids = set(r.get("variant_id", "") for r in target_controls if r.get("variant_id"))

    eligible = []
    rejection = Counter()
    for row in source_rows:
        if row.get("target") != target:
            continue
        if row.get("variant_id") in control_ids or row.get("canonical_unique_key") in control_keys or row.get("canonical_sequence_hash_full") in control_hashes:
            rejection["duplicate_with_control"] += 1
            continue
        if cfg.get("exclude_ag102h_av105h_main") and contains_ag102h_av105h(row):
            rejection["ag102h_av105h_main_excluded"] += 1
            continue
        failures = row_hard_failures(row, parents)
        if failures:
            rejection["hard_failure"] += 1
            continue
        eligible.append(row)

    eligible.sort(key=lambda r: (-selection_score(r, target), r.get("variant_id", "")))
    selected = []
    counts = make_counts()
    for row in target_controls:
        add_count(row, counts)
    for row in eligible:
        if len(selected) + len(target_controls) >= TARGET_SIZE:
            break
        ok, reason = can_add(row, counts, cfg)
        if not ok:
            rejection[reason] += 1
            continue
        selected.append(decorate_main(row, target))
        add_count(row, counts)
    final_rows = selected + target_controls
    final_rows.sort(key=lambda r: (0 if r.get("candidate_type") == "main_candidate" else 1, -selection_score(r, target), r.get("variant_id", "")))
    for idx, row in enumerate(final_rows, 1):
        row["per_target_15k_row_id"] = "%s_15K_%05d" % (cfg["short_name"], idx)
        row["draft_status"] = "draft_not_final_order_list"
    return final_rows, rejection, eligible


def audit_target(target, rows, parents):
    total = len(rows)
    main = [r for r in rows if r.get("candidate_type") == "main_candidate"]
    controls = [r for r in rows if r.get("candidate_type") == "control_anchor"]
    failures = []
    for row in rows:
        fail = row_hard_failures(row, parents)
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
    summary = {
        "target": target,
        "selected_count": total,
        "main_count": len(main),
        "control_anchor_count": len(controls),
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
        "four_mut_count": sum(1 for r in rows if is_four_mut(r)),
        "four_mut_fraction": "%.6f" % (float(sum(1 for r in rows if is_four_mut(r))) / total if total else 0),
        "local_expansion_count": sum(1 for r in rows if is_local(r)),
        "local_expansion_fraction": "%.6f" % (float(sum(1 for r in rows if is_local(r))) / total if total else 0),
        "seed_only_count": sum(1 for r in rows if is_seed_only(r)),
        "seed_only_fraction": "%.6f" % (float(sum(1 for r in rows if is_seed_only(r))) / total if total else 0),
        "neutral_boundary_or_high_risk_count": sum(1 for r in rows if is_neutral_boundary(r)),
        "neutral_boundary_or_high_risk_fraction": "%.6f" % (float(sum(1 for r in rows if is_neutral_boundary(r))) / total if total else 0),
        "glycan_low_risk_claim_count": sum(1 for r in rows if "low" in str(r.get("glycan_status", "")).lower()),
        "sdAb_complex_diagnostic_required_count": "",
        "AG102H_AV105H_main_count": "",
        "AY111H_containing_count": "",
        "AY111H_containing_fraction": "",
    }
    seed = Counter(r.get("his_seed_set", "") or "NA" for r in rows)
    cluster = Counter(r.get("near_duplicate_cluster_id", "") or "NA" for r in rows)
    if seed:
        summary["top_seed"], summary["top_seed_count"] = seed.most_common(1)[0]
        summary["top_seed_fraction"] = "%.6f" % (float(summary["top_seed_count"]) / total)
    if cluster:
        summary["top_cluster"], summary["top_cluster_count"] = cluster.most_common(1)[0]
        summary["top_cluster_fraction"] = "%.6f" % (float(summary["top_cluster_count"]) / total)
    if target == "Ab_sdAb":
        diag = sum(1 for r in rows if truthy(r.get("complex_diagnostic_required", "")))
        summary["sdAb_complex_diagnostic_required_count"] = diag
        summary["AG102H_AV105H_main_count"] = sum(1 for r in main if contains_ag102h_av105h(r))
        ay = sum(1 for r in rows if contains_ay111h(r))
        summary["AY111H_containing_count"] = ay
        summary["AY111H_containing_fraction"] = "%.6f" % (float(ay) / total if total else 0)
    verdict = "PASS"
    reasons = []
    if total != TARGET_SIZE:
        verdict = "FAIL"; reasons.append("selected_count_not_15000")
    if summary["canonical_duplicate_count"] or summary["sequence_duplicate_count"] or summary["variant_id_duplicate_count"]:
        verdict = "FAIL"; reasons.append("duplicate_detected")
    if failures:
        verdict = "FAIL"; reasons.append("hard_failure_detected")
    if target == "Ab_sdAb":
        if summary["sdAb_complex_diagnostic_required_count"] != total:
            verdict = "FAIL"; reasons.append("sdAb_complex_diagnostic_not_100pct")
        if summary["AG102H_AV105H_main_count"]:
            verdict = "FAIL"; reasons.append("AG102H_AV105H_main_present")
    if summary["glycan_low_risk_claim_count"]:
        verdict = "FAIL"; reasons.append("glycan_low_risk_claim_present")
    summary["verdict"] = verdict
    summary["reasons"] = ";".join(reasons) if reasons else "all_per_target_15k_gates_pass"
    return summary, failures


def md_table(rows, columns):
    out = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        out.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(out)


def by_count(rows, fields):
    c = Counter(tuple(r.get(f, "") for f in fields) for r in rows)
    out = []
    for key, count in sorted(c.items(), key=lambda kv: (-kv[1], kv[0])):
        row = {field: value for field, value in zip(fields, key)}
        row["count"] = count
        out.append(row)
    return out


def top_rows(rows, field, limit=15):
    c = Counter(r.get(field, "") or "NA" for r in rows)
    total = len(rows)
    return [{field: k, "count": v, "fraction": "%.6f" % (float(v) / total if total else 0)} for k, v in c.most_common(limit)]


def write_target_audit(path, target, rows, summary, failures, rejection):
    with open(path, "w") as h:
        h.write("# %s per-target 15K candidate draft audit\n\n" % TARGET_CONFIG[target]["short_name"])
        h.write("Verdict: `%s`.\n\n" % summary["verdict"])
        h.write("This is a pre-synthesis candidate-pool draft, not a final order list.\n")
        h.write("No new variants, AF3, SimpleFold, PyRosetta, FoldX, MD, or glycan modeling were run.\n\n")
        h.write("## Summary\n\n")
        cols = ["target", "selected_count", "main_count", "control_anchor_count", "canonical_duplicate_count", "sequence_duplicate_count", "variant_id_duplicate_count", "hard_failure_count", "top_seed", "top_seed_count", "top_seed_fraction", "top_cluster", "top_cluster_count", "top_cluster_fraction", "four_mut_count", "four_mut_fraction", "local_expansion_count", "local_expansion_fraction", "seed_only_count", "seed_only_fraction", "neutral_boundary_or_high_risk_count", "neutral_boundary_or_high_risk_fraction", "glycan_low_risk_claim_count", "sdAb_complex_diagnostic_required_count", "AG102H_AV105H_main_count", "AY111H_containing_count", "AY111H_containing_fraction", "reasons"]
        h.write(md_table([summary], cols))
        h.write("\n\n## Evidence class distribution\n\n")
        h.write(md_table(by_count(rows, ["candidate_type", "evidence_class", "quality_tier", "template_match_tier"])[:40], ["candidate_type", "evidence_class", "quality_tier", "template_match_tier", "count"]))
        h.write("\n\n## Source distribution\n\n")
        h.write(md_table(by_count(rows, ["candidate_type", "source_pool", "quality_tier", "template_match_tier"])[:40], ["candidate_type", "source_pool", "quality_tier", "template_match_tier", "count"]))
        h.write("\n\n## Mutation count distribution\n\n")
        h.write(md_table(by_count(rows, ["mutation_count"])[:20], ["mutation_count", "count"]))
        h.write("\n\n## Seed distribution\n\n")
        h.write(md_table(top_rows(rows, "his_seed_set", 20), ["his_seed_set", "count", "fraction"]))
        h.write("\n\n## Cluster distribution\n\n")
        h.write(md_table(top_rows(rows, "near_duplicate_cluster_id", 20), ["near_duplicate_cluster_id", "count", "fraction"]))
        h.write("\n\n## Rejection / cap pressure summary\n\n")
        rej = [{"reject_reason": k, "count": v} for k, v in rejection.most_common()]
        h.write(md_table(rej, ["reject_reason", "count"]) if rej else "No rejected rows under target construction.\n")
        h.write("\n\n## Conservative labels\n\n")
        h.write("- Glycan remains unchecked; no low-risk glycan claim is made.\n")
        if target == "Ab_sdAb":
            h.write("- sdAb remains secondary and complex-diagnostic-required for all rows.\n")
            h.write("- AG102H / AV105H are excluded from main candidates; if present, they are controls only.\n")
        h.write("- Manual review is required before final wet-lab order generation.\n")


def write_comparison(path, summaries):
    with open(path, "w") as h:
        h.write("# Per-target 15K candidate draft comparison\n\n")
        overall = "PASS_PER_TARGET_15K_DRAFTS_BUILT__FINAL_ORDER_LOCKED"
        if any(s["verdict"] != "PASS" for s in summaries):
            overall = "PATCH_OR_FAIL_REVIEW_REQUIRED__FINAL_ORDER_LOCKED"
        h.write("Verdict: `%s`.\n\n" % overall)
        h.write("Each antibody now has its own 15K pre-synthesis candidate draft. These are candidate pools, not final order lists.\n\n")
        cols = ["target", "selected_count", "main_count", "control_anchor_count", "top_seed", "top_seed_fraction", "top_cluster_fraction", "four_mut_fraction", "local_expansion_fraction", "seed_only_fraction", "neutral_boundary_or_high_risk_fraction", "sdAb_complex_diagnostic_required_count", "AG102H_AV105H_main_count", "glycan_low_risk_claim_count", "verdict", "reasons"]
        h.write(md_table(summaries, cols))
        h.write("\n\n## Output files\n\n")
        h.write("- `final_candidate_pool_1E62_15k_draft.csv`\n")
        h.write("- `final_candidate_pool_sdAb_15k_draft.csv`\n")
        h.write("- `final_candidate_pool_1E62_15k_audit.md`\n")
        h.write("- `final_candidate_pool_sdAb_15k_audit.md`\n")
        h.write("- `control_anchor_15k_plan.md`\n")
        h.write("\n## Locked items\n\n")
        h.write("- Direct synthesis / wet-lab order generation remains locked.\n")
        h.write("- sdAb remains secondary / complex-diagnostic-required, not structure-confirmed.\n")
        h.write("- Glycan remains unchecked, not low risk.\n")


def write_control_plan(path, target_rows):
    controls = []
    for rows in target_rows.values():
        controls.extend([r for r in rows if r.get("candidate_type") == "control_anchor"])
    with open(path, "w") as h:
        h.write("# Control anchor 15K plan\n\n")
        h.write("Controls / anchors are counted inside each target's 15K candidate pool.\n\n")
        h.write(md_table(by_count(controls, ["target", "his_seed_set", "mutation_count", "sequence_reconstruction_status"]), ["target", "his_seed_set", "mutation_count", "sequence_reconstruction_status", "count"]))
        h.write("\n\nCurrent explicit control panel contains 20 rows per target. No additional synthetic controls were invented in this step.\n")


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)
    source_fields, source_rows = read_csv(EXPANDED_V2)
    control_fields, control_rows = read_csv(CONTROL_PANEL)
    parents = derive_parents(source_rows)
    lookup = build_sequence_lookup(LOOKUP_SOURCES)

    target_rows = {}
    summaries = []
    all_failures = []
    all_rejections = {}
    output_fields = []
    required = ["variant_id", "target", "sequence", "mutation_list", "his_seed_set", "near_duplicate_cluster_id", "candidate_type", "evidence_class", "quality_tier", "template_match_tier", "source_pool", "local_expansion_flag", "seed_only_flag", "complex_diagnostic_required", "glycan_status", "buildability_status", "selection_reason"]
    for field in required + ["per_target_15k_row_id", "draft_status", "not_final_order_list", "manual_review_required"] + source_fields + control_fields:
        if field not in output_fields:
            output_fields.append(field)

    for target in TARGETS:
        rows, rejection, eligible = select_target(target, source_rows, control_rows, lookup, parents)
        target_rows[target] = rows
        summary, failures = audit_target(target, rows, parents)
        summaries.append(summary)
        all_failures.extend(failures)
        all_rejections[target] = rejection
        short = TARGET_CONFIG[target]["short_name"]
        write_csv(os.path.join(OUT_DIR, "final_candidate_pool_%s_15k_draft.csv" % short), output_fields, rows)
        write_target_audit(os.path.join(OUT_DIR, "final_candidate_pool_%s_15k_audit.md" % short), target, rows, summary, failures, rejection)

    write_comparison(os.path.join(OUT_DIR, "final_candidate_pool_15k_by_target_comparison.md"), summaries)
    write_control_plan(os.path.join(OUT_DIR, "control_anchor_15k_plan.md"), target_rows)
    write_csv(
        os.path.join(OUT_DIR, "per_target_15k_audit_summary.csv"),
        list(summaries[0].keys()),
        summaries,
    )
    write_csv(
        os.path.join(OUT_DIR, "per_target_15k_hard_failure_rows.csv"),
        ["target", "variant_id", "candidate_type", "failure"],
        all_failures,
    )
    manifest = [
        {"input": "expanded_v2", "path": os.path.relpath(EXPANDED_V2, ROOT), "rows": len(source_rows), "sha256": sha256_file(EXPANDED_V2)},
        {"input": "control_panel", "path": os.path.relpath(CONTROL_PANEL, ROOT), "rows": len(control_rows), "sha256": sha256_file(CONTROL_PANEL)},
    ]
    write_csv(os.path.join(OUT_DIR, "per_target_15k_input_manifest.csv"), ["input", "path", "rows", "sha256"], manifest)

    with open(os.path.join(OUT_DIR, "final_candidate_pool_15k_by_target_comparison.md")) as src, open(TASK_REPORT, "w") as dst:
        dst.write(src.read())

    verdict = "PASS_PER_TARGET_15K_DRAFTS_BUILT__FINAL_ORDER_LOCKED"
    if any(s["verdict"] != "PASS" for s in summaries):
        verdict = "PATCH_OR_FAIL_REVIEW_REQUIRED__FINAL_ORDER_LOCKED"
    print(verdict)
    print(os.path.relpath(OUT_DIR, ROOT))


if __name__ == "__main__":
    main()
