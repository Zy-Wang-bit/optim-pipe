#!/usr/bin/env python3
"""Build and audit a synthesis-ready draft candidate list.

This script executes the current planning policy:

- use the 15K refill-v3 draft as the primary review object;
- retain the 40 control / anchor rows inside the realistic 14,220 capacity;
- remove low-priority main rows per target to make room for controls;
- keep all outputs explicitly marked as draft / not final order list.

The script intentionally uses only the Python standard library so it can run in
the legacy Python 3.6 environment available on this host.
"""

import csv
import hashlib
import os
import re
from collections import Counter, defaultdict


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

PRIMARY_DRAFT = os.path.join(
    ROOT,
    "results/initial_design_generation/final_candidate_pool_refill_v3/"
    "final_candidate_pool_15k_refill_v3_draft.csv",
)
BACKUP_DRAFT = os.path.join(
    ROOT,
    "results/initial_design_generation/final_candidate_pool_refill_v3/"
    "final_candidate_pool_10k_refill_v3_draft.csv",
)
CONTROL_PANEL = os.path.join(
    ROOT,
    "results/initial_design_generation/final_candidate_pool_planning/"
    "final_candidate_pool_control_anchor_panel.csv",
)
OUT_DIR = os.path.join(
    ROOT,
    "results/initial_design_generation/final_synthesis_ready_selection_planning",
)
TASK_REPORT = os.path.join(
    ROOT,
    ".tasks/active/initial-design-generation/current_stage_report.md",
)

SEQUENCE_LOOKUP_SOURCES = [
    PRIMARY_DRAFT,
    BACKUP_DRAFT,
    os.path.join(
        ROOT,
        "results/initial_design_generation/expanded_tier2_candidate_pool_v2/"
        "expanded_tier2_candidate_pool_v2.csv",
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
ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")
CAPACITY_POLICY = "Option_A_realistic_14220_controls_counted_inside"
PRIMARY_CAPACITY = 14220
CONTROL_ROLE = "control_anchor"
MAIN_ROLE = "main_candidate"


def read_csv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [dict(row) for row in reader]
        return reader.fieldnames or [], rows


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
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
    mutations = []
    if not mutation_list:
        return mutations
    for token in str(mutation_list).split(";"):
        token = token.strip()
        if not token:
            continue
        match = re.match(r"^([A-Z])([A-Z])(\d+)([A-Z])$", token)
        if not match:
            mutations.append((token, None, None, None, None, "parse_error"))
            continue
        chain, old, pos, new = match.groups()
        mutations.append((token, chain, old, int(pos), new, "ok"))
    return mutations


def reverse_to_parent(row):
    seq = row.get("sequence", "")
    if not seq:
        return None
    parent = list(seq)
    for token, chain, old, pos, new, status in parse_mutations(row.get("mutation_list", "")):
        if status != "ok":
            return None
        if pos < 1 or pos > len(parent):
            return None
        if parent[pos - 1] != new:
            return None
        parent[pos - 1] = old
    return "".join(parent)


def derive_parent_sequences(rows):
    parent = {}
    for target in TARGETS:
        for row in rows:
            if row.get("target") != target:
                continue
            pseq = reverse_to_parent(row)
            if pseq:
                parent[target] = pseq
                break
    return parent


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
            problems.append(
                "%s:parent_mismatch_expected_%s_observed_%s"
                % (token, old, seq[pos - 1])
            )
            continue
        seq[pos - 1] = new
    if problems:
        return "".join(seq), "mutation_reconstruction_warning:" + "|".join(problems)
    return "".join(seq), "reconstructed_from_parent_mutation_list"


def build_sequence_lookup(sources):
    lookup = {}
    source_hits = Counter()
    for path in sources:
        if not os.path.exists(path):
            continue
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle)
            if "variant_id" not in (reader.fieldnames or []):
                continue
            has_sequence = "sequence" in (reader.fieldnames or [])
            if not has_sequence:
                continue
            for row in reader:
                vid = row.get("variant_id", "")
                seq = row.get("sequence", "")
                if vid and seq and vid not in lookup:
                    lookup[vid] = (seq, os.path.relpath(path, ROOT))
                    source_hits[os.path.relpath(path, ROOT)] += 1
    return lookup, source_hits


def is_bad_status(value):
    value = str(value or "").strip().lower()
    if value == "":
        return False
    good = set(["pass", "passed", "ok", "true", "eligible", "accepted"])
    return value not in good


def motif_positions(seq):
    positions = set()
    for i in range(0, max(0, len(seq) - 2)):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in ("S", "T"):
            positions.add(i + 1)
    return positions


def classify_hard_failures(row, parents):
    failures = []
    seq = row.get("sequence", "")
    target = row.get("target", "")
    parent_seq = parents.get(target, "")
    if not seq:
        failures.append("missing_sequence")
        return failures
    illegal = sorted(set(seq) - ALLOWED_AA)
    if illegal:
        failures.append("illegal_aa:" + "".join(illegal))
    if parent_seq and len(seq) != len(parent_seq):
        failures.append("sequence_length_mismatch")
    if parent_seq and seq.count("C") > parent_seq.count("C"):
        failures.append("new_cys_by_sequence_count")
    if "C" in [m[4] for m in parse_mutations(row.get("mutation_list", "")) if m[5] == "ok"]:
        failures.append("new_cys_by_mutation")
    if parent_seq:
        new_motifs = motif_positions(seq) - motif_positions(parent_seq)
        if new_motifs:
            failures.append("new_canonical_NXS_T_motif:" + ",".join(map(str, sorted(new_motifs))))
    if target == "Ab_sdAb":
        muts = set(row.get("mutation_list", "").split(";"))
        if "AV105H" in muts and "AD110H" in muts:
            failures.append("forbidden_pair_AV105H_AD110H")
    for field in ("hard_filter_status", "buildability_light_status", "forbidden_pair_status"):
        if is_bad_status(row.get(field, "")):
            failures.append("bad_%s:%s" % (field, row.get(field, "")))
    class_text = " ".join(
        str(row.get(field, ""))
        for field in (
            "candidate_quality_tier",
            "stage2a_final_class",
            "tier2_class",
            "t2_class_current",
            "tier1_review_class",
            "selection_eligibility",
        )
    ).lower()
    severe_terms = ["severe", "unsupported", "reject", "structure_risk", "mechanism_weak"]
    if any(term in class_text for term in severe_terms):
        failures.append("severe_or_unsupported_class")
    return failures


def sort_key(row):
    rank = as_int(row.get("selection_rank", ""), 10 ** 9)
    score = as_float(row.get("final_planning_score", ""), -10 ** 9)
    return (rank, -score, row.get("variant_id", ""))


def normalize_control_row(row, sequence_lookup, parents):
    out = dict(row)
    out.setdefault("source_variant_id", row.get("variant_id", ""))
    seq_source = ""
    seq_status = ""
    if row.get("variant_id", "") in sequence_lookup:
        seq, source_path = sequence_lookup[row.get("variant_id", "")]
        seq_source = source_path
        seq_status = "sequence_lookup_by_variant_id"
    else:
        seq, seq_status = apply_mutations(parents.get(row.get("target", ""), ""), row.get("mutation_list", ""))
        seq_source = "parent_sequence_plus_mutation_list"
    out["sequence"] = seq
    full_hash = sha256_text(seq) if seq else ""
    out["sequence_hash"] = full_hash[:10] if full_hash else ""
    out["canonical_sequence_hash_full"] = full_hash
    out["canonical_unique_key"] = "%s|%s" % (row.get("target", ""), full_hash) if full_hash else ""
    out["source_pool"] = "control_anchor_panel"
    out["source_pool_list"] = "control_anchor_panel"
    out["candidate_quality_tier"] = "control_anchor"
    out["expanded_pool_role"] = "control_anchor_panel"
    out["template_match_tier"] = "control_anchor"
    out["template_match_level"] = "control_anchor"
    out["complex_diagnostic_required"] = "True" if row.get("target") == "Ab_sdAb" else ""
    out["is_local_expansion"] = "False"
    out["is_seed_only"] = "False"
    out["is_four_mut"] = "True" if as_int(row.get("mutation_count", ""), 0) >= 4 else "False"
    out["is_neutral_boundary_or_high_risk"] = "False"
    out["contains_AY111H"] = "True" if "AY111H" in row.get("mutation_list", "") else "False"
    out["contains_AG102H_or_AV105H"] = (
        "True"
        if ("AG102H" in row.get("mutation_list", "") or "AV105H" in row.get("mutation_list", ""))
        else "False"
    )
    out["final_selection_role"] = "control_anchor_panel"
    out["selection_scenario"] = "synthesis_ready_draft_option_A_14220"
    out["sequence_source"] = seq_source
    out["sequence_reconstruction_status"] = seq_status
    return out


def prepare_final_rows(main_rows, control_rows, sequence_lookup, parents):
    controls_by_target = defaultdict(list)
    control_ids = set()
    control_unique_keys = set()
    control_sequence_hashes = set()
    for row in control_rows:
        norm = normalize_control_row(row, sequence_lookup, parents)
        controls_by_target[norm.get("target", "")].append(norm)
        control_ids.add(norm.get("variant_id", ""))
        if norm.get("canonical_unique_key", ""):
            control_unique_keys.add(norm.get("canonical_unique_key", ""))
        if norm.get("canonical_sequence_hash_full", ""):
            control_sequence_hashes.add(norm.get("canonical_sequence_hash_full", ""))

    main_by_target = defaultdict(list)
    for row in main_rows:
        if row.get("variant_id", "") in control_ids:
            continue
        if row.get("canonical_unique_key", "") in control_unique_keys:
            continue
        if row.get("canonical_sequence_hash_full", "") in control_sequence_hashes:
            continue
        main_by_target[row.get("target", "")].append(dict(row))

    original_counts = Counter(row.get("target", "") for row in main_rows)
    final_rows = []
    dropped_rows = []

    for target in TARGETS:
        target_total = original_counts[target]
        target_controls = controls_by_target[target]
        main_needed = target_total - len(target_controls)
        candidates = sorted(main_by_target[target], key=sort_key)
        selected = candidates[:main_needed]
        dropped = candidates[main_needed:]
        for row in selected:
            row["final_candidate_type"] = MAIN_ROLE
            row["final_synthesis_status"] = "draft_not_final_order_list"
            row["manual_review_required"] = "True"
            row["not_final_order_list"] = "True"
            row["final_capacity_policy"] = CAPACITY_POLICY
            row["sequence_source"] = "primary_15k_refill_v3_draft"
            row["sequence_reconstruction_status"] = "not_required"
            final_rows.append(row)
        for row in target_controls:
            row["final_candidate_type"] = CONTROL_ROLE
            row["final_synthesis_status"] = "draft_not_final_order_list"
            row["manual_review_required"] = "True"
            row["not_final_order_list"] = "True"
            row["final_capacity_policy"] = CAPACITY_POLICY
            final_rows.append(row)
        for row in dropped:
            row["drop_reason"] = "control_anchor_counted_inside_capacity"
            dropped_rows.append(row)

    final_rows = sorted(
        final_rows,
        key=lambda row: (
            row.get("target", ""),
            1 if row.get("final_candidate_type") == CONTROL_ROLE else 0,
            as_int(row.get("selection_rank", ""), 10 ** 9),
            row.get("variant_id", ""),
        ),
    )
    for i, row in enumerate(final_rows, 1):
        row["synthesis_draft_row_id"] = "SYN_DRAFT_%05d" % i
    return final_rows, dropped_rows


def count_by(rows, fields):
    counter = Counter()
    for row in rows:
        key = tuple(row.get(field, "") for field in fields)
        counter[key] += 1
    out = []
    for key, count in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
        row = {field: value for field, value in zip(fields, key)}
        row["count"] = count
        out.append(row)
    return out


def top_by(rows, target, field, limit=15):
    target_rows = [r for r in rows if r.get("target") == target]
    total = len(target_rows)
    counter = Counter(r.get(field, "") or "NA" for r in target_rows)
    out = []
    for value, count in counter.most_common(limit):
        out.append(
            {
                "target": target,
                field: value,
                "count": count,
                "fraction": "%.6f" % (float(count) / total if total else 0.0),
            }
        )
    return out


def audit_rows(rows, parents):
    target_summary = []
    failures_by_row = []
    duplicate_variant_count = 0
    duplicate_key_count = 0
    duplicate_sequence_count = 0
    seen_variant = Counter(r.get("variant_id", "") for r in rows)
    seen_key = Counter(r.get("canonical_unique_key", "") for r in rows)
    seen_seq = Counter(r.get("canonical_sequence_hash_full", "") for r in rows)
    duplicate_variant_count = sum(1 for k, c in seen_variant.items() if k and c > 1)
    duplicate_key_count = sum(1 for k, c in seen_key.items() if k and c > 1)
    duplicate_sequence_count = sum(1 for k, c in seen_seq.items() if k and c > 1)

    for row in rows:
        failures = classify_hard_failures(row, parents)
        if failures:
            failures_by_row.append(
                {
                    "variant_id": row.get("variant_id", ""),
                    "target": row.get("target", ""),
                    "final_candidate_type": row.get("final_candidate_type", ""),
                    "hard_failure_reasons": ";".join(failures),
                }
            )

    for target in TARGETS:
        trs = [r for r in rows if r.get("target") == target]
        total = len(trs)
        controls = [r for r in trs if r.get("final_candidate_type") == CONTROL_ROLE]
        main = [r for r in trs if r.get("final_candidate_type") == MAIN_ROLE]
        top_seed_count = Counter(r.get("his_seed_set", "") or "NA" for r in trs).most_common(1)
        top_cluster_count = Counter(r.get("near_duplicate_cluster_id", "") or "NA" for r in trs).most_common(1)
        sdab_diag = sum(1 for r in trs if r.get("target") == "Ab_sdAb" and truthy(r.get("complex_diagnostic_required", "")))
        target_summary.append(
            {
                "target": target,
                "total_count": total,
                "main_count": len(main),
                "control_anchor_count": len(controls),
                "top_seed": top_seed_count[0][0] if top_seed_count else "",
                "top_seed_count": top_seed_count[0][1] if top_seed_count else 0,
                "top_seed_fraction": "%.6f" % (float(top_seed_count[0][1]) / total if total and top_seed_count else 0.0),
                "top_cluster": top_cluster_count[0][0] if top_cluster_count else "",
                "top_cluster_count": top_cluster_count[0][1] if top_cluster_count else 0,
                "top_cluster_fraction": "%.6f" % (float(top_cluster_count[0][1]) / total if total and top_cluster_count else 0.0),
                "four_mut_count": sum(1 for r in trs if truthy(r.get("is_four_mut", ""))),
                "four_mut_fraction": "%.6f" % (float(sum(1 for r in trs if truthy(r.get("is_four_mut", "")))) / total if total else 0.0),
                "local_expansion_count": sum(1 for r in trs if truthy(r.get("is_local_expansion", ""))),
                "local_expansion_fraction": "%.6f" % (float(sum(1 for r in trs if truthy(r.get("is_local_expansion", "")))) / total if total else 0.0),
                "seed_only_count": sum(1 for r in trs if truthy(r.get("is_seed_only", ""))),
                "seed_only_fraction": "%.6f" % (float(sum(1 for r in trs if truthy(r.get("is_seed_only", "")))) / total if total else 0.0),
                "neutral_boundary_or_high_risk_count": sum(1 for r in trs if truthy(r.get("is_neutral_boundary_or_high_risk", ""))),
                "neutral_boundary_or_high_risk_fraction": "%.6f" % (float(sum(1 for r in trs if truthy(r.get("is_neutral_boundary_or_high_risk", "")))) / total if total else 0.0),
                "AG102H_AV105H_main_count": sum(
                    1
                    for r in main
                    if target == "Ab_sdAb" and truthy(r.get("contains_AG102H_or_AV105H", ""))
                ),
                "sdAb_complex_diagnostic_required_count": sdab_diag if target == "Ab_sdAb" else "",
            }
        )
    hard_counts = {
        "duplicate_variant_id_count": duplicate_variant_count,
        "duplicate_canonical_unique_key_count": duplicate_key_count,
        "duplicate_canonical_sequence_hash_count": duplicate_sequence_count,
        "hard_failure_row_count": len(failures_by_row),
        "glycan_low_risk_claim_count": sum(
            1
            for row in rows
            if "low" in str(row.get("glycan_status", "") + row.get("glycan_risk_status", "") + row.get("glycan_claim", "")).lower()
        ),
    }
    return target_summary, failures_by_row, hard_counts


def md_table(rows, columns):
    lines = []
    lines.append("| " + " | ".join(columns) + " |")
    lines.append("| " + " | ".join(["---"] * len(columns)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def write_audit_report(path, manifest_rows, target_summary, hard_counts, failures, rows, dropped_rows):
    verdict = "FAIL" if failures or hard_counts["duplicate_canonical_unique_key_count"] or hard_counts["duplicate_canonical_sequence_hash_count"] or hard_counts["duplicate_variant_id_count"] or hard_counts["glycan_low_risk_claim_count"] else "PASS"
    if verdict == "PASS" and len(rows) < 15000:
        verdict = "PASS_OPTION_A_REALISTIC_CAPACITY"
    source_summary = count_by(rows, ["target", "source_pool", "candidate_quality_tier", "template_match_tier"])
    evidence_summary = count_by(rows, ["target", "final_candidate_type", "candidate_quality_tier", "assigned_template_status"])
    seed_rows = []
    cluster_rows = []
    for target in TARGETS:
        seed_rows.extend(top_by(rows, target, "his_seed_set", 12))
        cluster_rows.extend(top_by(rows, target, "near_duplicate_cluster_id", 12))
    control_summary = count_by(
        [r for r in rows if r.get("final_candidate_type") == CONTROL_ROLE],
        ["target", "his_seed_set", "mutation_count", "sequence_reconstruction_status"],
    )
    with open(path, "w") as handle:
        handle.write("# Synthesis-ready candidate list draft audit\n\n")
        handle.write("## Executive Summary\n\n")
        handle.write("Verdict: `%s`.\n\n" % verdict)
        handle.write(
            "This stage built a synthesis-ready review draft under Option A: "
            "14,220 total rows with controls counted inside capacity. "
            "It did not generate new variants, did not run AF3, SimpleFold, "
            "PyRosetta, FoldX, MD, or glycan modeling, and did not create a final order list.\n\n"
        )
        handle.write("All outputs remain `draft`, `not final order list`, and `manual review required`.\n\n")
        handle.write("## Input Manifest\n\n")
        handle.write(md_table(manifest_rows, ["input_name", "path", "rows", "sha256", "role"]))
        handle.write("\n\n## Capacity And Target Allocation\n\n")
        handle.write(md_table(target_summary, [
            "target", "total_count", "main_count", "control_anchor_count",
            "top_seed", "top_seed_count", "top_seed_fraction",
            "top_cluster", "top_cluster_count", "top_cluster_fraction",
            "four_mut_count", "four_mut_fraction",
            "local_expansion_count", "local_expansion_fraction",
            "seed_only_count", "seed_only_fraction",
            "neutral_boundary_or_high_risk_count", "neutral_boundary_or_high_risk_fraction",
            "AG102H_AV105H_main_count", "sdAb_complex_diagnostic_required_count",
        ]))
        handle.write("\n\n## Hard Buildability / Identity Checks\n\n")
        hard_rows = [{"metric": key, "value": value} for key, value in sorted(hard_counts.items())]
        handle.write(md_table(hard_rows, ["metric", "value"]))
        handle.write("\n\n")
        if failures:
            handle.write("### Hard failure rows\n\n")
            handle.write(md_table(failures[:50], ["variant_id", "target", "final_candidate_type", "hard_failure_reasons"]))
            handle.write("\n\n")
        else:
            handle.write("No hard buildability / identity failure rows were detected.\n\n")
        handle.write("## Evidence / Quality Summary\n\n")
        handle.write(md_table(evidence_summary[:40], ["target", "final_candidate_type", "candidate_quality_tier", "assigned_template_status", "count"]))
        handle.write("\n\n## Source Summary\n\n")
        handle.write(md_table(source_summary[:50], ["target", "source_pool", "candidate_quality_tier", "template_match_tier", "count"]))
        handle.write("\n\n## Seed Summary\n\n")
        handle.write(md_table(seed_rows, ["target", "his_seed_set", "count", "fraction"]))
        handle.write("\n\n## Cluster Summary\n\n")
        handle.write(md_table(cluster_rows, ["target", "near_duplicate_cluster_id", "count", "fraction"]))
        handle.write("\n\n## Control / Anchor Summary\n\n")
        handle.write(md_table(control_summary, ["target", "his_seed_set", "mutation_count", "sequence_reconstruction_status", "count"]))
        handle.write("\n\n## Dropped Main Rows\n\n")
        handle.write(
            "Dropped main rows to count controls inside capacity: `%d`.\n\n" % len(dropped_rows)
        )
        drop_summary = count_by(dropped_rows, ["target", "drop_reason"])
        handle.write(md_table(drop_summary, ["target", "drop_reason", "count"]))
        handle.write("\n\n## Conservative Labels\n\n")
        handle.write("- sdAb remains secondary / complex-diagnostic-required.\n")
        handle.write("- Glycan remains unchecked; glycan low-risk claim count must remain 0.\n")
        handle.write("- This report is not a wet-lab order list.\n")
        handle.write("- Manual review is required before any synthesis-ready final list.\n")
    return verdict


def write_control_plan(path, controls):
    summary = count_by(controls, ["target", "his_seed_set", "mutation_count", "sequence_reconstruction_status"])
    with open(path, "w") as handle:
        handle.write("# Control / anchor final plan\n\n")
        handle.write("Status: `draft_not_final_order_list`; manual review required.\n\n")
        handle.write("Controls / anchors are counted inside the 14,220-row Option A draft capacity.\n\n")
        handle.write(md_table(summary, ["target", "his_seed_set", "mutation_count", "sequence_reconstruction_status", "count"]))
        handle.write("\n\nInterpretation:\n\n")
        handle.write("- 1E62 controls calibrate light-chain His seed and prior risk/anchor behavior.\n")
        handle.write("- sdAb controls retain AD110H-rich calibration plus AG102H / AV105H risk controls.\n")
        handle.write("- AG102H / AV105H remain control/risk rows and are not main sdAb candidates.\n")
        handle.write("- Controls are not evidence of final wet-lab success; they are calibration rows.\n")


def write_wetlab_plan(path, rows):
    target_counts = Counter(r.get("target", "") for r in rows)
    control_counts = Counter(r.get("target", "") for r in rows if r.get("final_candidate_type") == CONTROL_ROLE)
    with open(path, "w") as handle:
        handle.write("# Wet-lab interpretation plan\n\n")
        handle.write("Status: `draft`; this is not a final order list.\n\n")
        handle.write("The draft should be interpreted as a synthesis-ready review object, not as a proven final library.\n\n")
        handle.write("## Library composition\n\n")
        handle.write("- Total draft rows: `%d`.\n" % len(rows))
        for target in TARGETS:
            handle.write(
                "- %s: `%d` rows, including `%d` controls / anchors.\n"
                % (target, target_counts[target], control_counts[target])
            )
        handle.write("\n## Interpretation rules\n\n")
        handle.write("- 1E62 remains the primary branch.\n")
        handle.write("- sdAb remains secondary and complex-diagnostic-required; do not describe it as structure-confirmed.\n")
        handle.write("- Glycan status remains unchecked; do not describe any row as low glycan risk.\n")
        handle.write("- Controls / anchors should be interpreted separately from main candidates.\n")
        handle.write("- Final wet-lab ordering requires manual approval after reviewing this draft audit.\n")


def build_manifest(main_rows, control_rows):
    return [
        {
            "input_name": "primary_15k_refill_v3_draft",
            "path": os.path.relpath(PRIMARY_DRAFT, ROOT),
            "rows": len(main_rows),
            "sha256": sha256_file(PRIMARY_DRAFT),
            "role": "primary_review_input",
        },
        {
            "input_name": "backup_10k_refill_v3_draft",
            "path": os.path.relpath(BACKUP_DRAFT, ROOT),
            "rows": sum(1 for _ in open(BACKUP_DRAFT)) - 1,
            "sha256": sha256_file(BACKUP_DRAFT),
            "role": "backup_lower_capacity_input",
        },
        {
            "input_name": "control_anchor_panel",
            "path": os.path.relpath(CONTROL_PANEL, ROOT),
            "rows": len(control_rows),
            "sha256": sha256_file(CONTROL_PANEL),
            "role": "controls_counted_inside_capacity",
        },
    ]


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)
    main_fields, main_rows = read_csv(PRIMARY_DRAFT)
    control_fields, control_rows = read_csv(CONTROL_PANEL)
    backup_fields, backup_rows = read_csv(BACKUP_DRAFT)
    parents = derive_parent_sequences(main_rows)
    sequence_lookup, lookup_source_hits = build_sequence_lookup(SEQUENCE_LOOKUP_SOURCES)
    final_rows, dropped_rows = prepare_final_rows(main_rows, control_rows, sequence_lookup, parents)
    target_summary, failures, hard_counts = audit_rows(final_rows, parents)
    manifest_rows = build_manifest(main_rows, control_rows)

    meta_fields = [
        "synthesis_draft_row_id",
        "final_synthesis_status",
        "final_candidate_type",
        "final_capacity_policy",
        "manual_review_required",
        "not_final_order_list",
        "sequence_source",
        "sequence_reconstruction_status",
    ]
    all_fields = []
    for field in meta_fields + main_fields + control_fields:
        if field not in all_fields:
            all_fields.append(field)
    for field in [
        "sequence",
        "sequence_hash",
        "canonical_sequence_hash_full",
        "canonical_unique_key",
        "source_pool",
        "candidate_quality_tier",
        "complex_diagnostic_required",
    ]:
        if field not in all_fields:
            all_fields.append(field)

    draft_csv = os.path.join(OUT_DIR, "synthesis_ready_candidate_list_draft.csv")
    write_csv(draft_csv, all_fields, final_rows)

    write_csv(
        os.path.join(OUT_DIR, "final_selection_input_manifest.csv"),
        ["input_name", "path", "rows", "sha256", "role"],
        manifest_rows,
    )
    write_csv(
        os.path.join(OUT_DIR, "synthesis_ready_target_audit.csv"),
        list(target_summary[0].keys()),
        target_summary,
    )
    write_csv(
        os.path.join(OUT_DIR, "synthesis_ready_hard_failure_rows.csv"),
        ["variant_id", "target", "final_candidate_type", "hard_failure_reasons"],
        failures,
    )
    write_csv(
        os.path.join(OUT_DIR, "synthesis_ready_dropped_main_rows.csv"),
        sorted(set(k for row in dropped_rows for k in row.keys())),
        dropped_rows,
    )

    audit_path = os.path.join(OUT_DIR, "synthesis_ready_candidate_list_audit.md")
    verdict = write_audit_report(
        audit_path, manifest_rows, target_summary, hard_counts, failures, final_rows, dropped_rows
    )
    controls_final = [r for r in final_rows if r.get("final_candidate_type") == CONTROL_ROLE]
    write_control_plan(os.path.join(OUT_DIR, "control_anchor_final_plan.md"), controls_final)
    write_wetlab_plan(os.path.join(OUT_DIR, "wetlab_interpretation_plan.md"), final_rows)

    with open(audit_path) as src, open(TASK_REPORT, "w") as dst:
        dst.write(src.read())

    print(verdict)
    print(os.path.relpath(OUT_DIR, ROOT))


if __name__ == "__main__":
    main()
