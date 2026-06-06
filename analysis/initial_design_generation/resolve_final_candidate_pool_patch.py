#!/usr/bin/env python3
"""Build PATCH-resolution final candidate-pool drafts.

This keeps the clean 1E62 subset obtained under the original default quotas and
fills the remaining draft capacity with sdAb. It is a review draft only, not a
synthesis unlock.
"""

import csv
import hashlib
import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BUILDER_PATH = ROOT / "analysis/initial_design_generation/build_final_candidate_pool_drafts.py"
OUT_DIR = ROOT / "results/initial_design_generation/final_candidate_pool_patch_resolution"
CURRENT_STAGE_REPORT = ROOT / ".tasks/active/initial-design-generation/current_stage_report.md"


spec = importlib.util.spec_from_file_location("final_builder", str(BUILDER_PATH))
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


SCENARIOS = {
    "10k_patch_resolved": {
        "total": 10000,
        "base_1e62_quota": 6000,
        "target": "10K",
    },
    "15k_patch_resolved": {
        "total": 15000,
        "base_1e62_quota": 8500,
        "target": "15K",
    },
}


def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


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


def md_table(rows, fields):
    lines = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(f, "")) for f in fields) + " |")
    return "\n".join(lines)


def build_resolution(rows, name, config):
    one_e62, excl1, audit1 = builder.select_target(rows, "Ab_1E62", config["base_1e62_quota"], name + "_1e62_clean_basis")
    remaining = config["total"] - len(one_e62)
    sdab, excl2, audit2 = builder.select_target(rows, "Ab_sdAb", remaining, name + "_sdab_topup")
    selected = []
    for row in one_e62:
        out = dict(row)
        out["selection_scenario"] = name
        out["resolution_role"] = "accepted_1E62_clean_undersized_subset"
        out["basis_quota"] = config["base_1e62_quota"]
        out["effective_target_quota"] = len(one_e62)
        selected.append(out)
    for row in sdab:
        out = dict(row)
        out["selection_scenario"] = name
        out["resolution_role"] = "sdAb_topup_to_total_size"
        out["basis_quota"] = remaining
        out["effective_target_quota"] = remaining
        selected.append(out)
    for idx, row in enumerate(selected, start=1):
        row["combined_selection_rank"] = idx

    # Report 1E62 against the original basis quota, because this resolution
    # explicitly accepts the clean undersized subset produced under that quota
    # instead of rebasing caps downward and incorrectly turning the subset into
    # a hard-cap failure.
    audit1_actual = dict(audit1)
    audit1_actual["selection_scenario"] = name
    audit2_actual = builder.audit_target([r for r in selected if r["target"] == "Ab_sdAb"], "Ab_sdAb", remaining, name)
    for audit, role, original in [
        (audit1_actual, "accepted_1E62_clean_undersized_subset", audit1),
        (audit2_actual, "sdAb_topup_to_total_size", audit2),
    ]:
        audit["resolution_role"] = role
        audit["basis_quota"] = config["base_1e62_quota"] if audit["target"] == "Ab_1E62" else remaining
        audit["original_basis_verdict"] = original["verdict"]
        audit["original_basis_reasons"] = original["reasons"]
        if audit["target"] == "Ab_1E62":
            audit["manual_waiver_required"] = True
            audit["waiver_reason"] = "accepted_undersized_clean_1E62_subset_instead_of_relaxing_hard_caps"
        else:
            audit["manual_waiver_required"] = audit["verdict"] != "PASS"
            audit["waiver_reason"] = "top_seed_soft_cap_or_hard_phase_review" if audit["verdict"] != "PASS" else ""

    excluded = []
    for row in excl1 + excl2:
        row = dict(row)
        row["selection_scenario"] = name
        excluded.append(row)
    return selected, excluded, [audit1_actual, audit2_actual]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = builder.add_features(builder.read_csv(builder.INPUT_POOL))
    original_fields = list(rows[0].keys())
    extra = [
        "selection_scenario",
        "resolution_role",
        "basis_quota",
        "effective_target_quota",
        "combined_selection_rank",
    ]
    fieldnames = original_fields + [x for x in extra if x not in original_fields]

    all_audits = []
    all_excluded = []
    selected_by_name = {}
    for name, config in SCENARIOS.items():
        selected, excluded, audits = build_resolution(rows, name, config)
        selected_by_name[name] = selected
        all_excluded.extend(excluded)
        all_audits.extend(audits)
        write_csv(OUT_DIR / ("%s_draft.csv" % name), selected, fieldnames)

    write_csv(OUT_DIR / "patch_resolution_audit_checks.csv", all_audits)
    write_csv(OUT_DIR / "patch_resolution_excluded_candidates.csv", all_excluded[:60000])

    target_rows = []
    for name, selected in selected_by_name.items():
        counts = {}
        for row in selected:
            counts[row["target"]] = counts.get(row["target"], 0) + 1
        for target, count in sorted(counts.items()):
            target_rows.append({"selection_scenario": name, "target": target, "count": count})

    audit_fields = [
        "selection_scenario",
        "target",
        "resolution_role",
        "selected_count",
        "target_quota",
        "basis_quota",
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
        "AY111H_containing_count",
        "AY111H_containing_fraction",
        "AG102H_AV105H_main_count",
        "sdAb_complex_diagnostic_required_count",
        "canonical_duplicate_count",
        "verdict",
        "manual_waiver_required",
        "waiver_reason",
        "reasons",
    ]
    overall = "PATCH_RESOLUTION_DRAFTS_BUILT__MANUAL_REVIEW_REQUIRED__FINAL_SELECTION_LOCKED"
    report = """# Final Candidate Pool PATCH Resolution

## Executive Summary

Verdict: `%s`.

This stage resolves the previous size shortfall by accepting the clean 1E62 subset selected under the original default quota basis, then topping up the remaining draft capacity with sdAb. It does not relax hard filters, does not generate new variants, and does not unlock final synthesis selection.

## Source

| input | sha256 |
| --- | --- |
| final_candidate_pool_planning_audit_report.md | %s |
| expanded_tier2_candidate_pool_v2.csv | %s |

## Target Counts

%s

## Audit

%s

## Interpretation

- 1E62 remains the primary branch, but the current 20K mother pool cannot satisfy the original 1E62 count targets under strict hard caps.
- The resolution draft accepts the clean undersized 1E62 subset instead of filling with extra 4-mut, local-expansion, seed-only, or dominant-cluster rows.
- sdAb fills the remaining slots and remains secondary / complex-diagnostic-required.
- This is a manual-review draft only; final synthesis-ready selection remains locked.

## Next Gate

Allowed now:

- Review whether the adjusted target allocation is acceptable.
- Decide whether to accept these PATCH-resolution drafts or request targeted 1E62 refill / Backfill v3.

Still locked:

- Final synthesis-ready 10K / 15K selection.
- Wet-lab order list generation.
- New AF3 / SimpleFold / PyRosetta / FoldX compute.
- MD / constant-pH MD.
- sdAb structure-confirmed upgrade.
- Glycan low-risk claim.
""" % (
        overall,
        sha256_file(ROOT / "results/initial_design_generation/final_candidate_pool_planning/final_candidate_pool_planning_audit_report.md"),
        sha256_file(builder.INPUT_POOL),
        md_table(target_rows, ["selection_scenario", "target", "count"]),
        md_table(all_audits, audit_fields),
    )
    (OUT_DIR / "patch_resolution_audit_report.md").write_text(report, encoding="utf-8")
    (OUT_DIR / "patch_resolution_next_gate_report.md").write_text(report, encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(report, encoding="utf-8")
    print(overall)
    print(OUT_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
