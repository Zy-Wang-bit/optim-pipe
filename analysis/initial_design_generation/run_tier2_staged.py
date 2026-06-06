#!/usr/bin/env python3
"""Build staged Tier 2 inputs before structural compute.

This script is intentionally table-only. It does not run PyRosetta, pKa,
FoldX, AF3, SimpleFold, or MD.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
TIER1_DIR = ROOT / "results/initial_design_generation/tier1_filtering"
DEFAULT_OUT = ROOT / "results/initial_design_generation/tier2_staged"
CONTROL_PANEL = ROOT / "results/initial_design_generation/dry_run/control_anchor_panel_draft.csv"

PRIMARY = TIER1_DIR / "tier2_core_primary_proposal.csv"
RESERVE = TIER1_DIR / "tier2_core_reserve_pool.csv"
FEATURE_TABLE = TIER1_DIR / "tier1_feature_table.csv"
CAP_POLICY = TIER1_DIR / "tier1_cap_policy.yaml"
TIER1_MANIFEST = TIER1_DIR / "tier1_filtering_manifest.yaml"

BACKBONES = {
    "Ab_1E62": {
        "display_name": "1E62",
        "window_id": "1E62_VL_001_040",
        "pdb": ROOT / "results/initial_design_generation/p0_mpnn_backbones/af3_1E62_seed7_model5_rank1.pdb",
        "design_chain": "L",
        "mutation_chain_prefixes": {"L"},
        "protected": {23},
        "stage1_size": 3000,
        "f_cap": 600,
        "seed_cap_fraction": 0.12,
        "cluster_cap_fraction": 0.02,
        "foldx_subset_cap": 1000,
    },
    "Ab_sdAb": {
        "display_name": "sdAb",
        "window_id": "sdAb_VHH_072_111",
        "pdb": ROOT / "results/initial_design_generation/p0_mpnn_backbones/af3_sdAb_seed16_model5_rank1.pdb",
        "design_chain": "A",
        "mutation_chain_prefixes": {"A"},
        "protected": {96},
        "stage1_size": 4000,
        "f_cap": 600,
        "seed_cap_fraction": 0.10,
        "cluster_cap_fraction": 0.02,
        "foldx_subset_cap": 1200,
    },
}

EXPECTED_PRIMARY = {"Ab_1E62": 12000, "Ab_sdAb": 15000}
EXPECTED_RESERVE = {"Ab_1E62": 8000, "Ab_sdAb": 10000}
EXPECTED_FEATURE = {"Ab_1E62": 120165, "Ab_sdAb": 150041}

MUT_RE = re.compile(r"^([A-Za-z])([A-Za-z])(\d+)([A-Za-z])$")
AA3 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy() if max_rows else df.copy()
    lines = [
        "| " + " | ".join(map(str, view.columns)) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]).replace("\n", " ") for c in view.columns) + " |")
    return "\n".join(lines)


def parse_mutations(value: object) -> list[tuple[str, str, int, str]]:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return []
    out: list[tuple[str, str, int, str]] = []
    for token in str(value).replace(",", ";").split(";"):
        token = token.strip()
        if not token or token.lower() == "nan":
            continue
        m = MUT_RE.match(token)
        if not m:
            raise ValueError(f"Cannot parse mutation token: {token}")
        chain, old, pos, new = m.groups()
        out.append((chain, old, int(pos), new))
    return out


def pdb_chain_sequence(path: Path, chain_id: str) -> str:
    residues: list[tuple[int, str]] = []
    seen: set[tuple[str, int, str]] = set()
    with path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip()
            if chain != chain_id:
                continue
            resseq = int(line[22:26])
            icode = line[26].strip()
            key = (chain, resseq, icode)
            if key in seen:
                continue
            seen.add(key)
            residues.append((resseq, AA3.get(line[17:20].strip(), "X")))
    residues.sort(key=lambda x: x[0])
    return "".join(aa for _, aa in residues)


def validate_rows(df: pd.DataFrame, source_name: str) -> list[dict[str, object]]:
    checks: list[dict[str, object]] = []
    required = [
        "target",
        "variant_id",
        "sequence",
        "canonical_sequence_hash_full",
        "canonical_mutated_window_sequence",
        "normalized_mutation_list",
        "tier1_review_class",
        "tier1_rank_score",
        "near_duplicate_cluster_id",
        "his_seed_set",
    ]
    for col in required:
        checks.append(
            {
                "scope": source_name,
                "check": f"required_column:{col}",
                "status": "PASS" if col in df.columns else "FAIL",
                "details": "" if col in df.columns else "missing",
            }
        )
    if "variant_id" in df.columns:
        dup = int(df["variant_id"].duplicated().sum())
        checks.append({"scope": source_name, "check": "variant_id_unique", "status": "PASS" if dup == 0 else "FAIL", "details": dup})
    if "canonical_sequence_hash_full" in df.columns:
        dup = int(df["canonical_sequence_hash_full"].duplicated().sum())
        checks.append(
            {
                "scope": source_name,
                "check": "canonical_sequence_hash_full_unique",
                "status": "PASS" if dup == 0 else "FAIL",
                "details": dup,
            }
        )
    return checks


def target_from_control(value: str) -> str:
    if value == "1E62":
        return "Ab_1E62"
    if value == "sdAb":
        return "Ab_sdAb"
    return value


def build_control_rows(reference_cols: list[str]) -> pd.DataFrame:
    if not CONTROL_PANEL.exists():
        return pd.DataFrame(columns=reference_cols)
    controls = pd.read_csv(CONTROL_PANEL)
    rows: list[dict[str, object]] = []
    for _, row in controls.iterrows():
        target = target_from_control(str(row["target"]))
        if target not in BACKBONES:
            continue
        info = BACKBONES[target]
        seq = pdb_chain_sequence(info["pdb"], info["design_chain"])
        muts = "" if pd.isna(row.get("mutation_list")) else str(row.get("mutation_list"))
        try:
            for chain, old, pos, new in parse_mutations(muts):
                if chain not in info["mutation_chain_prefixes"]:
                    raise ValueError(f"control chain {chain} incompatible with {target}")
                if pos < 1 or pos > len(seq):
                    raise ValueError(f"control position {pos} out of range for {target}")
                if seq[pos - 1] != old:
                    raise ValueError(f"control parent mismatch {target} {chain}{pos}: expected {old}, found {seq[pos - 1]}")
                seq = seq[: pos - 1] + new + seq[pos:]
        except Exception as exc:
            rows.append({"target": target, "variant_id": row["variant_id"], "tier2_control_error": str(exc)})
            continue
        out = {col: pd.NA for col in reference_cols}
        out.update(
            {
                "target": target,
                "variant_id": str(row["variant_id"]),
                "sequence_hash": str(row.get("sequence_hash", "")),
                "mutation_list": muts,
                "normalized_mutation_list": muts,
                "mutation_count": len(parse_mutations(muts)),
                "His_count": sum(1 for _, _, _, new in parse_mutations(muts) if new == "H"),
                "his_seed_set": muts,
                "primary_generation_route": "control_anchor",
                "all_source_routes": "control_anchor",
                "rescue_mutation_list": "",
                "rescue_signature": "",
                "rescue_count": 0,
                "hard_filter_status": "pass",
                "liability_flags": "",
                "near_duplicate_cluster_id": f"{target}_control_{row['variant_id']}",
                "cluster_size": 1,
                "hit_likelihood_score_v0": 0.0,
                "neutral_retention_score": 1.0 if muts == "" else 0.75,
                "acidic_release_support_score": 0.0 if muts == "" else 0.55,
                "global_weakening_risk_score": 0.0 if muts == "" else 0.25,
                "glycan_or_epitope_risk_score": 0.0,
                "foldx_proxy_feature": 0.0,
                "interface_hotspot_proxy": "control_anchor",
                "buildability_light_status": row.get("buildability_light_status", "pass"),
                "mpnn_score_status": "not_applicable_control",
                "source_universe": "control_anchor_panel",
                "source_branch": "control_anchor",
                "production_pool_overlap_status": "control_anchor",
                "selection_eligibility": "control_only",
                "mpnn_score_type": "not_applicable",
                "p0_evidence_class": "control_anchor",
                "p0_evidence_flag": "control_anchor",
                "tier1_review_status": "control_anchor",
                "generation_record_id": str(row["variant_id"]),
                "sequence": seq,
                "window": info["window_id"],
                "window_id": info["window_id"],
                "canonical_sequence_hash_full": hashlib.sha256(seq.encode()).hexdigest(),
                "canonical_sequence_hash_short": hashlib.sha256(seq.encode()).hexdigest()[:12],
                "sequence_hash_length": 64,
                "input_sequence_hash_length": len(str(row.get("sequence_hash", ""))),
                "canonical_mutated_window_sequence": seq,
                "generation_record_ids": str(row["variant_id"]),
                "mpnn_sparse_overlap": False,
                "mpnn_sparse_novel": False,
                "source_admission_status": "control_anchor",
                "tier2_eligibility_status": "control_anchor",
                "sequence_hamming_cluster_id": f"{target}_control_{row['variant_id']}",
                "mutation_set_jaccard_cluster_id": f"{target}_control_{row['variant_id']}",
                "his_seed_cluster_id": f"{target}_control_{muts or 'WT'}",
                "mpnn_score_comparable": False,
                "mpnn_score_direction": "not_applicable",
                "mpnn_score_lower_is_better": False,
                "mpnn_score_comparable_group": "control_anchor",
                "acidic_release_support_score_v0": 0.0 if muts == "" else 0.55,
                "neutral_retention_score_v0": 1.0 if muts == "" else 0.75,
                "global_weakening_risk_score_v0": 0.0 if muts == "" else 0.25,
                "liability_risk_score_v0": 0.0,
                "tier1_acidic_release_support_flag": "control",
                "tier1_neutral_retention_flag": "control",
                "tier1_global_weakening_risk_flag": "control",
                "tier1_liability_flag": "pass",
                "tier1_mpnn_support_flag": "not_applicable",
                "tier1_discordance_flag": "control",
                "feature_missing_reason": "control_anchor_not_production_candidate",
                "feature_confidence": "control_anchor",
                "tier1_feature_completeness_class": "control_anchor",
                "tier1_review_class": "control_anchor",
                "tier1_review_reason": str(row.get("control_type", "control_anchor")),
                "tier1_review_confidence": "control",
                "tier1_review_class_priority": -1,
                "tier1_rank_score": 1.0,
                "tier1_target_rank": 0,
                "tier2_proposal_role": "control_anchor",
                "tier2_source_pool": "control_anchor",
                "control_type": str(row.get("control_type", "control_anchor")),
                "tier2_control_error": "",
            }
        )
        rows.append(out)
    return pd.DataFrame(rows, columns=sorted(set(reference_cols) | {"control_type", "tier2_control_error"}))


def add_mapping_checks(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = df.copy()
    records: list[dict[str, object]] = []
    parent_sequences = {
        target: pdb_chain_sequence(info["pdb"], info["design_chain"]) for target, info in BACKBONES.items()
    }
    statuses: list[str] = []
    reasons: list[str] = []
    for row in out.itertuples(index=False):
        target = str(getattr(row, "target"))
        vid = str(getattr(row, "variant_id"))
        seq = str(getattr(row, "sequence", ""))
        info = BACKBONES.get(target)
        row_reasons: list[str] = []
        if info is None:
            row_reasons.append("unknown_target")
        else:
            if not info["pdb"].exists():
                row_reasons.append("backbone_missing")
            parent_seq = parent_sequences[target]
            if not seq or seq == "nan":
                row_reasons.append("sequence_missing")
            muts_value = getattr(row, "normalized_mutation_list", getattr(row, "mutation_list", ""))
            try:
                muts = parse_mutations(muts_value)
                for chain, old, pos, new in muts:
                    if chain not in info["mutation_chain_prefixes"]:
                        row_reasons.append(f"invalid_chain:{chain}")
                    if pos in info["protected"]:
                        row_reasons.append(f"protected_residue_mutated:{chain}{pos}")
                    if pos < 1 or pos > len(parent_seq):
                        row_reasons.append(f"position_out_of_range:{chain}{pos}")
                    elif parent_seq[pos - 1] != old:
                        row_reasons.append(f"parent_aa_mismatch:{chain}{old}{pos}{new}:pdb={parent_seq[pos - 1]}")
                    if seq and seq != "nan" and pos <= len(seq) and seq[pos - 1] != new:
                        row_reasons.append(f"sequence_mutation_mismatch:{chain}{old}{pos}{new}:seq={seq[pos - 1]}")
            except Exception as exc:
                row_reasons.append(f"mutation_parse_error:{exc}")
        status = "PASS" if not row_reasons else "FAIL"
        statuses.append(status)
        reasons.append(";".join(row_reasons))
        if row_reasons:
            records.append({"variant_id": vid, "target": target, "status": status, "reason": ";".join(row_reasons)})
    out["tier2_mapping_status"] = statuses
    out["tier2_mapping_reason"] = reasons
    return out, pd.DataFrame(records)


def freeze_inputs(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    primary = pd.read_csv(PRIMARY)
    reserve = pd.read_csv(RESERVE)
    feature = pd.read_csv(FEATURE_TABLE, usecols=["target", "variant_id", "canonical_sequence_hash_full"])

    checks: list[dict[str, object]] = []
    checks.extend(validate_rows(primary, "primary"))
    checks.extend(validate_rows(reserve, "reserve"))
    for label, df, expected in [("primary", primary, EXPECTED_PRIMARY), ("reserve", reserve, EXPECTED_RESERVE), ("feature", feature, EXPECTED_FEATURE)]:
        counts = df["target"].value_counts().to_dict()
        for target, n in expected.items():
            got = int(counts.get(target, 0))
            checks.append({"scope": label, "check": f"expected_count:{target}", "status": "PASS" if got == n else "FAIL", "details": f"{got}/{n}"})

    primary["tier2_source_pool"] = "primary"
    reserve["tier2_source_pool"] = "reserve"
    snapshot = pd.concat([primary, reserve], ignore_index=True, sort=False)
    controls = build_control_rows(snapshot.columns.tolist())
    if not controls.empty:
        snapshot = pd.concat([snapshot, controls], ignore_index=True, sort=False)
    snapshot, mapping_failures = add_mapping_checks(snapshot)
    if not mapping_failures.empty:
        checks.append({"scope": "snapshot", "check": "chain_coordinate_residue_mapping", "status": "FAIL", "details": len(mapping_failures)})
    else:
        checks.append({"scope": "snapshot", "check": "chain_coordinate_residue_mapping", "status": "PASS", "details": len(snapshot)})

    duplicate_hashes = int(snapshot["canonical_sequence_hash_full"].duplicated().sum())
    checks.append(
        {
            "scope": "snapshot",
            "check": "canonical_sequence_hash_full_unique_all_inputs",
            "status": "PASS" if duplicate_hashes == 0 else "WARN",
            "details": duplicate_hashes,
        }
    )
    validation = pd.DataFrame(checks)

    write_csv(snapshot, out_dir / "tier2_candidate_snapshot.csv")
    write_csv(validation, out_dir / "tier2_input_validation.csv")
    write_csv(mapping_failures, out_dir / "tier2_mapping_failures.csv")

    verdict = "FAIL" if validation["status"].eq("FAIL").any() else ("WARN" if validation["status"].eq("WARN").any() else "PASS")
    manifest = {
        "version": 1,
        "status": verdict,
        "tier2_compute_started": False,
        "final_10k_selection_started": False,
        "inputs": {
            str(PRIMARY): file_sha(PRIMARY),
            str(RESERVE): file_sha(RESERVE),
            str(FEATURE_TABLE): file_sha(FEATURE_TABLE),
            str(CAP_POLICY): file_sha(CAP_POLICY) if CAP_POLICY.exists() else None,
            str(TIER1_MANIFEST): file_sha(TIER1_MANIFEST) if TIER1_MANIFEST.exists() else None,
        },
        "row_counts": {
            "snapshot": int(len(snapshot)),
            "primary": int((snapshot["tier2_source_pool"] == "primary").sum()),
            "reserve": int((snapshot["tier2_source_pool"] == "reserve").sum()),
            "control_anchor": int((snapshot["tier2_source_pool"] == "control_anchor").sum()) if "tier2_source_pool" in snapshot else 0,
        },
    }
    (out_dir / "tier2_input_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = [
        "# Tier2 Input Freeze Report",
        "",
        f"Verdict: `{verdict}`",
        "",
        "## Validation",
        markdown_table(validation),
        "",
        "## Counts",
        markdown_table(snapshot.groupby(["target", "tier2_source_pool"], dropna=False).size().reset_index(name="count")),
    ]
    if not mapping_failures.empty:
        report += ["", "## Mapping Failures", markdown_table(mapping_failures, max_rows=50)]
    (out_dir / "tier2_input_freeze_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    if verdict == "FAIL":
        raise SystemExit("Tier2 input freeze failed; inspect tier2_input_freeze_report.md")


def add_precheck_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    numeric_cols = [
        "neutral_retention_score_v0",
        "acidic_release_support_score_v0",
        "global_weakening_risk_score_v0",
        "mpnn_score_percentile_within_route",
        "rescue_count",
        "His_count",
    ]
    for col in numeric_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    out["t2a_neutral_retention_proxy"] = out["neutral_retention_score_v0"].fillna(out.get("neutral_retention_score", 0.0)).fillna(0.0)
    out["t2a_acidic_release_proxy"] = out["acidic_release_support_score_v0"].fillna(out.get("acidic_release_support_score", 0.0)).fillna(0.0)
    out["t2a_global_weakening_risk_proxy"] = out["global_weakening_risk_score_v0"].fillna(out.get("global_weakening_risk_score", 1.0)).fillna(1.0)
    out["t2a_mpnn_compatibility_proxy"] = (1.0 - out["mpnn_score_percentile_within_route"]).clip(lower=0, upper=1).fillna(0.5)
    out["t2a_has_rescue"] = out["rescue_count"].fillna(0).astype(float) > 0
    out["t2a_is_f_class"] = out["tier1_review_class"].eq("F_neutral_boundary_or_high_risk")
    out["t2a_is_ab_priority"] = out["tier1_review_class"].isin(["A_tier1_priority_candidate", "B_rescue_enriched_candidate"])
    out["t2a_neutral_boundary"] = out["tier1_neutral_retention_flag"].isin(["boundary", "risk"])
    # Tier1 near-duplicate IDs collapse many A/B rows to His seed clusters.
    # Tier2 keeps the His seed cap separately, so use a finer cluster key for
    # A/B rows to avoid applying the same cap twice under different names.
    out["tier2_diversity_cluster_id"] = out["near_duplicate_cluster_id"].astype(str)
    fine_mask = out["t2a_is_ab_priority"] | out["tier1_review_class"].eq("control_anchor")
    out.loc[fine_mask, "tier2_diversity_cluster_id"] = (
        out.loc[fine_mask, "near_duplicate_cluster_id"].astype(str)
        + "|"
        + out.loc[fine_mask, "rescue_signature"].fillna("").astype(str)
        + "|"
        + out.loc[fine_mask, "primary_generation_route"].fillna("").astype(str)
    )
    out["t2a_priority_score"] = (
        0.34 * out["t2a_neutral_retention_proxy"]
        + 0.24 * (1.0 - out["t2a_global_weakening_risk_proxy"]).clip(0, 1)
        + 0.22 * out["t2a_acidic_release_proxy"]
        + 0.12 * out["t2a_mpnn_compatibility_proxy"]
        + 0.08 * out["t2a_has_rescue"].astype(float)
    ).round(6)
    out["feature_source_mutation_region"] = "metadata"
    out["feature_source_his_seed_environment"] = "cheap_proxy"
    out["feature_source_rescue_objective"] = "metadata"
    out["feature_source_neutral_retention"] = "tier1_existing"
    out["feature_source_global_weakening"] = "tier1_existing"
    out["feature_source_mpnn"] = "p0_mpnn_existing"
    out["feature_source_glycan_or_epitope"] = "cheap_proxy"
    out["feature_source_structure_compute"] = "not_run"
    out["t2a_precheck_class"] = "T2A_deprioritized"
    out.loc[out["tier1_review_class"].eq("control_anchor"), "t2a_precheck_class"] = "T2A_control_anchor"
    out.loc[out["tier1_review_class"].eq("E_sparse_mpnn_supplement_review"), "t2a_precheck_class"] = "T2A_sparse_supplement_probe"
    out.loc[out["t2a_is_f_class"], "t2a_precheck_class"] = "T2A_pass_boundary_probe"
    out.loc[out["t2a_is_ab_priority"], "t2a_precheck_class"] = "T2A_pass_high_priority"
    out.loc[
        out["t2a_is_ab_priority"] & out["primary_generation_route"].isin(["ProteinMPNN_seeded_rescue", "structure_or_interface_guided"]),
        "t2a_precheck_class",
    ] = "T2A_pass_route_representative"
    out.loc[out["tier2_mapping_status"].ne("PASS"), "t2a_precheck_class"] = "T2A_fail_obvious_risk"
    return out


def run_precheck(out_dir: Path) -> None:
    snapshot = pd.read_csv(out_dir / "tier2_candidate_snapshot.csv")
    features = add_precheck_features(snapshot)
    write_csv(features, out_dir / "tier2a_precheck_features.csv")
    write_csv(features.groupby(["target", "t2a_precheck_class"], dropna=False).size().reset_index(name="count"), out_dir / "tier2a_precheck_summary_by_class.csv")
    write_csv(features.groupby(["target", "tier2_source_pool", "tier1_review_class"], dropna=False).size().reset_index(name="count"), out_dir / "tier2a_precheck_summary_by_target.csv")
    report = [
        "# Tier2-A Cheap Precheck Audit",
        "",
        "No structural compute was run in this stage.",
        "",
        "## Precheck Classes",
        markdown_table(features.groupby(["target", "t2a_precheck_class"], dropna=False).size().reset_index(name="count")),
        "",
        "## Feature Sources",
        markdown_table(pd.DataFrame({"feature_source_column": [c for c in features.columns if c.startswith("feature_source_")], "values": [features[c].dropna().astype(str).unique().tolist() for c in features.columns if c.startswith("feature_source_")]})),
    ]
    (out_dir / "tier2a_precheck_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def add_with_caps(
    pool: pd.DataFrame,
    selected_ids: set[str],
    selected_rows: list[pd.Series],
    count: int,
    target: str,
    bucket: str,
    counters: dict[str, defaultdict[str, int]],
) -> None:
    info = BACKBONES[target]
    cluster_limit = max(1, int(info["stage1_size"] * info["cluster_cap_fraction"]))
    seed_limit = max(1, int(info["stage1_size"] * info["seed_cap_fraction"]))
    f_limit = info["f_cap"]
    bucket_key = f"{target}:{bucket}"
    for _, row in pool.iterrows():
        if counters["bucket"][bucket_key] >= count:
            return
        vid = str(row["variant_id"])
        if vid in selected_ids:
            continue
        cluster = str(row.get("tier2_diversity_cluster_id", row.get("near_duplicate_cluster_id", "")))
        seed = str(row.get("his_seed_set", ""))
        is_f = row.get("tier1_review_class") == "F_neutral_boundary_or_high_risk"
        if counters["cluster"][cluster] >= cluster_limit:
            continue
        if counters["seed"][seed] >= seed_limit:
            continue
        if is_f and counters["f_class"][target] >= f_limit:
            continue
        row = row.copy()
        row["tier2_stage1_bucket"] = bucket
        selected_rows.append(row)
        selected_ids.add(vid)
        counters["bucket"][bucket_key] += 1
        counters["cluster"][cluster] += 1
        counters["seed"][seed] += 1
        if is_f:
            counters["f_class"][target] += 1


def build_stage1(out_dir: Path) -> None:
    df = pd.read_csv(out_dir / "tier2a_precheck_features.csv")
    df = df[df["tier2_mapping_status"].eq("PASS")].copy()
    df["sort_rank"] = df["t2a_priority_score"].rank(method="first", ascending=False)
    df = df.sort_values(["target", "t2a_priority_score", "tier1_rank_score"], ascending=[True, False, False])
    selected_rows: list[pd.Series] = []
    selected_ids: set[str] = set()
    counters: dict[str, defaultdict[str, int]] = {
        "bucket": defaultdict(int),
        "cluster": defaultdict(int),
        "seed": defaultdict(int),
        "f_class": defaultdict(int),
    }

    quotas = {
        "Ab_1E62": [
            ("control_anchor", 12),
            ("A_B_priority", 1500),
            ("His_plus_rescue", 600),
            ("F_boundary_representatives", 500),
            ("route_seed_cluster_representatives", 300),
        ],
        "Ab_sdAb": [
            ("control_anchor", 8),
            ("A_B_priority", 1600),
            ("neutral_retention_prioritized", 1000),
            ("His_plus_rescue", 600),
            ("F_boundary_representatives", 450),
            ("route_seed_cluster_representatives", 300),
        ],
    }
    for target, target_quotas in quotas.items():
        t = df[df["target"].eq(target)].copy()
        for bucket, n in target_quotas:
            if bucket == "control_anchor":
                pool = t[t["tier1_review_class"].eq("control_anchor")]
            elif bucket == "A_B_priority":
                pool = t[t["tier1_review_class"].isin(["A_tier1_priority_candidate", "B_rescue_enriched_candidate"])]
            elif bucket == "neutral_retention_prioritized":
                pool = t[t["tier1_neutral_retention_flag"].isin(["pass", "control"]) & ~t["tier1_review_class"].eq("F_neutral_boundary_or_high_risk")]
            elif bucket == "His_plus_rescue":
                pool = t[(t["rescue_count"].fillna(0).astype(float) > 0) & t["His_count"].fillna(0).astype(float).ge(1)]
            elif bucket == "F_boundary_representatives":
                pool = t[t["tier1_review_class"].eq("F_neutral_boundary_or_high_risk")]
            else:
                route_pieces = []
                for route, group in t.groupby("primary_generation_route", dropna=False):
                    route_pieces.append(group.head(max(10, n // max(1, t["primary_generation_route"].nunique()))))
                pool = pd.concat(route_pieces, ignore_index=False) if route_pieces else t.iloc[0:0]
            add_with_caps(pool, selected_ids, selected_rows, n, target, bucket, counters)

        target_size = BACKBONES[target]["stage1_size"]
        remainder_pool = t[~t["tier1_review_class"].isin(["F_neutral_boundary_or_high_risk", "G_audit_or_demoted"])].copy()
        add_with_caps(remainder_pool, selected_ids, selected_rows, target_size, target, "balanced_fill", counters)
        if sum(1 for r in selected_rows if r["target"] == target) < target_size:
            add_with_caps(t, selected_ids, selected_rows, target_size, target, "last_resort_fill", counters)

    stage1 = pd.DataFrame(selected_rows).reset_index(drop=True)
    # Trim exactly per target in case a bucket overfilled through fill logic.
    stage1 = (
        stage1.sort_values(["target", "tier2_stage1_bucket", "t2a_priority_score"], ascending=[True, True, False])
        .groupby("target", group_keys=False)
        .head(0)
    ) if stage1.empty else stage1
    parts = []
    bucket_priority = {
        "control_anchor": 0,
        "A_B_priority": 1,
        "neutral_retention_prioritized": 2,
        "His_plus_rescue": 3,
        "F_boundary_representatives": 4,
        "route_seed_cluster_representatives": 5,
        "balanced_fill": 6,
        "last_resort_fill": 7,
    }
    for target, info in BACKBONES.items():
        target_df = pd.DataFrame(selected_rows)
        target_df = target_df[target_df["target"].eq(target)].drop_duplicates("variant_id", keep="first")
        target_df["tier2_stage1_bucket_priority"] = target_df["tier2_stage1_bucket"].map(bucket_priority).fillna(99).astype(int)
        target_df = target_df.sort_values(["tier2_stage1_bucket_priority", "t2a_priority_score", "tier1_rank_score"], ascending=[True, False, False]).head(info["stage1_size"])
        parts.append(target_df)
    stage1 = pd.concat(parts, ignore_index=True, sort=False)

    stage1["tier2_stage1_role"] = "full_stage1"
    stage1["tier2_smoke_selected"] = False
    stage1["tier2_calibration_selected"] = False
    for target in BACKBONES:
        target_idx = stage1[stage1["target"].eq(target)].index
        smoke_idx = []
        for bucket, group in stage1.loc[target_idx].groupby("tier2_stage1_bucket", dropna=False):
            smoke_idx.extend(group.head(max(1, 100 // max(1, stage1.loc[target_idx, "tier2_stage1_bucket"].nunique()))).index.tolist())
        smoke_idx = smoke_idx[: 100 if target == "Ab_1E62" else 100]
        stage1.loc[smoke_idx, "tier2_smoke_selected"] = True
        calib_idx = []
        for bucket, group in stage1.loc[target_idx].groupby("tier2_stage1_bucket", dropna=False):
            calib_idx.extend(group.head(max(1, (500 if target == "Ab_1E62" else 600) // max(1, stage1.loc[target_idx, "tier2_stage1_bucket"].nunique()))).index.tolist())
        calib_idx = calib_idx[: (500 if target == "Ab_1E62" else 600)]
        stage1.loc[calib_idx, "tier2_calibration_selected"] = True

    stage1["tier2_foldx_subset_selected"] = False
    for target, info in BACKBONES.items():
        foldx_idx = []
        target_df = stage1[stage1["target"].eq(target)]
        for bucket, group in target_df.groupby("tier2_stage1_bucket", dropna=False):
            foldx_idx.extend(group.head(max(1, info["foldx_subset_cap"] // max(1, target_df["tier2_stage1_bucket"].nunique()))).index.tolist())
        stage1.loc[foldx_idx[: info["foldx_subset_cap"]], "tier2_foldx_subset_selected"] = True

    smoke = stage1[stage1["tier2_smoke_selected"]].copy()
    calibration = stage1[stage1["tier2_calibration_selected"]].copy()
    write_csv(stage1, out_dir / "tier2b_stage1_candidate_list.csv")
    write_csv(smoke, out_dir / "tier2b_smoke_candidate_list.csv")
    write_csv(calibration, out_dir / "tier2b_calibration_candidate_list.csv")
    summary = stage1.groupby(["target", "tier2_stage1_bucket", "tier1_review_class"], dropna=False).size().reset_index(name="count")
    write_csv(summary, out_dir / "tier2b_stage1_composition_summary.csv")
    manifest = {
        "version": 1,
        "status": "ready_for_smoke_compute",
        "tier2_compute_started": False,
        "full_stage1_compute_started": False,
        "counts": {
            "stage1": stage1["target"].value_counts().to_dict(),
            "smoke": smoke["target"].value_counts().to_dict(),
            "calibration": calibration["target"].value_counts().to_dict(),
            "foldx_subset": stage1[stage1["tier2_foldx_subset_selected"]]["target"].value_counts().to_dict(),
        },
    }
    (out_dir / "tier2b_stage1_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = [
        "# Tier2-B Stage-1 Selection Audit",
        "",
        "## Counts",
        markdown_table(stage1.groupby(["target"], dropna=False).size().reset_index(name="count")),
        "",
        "## Smoke Counts",
        markdown_table(smoke.groupby(["target", "tier2_stage1_bucket"], dropna=False).size().reset_index(name="count")),
        "",
        "## Calibration Counts",
        markdown_table(calibration.groupby(["target", "tier2_stage1_bucket"], dropna=False).size().reset_index(name="count")),
        "",
        "## Full Stage-1 Composition",
        markdown_table(summary),
    ]
    (out_dir / "tier2b_stage1_selection_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["freeze", "precheck", "stage1", "all"], default="all")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    if args.stage in {"freeze", "all"}:
        freeze_inputs(out_dir)
    if args.stage in {"precheck", "all"}:
        run_precheck(out_dir)
    if args.stage in {"stage1", "all"}:
        build_stage1(out_dir)


if __name__ == "__main__":
    main()
