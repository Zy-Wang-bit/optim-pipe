#!/usr/bin/env python3
"""Collect P0 ProteinMPNN score-only and constrained-generation outputs.

The collector is safe to run before ProteinMPNN jobs finish. Missing outputs are
kept as explicit pending/missing statuses; no score is fabricated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from analysis.initial_design_generation import run_dry_run as dry
from analysis.initial_design_generation import run_p0_mpnn_integration as p0


ROOT = dry.ROOT
P0 = ROOT / "results/initial_design_generation/p0_mpnn"
RUNNER = ROOT / "results/initial_design_generation/p0_mpnn_runner_inputs"
POOL = ROOT / "results/initial_design_generation/production_initial_pool/production_initial_pool_candidates_all.csv"
REFERENCE = ROOT / "results/ph_sensitive_40aa_window/tables/reference_sequence_map.csv"
CONFIG = dry.CONFIG
TARGET_REVERSE = {"Ab_1E62": "1E62", "Ab_sdAb": "sdAb"}
TARGET_SAFE = {"1E62": "Ab_1E62", "sdAb": "Ab_sdAb"}


def normalize_target(value: str) -> str:
    return TARGET_REVERSE.get(str(value), str(value))


def safe_target(value: str) -> str:
    value = normalize_target(value)
    return TARGET_SAFE.get(value, value)


def split_positions(value: object) -> list[int]:
    if value is None or pd.isna(value):
        return []
    out: list[int] = []
    for item in str(value).split(";"):
        item = item.strip()
        if item and item.lower() != "nan":
            out.append(int(float(item)))
    return out


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


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


def load_npz_score(path: Path) -> tuple[float, float]:
    data = np.load(path, allow_pickle=True)
    score = float(np.asarray(data["score"]).mean())
    global_score = float(np.asarray(data["global_score"]).mean())
    return score, global_score


def reference_sequences() -> dict[tuple[str, str], str]:
    ref = pd.read_csv(REFERENCE)
    seqs: dict[tuple[str, str], str] = {}
    for (target, chain), sub in ref.groupby([ref["target"].map(normalize_target), "chain"]):
        seqs[(target, chain)] = "".join(sub.sort_values("local_pos")["aa"].astype(str).tolist())
    return seqs


def mutation_list(parent: str, seq: str, chain: str) -> list[str]:
    muts = []
    for i, (a, b) in enumerate(zip(parent, seq), start=1):
        if a != b:
            muts.append(f"{chain}{a}{i}{b}")
    return muts


def nxs_t_sites(seq: str) -> set[int]:
    sites = set()
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in {"S", "T"}:
            sites.add(i + 1)
    return sites


def parse_mpnn_header(header: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in header.split(","):
        if "=" in part:
            k, v = part.strip().split("=", 1)
            out[k.strip()] = v.strip()
    return out


def summarize_pending(score: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    summary = score.groupby(group_cols, dropna=False).size().reset_index(name="candidate_count")
    summary["scored_count"] = score.groupby(group_cols, dropna=False)["mpnn_score_status"].apply(
        lambda s: int((s == "scored_by_mpnn").sum())
    ).reset_index(drop=True)
    summary["pending_count"] = summary["candidate_count"] - summary["scored_count"]
    scored = score[score["mpnn_score_status"] == "scored_by_mpnn"]
    if scored.empty:
        summary["mpnn_score_mean"] = pd.NA
        summary["mpnn_parent_delta_mean"] = pd.NA
        return summary
    means = scored.groupby(group_cols, dropna=False).agg(
        mpnn_score_mean=("mpnn_total_score_per_residue", "mean"),
        mpnn_parent_delta_mean=("mpnn_parent_delta", "mean"),
    ).reset_index()
    return summary.merge(means, on=group_cols, how="left")


def collect_score_only(runner_dir: Path, p0_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest_path = runner_dir / "mpnn_score_only_input_manifest.csv"
    pool = pd.read_csv(POOL)
    pool_by_variant = pool.set_index("variant_id", drop=False)
    manifest = pd.read_csv(manifest_path)
    rows: list[dict] = []
    shard_rows: list[dict] = []
    for shard in manifest.itertuples(index=False):
        fasta = Path(shard.fasta_path)
        records = parse_fasta(fasta)
        score_dir = Path(shard.score_output_dir) / "score_only"
        pdb_stem = Path(shard.structure_path).stem
        parent_file = score_dir / f"{pdb_stem}_pdb.npz"
        parent_score = pd.NA
        if parent_file.exists():
            parent_score, _ = load_npz_score(parent_file)
        shard_scored = 0
        shard_missing = 0
        for i, (variant_id, seq) in enumerate(records, start=1):
            if variant_id not in pool_by_variant.index:
                base = {
                    "target": shard.target,
                    "variant_id": variant_id,
                    "sequence_hash": "",
                    "mutation_list": "",
                    "primary_generation_route": "",
                    "his_seed_set": "",
                    "rescue_signature": "",
                }
            else:
                p = pool_by_variant.loc[variant_id]
                base = {
                    "target": p["target"],
                    "variant_id": variant_id,
                    "sequence_hash": p.get("sequence_hash", ""),
                    "mutation_list": p.get("mutation_list", ""),
                    "primary_generation_route": p.get("primary_generation_route", ""),
                    "his_seed_set": p.get("his_seed_set", ""),
                    "rescue_signature": p.get("rescue_signature", ""),
                }
            npz = score_dir / f"{pdb_stem}_fasta_{i}.npz"
            row = {
                **base,
                "mpnn_generation_status": "not_generated_by_mpnn",
                "mpnn_score_status": "pending_mpnn_run",
                "mpnn_failure_reason": "score_output_missing",
                "mpnn_model_version": "",
                "mpnn_input_backbone_id": shard.backbone_id,
                "mpnn_fixed_positions": "",
                "mpnn_designed_positions": shard.mpnn_design_chain_id,
                "mpnn_temperature": "",
                "mpnn_sample_id": f"{Path(fasta).stem}:{i}",
                "mpnn_random_seed": "",
                "mpnn_total_score": pd.NA,
                "mpnn_total_score_per_residue": pd.NA,
                "mpnn_window_score": pd.NA,
                "mpnn_mutated_position_score": pd.NA,
                "mpnn_his_seed_score": pd.NA,
                "mpnn_rescue_score": pd.NA,
                "mpnn_parent_delta": pd.NA,
                "mpnn_route_percentile": pd.NA,
                "mpnn_his_seed_percentile": pd.NA,
                "mpnn_rescue_type_percentile": pd.NA,
                "mpnn_score_confidence": "not_scored",
            }
            if npz.exists():
                score, global_score = load_npz_score(npz)
                row.update(
                    {
                        "mpnn_score_status": "scored_by_mpnn",
                        "mpnn_failure_reason": "",
                        "mpnn_total_score_per_residue": score,
                        "mpnn_total_score": score * len(seq),
                        "mpnn_window_score": pd.NA,
                        "mpnn_parent_delta": score - parent_score if not pd.isna(parent_score) else pd.NA,
                        "mpnn_score_confidence": "score_only_full_design_chain",
                    }
                )
                shard_scored += 1
            else:
                shard_missing += 1
            rows.append(row)
        shard_rows.append(
            {
                "target": shard.target,
                "backbone_id": shard.backbone_id,
                "shard_index": shard.shard_index,
                "candidate_count": len(records),
                "scored_count": shard_scored,
                "missing_count": shard_missing,
                "score_output_dir": shard.score_output_dir,
                "score_status": "complete" if shard_scored == len(records) else "pending_or_partial",
            }
        )
    score = pd.DataFrame(rows)
    for col in ["primary_generation_route", "his_seed_set", "rescue_signature"]:
        if col not in score:
            score[col] = ""
    for group_col, out_col in [
        ("primary_generation_route", "mpnn_route_percentile"),
        ("his_seed_set", "mpnn_his_seed_percentile"),
        ("rescue_signature", "mpnn_rescue_type_percentile"),
    ]:
        scored = score["mpnn_score_status"] == "scored_by_mpnn"
        if scored.any():
            score.loc[scored, out_col] = score[scored].groupby(["target", group_col])[
                "mpnn_total_score_per_residue"
            ].rank(pct=True, ascending=True)
    score = score[p0.MPNN_SCORE_COLUMNS]
    dry.write_csv(score, p0_dir / "mpnn_scores_current_pool.csv")
    dry.write_csv(summarize_pending(score, ["target", "primary_generation_route"]), p0_dir / "mpnn_score_distribution_by_route.csv")
    dry.write_csv(summarize_pending(score, ["target", "his_seed_set"]), p0_dir / "mpnn_score_distribution_by_his_seed.csv")
    dry.write_csv(summarize_pending(score, ["target", "primary_generation_route", "rescue_signature"]), p0_dir / "mpnn_score_distribution_by_rescue_type.csv")
    dry.write_csv(summarize_pending(score, ["target", "primary_generation_route"]), p0_dir / "mpnn_parent_delta_by_route.csv")
    audit = pd.DataFrame(shard_rows)
    dry.write_csv(audit, p0_dir / "mpnn_score_only_collection_audit.csv")
    return score, audit


def collect_generation(runner_dir: Path, p0_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    config = yaml.safe_load(CONFIG.read_text())
    jobs = pd.read_csv(runner_dir / "mpnn_runner_jobs_constrained.csv")
    refs = reference_sequences()
    rows: list[dict] = []
    audit_rows: list[dict] = []
    rescue_rows: list[dict] = []
    seq_root = runner_dir / "mpnn_out"
    for job in jobs.itertuples(index=False):
        target = normalize_target(job.target)
        safe = safe_target(target)
        chain = job.mpnn_design_chain_id
        parent = refs[(target, chain)]
        hard_positions = set(config["targets"][target]["hard_protect_positions"])
        parent_nxs = nxs_t_sites(parent)
        seed_positions = split_positions(job.seed_positions_forced_to_H)
        fasta_candidates = list(seq_root.glob(f"**/{job.job_name}.fa")) + list(seq_root.glob(f"**/{job.job_name}*.fasta"))
        if not fasta_candidates:
            audit_rows.append(
                {
                    "target": safe,
                    "job_name": job.job_name,
                    "generated_count": 0,
                    "passed_hard_filter_count": 0,
                    "status": "pending_mpnn_generation_run",
                    "failure_reason": "generation_fasta_missing",
                }
            )
            continue
        generated = 0
        passed = 0
        for fasta in fasta_candidates:
            for header, seq in parse_fasta(fasta):
                if not header.startswith("T="):
                    continue
                generated += 1
                meta = parse_mpnn_header(header)
                seq_clean = seq.replace("/", "")
                muts = mutation_list(parent, seq_clean, chain)
                hard_fail_reasons = []
                if len(seq_clean) != len(parent):
                    hard_fail_reasons.append("sequence_length_mismatch")
                if set(seq_clean) - set("ACDEFGHIKLMNPQRSTVWY"):
                    hard_fail_reasons.append("nonstandard_amino_acid")
                allowed_positions = set(seed_positions) | set(split_positions(job.rescue_design_positions))
                mutated_positions = {
                    int(re.match(rf"{re.escape(chain)}[A-Z](\d+)[A-Z]", m).group(1))
                    for m in muts
                    if re.match(rf"{re.escape(chain)}[A-Z](\d+)[A-Z]", m)
                }
                unexpected_positions = mutated_positions - allowed_positions
                if unexpected_positions:
                    hard_fail_reasons.append(
                        "unexpected_mutation_positions:" + ";".join(str(x) for x in sorted(unexpected_positions))
                    )
                for pos in seed_positions:
                    if len(seq_clean) < pos or seq_clean[pos - 1] != "H":
                        hard_fail_reasons.append(f"seed_{pos}_not_H")
                for pos in hard_positions:
                    if len(seq_clean) >= pos and seq_clean[pos - 1] != parent[pos - 1]:
                        hard_fail_reasons.append(f"hard_protect_{pos}_mutated")
                rescue_positions = set(split_positions(job.rescue_design_positions))
                for pos in rescue_positions:
                    if len(seq_clean) >= pos and seq_clean[pos - 1] in {"C", "H", "X"}:
                        hard_fail_reasons.append(f"rescue_{pos}_disallowed_{seq_clean[pos - 1]}")
                if any(m.endswith("C") for m in muts):
                    hard_fail_reasons.append("new_cys")
                if any(m.endswith("H") and int(re.match(rf"{re.escape(chain)}[A-Z](\d+)H", m).group(1)) not in seed_positions for m in muts if re.match(rf"{re.escape(chain)}[A-Z](\d+)H", m)):
                    hard_fail_reasons.append("non_seed_his")
                if target == "sdAb" and len(seq_clean) >= 110 and seq_clean[104] == "H" and seq_clean[109] == "H":
                    hard_fail_reasons.append("sdAb_V105H_D110H")
                if nxs_t_sites(seq_clean) - parent_nxs:
                    hard_fail_reasons.append("new_canonical_nxs_t")
                his_mut_count = sum(1 for m in muts if m.endswith("H"))
                if his_mut_count > int(config["targets"][target]["main_max_his"]):
                    hard_fail_reasons.append("his_count_exceeds_main_max")
                non_his_rescue_count = sum(
                    1
                    for m in muts
                    if (
                        re.match(rf"{re.escape(chain)}[A-Z](\d+)[A-Z]", m)
                        and int(re.match(rf"{re.escape(chain)}[A-Z](\d+)[A-Z]", m).group(1)) in rescue_positions
                        and not m.endswith("H")
                    )
                )
                if non_his_rescue_count > int(config["targets"][target]["max_non_his_rescue_count"]):
                    hard_fail_reasons.append("non_his_rescue_count_exceeds_max")
                actual_temp = str(meta.get("T", ""))
                planned_temp = str(job.mpnn_temperature)
                if actual_temp and actual_temp != planned_temp:
                    hard_fail_reasons.append(f"temperature_mismatch_planned_{planned_temp}_actual_{actual_temp}")
                status = "pass" if not hard_fail_reasons else "fail"
                if status == "pass":
                    passed += 1
                variant_hash = hashlib.sha256(f"{job.job_name}:{header}:{seq}".encode()).hexdigest()[:12]
                variant_id = f"mpnn_{target}_{variant_hash}"
                rows.append(
                    {
                        "target": safe,
                        "variant_id": variant_id,
                        "generation_record_id": job.job_name,
                        "sequence": seq_clean,
                        "sequence_hash": hashlib.sha256(seq_clean.encode()).hexdigest()[:12],
                        "mutation_list": ";".join(muts),
                        "mutation_count": len(muts),
                        "His_count": his_mut_count,
                        "his_seed_set": job.his_seed_set,
                        "primary_generation_route": "constrained_MPNN_rescue",
                        "rescue_mutation_list": "",
                        "rescue_signature": "",
                        "rescue_count": max(0, len(muts) - len(seed_positions)),
                        "mpnn_generation_status": "generated_by_constrained_mpnn",
                        "mpnn_model_version": "",
                        "mpnn_input_backbone_id": job.backbone_id,
                        "mpnn_fixed_positions": "",
                        "mpnn_designed_positions": job.rescue_design_positions,
                        "mpnn_temperature": meta.get("T", job.mpnn_temperature),
                        "mpnn_sample_id": meta.get("sample", ""),
                        "mpnn_random_seed": "",
                        "hard_filter_status": status,
                        "liability_flags": ";".join(hard_fail_reasons),
                        "selection_reason": "p0_constrained_mpnn_generated_candidate",
                    }
                )
        audit_rows.append(
            {
                "target": safe,
                "job_name": job.job_name,
                "generated_count": generated,
                "passed_hard_filter_count": passed,
                "status": "complete_or_partial" if generated else "pending_mpnn_generation_run",
                "failure_reason": "" if generated else "generation_fasta_missing",
            }
        )
    pool = pd.DataFrame(rows, columns=p0.POOL_COLUMNS)
    audit = pd.DataFrame(audit_rows)
    for target in ["1E62", "sdAb"]:
        safe = safe_target(target)
        sub = pool[pool["target"] == safe]
        if sub.empty:
            sub = p0.empty_pool(target, "pending_constrained_mpnn_run")
        dry.write_csv(sub, p0_dir / f"constrained_mpnn_rescue_pool_{target}.csv")
    dry.write_csv(audit, p0_dir / "constrained_mpnn_generation_collection_audit.csv")
    dry.write_csv(audit.groupby("target", dropna=False).agg(
        generated_count=("generated_count", "sum"),
        passed_hard_filter_count=("passed_hard_filter_count", "sum"),
        pending_job_count=("status", lambda s: int((s == "pending_mpnn_generation_run").sum())),
    ).reset_index(), p0_dir / "constrained_mpnn_rescue_audit.csv")
    return pool, audit


def collect_relaxed_generation(runner_dir: Path, p0_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    config = yaml.safe_load(CONFIG.read_text())
    jobs_path = runner_dir / "mpnn_runner_jobs_relaxed.csv"
    if not jobs_path.exists():
        audit = pd.DataFrame(columns=["target", "job_name", "generated_count", "passed_hard_filter_count", "status", "failure_reason"])
        for target in ["1E62", "sdAb"]:
            dry.write_csv(p0.empty_pool(target, "pending_relaxed_mpnn_run"), p0_dir / f"relaxed_mpnn_counterfactual_pool_{target}.csv")
        dry.write_csv(audit, p0_dir / "relaxed_mpnn_generation_collection_audit.csv")
        return pd.DataFrame(columns=p0.POOL_COLUMNS), audit
    jobs = pd.read_csv(jobs_path)
    refs = reference_sequences()
    rows: list[dict] = []
    audit_rows: list[dict] = []
    seq_root = runner_dir / "relaxed_mpnn_out"
    for job in jobs.itertuples(index=False):
        target = normalize_target(job.target)
        safe = safe_target(target)
        chain = job.mpnn_design_chain_id
        parent = refs[(target, chain)]
        hard_positions = set(config["targets"][target]["hard_protect_positions"])
        parent_nxs = nxs_t_sites(parent)
        seed_positions = split_positions(job.seed_positions_forced_to_H)
        design_positions = set(split_positions(job.relaxed_design_positions))
        fasta_candidates = list(seq_root.glob(f"**/{job.job_name}.fa")) + list(seq_root.glob(f"**/{job.job_name}*.fasta"))
        if not fasta_candidates:
            audit_rows.append(
                {
                    "target": safe,
                    "job_name": job.job_name,
                    "branch": job.branch,
                    "generated_count": 0,
                    "passed_hard_filter_count": 0,
                    "status": "pending_relaxed_mpnn_generation_run",
                    "failure_reason": "generation_fasta_missing",
                }
            )
            continue
        generated = 0
        passed = 0
        for fasta in fasta_candidates:
            for header, seq in parse_fasta(fasta):
                if not header.startswith("T="):
                    continue
                generated += 1
                meta = parse_mpnn_header(header)
                seq_clean = seq.replace("/", "")
                muts = mutation_list(parent, seq_clean, chain)
                hard_fail_reasons = []
                if len(seq_clean) != len(parent):
                    hard_fail_reasons.append("sequence_length_mismatch")
                if set(seq_clean) - set("ACDEFGHIKLMNPQRSTVWY"):
                    hard_fail_reasons.append("nonstandard_amino_acid")
                mutated_positions = {
                    int(re.match(rf"{re.escape(chain)}[A-Z](\d+)[A-Z]", m).group(1))
                    for m in muts
                    if re.match(rf"{re.escape(chain)}[A-Z](\d+)[A-Z]", m)
                }
                unexpected_positions = mutated_positions - design_positions
                if unexpected_positions:
                    hard_fail_reasons.append(
                        "unexpected_mutation_positions:" + ";".join(str(x) for x in sorted(unexpected_positions))
                    )
                for pos in seed_positions:
                    if len(seq_clean) < pos or seq_clean[pos - 1] != "H":
                        hard_fail_reasons.append(f"seed_{pos}_not_H")
                for pos in hard_positions:
                    if len(seq_clean) >= pos and seq_clean[pos - 1] != parent[pos - 1]:
                        hard_fail_reasons.append(f"hard_protect_{pos}_mutated")
                if any(m.endswith("C") for m in muts):
                    hard_fail_reasons.append("new_cys")
                if target == "sdAb" and len(seq_clean) >= 110 and seq_clean[104] == "H" and seq_clean[109] == "H":
                    hard_fail_reasons.append("sdAb_V105H_D110H")
                if nxs_t_sites(seq_clean) - parent_nxs:
                    hard_fail_reasons.append("new_canonical_nxs_t")
                his_mut_count = sum(1 for m in muts if m.endswith("H"))
                if his_mut_count > int(config["targets"][target]["main_max_his"]):
                    hard_fail_reasons.append("his_count_exceeds_main_max")
                actual_temp = str(meta.get("T", ""))
                planned_temp = str(job.mpnn_temperature)
                if actual_temp and actual_temp != planned_temp:
                    hard_fail_reasons.append(f"temperature_mismatch_planned_{planned_temp}_actual_{actual_temp}")
                status = "pass" if not hard_fail_reasons else "fail"
                if status == "pass":
                    passed += 1
                variant_hash = hashlib.sha256(f"{job.job_name}:{header}:{seq}".encode()).hexdigest()[:12]
                variant_id = f"relaxed_{target}_{variant_hash}"
                rows.append(
                    {
                        "target": safe,
                        "variant_id": variant_id,
                        "generation_record_id": job.job_name,
                        "sequence": seq_clean,
                        "sequence_hash": hashlib.sha256(seq_clean.encode()).hexdigest()[:12],
                        "mutation_list": ";".join(muts),
                        "mutation_count": len(muts),
                        "His_count": his_mut_count,
                        "his_seed_set": job.his_seed_set,
                        "primary_generation_route": str(job.branch),
                        "rescue_mutation_list": "",
                        "rescue_signature": "",
                        "rescue_count": 0,
                        "mpnn_generation_status": "generated_by_relaxed_mpnn",
                        "mpnn_model_version": "",
                        "mpnn_input_backbone_id": job.backbone_id,
                        "mpnn_fixed_positions": "",
                        "mpnn_designed_positions": job.relaxed_design_positions,
                        "mpnn_temperature": meta.get("T", job.mpnn_temperature),
                        "mpnn_sample_id": meta.get("sample", ""),
                        "mpnn_random_seed": "",
                        "hard_filter_status": status,
                        "liability_flags": ";".join(hard_fail_reasons),
                        "selection_reason": "p0_relaxed_counterfactual_audit_only_no_preallocated_quota",
                    }
                )
        audit_rows.append(
            {
                "target": safe,
                "job_name": job.job_name,
                "branch": job.branch,
                "generated_count": generated,
                "passed_hard_filter_count": passed,
                "status": "complete_or_partial" if generated else "pending_relaxed_mpnn_generation_run",
                "failure_reason": "" if generated else "generation_fasta_missing",
            }
        )
    pool = pd.DataFrame(rows, columns=p0.POOL_COLUMNS)
    audit = pd.DataFrame(audit_rows)
    for target in ["1E62", "sdAb"]:
        safe = safe_target(target)
        sub = pool[pool["target"] == safe]
        if sub.empty:
            sub = p0.empty_pool(target, "pending_relaxed_mpnn_run")
        dry.write_csv(sub, p0_dir / f"relaxed_mpnn_counterfactual_pool_{target}.csv")
    dry.write_csv(audit, p0_dir / "relaxed_mpnn_generation_collection_audit.csv")
    return pool, audit


def write_report(score: pd.DataFrame, score_audit: pd.DataFrame, gen_pool: pd.DataFrame, gen_audit: pd.DataFrame, relaxed_pool: pd.DataFrame, relaxed_audit: pd.DataFrame, p0_dir: Path) -> None:
    scored_count = int((score["mpnn_score_status"] == "scored_by_mpnn").sum())
    generated_count = int(gen_audit["generated_count"].sum()) if not gen_audit.empty else 0
    relaxed_generated_count = int(relaxed_audit["generated_count"].sum()) if not relaxed_audit.empty else 0
    lines = [
        "# P0 ProteinMPNN Result Collection Report",
        "",
        f"Score rows: {len(score)}",
        f"Scored rows: {scored_count}",
        f"Score shards: {len(score_audit)}",
        f"Constrained generated candidates: {generated_count}",
        f"Relaxed counterfactual generated candidates: {relaxed_generated_count}",
        "",
        "Missing ProteinMPNN outputs are retained as pending/missing statuses; no scores are fabricated.",
    ]
    (p0_dir / "p0_mpnn_result_collection_report.md").write_text("\n".join(lines) + "\n")


def write_manifest(p0_dir: Path, runner_dir: Path) -> None:
    output_files = [
        p0_dir / "mpnn_scores_current_pool.csv",
        p0_dir / "mpnn_score_only_collection_audit.csv",
        p0_dir / "constrained_mpnn_generation_collection_audit.csv",
        p0_dir / "relaxed_mpnn_generation_collection_audit.csv",
        p0_dir / "p0_mpnn_result_collection_report.md",
    ]
    manifest = {
        "run_mode": "p0_mpnn_result_collection",
        "proteinmpnn_outputs_required": True,
        "runner_input_dir": str(runner_dir.relative_to(ROOT)),
        "outputs": {
            str(p.relative_to(ROOT)): sha_file(p)
            for p in output_files
            if p.exists()
        },
    }
    with (p0_dir / "p0_mpnn_result_collection_manifest.yaml").open("w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, allow_unicode=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect P0 ProteinMPNN results.")
    parser.add_argument("--runner-dir", default=str(RUNNER.relative_to(ROOT)))
    parser.add_argument("--p0-dir", default=str(P0.relative_to(ROOT)))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    runner_dir = ROOT / args.runner_dir
    p0_dir = ROOT / args.p0_dir
    score, score_audit = collect_score_only(runner_dir, p0_dir)
    gen_pool, gen_audit = collect_generation(runner_dir, p0_dir)
    relaxed_pool, relaxed_audit = collect_relaxed_generation(runner_dir, p0_dir)
    write_report(score, score_audit, gen_pool, gen_audit, relaxed_pool, relaxed_audit, p0_dir)
    write_manifest(p0_dir, runner_dir)
    print(
        {
            "status": "p0_mpnn_results_collected",
            "score_rows": len(score),
            "scored_rows": int((score["mpnn_score_status"] == "scored_by_mpnn").sum()),
            "generation_jobs": len(gen_audit),
            "generated_candidates": int(gen_audit["generated_count"].sum()) if not gen_audit.empty else 0,
            "relaxed_generation_jobs": len(relaxed_audit),
            "relaxed_generated_candidates": int(relaxed_audit["generated_count"].sum()) if not relaxed_audit.empty else 0,
        }
    )


if __name__ == "__main__":
    main()
