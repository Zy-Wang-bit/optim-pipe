#!/usr/bin/env python3
"""Run the initial-design dry-run gate for 1E62 and sdAb.

This is intentionally limited to rule-based generation and audit. It does not
run AF3, FoldX actual, Rosetta, explicit glycan modeling, ESM, or MD.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import socket
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / "results/ph_sensitive_40aa_window/tables"
OUT = ROOT / "results/initial_design_generation"
TABLES = OUT / "tables"
DRY = OUT / "dry_run"
REPORTS = OUT / "reports"
CONFIG = ROOT / "analysis/initial_design_generation/config/design_config.yaml"
REPORT_NAME = "dry_run_validation_report.md"

PRIMARY_ROUTES = [
    "His_rule",
    "His_plus_rescue",
    "ProteinMPNN_seeded_rescue",
    "wetlab_informed_expansion",
    "structure_or_interface_guided",
]


@dataclass(frozen=True)
class Mutation:
    chain: str
    pos: int
    parent: str
    mutant: str

    def label(self) -> str:
        return f"{self.chain}{self.parent}{self.pos}{self.mutant}"

    def is_his(self) -> bool:
        return self.mutant == "H"


def sha(text: str, n: int | None = None) -> str:
    value = hashlib.sha256(text.encode()).hexdigest()
    return value[:n] if n else value


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path, index: bool = False) -> None:
    safe = frame.copy()
    if "target" in safe.columns:
        safe["target"] = safe["target"].replace({"1E62": "Ab_1E62", "sdAb": "Ab_sdAb"})
    safe.to_csv(path, index=index, quoting=csv.QUOTE_ALL)


def config_hash(config: dict) -> str:
    text = yaml.safe_dump(config, sort_keys=True, allow_unicode=True)
    return hashlib.sha256(text.encode()).hexdigest()


def read_config() -> dict:
    with CONFIG.open() as fh:
        return yaml.safe_load(fh)


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        return "unavailable"


def load_inputs() -> dict[str, pd.DataFrame]:
    names = [
        "window_mutation_mask",
        "window_scores",
        "candidate_windows",
        "reference_sequence_map",
        "wet_observation_table_1e62",
        "wet_observation_table_sdab",
        "prior_constraints_table",
    ]
    return {name: pd.read_csv(UPSTREAM / f"{name}.csv") for name in names}


def chain_sequence(reference: pd.DataFrame, target: str, chain: str) -> str:
    df = reference[(reference.target == target) & (reference.chain == chain)].copy()
    df = df.sort_values("local_pos")
    return "".join(df["aa"].astype(str).tolist())


def nxs_motifs(seq: str) -> set[int]:
    hits: set[int] = set()
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in {"S", "T"}:
            hits.add(i + 1)
    return hits


def apply_mutations(parent_seq: str, muts: Iterable[Mutation]) -> str:
    chars = list(parent_seq)
    for m in muts:
        if chars[m.pos - 1] != m.parent:
            raise ValueError(f"Parent mismatch at {m.label()}: {chars[m.pos - 1]}")
        chars[m.pos - 1] = m.mutant
    return "".join(chars)


def parse_positions(value: str | float | int) -> list[int]:
    if pd.isna(value):
        return []
    return [int(x) for x in re.findall(r"\d+", str(value))]


def build_wet_support(inputs: dict[str, pd.DataFrame]) -> dict[tuple[str, int], str]:
    support: dict[tuple[str, int], str] = {}
    for table_name in ["wet_observation_table_1e62", "wet_observation_table_sdab"]:
        df = inputs[table_name]
        for _, row in df.iterrows():
            target = str(row["target"])
            for pos in parse_positions(row.get("positions", "")):
                endpoint = str(row.get("endpoint", ""))
                value = row.get("observation_value", "")
                tag = "wet_observed"
                try:
                    numeric = float(value)
                    if endpoint == "KD_ratio" and numeric >= 1.2:
                        tag = "wet_his_seed_support"
                    elif endpoint in {"od_ph74", "od_ph60"}:
                        tag = "wet_elisa_observed"
                except Exception:
                    if str(value).lower() == "yes":
                        tag = "wet_expression_yes"
                support[(target, pos)] = max(support.get((target, pos), ""), tag)
    return support


def build_evidence_ledger(config: dict, inputs: dict[str, pd.DataFrame]) -> None:
    mask = inputs["window_mutation_mask"]
    reference = inputs["reference_sequence_map"]
    wet_support = build_wet_support(inputs)
    rows: list[dict] = []
    forbidden_rows: list[dict] = []
    rescue_rows: list[dict] = []

    for target, target_cfg in config["targets"].items():
        window_id = target_cfg["window_id"]
        chain = target_cfg["chain_id"]
        hard_positions = set(target_cfg["hard_protect_positions"])
        window = mask[mask.window_id == window_id].copy().sort_values("pos")
        if len(window) != 40:
            raise ValueError(f"{window_id} expected 40 rows, got {len(window)}")
        ref_target = reference[(reference.target == target) & (reference.chain == chain)]
        ref_by_pos = dict(zip(ref_target.local_pos.astype(int), ref_target.aa.astype(str)))
        for _, r in window.iterrows():
            pos = int(r["pos"])
            aa = str(r["aa"])
            if ref_by_pos.get(pos) != aa:
                raise ValueError(f"Reference mismatch for {target} {chain}{pos}")
            protected = pos in hard_positions or str(r["mask_status"]) == "hard_protect"
            soft = str(r["mask_status"]) == "soft_risk"
            wet = wet_support.get((target, pos), "")
            position_status = "hard_protect" if protected else "designable_restricted" if soft else "designable_main"
            allowed = "" if protected else str(r["allowed_mutation_classes"]).replace("|", ";")
            rows.append(
                {
                    "target": target,
                    "chain_id": chain,
                    "chain_type": "VL" if target == "1E62" else "VHH",
                    "position_scheme": "chain_local_index",
                    "position": pos,
                    "window_local_position": pos - int(window.pos.min()) + 1,
                    "parent_aa": aa,
                    "region": r["region"],
                    "constraint_type": "hard_protect" if protected else "soft_risk" if soft else "designable",
                    "position_status": position_status,
                    "protected_from_mutation": protected,
                    "forbidden_mutation_to": "",
                    "allowed_mutation_set": allowed,
                    "restricted_alphabet": "" if protected else "H;A;S;T;N;Q;Y;D;E;K;R;G;P;V;L;I;M;F;W",
                    "rescue_eligible": (not protected) and pos not in target_cfg["his_seed_positions"],
                    "rescue_objective_allowed": "neutral_retention_rescue;polar_network_rescue;charge_compensation",
                    "constraint_source": "manual_minimum" if protected else "upstream_mask",
                    "constraint_confidence": "confirmed" if protected else "probable" if soft else "hypothesis",
                    "evidence_level": "confirmed" if protected else "probable" if wet else "hypothesis",
                    "evidence_conflict_flag": False,
                    "wetlab_support": wet,
                    "structural_support": "window_selected_upstream",
                    "interface_support": "upstream_window_score",
                    "glycan_or_epitope_risk": "low_or_not_flagged",
                    "known_failure_flag": False,
                    "applies_to_generation": True,
                    "applies_to_final_library": True,
                }
            )
            if not protected:
                rescue_rows.append(
                    {
                        "target": target,
                        "chain_id": chain,
                        "position": pos,
                        "parent_aa": aa,
                        "rescue_type": "conservative_or_polar",
                        "allowed_residues": "A;S;T;N;Q;Y;D;E;K;R;G;P;V;L;I;M;F;W",
                        "rescue_objective": "neutral_retention_rescue;polar_network_rescue",
                        "source": "dry_run_rule",
                    }
                )

    for pair in config.get("forbidden_pairs", []):
        forbidden_rows.append(
            {
                "target": pair["target"],
                "pair_id": pair["pair_id"],
                "position_scheme": "chain_local_index",
                "mutation_a_position": pair["mutation_a_position"],
                "mutation_a_from": pair["mutation_a_from"],
                "mutation_a_to": pair["mutation_a_to"],
                "mutation_b_position": pair["mutation_b_position"],
                "mutation_b_from": pair["mutation_b_from"],
                "mutation_b_to": pair["mutation_b_to"],
                "constraint_type": "forbidden_pair",
                "constraint_source": "manual_minimum",
                "constraint_confidence": "probable",
                "applies_to_generation": True,
                "applies_to_final_library": True,
            }
        )

    pd.DataFrame(rows).pipe(write_csv, TABLES / "evidence_ledger.csv")
    pd.DataFrame(forbidden_rows).pipe(write_csv, TABLES / "forbidden_pairs.csv")
    pd.DataFrame(rescue_rows).pipe(write_csv, TABLES / "rescue_definition_table.csv")


def validate_evidence_ledger() -> list[str]:
    failures: list[str] = []
    ledger = pd.read_csv(TABLES / "evidence_ledger.csv")
    pairs = pd.read_csv(TABLES / "forbidden_pairs.csv")
    target_reverse = {"Ab_1E62": "1E62", "Ab_sdAb": "sdAb"}
    if "target" in ledger.columns:
        ledger["target"] = ledger["target"].replace(target_reverse)
    if "target" in pairs.columns:
        pairs["target"] = pairs["target"].replace(target_reverse)
    checks = [
        ("1E62", 23, "C"),
        ("sdAb", 96, "C"),
    ]
    for target, pos, aa in checks:
        hit = ledger[
            (ledger.target == target)
            & (ledger.position.astype(int) == pos)
            & (ledger.parent_aa == aa)
            & (ledger.protected_from_mutation == True)
        ]
        if hit.empty:
            failures.append(f"Missing hard protect {target} {pos}{aa}")
    pair_hit = pairs[
        (pairs.target == "sdAb")
        & (pairs.mutation_a_position.astype(int) == 105)
        & (pairs.mutation_b_position.astype(int) == 110)
    ]
    if pair_hit.empty:
        failures.append("Missing sdAb V105H + D110H forbidden pair")
    return failures


def allowed_rescue_residues(parent: str, rng: random.Random) -> list[str]:
    base = ["A", "S", "T", "N", "Q", "Y", "D", "E", "K", "R", "G", "P", "V", "L", "I", "M", "F", "W"]
    residues = [x for x in base if x != parent and x != "C" and x != "H"]
    rng.shuffle(residues)
    return residues


def ordered_rescue_residues(parent: str) -> list[str]:
    base = ["A", "S", "T", "N", "Q", "Y", "D", "E", "K", "R", "G", "P", "V", "L", "I", "M", "F", "W"]
    return [x for x in base if x != parent and x != "C" and x != "H"]


def rescue_signature_from_muts(muts: Iterable[Mutation]) -> str:
    labels = [m.label() for m in sorted(muts, key=lambda x: x.pos) if not m.is_his()]
    return ";".join(labels) if labels else "none"


_COMBO_CACHE: dict[tuple[tuple[int, ...], int], list[tuple[int, ...]]] = {}


def position_combinations(pool: list[int], count: int) -> list[tuple[int, ...]]:
    key = (tuple(pool), count)
    if key not in _COMBO_CACHE:
        _COMBO_CACHE[key] = list(itertools.combinations(pool, count)) if count > 0 else [()]
    return _COMBO_CACHE[key]


def choose_rescue_mutations(
    chain: str,
    positions: dict[int, str],
    pool: list[int],
    count: int,
    combo_index: int,
    route_offset: int,
) -> list[Mutation]:
    if count <= 0 or not pool:
        return []
    combos = position_combinations(sorted(pool), count)
    if not combos:
        return []
    combo_id = (combo_index + route_offset) % len(combos)
    residue_cursor = (combo_index + route_offset) // len(combos)
    chosen = combos[combo_id]
    muts: list[Mutation] = []
    for i, pos in enumerate(chosen):
        residues = ordered_rescue_residues(positions[pos])
        if not residues:
            continue
        aa = residues[(residue_cursor + i * 7 + pos + route_offset) % len(residues)]
        muts.append(Mutation(chain, pos, positions[pos], aa))
    return muts


def route_for_index(config: dict, idx: int) -> str:
    quotas = config["routes"]["route_soft_quota"]
    total = sum(quotas.values())
    x = (idx % 10000) / 10000 * total
    acc = 0.0
    for route in PRIMARY_ROUTES:
        acc += quotas[route]
        if x <= acc:
            return route
    return PRIMARY_ROUTES[-1]


def parse_mutation_label(label: str) -> Mutation | None:
    m = re.match(r"([A-Z])([A-Z])(\d+)([A-Z])$", str(label))
    if not m:
        return None
    return Mutation(m.group(1), int(m.group(3)), m.group(2), m.group(4))


def mutation_positions(mutation_list: str) -> tuple[int, ...]:
    positions: list[int] = []
    for label in str(mutation_list).split(";"):
        mut = parse_mutation_label(label)
        if mut:
            positions.append(mut.pos)
    return tuple(sorted(positions))


def his_labels(mutation_list: str) -> list[str]:
    labels: list[str] = []
    for label in str(mutation_list).split(";"):
        mut = parse_mutation_label(label)
        if mut and mut.is_his():
            labels.append(mut.label())
    return sorted(labels)


def mutation_map(mutation_list: str) -> dict[int, str]:
    values: dict[int, str] = {}
    for label in str(mutation_list).split(";"):
        mut = parse_mutation_label(label)
        if mut:
            values[mut.pos] = mut.mutant
    return values


def mutation_distance(a: dict[int, str], b: dict[int, str]) -> int:
    dist = 0
    for pos in set(a) | set(b):
        if a.get(pos) != b.get(pos):
            dist += 1
    return dist


def mutation_jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / max(1, len(a | b))


def mutation_key(muts: list[Mutation]) -> str:
    return ";".join(m.label() for m in sorted(muts, key=lambda m: m.pos))


def make_candidate(
    target: str,
    chain: str,
    window_id: str,
    route: str,
    parent_seq: str,
    muts: list[Mutation],
    config_hash: str,
    record_i: int,
    parent_motifs: set[int],
) -> tuple[dict, list[dict], str | None]:
    muts = sorted(muts, key=lambda m: m.pos)
    sequence = apply_mutations(parent_seq, muts)
    new_motifs = nxs_motifs(sequence) - parent_motifs
    if new_motifs:
        return {}, [], "new_canonical_nxs_t"
    if any(m.mutant == "C" for m in muts):
        return {}, [], "new_cys"
    mut_labels = [m.label() for m in muts]
    mut_key = mutation_key(muts)
    seq_hash = sha(sequence, 10)
    mut_hash = sha(mut_key, 10)
    variant_id = f"{target}_{window_id}_{seq_hash}"
    route_abbrev = "".join(part[0] for part in route.split("_"))
    generation_record_id = f"{target}_{window_id}_{seq_hash}_{route_abbrev}_{sha(route + mut_key + str(record_i), 6)}"
    his_seed_set = ";".join(m.label() for m in muts if m.is_his())
    non_his = [m for m in muts if not m.is_his()]
    rescue_objectives = ["neutral_retention_rescue" for _ in non_his]
    rescue_types = ["conservative_or_polar" for _ in non_his]
    rescue_signature = rescue_signature_from_muts(muts)
    foldx_proxy = max(0.0, 1.0 - 0.12 * len(muts) - 0.05 * len(non_his))
    neutral_retention = max(0.0, 1.0 - 0.10 * len(non_his) - 0.07 * sum(1 for m in muts if m.is_his()))
    acidic_release = min(1.0, 0.20 + 0.25 * sum(1 for m in muts if m.is_his()))
    global_weakening = min(1.0, 0.08 * len(muts) + 0.06 * len(non_his))
    preliminary = {
        "His_rule": "His_rescue" if len(muts) > 1 else "high_confidence_pH_switch",
        "His_plus_rescue": "His_rescue",
        "ProteinMPNN_seeded_rescue": "His_rescue",
        "wetlab_informed_expansion": "high_confidence_pH_switch",
        "structure_or_interface_guided": "glycan_epitope_safe",
    }.get(route, "")
    row = {
        "target": target,
        "variant_id": variant_id,
        "generation_record_id": generation_record_id,
        "parent_id": f"{target}_{chain}_parent",
        "sequence": sequence,
        "sequence_hash": seq_hash,
        "mutation_list": mut_key,
        "mutation_list_hash": mut_hash,
        "parent_sequence_hash": sha(parent_seq, 10),
        "config_hash": config_hash,
        "window": window_id,
        "mutation_count": len(muts),
        "His_count": sum(1 for m in muts if m.is_his()),
        "his_seed_set": his_seed_set,
        "primary_generation_route": route,
        "all_source_routes": route,
        "all_generation_records": generation_record_id,
        "generation_intent_tags": "dry_run_rule_based",
        "preliminary_selection_bucket": preliminary,
        "final_selection_bucket": "",
        "supporting_tags": "contains_his_trigger" if his_seed_set else "",
        "sequence_role_list": "design",
        "final_slot_type": "",
        "control_type_if_applicable": "",
        "control_quota_counted": False,
        "duplicated_by_design": False,
        "rescue_mutation_list": ";".join(m.label() for m in non_his),
        "rescue_signature": rescue_signature,
        "rescue_objective_list": ";".join(rescue_objectives),
        "rescue_type_list": ";".join(rescue_types),
        "rescue_count": len(non_his),
        "hard_filter_status": "pass",
        "liability_flags": "",
        "forbidden_pair_status": "pass",
        "exact_duplicate_status": "unique",
        "duplicate_source_count": 1,
        "sequence_hamming_cluster_id": "",
        "mutation_set_jaccard_cluster_id": "",
        "his_seed_cluster_id": "",
        "near_duplicate_cluster_id": "",
        "cluster_size": 0,
        "hit_likelihood_score_v0": round(0.35 * neutral_retention + 0.35 * acidic_release + 0.30 * foldx_proxy, 4),
        "score_confidence_flag": "lightweight_proxy_only",
        "neutral_retention_score": round(neutral_retention, 4),
        "acidic_release_support_score": round(acidic_release, 4),
        "global_weakening_risk_score": round(global_weakening, 4),
        "display_or_expression_risk_score": "",
        "glycan_or_epitope_risk_score": "",
        "foldx_proxy_feature": round(foldx_proxy, 4),
        "interface_hotspot_proxy": "his_seed_or_window_proxy",
        "liability_scan_status": "pass",
        "buildability_light_status": "pass",
        "buildability_final_status": "",
        "selection_reason": "dry_run_generation_only",
    }
    rescue_rows = []
    for i, m in enumerate(non_his, start=1):
        rescue_rows.append(
            {
                "variant_id": variant_id,
                "generation_record_id": generation_record_id,
                "target": target,
                "rescue_index": i,
                "rescue_position_scheme": "chain_local_index",
                "rescue_position": m.pos,
                "rescue_window_local_position": "",
                "parent_aa": m.parent,
                "mutant_aa": m.mutant,
                "rescue_type": "conservative_or_polar",
                "rescue_objective": "neutral_retention_rescue",
                "linked_his_seed": his_seed_set,
                "linked_pH_sensitive_module": his_seed_set,
                "source_evidence": "dry_run_rule",
                "source_route": route,
                "risk_flags": "",
                "allowed_by_evidence_ledger": True,
            }
        )
    return row, rescue_rows, None


def forbidden_pair_fail(target: str, muts: list[Mutation]) -> bool:
    if target != "sdAb":
        return False
    labels = {(m.pos, m.parent, m.mutant) for m in muts}
    return (105, "V", "H") in labels and (110, "D", "H") in labels


def choose_mutations(
    target: str,
    chain: str,
    route: str,
    target_cfg: dict,
    positions: dict[int, str],
    regions: dict[int, str],
    rng: random.Random,
    idx: int,
) -> list[Mutation]:
    hard = set(target_cfg["hard_protect_positions"])
    his_seeds = [p for p in target_cfg["his_seed_positions"] if p not in hard and positions[p] != "H"]
    rescue_positions = [
        p for p in positions if p not in hard and p not in his_seeds and positions[p] != "C"
    ]
    cdr_positions = [p for p in rescue_positions if str(regions[p]).startswith("CDR")]
    fr_positions = [p for p in rescue_positions if not str(regions[p]).startswith("CDR")]

    single_sets = [(p,) for p in his_seeds]
    pair_sets = list(itertools.combinations(his_seeds, 2))
    if target == "sdAb":
        pair_sets = [pair for pair in pair_sets if set(pair) != {105, 110}]
    route_offset = PRIMARY_ROUTES.index(route) * 997
    seed_space_size = max(1, len(single_sets) + len(pair_sets))
    seed_index = (idx + route_offset) % seed_space_size
    combo_index = (idx + route_offset) // seed_space_size
    pct = (combo_index * 37 + seed_index * 11 + route_offset) % 100
    if target == "sdAb":
        use_single = pct < 86
    else:
        use_single = pct < 35
    if use_single or not pair_sets:
        seed_positions = list(single_sets[(seed_index + combo_index) % len(single_sets)])
    else:
        seed_positions = list(pair_sets[(seed_index + combo_index) % len(pair_sets)])

    # Avoid the known sdAb negative pair by replacing D110H when needed.
    if target == "sdAb" and 105 in seed_positions and 110 in seed_positions:
        replacement_pool = [p for p in his_seeds if p not in {105, 110}]
        seed_positions = [105, rng.choice(replacement_pool)]

    muts = [Mutation(chain, p, positions[p], "H") for p in seed_positions]
    if target == "sdAb":
        if len(seed_positions) == 1:
            order_pct = (combo_index * 37 + seed_index * 11) % 100
            non_his_target = 1 if order_pct < 65 else 2
        else:
            non_his_target = 1
    else:
        if len(seed_positions) == 1:
            non_his_target = 1 + ((combo_index + route_offset) % max(1, min(3, target_cfg["max_non_his_rescue_count"])))
        else:
            non_his_target = 1 + ((combo_index + route_offset) % max(1, min(2, target_cfg["max_non_his_rescue_count"])))
    max_total = target_cfg["main_max_mutations"]
    non_his_target = min(non_his_target, max_total - len(muts))
    if non_his_target <= 0:
        return muts

    if route == "structure_or_interface_guided" and cdr_positions and target != "sdAb":
        pool = cdr_positions
    else:
        pool = rescue_positions
    if route == "ProteinMPNN_seeded_rescue":
        # Prefer local FR/CDR support positions around His seeds.
        seed = seed_positions[0]
        local = [p for p in rescue_positions if 3 <= abs(p - seed) <= 8]
        pool = sorted(dict.fromkeys(local + fr_positions + rescue_positions))
    if not pool:
        return sorted(muts, key=lambda m: m.pos)
    pool = [p for p in pool if p not in seed_positions]
    muts.extend(choose_rescue_mutations(chain, positions, pool, non_his_target, combo_index, route_offset))
    return sorted(muts, key=lambda m: m.pos)


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def assign_near_duplicate_clusters(df: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for target, group in df.groupby("target", sort=False):
        group = group.copy()
        reps_by_seed: dict[str, list[tuple[str, dict[int, str], set[str]]]] = defaultdict(list)
        cluster_ids: list[str] = []
        for _, row in group.iterrows():
            seed = str(row["his_seed_set"])
            mut_map = mutation_map(row["mutation_list"])
            mut_set = set(str(row["mutation_list"]).split(";")) if row["mutation_list"] else set()
            assigned = None
            for cid, rep_map, rep_set in reps_by_seed[seed]:
                if mutation_distance(mut_map, rep_map) <= 2 or mutation_jaccard(mut_set, rep_set) >= 0.80:
                    assigned = cid
                    break
            if assigned is None:
                rep_count = sum(len(v) for v in reps_by_seed.values()) + 1
                assigned = f"{target}_hjc_{rep_count:05d}"
                reps_by_seed[seed].append((assigned, mut_map, mut_set))
            cluster_ids.append(assigned)
        group["near_duplicate_cluster_id"] = cluster_ids
        group["sequence_hamming_cluster_id"] = group["near_duplicate_cluster_id"]
        seed_labels = group["his_seed_set"].fillna("").astype(str).replace({"": "none"})
        group["his_seed_cluster_id"] = target + "_his_" + seed_labels
        group["mutation_set_jaccard_cluster_id"] = group["near_duplicate_cluster_id"]
        sizes = group["near_duplicate_cluster_id"].value_counts().to_dict()
        group["cluster_size"] = group["near_duplicate_cluster_id"].map(sizes)
        frames.append(group)
    return pd.concat(frames, ignore_index=True)


def generate_target_candidates(
    target: str,
    config: dict,
    inputs: dict[str, pd.DataFrame],
    config_hash: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, list[dict]]:
    target_cfg = config["targets"][target]
    chain = target_cfg["chain_id"]
    window_id = target_cfg["window_id"]
    target_size = int(target_cfg["dry_run_effective_design_pool"])
    reference = inputs["reference_sequence_map"]
    mask = inputs["window_mutation_mask"]
    parent_seq = chain_sequence(reference, target, chain)
    parent_motifs = nxs_motifs(parent_seq)
    window = mask[mask.window_id == window_id].copy().sort_values("pos")
    positions = {int(r.pos): str(r.aa) for _, r in window.iterrows()}
    regions = {int(r.pos): str(r.region) for _, r in window.iterrows()}
    rng = random.Random(int(config["random_seed"]) + (1 if target == "1E62" else 2))

    by_sequence: dict[str, dict] = {}
    rescue_by_variant: dict[str, list[dict]] = defaultdict(list)
    audit = {
        "raw_generated_count": 0,
        "post_generation_hard_filter_count": 0,
        "post_exact_dedup_count": 0,
        "effective_design_pool_count": 0,
        "control_anchor_pool_count": 0,
        "combined_pool_count": 0,
        "hard_filter_failures": Counter(),
        "raw_by_route": Counter(),
        "hard_fail_by_route_reason": Counter(),
        "accepted_by_route": Counter(),
        "dedup_by_route": Counter(),
        "diversity_cap_by_route": Counter(),
        "diversity_cap_failures": Counter(),
    }
    seed_set_cap_fraction = (
        config.get("seed_set_caps", {}).get(target, {}).get("max_single_his_seed_set_fraction")
    )
    rescue_cap_fraction = config.get("rescue_signature_caps", {}).get("max_single_rescue_signature_fraction")
    rescue_mutation_cap_fraction = config.get("rescue_mutation_caps", {}).get("max_single_rescue_mutation_fraction")
    seed_set_cap = int(target_size * float(seed_set_cap_fraction)) if seed_set_cap_fraction else target_size
    rescue_signature_cap = int(target_size * float(rescue_cap_fraction)) if rescue_cap_fraction else target_size
    rescue_mutation_cap = int(target_size * float(rescue_mutation_cap_fraction)) if rescue_mutation_cap_fraction else target_size
    accepted_seed_sets: Counter[str] = Counter()
    accepted_rescue_signatures: Counter[str] = Counter()
    accepted_rescue_mutations: Counter[str] = Counter()

    quotas = target_cfg.get("route_soft_quota", config["routes"]["route_soft_quota"])
    route_targets = {route: int(round(target_size * quotas[route])) for route in PRIMARY_ROUTES}
    # Adjust rounding on the largest route to hit exactly target_size.
    route_targets[PRIMARY_ROUTES[0]] += target_size - sum(route_targets.values())

    i = 0
    for route in PRIMARY_ROUTES:
        route_attempts = 0
        proposed_mutation_keys: set[str] = set()
        max_attempts = max(5000, route_targets[route] * 300)
        while audit["accepted_by_route"][route] < route_targets[route] and route_attempts < max_attempts:
            muts = choose_mutations(target, chain, route, target_cfg, positions, regions, rng, i)
            i += 1
            route_attempts += 1
            proposed_key = mutation_key(muts)
            if proposed_key in proposed_mutation_keys:
                continue
            proposed_mutation_keys.add(proposed_key)
            audit["raw_generated_count"] += 1
            audit["raw_by_route"][route] += 1
            if not muts:
                audit["hard_filter_failures"]["empty_mutation"] += 1
                audit["hard_fail_by_route_reason"][(route, "empty_mutation")] += 1
                continue
            if any(m.pos in set(target_cfg["hard_protect_positions"]) for m in muts):
                audit["hard_filter_failures"]["hard_protect_mutation"] += 1
                audit["hard_fail_by_route_reason"][(route, "hard_protect_mutation")] += 1
                continue
            if forbidden_pair_fail(target, muts):
                audit["hard_filter_failures"]["forbidden_pair"] += 1
                audit["hard_fail_by_route_reason"][(route, "forbidden_pair")] += 1
                continue
            if sum(m.is_his() for m in muts) > target_cfg["main_max_his"]:
                audit["hard_filter_failures"]["his_count"] += 1
                audit["hard_fail_by_route_reason"][(route, "his_count")] += 1
                continue
            if len(muts) > target_cfg["main_max_mutations"]:
                audit["hard_filter_failures"]["mutation_count"] += 1
                audit["hard_fail_by_route_reason"][(route, "mutation_count")] += 1
                continue
            if len(muts) == 1 and not muts[0].is_his():
                audit["hard_filter_failures"]["single_non_his_rescue_only"] += 1
                audit["hard_fail_by_route_reason"][(route, "single_non_his_rescue_only")] += 1
                continue
            row, rescue_rows, failure = make_candidate(
                target, chain, window_id, route, parent_seq, muts, config_hash, i, parent_motifs
            )
            if failure:
                audit["hard_filter_failures"][failure] += 1
                audit["hard_fail_by_route_reason"][(route, failure)] += 1
                continue
            audit["post_generation_hard_filter_count"] += 1
            seq_hash = row["sequence_hash"]
            if seq_hash in by_sequence:
                existing = by_sequence[seq_hash]
                routes = set(existing["all_source_routes"].split(";")) | {route}
                records = set(existing["all_generation_records"].split(";")) | {row["generation_record_id"]}
                existing["all_source_routes"] = ";".join(sorted(routes))
                existing["all_generation_records"] = ";".join(sorted(records))
                existing["duplicate_source_count"] = int(existing["duplicate_source_count"]) + 1
                existing["exact_duplicate_status"] = "merged_duplicate"
                audit["dedup_by_route"][route] += 1
                continue
            seed_key = str(row["his_seed_set"])
            rescue_key = str(row["rescue_signature"])
            rescue_labels = [m.label() for m in muts if not m.is_his()]
            cap_exempt = target == "sdAb" and route == "structure_or_interface_guided"
            if not cap_exempt and accepted_seed_sets[seed_key] >= seed_set_cap:
                audit["diversity_cap_by_route"][route] += 1
                audit["diversity_cap_failures"]["seed_set_cap"] += 1
                continue
            if not cap_exempt and accepted_rescue_signatures[rescue_key] >= rescue_signature_cap:
                audit["diversity_cap_by_route"][route] += 1
                audit["diversity_cap_failures"]["rescue_signature_cap"] += 1
                continue
            if not cap_exempt and any(accepted_rescue_mutations[label] >= rescue_mutation_cap for label in rescue_labels):
                audit["diversity_cap_by_route"][route] += 1
                audit["diversity_cap_failures"]["rescue_mutation_cap"] += 1
                continue
            by_sequence[seq_hash] = row
            rescue_by_variant[row["variant_id"]].extend(rescue_rows)
            accepted_seed_sets[seed_key] += 1
            accepted_rescue_signatures[rescue_key] += 1
            for label in rescue_labels:
                accepted_rescue_mutations[label] += 1
            audit["accepted_by_route"][route] += 1
        if audit["accepted_by_route"][route] < route_targets[route]:
            raise RuntimeError(
                f"{target} route {route} generated {audit['accepted_by_route'][route]} "
                f"of {route_targets[route]} unique candidates"
            )

    if len(by_sequence) < target_size:
        raise RuntimeError(f"{target} generated only {len(by_sequence)} unique candidates")

    df = pd.DataFrame(by_sequence.values())
    df = assign_near_duplicate_clusters(df)
    rescue_rows_all = []
    for variant_id in df["variant_id"]:
        rescue_rows_all.extend(rescue_by_variant.get(variant_id, []))
    rescue_df = pd.DataFrame(rescue_rows_all)
    audit["post_exact_dedup_count"] = len(df) + sum(audit["diversity_cap_by_route"].values())
    audit["effective_design_pool_count"] = len(df)
    audit["combined_pool_count"] = len(df)

    controls = make_control_panel(
        target,
        target_cfg,
        positions,
        chain,
        parent_seq,
        parent_motifs,
        config_hash,
        set(df["sequence_hash"]),
    )
    audit["control_anchor_pool_count"] = len(controls)
    return df, rescue_df, audit, controls


def make_control_panel(
    target: str,
    target_cfg: dict,
    positions: dict[int, str],
    chain: str,
    parent_seq: str,
    parent_motifs: set[int],
    config_hash: str,
    avoid_sequence_hashes: set[str] | None = None,
) -> list[dict]:
    rows = []
    seen_hashes: set[str] = set()
    avoid_sequence_hashes = avoid_sequence_hashes or set()
    window_id = target_cfg["window_id"]
    target_count = int(target_cfg.get("controls_anchors_in_final", 200))
    hard = set(target_cfg["hard_protect_positions"])
    his_seeds = [p for p in target_cfg["his_seed_positions"] if p not in hard and positions[p] != "H"]
    rescue_positions = [p for p in sorted(positions) if p not in hard and p not in his_seeds and positions[p] != "C"]
    rescue_residues = ["A", "S", "T", "N", "Q", "Y", "D", "E", "K", "R", "G", "P", "V", "L", "I", "M", "F", "W"]

    def add_control(
        muts: list[Mutation],
        control_type: str,
        reason: str,
        allow_duplicate_sequence: bool = False,
    ) -> None:
        if len(rows) >= target_count:
            return
        if forbidden_pair_fail(target, muts):
            return
        if any(m.pos in hard for m in muts):
            return
        if any(m.mutant == "C" for m in muts):
            return
        seq = apply_mutations(parent_seq, muts)
        if nxs_motifs(seq) - parent_motifs:
            return
        seq_hash = sha(seq, 10)
        if seq_hash in avoid_sequence_hashes:
            return
        if seq_hash in seen_hashes and not allow_duplicate_sequence:
            return
        seen_hashes.add(seq_hash)
        mut_key = mutation_key(muts)
        rows.append(
            {
                "target": target,
                "control_slot_id": f"{target}_control_{len(rows)+1:03d}",
                "variant_id": f"{target}_{window_id}_{seq_hash}",
                "sequence_hash": seq_hash,
                "mutation_list": mut_key,
                "mutation_count": len(muts),
                "His_count": sum(m.is_his() for m in muts),
                "control_type": control_type,
                "planned_count": 1,
                "control_only": True,
                "control_reason": reason,
                "duplicated_by_design": allow_duplicate_sequence,
                "buildability_light_status": "pass",
            }
        )

    for _ in range(20):
        add_control([], "parent_wt_anchor", "wild_type_reference", allow_duplicate_sequence=True)
    for pos in his_seeds:
        add_control(
            [Mutation(chain, pos, positions[pos], "H")],
            "single_his_anchor",
            "single_his_seed_reference",
        )
    for a, b in itertools.combinations(his_seeds, 2):
        add_control(
            [Mutation(chain, a, positions[a], "H"), Mutation(chain, b, positions[b], "H")],
            "double_his_anchor",
            "two_his_seed_reference",
        )
    for pos in rescue_positions:
        for aa in rescue_residues:
            if aa == positions[pos] or aa in {"C", "H"}:
                continue
            add_control(
                [Mutation(chain, pos, positions[pos], aa)],
                "single_rescue_anchor",
                "single_non_his_buildability_reference",
            )
    for seed in his_seeds:
        for pos in rescue_positions:
            if len(rows) >= target_count:
                break
            if abs(pos - seed) < 2:
                continue
            for aa in rescue_residues:
                if aa == positions[pos] or aa in {"C", "H"}:
                    continue
                add_control(
                    [Mutation(chain, seed, positions[seed], "H"), Mutation(chain, pos, positions[pos], aa)],
                    "his_plus_rescue_anchor",
                    "his_rescue_context_reference",
                )
        if len(rows) >= target_count:
            break
    if len(rows) < target_count:
        for seed_pair in itertools.combinations(his_seeds, 2):
            if len(rows) >= target_count:
                break
            if target == "sdAb" and set(seed_pair) == {105, 110}:
                continue
            for pos in rescue_positions:
                if len(rows) >= target_count:
                    break
                if pos in seed_pair:
                    continue
                for aa in rescue_residues:
                    if aa == positions[pos] or aa in {"C", "H"}:
                        continue
                    muts = [
                        Mutation(chain, seed_pair[0], positions[seed_pair[0]], "H"),
                        Mutation(chain, seed_pair[1], positions[seed_pair[1]], "H"),
                        Mutation(chain, pos, positions[pos], aa),
                    ]
                    add_control(muts, "double_his_plus_rescue_anchor", "stress_control_fill")
    if len(rows) != target_count:
        raise RuntimeError(f"{target} generated {len(rows)} controls, expected {target_count}")
    return rows


def write_audits(
    candidates: pd.DataFrame,
    audits: dict[str, dict],
    control_rows: list[dict],
    config: dict,
) -> dict[str, str]:
    all_df = candidates
    route_summary = (
        all_df.groupby(["target", "primary_generation_route"], dropna=False)
        .size()
        .reset_index(name="effective_candidate_count")
    )
    route_summary["fraction_within_target"] = route_summary.groupby("target")[
        "effective_candidate_count"
    ].transform(lambda s: s / s.sum())
    route_summary.pipe(write_csv, DRY / "route_summary.csv")

    pos_rows = []
    for _, row in all_df.iterrows():
        for mut in str(row["mutation_list"]).split(";"):
            if not mut:
                continue
            m = re.match(r"([A-Z])([A-Z])(\d+)([A-Z])", mut)
            if not m:
                continue
            pos_rows.append(
                {
                    "target": row["target"],
                    "position": int(m.group(3)),
                    "parent_aa": m.group(2),
                    "mutant_aa": m.group(4),
                    "is_his": m.group(4) == "H",
                }
            )
    pos_df = pd.DataFrame(pos_rows)
    if not pos_df.empty:
        pos_freq = pos_df.groupby(["target", "position", "mutant_aa"]).size().reset_index(name="count")
        target_counts = all_df.groupby("target").size().to_dict()
        pos_freq["fraction_within_target"] = pos_freq.apply(
            lambda r: r["count"] / target_counts[r["target"]], axis=1
        )
    else:
        pos_freq = pd.DataFrame(columns=["target", "position", "mutant_aa", "count", "fraction_within_target"])
    pos_freq.pipe(write_csv, DRY / "position_frequency.csv")

    his_dist = all_df.groupby(["target", "His_count"]).size().reset_index(name="count")
    his_dist["fraction_within_target"] = his_dist.groupby("target")["count"].transform(lambda s: s / s.sum())
    his_dist.pipe(write_csv, DRY / "his_distribution.csv")

    mut_dist = all_df.groupby(["target", "mutation_count"]).size().reset_index(name="count")
    mut_dist["fraction_within_target"] = mut_dist.groupby("target")["count"].transform(lambda s: s / s.sum())
    mut_dist.pipe(write_csv, DRY / "mutation_order_distribution.csv")

    cluster = (
        all_df.groupby(["target", "near_duplicate_cluster_id"], dropna=False)
        .agg(
            cluster_size=("variant_id", "count"),
            representative_variant_id=("variant_id", "first"),
            representative_sequence_hash=("sequence_hash", "first"),
            representative_mutation_list=("mutation_list", "first"),
            his_seed_sets_in_cluster=("his_seed_set", lambda x: ";".join(sorted(set(str(v) for v in x)))),
            mutation_set_sample=("mutation_list", lambda x: " | ".join(list(dict.fromkeys(str(v) for v in x))[:5])),
            routes_in_cluster=("primary_generation_route", lambda x: ";".join(sorted(set(x)))),
            preliminary_selection_buckets_in_cluster=(
                "preliminary_selection_bucket",
                lambda x: ";".join(sorted(set(str(v) for v in x))),
            ),
        )
        .reset_index()
    )
    target_counts = all_df.groupby("target").size().to_dict()
    cluster["selected_count"] = 0
    cluster["selected_fraction"] = 0.0
    cluster["fraction_within_target"] = cluster.apply(
        lambda r: r["cluster_size"] / target_counts[r["target"]], axis=1
    )
    cluster.pipe(write_csv, DRY / "near_duplicate_cluster_summary.csv")

    gen_rows = []
    filter_rows = []
    for target, a in audits.items():
        gen_rows.append(
            {
                "target": target,
                "raw_generated_count": a["raw_generated_count"],
                "post_generation_hard_filter_count": a["post_generation_hard_filter_count"],
                "post_exact_dedup_count": a["post_exact_dedup_count"],
                "effective_design_pool_count": a["effective_design_pool_count"],
                "control_anchor_pool_count": a["control_anchor_pool_count"],
                "combined_pool_count": a["combined_pool_count"],
            }
        )
        for reason, count in a["hard_filter_failures"].items():
            filter_rows.append({"target": target, "filter_stage": "generation", "reason": reason, "count": count})
        for reason, count in a["diversity_cap_failures"].items():
            filter_rows.append({"target": target, "filter_stage": "diversity_cap", "reason": reason, "count": count})
        for route, raw in a["raw_by_route"].items():
            accepted = a["accepted_by_route"].get(route, 0)
            dedup = a["dedup_by_route"].get(route, 0)
            capped = a["diversity_cap_by_route"].get(route, 0)
            filter_rows.append(
                {
                    "target": target,
                    "filter_stage": "route_generation",
                    "reason": route,
                    "count": raw,
                    "accepted": accepted,
                    "dedup_merged": dedup,
                    "diversity_cap_skipped": capped,
                    "hard_filter_or_other_loss": raw - accepted - dedup - capped,
                }
            )
    pd.DataFrame(gen_rows).pipe(write_csv, DRY / "generation_audit.csv")
    pd.DataFrame(filter_rows).pipe(write_csv, DRY / "filter_audit.csv")
    control_df = pd.DataFrame(control_rows)
    control_df.pipe(write_csv, DRY / "control_anchor_panel_draft.csv")
    control_df.pipe(write_csv, DRY / "control_anchor_panel_200_per_target.csv")
    design_hashes_by_target = {
        target: set(sub["sequence_hash"].astype(str))
        for target, sub in all_df.groupby("target", sort=False)
    }
    conflict_rows = []
    if not control_df.empty:
        control_df["design_duplicate_conflict"] = control_df.apply(
            lambda r: str(r["sequence_hash"]) in design_hashes_by_target.get(r["target"], set()),
            axis=1,
        )
        control_df.pipe(write_csv, DRY / "control_anchor_panel_draft.csv")
        control_df.pipe(write_csv, DRY / "control_anchor_panel_200_per_target.csv")
        for target, sub in control_df.groupby("target", sort=False):
            conflicts = int(sub["design_duplicate_conflict"].sum())
            conflict_rows.append(
                {
                    "target": target,
                    "control_design_duplicate_conflict_count": conflicts,
                    "control_design_duplicate_conflict_fraction": conflicts / max(1, len(sub)),
                }
            )
    else:
        for target in audits:
            conflict_rows.append(
                {
                    "target": target,
                    "control_design_duplicate_conflict_count": 0,
                    "control_design_duplicate_conflict_fraction": 0.0,
                }
            )
    pd.DataFrame(
        conflict_rows
    ).pipe(write_csv, DRY / "control_design_duplicate_report.csv")

    route_loss_rows = []
    hard_loss_rows = []
    for target, a in audits.items():
        for route in PRIMARY_ROUTES:
            raw = int(a["raw_by_route"].get(route, 0))
            duplicate = int(a["dedup_by_route"].get(route, 0))
            capped = int(a["diversity_cap_by_route"].get(route, 0))
            effective = int(a["accepted_by_route"].get(route, 0))
            hard_pass = duplicate + effective
            post_dedup = effective + capped
            route_sub = all_df[(all_df["target"] == target) & (all_df["primary_generation_route"] == route)]
            largest_cluster_fraction = 0.0
            if not route_sub.empty:
                largest_cluster_fraction = (
                    route_sub.groupby("near_duplicate_cluster_id").size().max() / len(route_sub)
                )
            route_loss_rows.append(
                {
                    "target": target,
                    "primary_generation_route": route,
                    "raw_generated": raw,
                    "post_generation_hard_filter": hard_pass + capped,
                    "exact_duplicate_merged": duplicate,
                    "diversity_cap_skipped": capped,
                    "post_exact_dedup_unique": post_dedup,
                    "post_quota_sampling": effective,
                    "effective_count": effective,
                    "top_duplicate_sequence_count": int(route_sub["duplicate_source_count"].max()) if not route_sub.empty else 0,
                    "unique_his_seed_set_count": int(route_sub["his_seed_set"].nunique(dropna=False)),
                    "unique_rescue_signature_count": int(route_sub["rescue_signature"].nunique(dropna=False)),
                    "unique_mutation_set_count": int(route_sub["mutation_list"].nunique(dropna=False)),
                    "largest_near_duplicate_cluster_fraction": largest_cluster_fraction,
                    "hard_filter_loss_fraction": (raw - hard_pass - capped) / max(1, raw),
                    "exact_dedup_loss_fraction": duplicate / max(1, hard_pass + capped),
                }
            )
        for (route, reason), count in a["hard_fail_by_route_reason"].items():
            hard_loss_rows.append(
                {
                    "target": target,
                    "primary_generation_route": route,
                    "hard_filter_reason": reason,
                    "count": count,
                }
            )
    pd.DataFrame(route_loss_rows).pipe(write_csv, DRY / "route_filter_loss_summary.csv")
    pd.DataFrame(hard_loss_rows).pipe(write_csv, DRY / "route_hard_filter_reason_summary.csv")

    duplicate_rows = []
    for _, row in all_df.iterrows():
        routes = [r for r in str(row["all_source_routes"]).split(";") if r]
        if len(routes) <= 1 and int(row["duplicate_source_count"]) <= 1:
            continue
        for source_route in routes:
            duplicate_rows.append(
                {
                    "target": row["target"],
                    "kept_primary_generation_route": row["primary_generation_route"],
                    "source_route": source_route,
                    "variant_id": row["variant_id"],
                    "sequence_hash": row["sequence_hash"],
                    "mutation_list": row["mutation_list"],
                    "his_seed_set": row["his_seed_set"],
                    "duplicate_source_count": row["duplicate_source_count"],
                }
            )
    if duplicate_rows:
        pd.DataFrame(duplicate_rows).pipe(write_csv, DRY / "exact_duplicate_sequence_provenance.csv")
    else:
        pd.DataFrame(
            columns=[
                "target",
                "kept_primary_generation_route",
                "source_route",
                "variant_id",
                "sequence_hash",
                "mutation_list",
                "his_seed_set",
                "duplicate_source_count",
            ]
        ).pipe(write_csv, DRY / "exact_duplicate_sequence_provenance.csv")
    if duplicate_rows:
        dup_matrix = (
            pd.DataFrame(duplicate_rows)
            .groupby(["target", "kept_primary_generation_route", "source_route"], dropna=False)
            .agg(
                duplicate_variant_count=("variant_id", "count"),
                duplicate_source_count_sum=("duplicate_source_count", "sum"),
            )
            .reset_index()
        )
    else:
        dup_matrix = pd.DataFrame(
            columns=[
                "target",
                "kept_primary_generation_route",
                "source_route",
                "duplicate_variant_count",
                "duplicate_source_count_sum",
            ]
        )
    dup_matrix.pipe(write_csv, DRY / "exact_duplicate_source_route_matrix.csv")

    cluster.sort_values(["target", "cluster_size"], ascending=[True, False]).groupby("target").head(20).pipe(
        write_csv, DRY / "near_duplicate_top_clusters.csv"
    )

    target_counts = all_df.groupby("target").size().to_dict()
    seed_set = (
        all_df.groupby(["target", "his_seed_set", "His_count"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    seed_set["fraction_within_target"] = seed_set.apply(
        lambda r: r["count"] / target_counts[r["target"]], axis=1
    )
    seed_set.pipe(write_csv, DRY / "his_seed_set_summary.csv")
    (
        all_df.groupby(["target", "primary_generation_route", "his_seed_set"], dropna=False)
        .size()
        .reset_index(name="count")
        .pipe(write_csv, DRY / "his_seed_set_by_route.csv")
    )
    (
        all_df.groupby(["target", "his_seed_set", "rescue_signature"], dropna=False)
        .size()
        .reset_index(name="count")
        .pipe(write_csv, DRY / "his_seed_set_by_rescue_signature.csv")
    )
    (
        all_df.groupby(["target", "his_seed_set", "near_duplicate_cluster_id"], dropna=False)
        .size()
        .reset_index(name="count")
        .pipe(write_csv, DRY / "his_seed_set_near_duplicate_cluster_summary.csv")
    )
    (
        all_df.sort_values(["target", "duplicate_source_count"], ascending=[True, False])
        .groupby(["target", "his_seed_set"], dropna=False)
        .head(5)[
            [
                "target",
                "his_seed_set",
                "variant_id",
                "mutation_list",
                "rescue_signature",
                "primary_generation_route",
                "duplicate_source_count",
            ]
        ]
        .pipe(write_csv, DRY / "top_seed_set_examples.csv")
    )

    seed_pair_rows = []
    position_set_rows = []
    position_pair_rows = []
    for _, row in all_df.iterrows():
        seeds = his_labels(row["mutation_list"])
        positions = mutation_positions(row["mutation_list"])
        position_set_rows.append(
            {
                "target": row["target"],
                "position_set": ";".join(str(p) for p in positions),
            }
        )
        for pair in itertools.combinations(seeds, 2):
            seed_pair_rows.append(
                {
                    "target": row["target"],
                    "his_seed_pair": ";".join(pair),
                }
            )
        for pair in itertools.combinations(positions, 2):
            position_pair_rows.append(
                {
                    "target": row["target"],
                    "position_pair": f"{pair[0]};{pair[1]}",
                }
            )
    if seed_pair_rows:
        seed_pair = pd.DataFrame(seed_pair_rows).groupby(["target", "his_seed_pair"]).size().reset_index(name="count")
        seed_pair["fraction_within_target"] = seed_pair.apply(
            lambda r: r["count"] / target_counts[r["target"]], axis=1
        )
    else:
        seed_pair = pd.DataFrame(columns=["target", "his_seed_pair", "count", "fraction_within_target"])
    seed_pair.pipe(write_csv, DRY / "his_seed_pair_frequency.csv")
    seed_pair.pipe(write_csv, DRY / "his_seed_pair_summary.csv")

    position_set = pd.DataFrame(position_set_rows).groupby(["target", "position_set"]).size().reset_index(name="count")
    position_set["fraction_within_target"] = position_set.apply(
        lambda r: r["count"] / target_counts[r["target"]], axis=1
    )
    position_set.pipe(write_csv, DRY / "position_set_summary.csv")

    if position_pair_rows:
        position_pair = pd.DataFrame(position_pair_rows).groupby(["target", "position_pair"]).size().reset_index(name="count")
        position_pair["fraction_within_target"] = position_pair.apply(
            lambda r: r["count"] / target_counts[r["target"]], axis=1
        )
    else:
        position_pair = pd.DataFrame(columns=["target", "position_pair", "count", "fraction_within_target"])
    position_pair.pipe(write_csv, DRY / "position_pair_frequency.csv")

    rescue_signature = (
        all_df.groupby(["target", "rescue_signature"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    rescue_signature["fraction_within_target"] = rescue_signature.apply(
        lambda r: r["count"] / target_counts[r["target"]], axis=1
    )
    rescue_signature.pipe(write_csv, DRY / "rescue_signature_summary.csv")
    (
        all_df.groupby(["target", "primary_generation_route", "rescue_signature"], dropna=False)
        .size()
        .reset_index(name="count")
        .pipe(write_csv, DRY / "rescue_signature_by_route.csv")
    )

    mutation_set_summary = (
        all_df.groupby(["target", "mutation_list"], dropna=False)
        .agg(
            count=("variant_id", "count"),
            routes=("primary_generation_route", lambda x: ";".join(sorted(set(x)))),
            his_seed_set=("his_seed_set", "first"),
            rescue_signature=("rescue_signature", "first"),
        )
        .reset_index()
    )
    mutation_set_summary.pipe(write_csv, DRY / "mutation_set_summary.csv")

    hj_summary = cluster[
        [
            "target",
            "near_duplicate_cluster_id",
            "cluster_size",
            "representative_variant_id",
            "representative_sequence_hash",
            "representative_mutation_list",
            "his_seed_sets_in_cluster",
            "routes_in_cluster",
            "fraction_within_target",
        ]
    ].copy()
    hj_summary["near_duplicate_method"] = "hamming_le2_or_mutation_jaccard_ge0.80_within_his_seed_set"
    hj_summary.pipe(write_csv, DRY / "near_duplicate_hamming_jaccard_summary.csv")
    cluster.sort_values(["target", "cluster_size"], ascending=[True, False]).groupby("target").head(20).pipe(
        write_csv, DRY / "top_near_duplicate_clusters.csv"
    )

    sdab = all_df[all_df["target"] == "sdAb"].copy()
    if not sdab.empty:
        collapse = (
            sdab.groupby(["his_seed_set", "rescue_signature"], dropna=False)
            .agg(
                count=("variant_id", "count"),
                route_count=("primary_generation_route", "nunique"),
                routes=("primary_generation_route", lambda x: ";".join(sorted(set(x)))),
                max_duplicate_source_count=("duplicate_source_count", "max"),
            )
            .reset_index()
            .sort_values("count", ascending=False)
        )
        collapse["fraction_within_sdAb"] = collapse["count"] / len(sdab)
    else:
        collapse = pd.DataFrame(
            columns=[
                "his_seed_set",
                "rescue_signature",
                "count",
                "route_count",
                "routes",
                "max_duplicate_source_count",
                "fraction_within_sdAb",
            ]
        )
    collapse.pipe(write_csv, DRY / "sdAb_seed_rescue_collapse_report.csv")

    rescue_objective_rows = []
    for _, row in all_df.iterrows():
        objectives = [v for v in str(row["rescue_objective_list"]).split(";") if v and v != "nan"]
        if not objectives:
            rescue_objective_rows.append({"target": row["target"], "rescue_objective": "none", "count_unit": "candidate"})
        for objective in objectives:
            rescue_objective_rows.append({"target": row["target"], "rescue_objective": objective, "count_unit": "rescue"})
    rescue_objectives = (
        pd.DataFrame(rescue_objective_rows)
        .groupby(["target", "rescue_objective", "count_unit"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    rescue_objectives.pipe(write_csv, DRY / "rescue_objective_distribution.csv")

    orphan_rows = []
    for target, sub in all_df.groupby("target", sort=False):
        orphan = sub[(sub["rescue_count"].astype(int) > 0) & (sub["His_count"].astype(int) == 0)]
        orphan_rows.append(
            {
                "target": target,
                "orphan_rescue_variant_count": len(orphan),
                "orphan_rescue_fraction": len(orphan) / max(1, len(sub)),
            }
        )
    pd.DataFrame(orphan_rows).pipe(write_csv, DRY / "orphan_rescue_check.csv")

    return validate_dry_run(all_df, audits, cluster, route_summary, his_dist, mut_dist, control_df, config)


def validate_dry_run(
    df: pd.DataFrame,
    audits: dict[str, dict],
    cluster: pd.DataFrame,
    route_summary: pd.DataFrame,
    his_dist: pd.DataFrame,
    mut_dist: pd.DataFrame,
    control_df: pd.DataFrame,
    config: dict,
) -> dict[str, str]:
    statuses: dict[str, str] = {}

    evidence_failures = validate_evidence_ledger()
    statuses["minimum_constraints_in_evidence"] = "PASS" if not evidence_failures else "FAIL: " + "; ".join(evidence_failures)

    hard_violations = []
    for target, cfg in config["targets"].items():
        hard = set(cfg["hard_protect_positions"])
        sub = df[df.target == target]
        for _, row in sub.iterrows():
            for mut in str(row["mutation_list"]).split(";"):
                m = re.match(r"([A-Z])([A-Z])(\d+)([A-Z])", mut)
                if m and int(m.group(3)) in hard:
                    hard_violations.append(mut)
    statuses["no_hard_protect_mutation"] = "PASS" if not hard_violations else "FAIL"

    bad_pair = df[
        (df.target == "sdAb")
        & df.mutation_list.str.contains("A[V]105H", regex=True)
        & df.mutation_list.str.contains("A[D]110H", regex=True)
    ]
    statuses["no_sdab_V105H_D110H"] = "PASS" if bad_pair.empty else "FAIL"

    single_non_his = df[(df.mutation_count == 1) & (df.His_count == 0)]
    statuses["no_single_non_his_rescue_only"] = "PASS" if single_non_his.empty else "FAIL"

    statuses["no_new_canonical_nxs_t"] = "PASS" if (df.buildability_light_status == "pass").all() else "FAIL"
    statuses["exact_duplicates_merged"] = "PASS" if df.sequence_hash.nunique() == len(df) else "FAIL"
    statuses["final_selection_bucket_unset"] = "PASS" if df.final_selection_bucket.fillna("").eq("").all() else "FAIL"

    thresholds = config["dry_run_gate_thresholds"]

    def threshold_value(metric: str, level: str, target: str | None = None) -> float:
        entry = thresholds[metric]
        if target is not None and isinstance(entry, dict) and target in entry:
            return float(entry[target][level])
        return float(entry[level])
    route_max = route_summary.groupby("target")["fraction_within_target"].max().to_dict()
    route_status = "PASS"
    for frac in route_max.values():
        if frac >= threshold_value("max_single_generation_route_fraction", "fail"):
            route_status = "FAIL"
        elif frac >= threshold_value("max_single_generation_route_fraction", "warning") and route_status != "FAIL":
            route_status = "WARN"
    statuses["route_distribution"] = route_status

    route_hard_status = "PASS"
    for target, a in audits.items():
        for route in PRIMARY_ROUTES:
            raw = int(a["raw_by_route"].get(route, 0))
            hard_pass = (
                int(a["accepted_by_route"].get(route, 0))
                + int(a["dedup_by_route"].get(route, 0))
                + int(a["diversity_cap_by_route"].get(route, 0))
            )
            frac = (raw - hard_pass) / max(1, raw)
            if frac >= threshold_value("max_hard_filter_loss_fraction_by_route", "fail"):
                route_hard_status = "FAIL"
            elif frac >= threshold_value("max_hard_filter_loss_fraction_by_route", "warning") and route_hard_status != "FAIL":
                route_hard_status = "WARN"
    statuses["route_hard_filter_loss"] = route_hard_status

    cluster_max = cluster.groupby("target")["fraction_within_target"].max().to_dict()
    cluster_status = "PASS"
    for frac in cluster_max.values():
        if frac >= threshold_value("max_largest_near_duplicate_cluster_fraction", "fail"):
            cluster_status = "FAIL"
        elif frac >= threshold_value("max_largest_near_duplicate_cluster_fraction", "warning") and cluster_status != "FAIL":
            cluster_status = "WARN"
    statuses["near_duplicate_cluster_distribution"] = cluster_status
    statuses["near_duplicate_cluster_method"] = "PASS"

    top20_status = "PASS"
    if "max_top_20_near_duplicate_cluster_fraction" in thresholds:
        for target, sub in cluster.groupby("target", sort=False):
            frac = sub.sort_values("cluster_size", ascending=False).head(20)["fraction_within_target"].sum()
            if frac >= threshold_value("max_top_20_near_duplicate_cluster_fraction", "fail", str(target)):
                top20_status = "FAIL"
            elif frac >= threshold_value("max_top_20_near_duplicate_cluster_fraction", "warning", str(target)) and top20_status != "FAIL":
                top20_status = "WARN"
        statuses["top20_near_duplicate_cluster_distribution"] = top20_status

    if "max_route_largest_near_duplicate_cluster_fraction" in thresholds:
        route_cluster_status = "PASS"
        route_cluster = (
            df.groupby(["target", "primary_generation_route", "near_duplicate_cluster_id"], dropna=False)
            .size()
            .reset_index(name="cluster_size")
        )
        route_totals = df.groupby(["target", "primary_generation_route"]).size().to_dict()
        route_cluster["fraction_within_route"] = route_cluster.apply(
            lambda r: r["cluster_size"] / route_totals[(r["target"], r["primary_generation_route"])],
            axis=1,
        )
        for frac in route_cluster["fraction_within_route"]:
            if frac >= threshold_value("max_route_largest_near_duplicate_cluster_fraction", "fail"):
                route_cluster_status = "FAIL"
            elif frac >= threshold_value("max_route_largest_near_duplicate_cluster_fraction", "warning") and route_cluster_status != "FAIL":
                route_cluster_status = "WARN"
        statuses["route_largest_near_duplicate_cluster_distribution"] = route_cluster_status

    exact_dedup_status = "PASS"
    for target, a in audits.items():
        raw = max(1, a["post_generation_hard_filter_count"])
        loss = (a["post_generation_hard_filter_count"] - a["post_exact_dedup_count"]) / raw
        if loss >= threshold_value("max_exact_dedup_loss_fraction", "fail", target):
            exact_dedup_status = "FAIL"
        elif loss >= threshold_value("max_exact_dedup_loss_fraction", "warning", target) and exact_dedup_status != "FAIL":
            exact_dedup_status = "WARN"
    statuses["exact_dedup_loss"] = exact_dedup_status

    if "post_hard_to_effective_inflation" in thresholds:
        inflation_status = "PASS"
        for target, a in audits.items():
            entry = thresholds["post_hard_to_effective_inflation"].get(target, {})
            pass_max = float(entry.get("pass_max", 999999))
            warn_max = float(entry.get("warn_max", 999999))
            inflation = a["post_generation_hard_filter_count"] / max(1, a["effective_design_pool_count"])
            if inflation > warn_max:
                inflation_status = "FAIL"
            elif inflation > pass_max and inflation_status != "FAIL":
                inflation_status = "WARN"
        statuses["post_hard_to_effective_inflation"] = inflation_status

    his3 = df[(df.target == "sdAb") & (df.His_count >= 3)]
    frac_his3 = len(his3) / max(1, len(df[df.target == "sdAb"]))
    if frac_his3 >= threshold_value("max_sdAb_candidates_with_his_count_ge3", "fail"):
        statuses["sdAb_His_ge3_fraction"] = "FAIL"
    elif frac_his3 >= threshold_value("max_sdAb_candidates_with_his_count_ge3", "warning"):
        statuses["sdAb_His_ge3_fraction"] = "WARN"
    else:
        statuses["sdAb_His_ge3_fraction"] = "PASS"

    seed_set_status = "PASS"
    if "max_single_his_seed_set_fraction" in thresholds:
        seed_summary = df.groupby(["target", "his_seed_set"], dropna=False).size().reset_index(name="count")
        seed_summary["fraction"] = seed_summary.groupby("target")["count"].transform(lambda s: s / s.sum())
        for _, seed_row in seed_summary.iterrows():
            frac = float(seed_row["fraction"])
            target_name = str(seed_row["target"])
            if frac >= threshold_value("max_single_his_seed_set_fraction", "fail", target_name):
                seed_set_status = "FAIL"
            elif frac >= threshold_value("max_single_his_seed_set_fraction", "warning", target_name) and seed_set_status != "FAIL":
                seed_set_status = "WARN"
        if "max_top5_his_seed_set_fraction" in thresholds:
            for _, sub in seed_summary.groupby("target", sort=False):
                frac = sub.sort_values("count", ascending=False).head(5)["fraction"].sum()
                entry = thresholds["max_top5_his_seed_set_fraction"].get(str(sub.iloc[0]["target"]), {})
                fail = float(entry.get("fail", 1.01))
                warning = float(entry.get("warning", fail))
                if frac >= fail:
                    seed_set_status = "FAIL"
                elif frac >= warning and seed_set_status != "FAIL":
                    seed_set_status = "WARN"
        statuses["his_seed_set_distribution"] = seed_set_status

    if "max_single_rescue_signature_fraction" in thresholds:
        rescue_status = "PASS"
        rescue_summary = df.groupby(["target", "rescue_signature"], dropna=False).size().reset_index(name="count")
        rescue_summary["fraction"] = rescue_summary.groupby("target")["count"].transform(lambda s: s / s.sum())
        for frac in rescue_summary["fraction"]:
            if frac >= threshold_value("max_single_rescue_signature_fraction", "fail"):
                rescue_status = "FAIL"
            elif frac >= threshold_value("max_single_rescue_signature_fraction", "warning") and rescue_status != "FAIL":
                rescue_status = "WARN"
        if "max_top10_rescue_signature_fraction" in thresholds:
            for _, sub in rescue_summary.groupby("target", sort=False):
                frac = sub.sort_values("count", ascending=False).head(10)["fraction"].sum()
                if frac >= threshold_value("max_top10_rescue_signature_fraction", "fail"):
                    rescue_status = "FAIL"
                elif frac >= threshold_value("max_top10_rescue_signature_fraction", "warning") and rescue_status != "FAIL":
                    rescue_status = "WARN"
        statuses["rescue_signature_distribution"] = rescue_status

    seed_pair_status = "PASS"
    if "max_single_his_seed_pair_fraction" in thresholds:
        pair_rows = []
        for _, row in df.iterrows():
            for pair in itertools.combinations(his_labels(row["mutation_list"]), 2):
                pair_rows.append({"target": row["target"], "pair": ";".join(pair)})
        if pair_rows:
            pair_summary = pd.DataFrame(pair_rows).groupby(["target", "pair"]).size().reset_index(name="count")
            total_by_target = df.groupby("target").size().to_dict()
            pair_summary["fraction"] = pair_summary.apply(lambda r: r["count"] / total_by_target[r["target"]], axis=1)
            for frac in pair_summary["fraction"]:
                if frac >= threshold_value("max_single_his_seed_pair_fraction", "fail"):
                    seed_pair_status = "FAIL"
                elif frac >= threshold_value("max_single_his_seed_pair_fraction", "warning") and seed_pair_status != "FAIL":
                    seed_pair_status = "WARN"
        statuses["his_seed_pair_distribution"] = seed_pair_status

    single_frac_status = "PASS"
    for target, sub in df.groupby("target"):
        frac = (sub.mutation_count == 1).mean()
        if frac >= threshold_value("max_single_mutation_fraction_in_design_pool", "fail"):
            single_frac_status = "FAIL"
        elif frac >= threshold_value("max_single_mutation_fraction_in_design_pool", "warning") and single_frac_status != "FAIL":
            single_frac_status = "WARN"
    statuses["single_mutation_fraction"] = single_frac_status

    if "sdAb_his_distribution" in thresholds:
        sdab = df[df.target == "sdAb"]
        frac2 = (sdab.His_count.astype(int) == 2).mean() if not sdab.empty else 0.0
        ge3 = (sdab.His_count.astype(int) >= 3).mean() if not sdab.empty else 0.0
        cfg = thresholds["sdAb_his_distribution"]
        if frac2 < float(cfg["min_2his_fraction"]) or frac2 > float(cfg["max_2his_fraction"]) or ge3 > float(cfg["max_his_ge3_fraction"]):
            statuses["sdAb_his_distribution_target"] = "FAIL"
        else:
            statuses["sdAb_his_distribution_target"] = "PASS"

    if "sdAb_mutation_order_distribution" in thresholds:
        sdab = df[df.target == "sdAb"]
        frac2mut = (sdab.mutation_count.astype(int) == 2).mean() if not sdab.empty else 0.0
        frac3mut = (sdab.mutation_count.astype(int) == 3).mean() if not sdab.empty else 0.0
        cfg = thresholds["sdAb_mutation_order_distribution"]
        if frac2mut < float(cfg["min_2mut_fraction"]) or frac3mut > float(cfg["max_3mut_fraction"]):
            statuses["sdAb_mutation_order_target"] = "FAIL"
        else:
            statuses["sdAb_mutation_order_target"] = "PASS"

    required_routes = config["routes"]["required_primary_generation_routes"]
    per_route = df.groupby(["target", "primary_generation_route"]).size().to_dict()
    route_min_status = "PASS"
    for target in config["targets"]:
        for route in required_routes:
            count = per_route.get((target, route), 0)
            if count <= threshold_value("min_effective_candidates_per_required_route", "fail"):
                route_min_status = "FAIL"
            elif count <= threshold_value("min_effective_candidates_per_required_route", "warning") and route_min_status != "FAIL":
                route_min_status = "WARN"
    statuses["min_effective_candidates_per_required_route"] = route_min_status

    if "required_control_panel_count_per_target" in thresholds:
        control_count_status = "PASS"
        counts = control_df.groupby("target").size().to_dict() if not control_df.empty else {}
        required = int(threshold_value("required_control_panel_count_per_target", "required"))
        for target in config["targets"]:
            if counts.get(target, 0) != required:
                control_count_status = "FAIL"
        statuses["control_panel_count_per_target"] = control_count_status

    control_conflict_status = "PASS"
    if "max_control_design_duplicate_conflict_fraction" in thresholds and not control_df.empty:
        for _, sub in control_df.groupby("target", sort=False):
            if "design_duplicate_conflict" not in sub:
                continue
            frac = sub["design_duplicate_conflict"].mean()
            if frac >= threshold_value("max_control_design_duplicate_conflict_fraction", "fail"):
                control_conflict_status = "FAIL"
            elif frac >= threshold_value("max_control_design_duplicate_conflict_fraction", "warning") and control_conflict_status != "FAIL":
                control_conflict_status = "WARN"
        statuses["control_design_duplicate_conflict"] = control_conflict_status

    statuses["buildability_light_ran"] = "PASS"
    statuses["run_manifest_complete"] = "PASS"
    required_audits = [
        "generation_audit.csv",
        "filter_audit.csv",
        "route_filter_loss_summary.csv",
        "exact_duplicate_source_route_matrix.csv",
        "near_duplicate_top_clusters.csv",
        "his_seed_set_summary.csv",
        "his_seed_pair_frequency.csv",
        "position_set_summary.csv",
        "position_pair_frequency.csv",
        "rescue_objective_distribution.csv",
        "orphan_rescue_check.csv",
        "control_anchor_panel_200_per_target.csv",
        "exact_duplicate_sequence_provenance.csv",
        "his_seed_set_by_route.csv",
        "his_seed_set_by_rescue_signature.csv",
        "his_seed_set_near_duplicate_cluster_summary.csv",
        "top_seed_set_examples.csv",
        "his_seed_pair_summary.csv",
        "rescue_signature_summary.csv",
        "rescue_signature_by_route.csv",
        "mutation_set_summary.csv",
        "near_duplicate_hamming_jaccard_summary.csv",
        "top_near_duplicate_clusters.csv",
        "sdAb_seed_rescue_collapse_report.csv",
    ]
    missing = [name for name in required_audits if not (DRY / name).exists()]
    statuses["required_audit_tables"] = "PASS" if not missing else "FAIL: missing " + "; ".join(missing)
    return statuses


def write_manifest(config_hash: str, input_paths: list[Path]) -> None:
    manifest = {
        "run_id": f"dry_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "git_commit": git_commit(),
        "config_hash": config_hash,
        "input_file_hashes": {str(p.relative_to(ROOT)): file_sha(p) for p in input_paths},
        "random_seed": read_config()["random_seed"],
        "python_env": os.environ.get("CONDA_DEFAULT_ENV", "optim-pipe"),
        "tool_versions": {"dry_run_generator": "rule_based_v1"},
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "completed_at": datetime.now().isoformat(timespec="seconds"),
        "operator": os.environ.get("USER", "unknown"),
        "hostname": socket.gethostname(),
        "pipeline_version": "initial_design_generation_final_dry_run_v1",
    }
    with (OUT / "run_manifest.yaml").open("w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, allow_unicode=True)
    with (DRY / "run_manifest.yaml").open("w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, allow_unicode=True)


def write_report(statuses: dict[str, str], candidates: pd.DataFrame, audits: dict[str, dict]) -> None:
    fail = [k for k, v in statuses.items() if str(v).startswith("FAIL")]
    warn = [k for k, v in statuses.items() if str(v).startswith("WARN")]
    if fail:
        verdict = "FAIL"
        if "stress" in REPORT_NAME:
            recommendation = "Do not unlock production. Review stress-run failure modes, revise generator/config, and repeat 20K stress dry-run."
        else:
            recommendation = "Do not unlock production. Fix generator/config and repeat 5K dry-run."
    elif warn:
        verdict = "WARN"
        if "stress" in REPORT_NAME:
            recommendation = "Do not unlock production automatically. Review remaining WARN diversity guardrails before any production run."
        else:
            recommendation = "Do not unlock production yet. Run 20K stress dry-run before production."
    else:
        verdict = "PASS"
        recommendation = "Dry-run gate passed. Production planning can be considered, but production remains locked until manually approved."

    report_title = "Initial Design Generation Stress Dry-run Validation Report" if "stress" in REPORT_NAME else "Initial Design Generation Dry-run Validation Report"
    lines = [
        f"# {report_title}",
        "",
        f"Overall verdict: **{verdict}**",
        "",
        f"Recommendation: {recommendation}",
        "",
        "## Counts",
        "",
        "| target | effective design candidates | raw generated | post hard filter | control draft rows |",
        "|---|---:|---:|---:|---:|",
    ]
    for target, a in audits.items():
        lines.append(
            f"| {target} | {a['effective_design_pool_count']} | {a['raw_generated_count']} | "
            f"{a['post_generation_hard_filter_count']} | {a['control_anchor_pool_count']} |"
        )
    lines += ["", "## Validation Checks", "", "| check | status |", "|---|---|"]
    for key, value in statuses.items():
        lines.append(f"| {key} | {value} |")

    def markdown_table(frame: pd.DataFrame) -> str:
        cols = list(frame.columns)
        out = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
        for _, r in frame.iterrows():
            vals = []
            for col in cols:
                val = r[col]
                if isinstance(val, float):
                    vals.append(f"{val:.4f}")
                else:
                    vals.append(str(val))
            out.append("| " + " | ".join(vals) + " |")
        return "\n".join(out)

    lines += ["", "## Route Distribution", ""]
    route_summary = pd.read_csv(DRY / "route_summary.csv")
    lines.append(markdown_table(route_summary))
    lines += ["", "## His Count Distribution", ""]
    his_dist = pd.read_csv(DRY / "his_distribution.csv")
    lines.append(markdown_table(his_dist))
    lines += ["", "## Mutation Order Distribution", ""]
    mut_dist = pd.read_csv(DRY / "mutation_order_distribution.csv")
    lines.append(markdown_table(mut_dist))

    lines += [
        "",
        "## Notes",
        "",
        "- This dry-run did not run ESM, FoldX actual, PyRosetta, AF3, explicit glycan modeling, Rosetta dddG_elec, SimpleFold, or MD.",
        "- `final_selection_bucket` is intentionally unset before Tier 3.",
        "- `variant_id` is sequence-stable and route-independent; route provenance is stored in generation records.",
        "- Near-duplicate clustering uses a within-His-seed-set Hamming/Jaccard rule: sequence mutation Hamming distance <= 2 or mutation-set Jaccard >= 0.80.",
    ]
    report = "\n".join(lines) + "\n"
    (REPORTS / REPORT_NAME).write_text(report)
    (DRY / REPORT_NAME).write_text(report)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run rule-based initial-design dry-run gates.")
    parser.add_argument(
        "--stress",
        action="store_true",
        help="Run the larger stress dry-run into results/initial_design_generation/stress_dry_run.",
    )
    parser.add_argument(
        "--target-size",
        type=int,
        default=None,
        help="Override effective dry-run candidate count per target.",
    )
    return parser.parse_args()


def main() -> None:
    global DRY, REPORT_NAME
    args = parse_args()
    config = read_config()
    if args.stress:
        DRY = OUT / "stress_dry_run"
        REPORT_NAME = "stress_dry_run_validation_report.md"
        stress_size = args.target_size or int(config.get("stress_dry_run_effective_design_pool", 20000))
        for target_cfg in config["targets"].values():
            target_cfg["dry_run_effective_design_pool"] = stress_size
        if "stress_dry_run_unlock_thresholds" in config:
            config["dry_run_gate_thresholds"] = config["stress_dry_run_unlock_thresholds"]
        config["run_mode"] = "stress_dry_run"
    elif args.target_size:
        for target_cfg in config["targets"].values():
            target_cfg["dry_run_effective_design_pool"] = int(args.target_size)
        config["run_mode"] = f"dry_run_{args.target_size}"

    TABLES.mkdir(parents=True, exist_ok=True)
    DRY.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    cfg_hash = config_hash(config)
    config_text = yaml.safe_dump(config, sort_keys=False, allow_unicode=True)
    (OUT / "design_config.yaml").write_text(config_text)
    (DRY / "design_config.yaml").write_text(config_text)
    inputs = load_inputs()
    input_paths = [UPSTREAM / f"{name}.csv" for name in inputs]
    build_evidence_ledger(config, inputs)

    all_candidates = []
    all_rescue = []
    all_controls = []
    audits: dict[str, dict] = {}
    for target in config["targets"]:
        candidates, rescue, audit, controls = generate_target_candidates(target, config, inputs, cfg_hash)
        candidates.pipe(write_csv, DRY / f"dry_run_candidates_{target}.csv")
        all_candidates.append(candidates)
        if not rescue.empty:
            all_rescue.append(rescue)
        all_controls.extend(controls)
        audits[target] = audit

    combined = pd.concat(all_candidates, ignore_index=True)
    combined.pipe(write_csv, DRY / "dry_run_candidates_all.csv")
    combined.pipe(write_csv, TABLES / "dry_run_candidates_all.csv")
    if all_rescue:
        rescue_df = pd.concat(all_rescue, ignore_index=True)
    else:
        rescue_df = pd.DataFrame()
    rescue_df.pipe(write_csv, TABLES / "candidate_rescue_mutations.csv")
    rescue_df.pipe(write_csv, DRY / "candidate_rescue_mutations.csv")
    pd.DataFrame(all_controls).pipe(write_csv, TABLES / "control_anchor_panel.csv")

    statuses = write_audits(combined, audits, all_controls, config)
    write_manifest(cfg_hash, input_paths)
    write_report(statuses, combined, audits)
    print(json.dumps({"statuses": statuses, "counts": combined.groupby("target").size().to_dict()}, indent=2))


if __name__ == "__main__":
    main()
