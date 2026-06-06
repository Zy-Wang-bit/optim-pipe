#!/usr/bin/env python
from __future__ import annotations

import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa

if __package__ is None or __package__ == "":
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from analysis.window_selection.common import (
    af3_input_dir,
    load_inputs,
    make_arg_parser,
    output_root,
    repo_path,
    write_csv,
)


CONTACT_CUTOFF_A = 5.0
N146_CUTOFF_A = 8.0
AG_LOOP = range(100, 165)
A_DETERMINANT = range(124, 148)
TM_SEGMENTS = [(1, 28), (80, 98), (170, 226)]

MANIFEST_COLUMNS = [
    "target",
    "seed",
    "model_index",
    "json_path",
    "structure_path",
    "ranking_score",
    "confidence_summary",
    "quality_status",
    "cluster_id",
    "blocker_flag",
    "ag_loop_contact_fraction",
    "tm_contact_fraction",
    "a_determinant_contact_fraction",
    "ag_loop_geometry_status",
]

FEATURE_COLUMNS = [
    "target",
    "seed",
    "model_index",
    "structure_path",
    "chain",
    "pos",
    "aa",
    "min_antigen_distance",
    "antigen_contact_count",
    "antigen_contact_positions",
    "contact_to_ag_100_164",
    "contact_to_a_determinant_124_147",
    "distance_to_N146",
    "n146_glycan_proximity",
    "seed_support_count",
    "model_support_count",
    "cluster_id",
    "cluster_support_count",
    "model_quality_status",
    "dominant_cluster_seed_count",
    "dominant_cluster_quality_status",
    "dominant_cluster_is_anomalous",
    "dominant_cluster_tm_interface_flag",
    "dominant_cluster_ag_loop_contact_fraction",
]


def rel_path(path: Path) -> str:
    try:
        return path.relative_to(repo_path(".")).as_posix()
    except ValueError:
        return path.as_posix()


def target_chain_ids(target_cfg: dict) -> tuple[list[str], str]:
    if target_cfg["antibody_type"] == "scFv":
        return ["H", "L"], target_cfg.get("antigen_chain", "C")
    return ["A"], target_cfg.get("antigen_chain", "B")


def is_tm_pos(pos: int) -> bool:
    return any(start <= pos <= end for start, end in TM_SEGMENTS)


def residue_pos(residue) -> int | None:
    hetflag, resseq, _icode = residue.id
    if hetflag.strip():
        return None
    return int(resseq)


def residue_aa(residue) -> str:
    name = residue.get_resname().strip()
    return protein_letters_3to1_extended.get(name.title(), name if name else "UNK")


def heavy_atom_coords(residue) -> list:
    coords = []
    for atom in residue.get_atoms():
        element = (atom.element or "").upper()
        if element == "H" or atom.get_name().upper().startswith("H"):
            continue
        coords.append(atom.coord)
    return coords


def min_distance(coords_a: list, coords_b: list) -> float:
    best = math.inf
    for coord_a in coords_a:
        for coord_b in coords_b:
            dist = float(((coord_a - coord_b) ** 2).sum() ** 0.5)
            if dist < best:
                best = dist
    return best


def parse_structure(path: Path):
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def chain_map(structure) -> dict[str, Any]:
    return {chain.id: chain for chain in structure.get_chains()}


def choose_antigen_chain(chains: dict[str, Any], antibody_ids: list[str], preferred_antigen_id: str) -> str | None:
    for candidate in [preferred_antigen_id, preferred_antigen_id.upper(), "C", "B"]:
        if candidate in chains and candidate not in antibody_ids:
            return candidate
    non_antibody = [chain_id for chain_id in chains if chain_id not in antibody_ids]
    if not non_antibody:
        return None
    return max(non_antibody, key=lambda chain_id: sum(1 for residue in chains[chain_id] if is_aa(residue, standard=False)))


def model_seed_and_index(path: Path, fallback_index: int) -> tuple[int | str, int | str]:
    for part in path.parts:
        sample_match = re.fullmatch(r"seed[-_](\d+)_sample[-_](\d+)", part, flags=re.IGNORECASE)
        if sample_match:
            return int(sample_match.group(1)), int(sample_match.group(2)) + 1

    text = path.as_posix()
    seed_match = re.search(r"(?:^|[/_-])seed[_-]?(\d+)(?:$|[/_-])", text, flags=re.IGNORECASE)
    model_match = re.search(r"model[_-]?(\d+)", text, flags=re.IGNORECASE)
    ranked_match = re.search(r"ranked[_-]?(\d+)", text, flags=re.IGNORECASE)
    seed: int | str = int(seed_match.group(1)) if seed_match else ""
    if model_match:
        model_index: int | str = int(model_match.group(1))
    elif ranked_match:
        model_index = int(ranked_match.group(1)) + 1
    else:
        model_index = fallback_index
    return seed, model_index


def expected_json_path(target: str, seed: int) -> str:
    matches = sorted(af3_input_dir(target).glob(f"*seed_{seed}.json"))
    if matches:
        return rel_path(matches[0])
    return rel_path(af3_input_dir(target) / f"{target.lower()}_aes_seed_{seed}.json")


def structure_files(target: str) -> list[Path]:
    root = output_root() / "structures" / "af3" / target
    if not root.exists():
        return []
    paths = []
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
            continue
        if "input" in path.relative_to(root).parts:
            continue
        paths.append(path)
    return sorted(paths)


def confidence_files_for(structure_path: Path) -> list[Path]:
    candidates = []
    for base in [structure_path.parent, *structure_path.parents[:3]]:
        if not base.exists():
            continue
        candidates.extend(base.glob("*summary_confidences.json"))
        candidates.extend(base.glob("*confidences.json"))
    return sorted(set(candidates))


def confidence_summary(structure_path: Path) -> tuple[str, str]:
    numeric_keys = {
        "ranking_score",
        "ranking_confidence",
        "confidence_score",
        "ptm",
        "iptm",
        "fraction_disordered",
    }
    for path in confidence_files_for(structure_path):
        try:
            data = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        if not isinstance(data, dict):
            continue
        summary = {key: data[key] for key in numeric_keys if key in data}
        ranking = next((summary[key] for key in ["ranking_score", "ranking_confidence", "confidence_score", "iptm", "ptm"] if key in summary), "")
        if summary:
            return str(ranking), json.dumps(summary, sort_keys=True)
    return "", ""


def contact_cluster_id(loop_positions: set[int], tm_positions: set[int], total_contacts: int) -> str:
    if loop_positions:
        bins = [
            ("100_111", range(100, 112)),
            ("112_123", range(112, 124)),
            ("124_135", range(124, 136)),
            ("136_147", range(136, 148)),
            ("148_164", range(148, 165)),
        ]
        active = [name for name, pos_range in bins if any(pos in pos_range for pos in loop_positions)]
        return "ag_loop_" + "_".join(active)
    if tm_positions:
        return "tm_contact"
    return "other_contact" if total_contacts else "no_contact"


def quality_status(total_contacts: int, ag_loop_fraction: float, tm_fraction: float, loop_residue_count: int) -> tuple[str, str]:
    if loop_residue_count == 0:
        return "hard_fail", "missing"
    if total_contacts == 0:
        return "hard_fail", "interpretable"
    if tm_fraction > 0.5 and ag_loop_fraction < 0.25:
        return "hard_fail", "interpretable"
    if ag_loop_fraction < 0.25:
        return "soft_fail", "interpretable"
    return "pass", "interpretable"


def analyze_structure(target: str, antibody_ids: list[str], antigen_id: str, path: Path, fallback_index: int) -> tuple[dict, list[dict]]:
    seed, model_index = model_seed_and_index(path, fallback_index)
    ranking_score, conf_summary = confidence_summary(path)
    try:
        structure = parse_structure(path)
        chains = chain_map(structure)
        actual_antigen_id = choose_antigen_chain(chains, antibody_ids, antigen_id)
        actual_antibody_ids = [chain_id for chain_id in antibody_ids if chain_id in chains]
        if actual_antigen_id is None or not actual_antibody_ids:
            raise ValueError("Expected antibody or antigen chains were not found.")

        antigen_residues = []
        antigen_n146_coords = []
        loop_residue_count = 0
        for residue in chains[actual_antigen_id]:
            if not is_aa(residue, standard=False):
                continue
            pos = residue_pos(residue)
            if pos is None:
                continue
            coords = heavy_atom_coords(residue)
            if not coords:
                continue
            antigen_residues.append((pos, coords))
            if pos in AG_LOOP:
                loop_residue_count += 1
            if pos == 146:
                antigen_n146_coords.extend(coords)

        features = []
        contact_pairs: set[tuple[str, int, int]] = set()
        loop_contact_positions: set[int] = set()
        tm_contact_positions: set[int] = set()
        a_det_contact_positions: set[int] = set()

        for chain_id in actual_antibody_ids:
            for residue in chains[chain_id]:
                if not is_aa(residue, standard=False):
                    continue
                pos = residue_pos(residue)
                if pos is None:
                    continue
                coords = heavy_atom_coords(residue)
                if not coords:
                    continue
                min_ag_dist = math.inf
                contacts = []
                for ag_pos, ag_coords in antigen_residues:
                    dist = min_distance(coords, ag_coords)
                    min_ag_dist = min(min_ag_dist, dist)
                    if dist <= CONTACT_CUTOFF_A:
                        contacts.append(ag_pos)
                        contact_pairs.add((chain_id, pos, ag_pos))
                        if ag_pos in AG_LOOP:
                            loop_contact_positions.add(ag_pos)
                        if ag_pos in A_DETERMINANT:
                            a_det_contact_positions.add(ag_pos)
                        if is_tm_pos(ag_pos):
                            tm_contact_positions.add(ag_pos)

                n146_dist = min_distance(coords, antigen_n146_coords) if antigen_n146_coords else math.inf
                features.append(
                    {
                        "target": target,
                        "seed": seed,
                        "model_index": model_index,
                        "structure_path": rel_path(path),
                        "chain": chain_id,
                        "pos": pos,
                        "aa": residue_aa(residue),
                        "min_antigen_distance": "" if math.isinf(min_ag_dist) else round(min_ag_dist, 3),
                        "antigen_contact_count": len(contacts),
                        "antigen_contact_positions": "|".join(str(value) for value in sorted(set(contacts))),
                        "contact_to_ag_100_164": any(value in AG_LOOP for value in contacts),
                        "contact_to_a_determinant_124_147": any(value in A_DETERMINANT for value in contacts),
                        "distance_to_N146": "" if math.isinf(n146_dist) else round(n146_dist, 3),
                        "n146_glycan_proximity": (not math.isinf(n146_dist)) and n146_dist <= N146_CUTOFF_A,
                    }
                )

        total_contacts = len(contact_pairs)
        ag_loop_fraction = len([pair for pair in contact_pairs if pair[2] in AG_LOOP]) / total_contacts if total_contacts else 0.0
        tm_fraction = len([pair for pair in contact_pairs if is_tm_pos(pair[2])]) / total_contacts if total_contacts else 0.0
        a_det_fraction = len([pair for pair in contact_pairs if pair[2] in A_DETERMINANT]) / total_contacts if total_contacts else 0.0
        status, loop_geometry = quality_status(total_contacts, ag_loop_fraction, tm_fraction, loop_residue_count)
        cluster_id = contact_cluster_id(loop_contact_positions, tm_contact_positions, total_contacts)
        for feature in features:
            feature.update(
                {
                    "cluster_id": cluster_id,
                    "model_quality_status": status,
                }
            )
        manifest = {
            "target": target,
            "seed": seed,
            "model_index": model_index,
            "json_path": expected_json_path(target, int(seed)) if isinstance(seed, int) else "",
            "structure_path": rel_path(path),
            "ranking_score": ranking_score,
            "confidence_summary": conf_summary,
            "quality_status": status,
            "cluster_id": cluster_id,
            "blocker_flag": False,
            "ag_loop_contact_fraction": round(ag_loop_fraction, 4),
            "tm_contact_fraction": round(tm_fraction, 4),
            "a_determinant_contact_fraction": round(a_det_fraction, 4),
            "ag_loop_geometry_status": loop_geometry,
        }
        return manifest, features
    except Exception as exc:
        manifest = {
            "target": target,
            "seed": seed,
            "model_index": model_index,
            "json_path": expected_json_path(target, int(seed)) if isinstance(seed, int) else "",
            "structure_path": rel_path(path),
            "ranking_score": ranking_score,
            "confidence_summary": f"parse_error={exc}",
            "quality_status": "hard_fail",
            "cluster_id": "parse_error",
            "blocker_flag": True,
            "ag_loop_contact_fraction": 0.0,
            "tm_contact_fraction": 0.0,
            "a_determinant_contact_fraction": 0.0,
            "ag_loop_geometry_status": "missing",
        }
        return manifest, []


def blocked_expected_rows(target: str, seeds: list[int], models_per_seed: int, seen: set[tuple[int, int]]) -> list[dict]:
    rows = []
    for seed in seeds:
        for model_index in range(1, models_per_seed + 1):
            if (seed, model_index) in seen:
                continue
            rows.append(
                {
                    "target": target,
                    "seed": seed,
                    "model_index": model_index,
                    "json_path": expected_json_path(target, seed),
                    "structure_path": "",
                    "ranking_score": "",
                    "confidence_summary": "no_AF3_structure_found",
                    "quality_status": "blocked_no_af3_structure",
                    "cluster_id": "",
                    "blocker_flag": True,
                    "ag_loop_contact_fraction": "",
                    "tm_contact_fraction": "",
                    "a_determinant_contact_fraction": "",
                    "ag_loop_geometry_status": "missing",
                }
            )
    return rows


def add_support_and_dominant_cluster(features: list[dict], manifest_rows: list[dict]) -> list[dict]:
    if not features:
        return features

    model_key = ["target", "chain", "pos"]
    feature_df = pd.DataFrame(features)
    support = feature_df[feature_df["antigen_contact_count"] > 0].groupby(model_key).agg(
        seed_support_count=("seed", "nunique"),
        model_support_count=("model_index", "count"),
    )
    cluster_counts = Counter(row["cluster_id"] for row in manifest_rows if row.get("cluster_id"))
    manifest_df = pd.DataFrame(manifest_rows)
    dominant_by_target: dict[str, dict[str, Any]] = {}
    for target, target_manifest in manifest_df.groupby("target"):
        usable = target_manifest[target_manifest["quality_status"].isin(["pass", "soft_fail"])]
        if usable.empty:
            dominant_by_target[target] = {
                "dominant_cluster_seed_count": 0,
                "dominant_cluster_quality_status": "blocked",
                "dominant_cluster_is_anomalous": True,
                "dominant_cluster_tm_interface_flag": False,
                "dominant_cluster_ag_loop_contact_fraction": 0.0,
            }
            continue
        dominant_cluster = usable["cluster_id"].value_counts().idxmax()
        members = usable[usable["cluster_id"] == dominant_cluster]
        mean_tm = pd.to_numeric(members["tm_contact_fraction"], errors="coerce").fillna(0).mean()
        mean_loop = pd.to_numeric(members["ag_loop_contact_fraction"], errors="coerce").fillna(0).mean()
        anomalous = dominant_cluster in {"no_contact", "tm_contact", "parse_error"} or mean_tm > 0.5
        dominant_by_target[target] = {
            "dominant_cluster_seed_count": members["seed"].nunique(),
            "dominant_cluster_quality_status": "pass" if (members["quality_status"] == "pass").any() else "soft_fail",
            "dominant_cluster_is_anomalous": bool(anomalous),
            "dominant_cluster_tm_interface_flag": bool(mean_tm > 0.5),
            "dominant_cluster_ag_loop_contact_fraction": round(float(mean_loop), 4),
        }

    enriched = []
    for row in features:
        key = tuple(row[col] for col in model_key)
        row.update(
            {
                "seed_support_count": int(support.loc[key, "seed_support_count"]) if key in support.index else 0,
                "model_support_count": int(support.loc[key, "model_support_count"]) if key in support.index else 0,
                "cluster_support_count": cluster_counts.get(row["cluster_id"], 0),
            }
        )
        row.update(dominant_by_target.get(row["target"], {}))
        enriched.append(row)
    return enriched


def parse_af3_outputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = load_inputs()
    seeds = [int(seed) for seed in cfg.get("af3", {}).get("seeds", [])]
    models_per_seed = int(cfg.get("af3", {}).get("models_per_seed", 5))
    manifest_rows: list[dict] = []
    feature_rows: list[dict] = []

    for target_cfg in cfg["targets"]:
        target = target_cfg["name"]
        antibody_ids, antigen_id = target_chain_ids(target_cfg)
        seen: set[tuple[int, int]] = set()
        for index, path in enumerate(structure_files(target), start=1):
            manifest, features = analyze_structure(target, antibody_ids, antigen_id, path, index)
            manifest_rows.append(manifest)
            feature_rows.extend(features)
            if isinstance(manifest["seed"], int) and isinstance(manifest["model_index"], int):
                seen.add((manifest["seed"], manifest["model_index"]))
        manifest_rows.extend(blocked_expected_rows(target, seeds, models_per_seed, seen))

    feature_rows = add_support_and_dominant_cluster(feature_rows, manifest_rows)
    return (
        pd.DataFrame(manifest_rows, columns=MANIFEST_COLUMNS),
        pd.DataFrame(feature_rows, columns=FEATURE_COLUMNS),
    )


def main() -> None:
    make_arg_parser("Parse AF3 structures into manifest and residue-level structure features.").parse_args()
    manifest, features = parse_af3_outputs()
    write_csv(manifest, "af3_model_manifest.csv")
    write_csv(features, "residue_structure_features.csv")
    blocked = int(manifest["blocker_flag"].fillna(False).astype(bool).sum()) if not manifest.empty else 0
    print(f"Wrote AF3 manifest rows={len(manifest)}, residue feature rows={len(features)}, blocked_rows={blocked}.")


if __name__ == "__main__":
    main()
