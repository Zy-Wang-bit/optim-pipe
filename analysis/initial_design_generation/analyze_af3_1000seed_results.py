#!/usr/bin/env python3
"""Analyze completed 64x1000 AF3 validation outputs.

This script is intentionally separate from the AF3 runner.  It reads existing
AF3 outputs, selects one best sample per random seed, computes interface-aware
metrics on those seed representatives, and writes static-review artifacts.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd

try:
    from scipy.spatial import cKDTree
except Exception:  # pragma: no cover
    cKDTree = None


ROOT = Path(__file__).resolve().parents[2]
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"
OUT_DIR = ROOT / "results/initial_design_generation/af3_1000seed_validation"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"
PROGRESS = TASK_DIR / "progress.md"

CONTACT_CUTOFF = 5.0
CLASH_CUTOFF = 2.0
BEST_MODEL_SELECTION_VERSION = "best_model_score_v1_seed_representative"

BEST_MODEL_WEIGHTS = {
    "normalized_complex_confidence": 0.30,
    "interface_plausibility_score": 0.25,
    "parent_contact_retention_score": 0.20,
    "his_interpretability_score": 0.15,
    "new_clash_or_bad_contact_penalty": -0.10,
}

TARGET_CFG = {
    "Ab_1E62": {
        "design_chain": "L",
        "antigen_chain": "C",
        "design_chain_index": 1,
        "antigen_chain_index": 2,
        "static_interface_threshold": 0.50,
        "branch": "primary_1E62",
    },
    "Ab_sdAb": {
        "design_chain": "A",
        "antigen_chain": "B",
        "design_chain_index": 0,
        "antigen_chain_index": 1,
        "static_interface_threshold": 0.35,
        "branch": "secondary_sdAb",
    },
}


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    text = df.copy()
    for col in text.columns:
        text[col] = text[col].map(lambda value: "" if pd.isna(value) else str(value))
    headers = list(text.columns)
    rows = text.values.tolist()
    widths = []
    for i, header in enumerate(headers):
        widths.append(max([len(header), *[len(row[i]) for row in rows]]))
    header_line = "| " + " | ".join(header.ljust(widths[i]) for i, header in enumerate(headers)) + " |"
    sep_line = "| " + " | ".join("-" * widths[i] for i in range(len(headers))) + " |"
    body = [
        "| " + " | ".join(row[i].ljust(widths[i]) for i in range(len(headers))) + " |"
        for row in rows
    ]
    return "\n".join([header_line, sep_line, *body])


def as_float(value: Any, default: float = math.nan) -> float:
    try:
        if value is None or value == "":
            return default
        out = float(value)
    except Exception:
        return default
    return out


def clamp01(value: float) -> float:
    if not math.isfinite(value):
        return 0.0
    return max(0.0, min(1.0, value))


def safe_get_matrix(matrix: Any, i: int, j: int, default: float = math.nan) -> float:
    try:
        return as_float(matrix[i][j], default)
    except Exception:
        return default


def parse_seed_sample(path: Path) -> Tuple[Optional[int], Optional[int]]:
    match = re.search(r"seed-(\d+)_sample-(\d+)", path.as_posix())
    if not match:
        return None, None
    return int(match.group(1)), int(match.group(2))


def iter_summary_files(output_dir: Path) -> Iterable[Path]:
    for root, _dirs, files in os.walk(output_dir, followlinks=True):
        if "summary_confidences.json" in files:
            yield Path(root) / "summary_confidences.json"


def load_panel_with_parents() -> pd.DataFrame:
    panel = pd.read_csv(OUT_DIR / "validation_panel_64.csv", low_memory=False)
    parents = pd.read_csv(OUT_DIR / "validation_panel_parent_baselines.csv", low_memory=False)
    for col in sorted(set(panel.columns).union(parents.columns)):
        if col not in panel.columns:
            panel[col] = ""
        if col not in parents.columns:
            parents[col] = ""
    return pd.concat([panel, parents[panel.columns]], ignore_index=True)


def parse_mutations(mutation_list: Any) -> List[Dict[str, Any]]:
    if not isinstance(mutation_list, str):
        return []
    out = []
    for item in re.split(r"[;,]\s*", mutation_list.strip()):
        if not item:
            continue
        match = re.match(r"([A-Z])([A-Z]?)(\d+)([A-Z])$", item)
        if not match:
            match = re.match(r"([A-Z])(\d+)([A-Z])$", item)
            if not match:
                continue
            old, pos, new = match.groups()
            chain = ""
        else:
            old, chain, pos, new = match.groups()
        out.append({"old": old, "chain": chain, "pos": int(pos), "new": new, "raw": item})
    return out


def confidence_proxy_score(row: Dict[str, Any]) -> float:
    return (
        0.35 * as_float(row.get("ranking_score"), 0.0)
        + 0.35 * as_float(row.get("design_antigen_iptm"), 0.0)
        + 0.20 * as_float(row.get("interface_confidence_proxy"), 0.0)
        + 0.10 * (1.0 - as_float(row.get("has_clash"), 0.0))
    )


def confidence_row(
    target: str,
    variant_id: str,
    panel_role: str,
    shard_id: int,
    summary_path: Path,
) -> Optional[Dict[str, Any]]:
    seed, sample = parse_seed_sample(summary_path)
    if seed is None:
        return None
    try:
        data = json.loads(summary_path.read_text())
    except Exception:
        return None
    cfg = TARGET_CFG[target]
    design_i = cfg["design_chain_index"]
    antigen_i = cfg["antigen_chain_index"]
    pair_iptm = safe_get_matrix(data.get("chain_pair_iptm"), design_i, antigen_i)
    pair_pae = safe_get_matrix(data.get("chain_pair_pae_min"), design_i, antigen_i)
    pae_score = 1.0 - min(max(pair_pae, 0.0), 30.0) / 30.0 if math.isfinite(pair_pae) else 0.0
    interface_proxy = clamp01(0.65 * pair_iptm + 0.35 * pae_score)
    row = {
        "target": target,
        "variant_id": variant_id,
        "panel_role": panel_role,
        "shard_id": shard_id,
        "seed": seed,
        "sample": sample,
        "summary_path": summary_path.relative_to(ROOT).as_posix(),
        "model_path": (summary_path.parent / "model.cif").relative_to(ROOT).as_posix(),
        "ptm": as_float(data.get("ptm")),
        "iptm": as_float(data.get("iptm")),
        "ranking_score": as_float(data.get("ranking_score")),
        "fraction_disordered": as_float(data.get("fraction_disordered")),
        "has_clash": as_float(data.get("has_clash"), 0.0),
        "design_antigen_iptm": pair_iptm,
        "design_antigen_pae_min": pair_pae,
        "interface_confidence_proxy": interface_proxy,
    }
    row["confidence_proxy_score"] = confidence_proxy_score(row)
    return row


def collect_confidence_rows(manifest: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for _, shard in manifest.iterrows():
        target = str(shard["target"])
        variant_id = str(shard["variant_id"])
        panel_role = str(shard.get("panel_role", ""))
        shard_id = int(shard["shard_id"])
        output_dir = ROOT / str(shard["output_dir"])
        for summary_path in iter_summary_files(output_dir):
            row = confidence_row(target, variant_id, panel_role, shard_id, summary_path)
            if row is not None:
                rows.append(row)
    return pd.DataFrame(rows)


def one_best_sample_per_seed(raw: pd.DataFrame) -> pd.DataFrame:
    sort_cols = ["target", "variant_id", "seed", "confidence_proxy_score", "ranking_score"]
    out = raw.sort_values(sort_cols, ascending=[True, True, True, False, False])
    return out.drop_duplicates(["target", "variant_id", "seed"], keep="first").copy()


def aa_atom_rows(cif_path: Path, target: str, mutation_list: str) -> Dict[str, Any]:
    cfg = TARGET_CFG[target]
    design_chain = cfg["design_chain"]
    antigen_chain = cfg["antigen_chain"]
    design_atoms: List[Tuple[float, float, float, int, str]] = []
    antigen_atoms: List[Tuple[float, float, float, int, str]] = []
    his_atoms: List[Tuple[float, float, float]] = []
    mut_his_positions = {
        item["pos"] for item in parse_mutations(mutation_list) if item.get("new") == "H"
    }
    try:
        handle = cif_path.open()
    except Exception as exc:
        return geometry_error(f"open_failed:{exc}")
    with handle:
        for line in handle:
            if not line.startswith("ATOM "):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            element = parts[2].upper()
            atom_name = parts[3].upper()
            if element == "H" or atom_name.startswith("H"):
                continue
            resname = parts[5].upper()
            label_chain = parts[6]
            try:
                label_pos = int(parts[8])
                x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
            except Exception:
                continue
            if label_chain == design_chain:
                design_atoms.append((x, y, z, label_pos, resname))
                if (label_pos in mut_his_positions) or resname == "HIS":
                    his_atoms.append((x, y, z))
            elif label_chain == antigen_chain:
                antigen_atoms.append((x, y, z, label_pos, resname))
    if not design_atoms or not antigen_atoms:
        return geometry_error("missing_design_or_antigen_atoms")
    return geometry_from_atoms(design_atoms, antigen_atoms, his_atoms)


def geometry_error(status: str) -> Dict[str, Any]:
    return {
        "structure_parse_status": status,
        "interface_contact_count": math.nan,
        "contact_signature": "",
        "antigen_contact_signature": "",
        "his_min_antigen_distance": math.nan,
        "clash_count": math.nan,
        "interface_geometry_score": 0.0,
        "his_interpretability_score": 0.0,
        "new_clash_or_bad_contact_penalty": 1.0,
    }


def geometry_from_atoms(
    design_atoms: Sequence[Tuple[float, float, float, int, str]],
    antigen_atoms: Sequence[Tuple[float, float, float, int, str]],
    his_atoms: Sequence[Tuple[float, float, float]],
) -> Dict[str, Any]:
    if cKDTree is None:
        return geometry_error("scipy_cKDTree_unavailable")
    antigen_coords = [(x, y, z) for x, y, z, _pos, _res in antigen_atoms]
    antigen_positions = [pos for _x, _y, _z, pos, _res in antigen_atoms]
    design_coords = [(x, y, z) for x, y, z, _pos, _res in design_atoms]
    design_positions = [pos for _x, _y, _z, pos, _res in design_atoms]
    tree = cKDTree(antigen_coords)

    contact_pairs: Set[Tuple[int, int]] = set()
    clash_count = 0
    for coord, dpos in zip(design_coords, design_positions):
        for idx in tree.query_ball_point(coord, CONTACT_CUTOFF):
            contact_pairs.add((dpos, antigen_positions[int(idx)]))
        clash_count += len(tree.query_ball_point(coord, CLASH_CUTOFF))

    his_min = math.nan
    if his_atoms:
        distances, _idxs = tree.query(his_atoms, k=1)
        try:
            his_min = float(min(distances))
        except Exception:
            his_min = math.nan
    if not math.isfinite(his_min):
        his_score = 0.50 if not his_atoms else 0.0
    elif his_min <= 6.0:
        his_score = 1.0
    elif his_min <= 10.0:
        his_score = 0.75
    elif his_min <= 14.0:
        his_score = 0.35
    else:
        his_score = 0.05

    contact_list = sorted(contact_pairs)
    antigen_signature = sorted({apos for _dpos, apos in contact_list})
    return {
        "structure_parse_status": "ok",
        "interface_contact_count": len(contact_list),
        "contact_signature": ";".join(f"{d}:{a}" for d, a in contact_list[:250]),
        "antigen_contact_signature": ";".join(map(str, antigen_signature[:250])),
        "his_min_antigen_distance": his_min,
        "clash_count": clash_count,
        "interface_geometry_score": clamp01(len(contact_list) / 25.0),
        "his_interpretability_score": his_score,
        "new_clash_or_bad_contact_penalty": clamp01(clash_count / 10.0),
    }


def geometry_worker(args: Tuple[str, str, str, str]) -> Dict[str, Any]:
    target, variant_id, mutation_list, model_path = args
    out = {
        "target": target,
        "variant_id": variant_id,
        "model_path": model_path,
    }
    out.update(aa_atom_rows(ROOT / model_path, target, mutation_list))
    return out


def add_geometry(seed_best: pd.DataFrame, panel: pd.DataFrame, workers: int) -> pd.DataFrame:
    meta = {
        (str(row["target"]), str(row["variant_id"])): str(row.get("mutation_list", ""))
        for _, row in panel.iterrows()
    }
    tasks = [
        (
            str(row["target"]),
            str(row["variant_id"]),
            meta.get((str(row["target"]), str(row["variant_id"])), ""),
            str(row["model_path"]),
        )
        for _, row in seed_best.iterrows()
    ]
    rows: List[Dict[str, Any]] = []
    if workers <= 1:
        for task in tasks:
            rows.append(geometry_worker(task))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(geometry_worker, task) for task in tasks]
            done = 0
            total = len(futures)
            for future in as_completed(futures):
                rows.append(future.result())
                done += 1
                if done % 5000 == 0:
                    print(f"geometry parsed {done}/{total}", flush=True)
    geom = pd.DataFrame(rows)
    return seed_best.merge(geom, on=["target", "variant_id", "model_path"], how="left")


def contact_set(signature: Any) -> Set[str]:
    if not isinstance(signature, str) or not signature:
        return set()
    return set(item for item in signature.split(";") if item)


def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / float(len(a | b))


def assign_interface_clusters(sub: pd.DataFrame) -> pd.DataFrame:
    ordered = sub.sort_values("best_model_score", ascending=False).copy()
    clusters: List[Dict[str, Any]] = []
    assignments: Dict[int, str] = {}
    for idx, row in ordered.iterrows():
        sig = contact_set(row.get("contact_signature", ""))
        assigned = None
        for cluster in clusters:
            if jaccard(sig, cluster["signature"]) >= 0.50:
                assigned = cluster["id"]
                break
        if assigned is None:
            assigned = f"ifc_{len(clusters):04d}"
            clusters.append({"id": assigned, "signature": sig})
        assignments[idx] = assigned
    ordered["cluster_id"] = pd.Series(assignments)
    return ordered.sort_index()


def parent_contact_signatures(seed_best: pd.DataFrame) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = {}
    parents = seed_best[seed_best["variant_id"].astype(str).str.contains("parent_WT", na=False)]
    for target, sub in parents.groupby("target"):
        best = sub.sort_values("best_model_score_no_parent", ascending=False).iloc[0]
        out[str(target)] = contact_set(best.get("contact_signature", ""))
    return out


def contact_retention(row: pd.Series, parent_signatures: Dict[str, Set[str]]) -> float:
    parent = parent_signatures.get(str(row["target"]), set())
    current = contact_set(row.get("contact_signature", ""))
    if not parent:
        return 0.5
    return len(parent & current) / float(max(len(parent), 1))


def compute_scores(seed_best: pd.DataFrame) -> pd.DataFrame:
    out = seed_best.copy()
    out["interface_plausibility_score"] = (
        0.50 * out["interface_confidence_proxy"].fillna(0)
        + 0.50 * out["interface_geometry_score"].fillna(0)
    ).clip(0, 1)
    out["normalized_complex_confidence"] = (
        0.50 * out["ranking_score"].fillna(0)
        + 0.50 * out["design_antigen_iptm"].fillna(0)
    ).clip(0, 1)
    out["parent_contact_retention_score"] = 0.5
    out["best_model_score_no_parent"] = (
        BEST_MODEL_WEIGHTS["normalized_complex_confidence"] * out["normalized_complex_confidence"]
        + BEST_MODEL_WEIGHTS["interface_plausibility_score"] * out["interface_plausibility_score"]
        + BEST_MODEL_WEIGHTS["parent_contact_retention_score"] * 0.5
        + BEST_MODEL_WEIGHTS["his_interpretability_score"] * out["his_interpretability_score"].fillna(0)
        + BEST_MODEL_WEIGHTS["new_clash_or_bad_contact_penalty"]
        * out["new_clash_or_bad_contact_penalty"].fillna(1)
    )
    parent_sigs = parent_contact_signatures(out)
    out["parent_contact_retention_score"] = out.apply(lambda row: contact_retention(row, parent_sigs), axis=1)
    parents = out["variant_id"].astype(str).str.contains("parent_WT", na=False)
    out.loc[parents, "parent_contact_retention_score"] = 1.0
    out["best_model_score"] = (
        BEST_MODEL_WEIGHTS["normalized_complex_confidence"] * out["normalized_complex_confidence"]
        + BEST_MODEL_WEIGHTS["interface_plausibility_score"] * out["interface_plausibility_score"]
        + BEST_MODEL_WEIGHTS["parent_contact_retention_score"] * out["parent_contact_retention_score"]
        + BEST_MODEL_WEIGHTS["his_interpretability_score"] * out["his_interpretability_score"].fillna(0)
        + BEST_MODEL_WEIGHTS["new_clash_or_bad_contact_penalty"]
        * out["new_clash_or_bad_contact_penalty"].fillna(1)
    )
    clustered = []
    for (_target, _variant), sub in out.groupby(["target", "variant_id"], dropna=False):
        clustered.append(assign_interface_clusters(sub))
    return pd.concat(clustered, ignore_index=True)


def parent_support_table(seed_best: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    parents = seed_best[seed_best["variant_id"].astype(str).str.contains("parent_WT", na=False)]
    for target, sub in parents.groupby("target"):
        threshold = TARGET_CFG[str(target)]["static_interface_threshold"]
        frac = float(sub["interface_plausibility_score"].ge(threshold).mean())
        best = sub.sort_values("best_model_score", ascending=False).iloc[0]
        out[str(target)] = {
            "parent_interface_plausible_fraction": frac,
            "parent_best_model_score": float(best["best_model_score"]),
            "parent_supported": bool(frac >= 0.20 and float(best["best_model_score"]) >= 0.35),
            "parent_best_model_path": best["model_path"],
        }
    return out


def parent_delta_class(target: str, variant_id: str, sub: pd.DataFrame, parent_stats: Dict[str, Dict[str, Any]]) -> str:
    threshold = TARGET_CFG[target]["static_interface_threshold"]
    support = float(sub["interface_plausibility_score"].ge(threshold).mean())
    retention = float(sub["parent_contact_retention_score"].mean())
    if "parent_WT" in variant_id:
        return "parent_supported" if support >= 0.20 else "parent_weak"
    parent = parent_stats.get(target)
    if not parent:
        return "parent_baseline_missing"
    if parent["parent_supported"] and support >= 0.20 and retention >= 0.30:
        return "parent_supported_mutant_supported"
    if parent["parent_supported"] and support >= 0.20 and retention >= 0.15:
        return "parent_supported_mutant_partial_or_shifted"
    if parent["parent_supported"] and (support < 0.05 or retention < 0.15):
        return "parent_supported_mutant_weak"
    if (not parent["parent_supported"]) and support >= 0.20:
        return "parent_weak_mutant_supported_manual_review"
    return "parent_weak_mutant_weak"


def static_class(best: pd.Series, support: float, best_cluster_fraction: float, parent_delta: str) -> str:
    good = (
        as_float(best.get("best_model_score"), 0.0) >= 0.35
        and as_float(best.get("new_clash_or_bad_contact_penalty"), 1.0) < 0.8
        and str(best.get("structure_parse_status")) == "ok"
    )
    if not good or parent_delta == "parent_supported_mutant_weak":
        return "unsupported"
    if best_cluster_fraction >= 0.20 and support >= 0.20:
        return "best-supported"
    if best_cluster_fraction >= 0.05 or support >= 0.05:
        return "best-plausible"
    return "best-rare"


def md_eligibility(cls: str) -> str:
    if cls in {"best-supported", "best-plausible"}:
        return "main_md_candidate"
    if cls == "best-rare":
        return "boundary_or_manual_review"
    return "control_only_or_reject"


def build_variant_outputs(seed_best: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    parent_stats = parent_support_table(seed_best)
    variant_rows = []
    cluster_rows = []
    parent_rows = []
    class_rows = []
    for (target, variant_id), sub in seed_best.groupby(["target", "variant_id"], dropna=False):
        target = str(target)
        variant_id = str(variant_id)
        threshold = TARGET_CFG[target]["static_interface_threshold"]
        best = sub.sort_values("best_model_score", ascending=False).iloc[0]
        support = float(sub["interface_plausibility_score"].ge(threshold).mean())
        best_cluster_id = str(best["cluster_id"])
        best_cluster_size = int(sub["cluster_id"].eq(best_cluster_id).sum())
        best_cluster_fraction = best_cluster_size / float(max(len(sub), 1))
        parent_delta = parent_delta_class(target, variant_id, sub, parent_stats)
        cls = static_class(best, support, best_cluster_fraction, parent_delta)
        if "parent_WT" in variant_id:
            cls = "parent_baseline"
        sdab_interp = ""
        if target == "Ab_sdAb" and "parent_WT" not in variant_id:
            if parent_delta == "parent_supported_mutant_weak":
                sdab_interp = "sdAb_complex_worse_than_parent"
            elif parent_delta == "parent_supported_mutant_partial_or_shifted":
                sdab_interp = "sdAb_complex_partially_supported"
            elif parent_delta == "parent_weak_mutant_weak":
                sdab_interp = "sdAb_parent_also_weak"
            elif support >= 0.20:
                sdab_interp = "sdAb_complex_partially_supported"
            else:
                sdab_interp = "sdAb_fold_supported_complex_weak"
        variant_rows.append(
            {
                "target": target,
                "variant_id": variant_id,
                "num_successful_seeds": int(sub["seed"].nunique()),
                "num_successful_models_seed_best": int(len(sub)),
                "variant_completion_status": "variant_full"
                if sub["seed"].nunique() >= 1000
                else ("variant_valid" if sub["seed"].nunique() >= 950 else "variant_incomplete"),
                "best_model_path": best["model_path"],
                "best_model_seed": int(best["seed"]),
                "best_model_sample": int(best["sample"]),
                "best_model_rank": 1,
                "best_model_score": float(best["best_model_score"]),
                "best_model_reason": "max_best_model_score_v1_over_seed_best_samples",
                "best_model_selection_version": BEST_MODEL_SELECTION_VERSION,
                "best_model_for_static_review": best["model_path"],
                "best_model_for_MD_start": best["model_path"],
                "best_cluster_id": best_cluster_id,
                "best_cluster_size": best_cluster_size,
                "best_cluster_fraction": best_cluster_fraction,
                "interface_plausible_fraction": support,
                "iptm_max": float(sub["iptm"].max()),
                "iptm_top10_mean": float(sub.sort_values("best_model_score", ascending=False)["iptm"].head(10).mean()),
                "iptm_top50_mean": float(sub.sort_values("best_model_score", ascending=False)["iptm"].head(50).mean()),
                "iptm_median": float(sub["iptm"].median()),
                "design_antigen_iptm_top10_mean": float(
                    sub.sort_values("best_model_score", ascending=False)["design_antigen_iptm"].head(10).mean()
                ),
                "his_min_antigen_distance_best": best["his_min_antigen_distance"],
                "interface_contact_count_best": best["interface_contact_count"],
                "clash_count_best": best["clash_count"],
                "static_classification": cls,
                "parent_comparison_class": parent_delta,
                "sdAb_specific_interpretation": sdab_interp,
            }
        )
        for cluster_id, csub in sub.groupby("cluster_id", dropna=False):
            rep = csub.sort_values("best_model_score", ascending=False).iloc[0]
            cluster_rows.append(
                {
                    "target": target,
                    "variant_id": variant_id,
                    "cluster_id": cluster_id,
                    "cluster_size": len(csub),
                    "cluster_fraction": len(csub) / float(max(len(sub), 1)),
                    "cluster_median_iptm": csub["iptm"].median(),
                    "cluster_interface_plausible_fraction": csub["interface_plausibility_score"].ge(threshold).mean(),
                    "cluster_contact_signature": rep.get("contact_signature", ""),
                    "cluster_representative_model_path": rep["model_path"],
                }
            )
        parent_rows.append(
            {
                "target": target,
                "variant_id": variant_id,
                "parent_signature_available": target in parent_stats,
                "parent_supported": parent_stats.get(target, {}).get("parent_supported", False),
                "parent_interface_plausible_fraction": parent_stats.get(target, {}).get(
                    "parent_interface_plausible_fraction", math.nan
                ),
                "mean_parent_contact_retention": float(sub["parent_contact_retention_score"].mean()),
                "delta_interpretation": parent_delta,
            }
        )
        class_rows.append(
            {
                "target": target,
                "variant_id": variant_id,
                "static_classification": cls,
                "parent_comparison_class": parent_delta,
                "best_cluster_fraction": best_cluster_fraction,
                "interface_plausible_fraction": support,
                "md_eligibility": "parent_baseline" if "parent_WT" in variant_id else md_eligibility(cls),
                "sdAb_specific_interpretation": sdab_interp,
            }
        )
    return (
        pd.DataFrame(variant_rows),
        pd.DataFrame(cluster_rows),
        pd.DataFrame(parent_rows),
        pd.DataFrame(class_rows),
    )


def top_cohort(raw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (target, variant_id), sub in raw.groupby(["target", "variant_id"], dropna=False):
        for rank, (_, row) in enumerate(sub.sort_values("confidence_proxy_score", ascending=False).head(50).iterrows(), start=1):
            rows.append(
                {
                    "target": target,
                    "variant_id": variant_id,
                    "top_rank": rank,
                    "seed": row["seed"],
                    "sample": row["sample"],
                    "model_path": row["model_path"],
                    "summary_path": row["summary_path"],
                    "ranking_score": row["ranking_score"],
                    "iptm": row["iptm"],
                    "ptm": row["ptm"],
                    "design_antigen_iptm": row["design_antigen_iptm"],
                    "confidence_proxy_score": row["confidence_proxy_score"],
                }
            )
    return pd.DataFrame(rows)


def write_reports(
    raw: pd.DataFrame,
    seed_best: pd.DataFrame,
    variant_summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    parent_compare: pd.DataFrame,
    classification: pd.DataFrame,
) -> None:
    cls_counts = (
        classification.groupby(["target", "static_classification"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    action_counts = (
        classification.groupby(["target", "md_eligibility"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    best_view = variant_summary[
        [
            "target",
            "variant_id",
            "num_successful_seeds",
            "best_model_seed",
            "best_model_score",
            "best_cluster_fraction",
            "interface_plausible_fraction",
            "static_classification",
            "parent_comparison_class",
        ]
    ].copy()
    lines = [
        "# AF3 1000-Seed Static Review",
        "",
        "AF3 1000 seeds are complex conformational sampling only; they are not pH 7.4 / pH 6.0 predictions.",
        "",
        "This review uses one best AF3 sample per random seed for sampling-support statistics. The full 5-sample-per-seed model table is retained for top-cohort review.",
        "",
        "## Scope",
        "",
        f"- Raw AF3 sample records: {len(raw):,}",
        f"- Seed-representative records: {len(seed_best):,}",
        f"- AF3 entries classified: {len(variant_summary):,}",
        "",
        "## Classification Counts",
        "",
        md_table(cls_counts),
        "",
        "## MD Eligibility Counts",
        "",
        md_table(action_counts),
        "",
        "## Best Model Summary",
        "",
        md_table(best_view),
        "",
        "## Parent Baseline Comparison",
        "",
        md_table(parent_compare),
    ]
    (OUT_DIR / "af3_1000seed_static_review.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    sdab = classification[classification["target"].eq("Ab_sdAb")].copy()
    sdab_lines = [
        "# sdAb Specific Interpretation Report",
        "",
        "sdAb remains a secondary / exploratory branch unless parent baseline and mutant comparison clearly support upgrade.",
        "",
        "## sdAb Interpretation Counts",
        "",
        md_table(
            sdab.groupby(["static_classification", "sdAb_specific_interpretation"], dropna=False)
            .size()
            .reset_index(name="count")
        )
        if not sdab.empty
        else "_No sdAb rows._",
        "",
        "## sdAb Rows",
        "",
        md_table(sdab),
    ]
    (OUT_DIR / "sdAb_specific_interpretation_report.md").write_text("\n".join(sdab_lines) + "\n", encoding="utf-8")

    parent_lines = [
        "# Parent Baseline Comparison Report",
        "",
        "Parent / WT 1000-seed baselines are required before strong mutant conclusions.",
        "",
        md_table(parent_compare),
    ]
    (OUT_DIR / "parent_baseline_comparison_report.md").write_text("\n".join(parent_lines) + "\n", encoding="utf-8")

    report_lines = [
        "# AF3 1000-Seed Static Analysis",
        "",
        "## Executive Summary",
        "",
        "Verdict: `AF3_STATIC_ANALYSIS_COMPLETE__MD_PANEL_PLANNING_ALLOWED__MD_RUNTIME_LOCKED`.",
        "",
        "This stage analyzed completed AF3 1000-seed complex sampling for 64 mutants and 2 parent / WT baselines. It did not run new AF3, MD, glycan modeling, or final library selection.",
        "",
        "AF3 is interpreted as complex-pose sampling only, not as pH-specific prediction.",
        "",
        "## Analysis Scope",
        "",
        f"- Raw AF3 sample records analyzed: {len(raw):,}",
        f"- Seed-level representatives analyzed structurally: {len(seed_best):,}",
        f"- Classified AF3 entries: {len(variant_summary):,}",
        "",
        "## Classification Counts",
        "",
        md_table(cls_counts),
        "",
        "## MD Eligibility Counts",
        "",
        md_table(action_counts),
        "",
        "## Parent Baseline Comparison",
        "",
        md_table(parent_compare),
        "",
        "## Interpretation",
        "",
        "- `best-supported` means the selected best model is good and the same interface hypothesis has substantial seed-level support.",
        "- `best-plausible` means the best model is reasonable but sampling support is moderate.",
        "- `best-rare` means a favorable model exists but is rarely sampled; it is mechanism-generating evidence, not strong structural evidence.",
        "- `unsupported` should not enter main MD except as explicit control / audit.",
        "- sdAb remains secondary / exploratory unless parent comparison and AF3 support justify upgrade.",
        "",
        "## Output Files",
        "",
        "- `results/initial_design_generation/af3_1000seed_validation/af3_1000seed_raw_model_summary.csv`",
        "- `results/initial_design_generation/af3_1000seed_validation/af3_1000seed_seed_best_model_summary.csv`",
        "- `results/initial_design_generation/af3_1000seed_validation/best_model_selection_audit.csv`",
        "- `results/initial_design_generation/af3_1000seed_validation/af3_1000seed_top_cohort_manifest.csv`",
        "- `results/initial_design_generation/af3_1000seed_validation/af3_1000seed_cluster_summary.csv`",
        "- `results/initial_design_generation/af3_1000seed_validation/af3_1000seed_interface_summary.csv`",
        "- `results/initial_design_generation/af3_1000seed_validation/af3_1000seed_parent_baseline_comparison.csv`",
        "- `results/initial_design_generation/af3_1000seed_validation/af3_1000seed_candidate_classification.csv`",
        "- `results/initial_design_generation/af3_1000seed_validation/af3_1000seed_static_review.md`",
        "- `results/initial_design_generation/af3_1000seed_validation/sdAb_specific_interpretation_report.md`",
        "- `results/initial_design_generation/af3_1000seed_validation/parent_baseline_comparison_report.md`",
        "",
        "## Still Locked",
        "",
        "- MD panel execution and MD runtime.",
        "- Glycan-risk modeling.",
        "- Final 10K / 15K synthesis-ready selection.",
        "- Wet-lab order list generation.",
    ]
    CURRENT_STAGE_REPORT.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    progress_block = (
        "## 2026-06-13 AF3 static analysis\n\n"
        "- Completed AF3 1000-seed static result analysis for 64 mutants plus 2 parent baselines.\n"
        "- Generated seed-level best-model, interface-cluster, parent-comparison, and candidate-classification outputs.\n"
        "- MD runtime remains locked pending review and MD panel construction.\n"
    )
    if PROGRESS.exists():
        text = PROGRESS.read_text(encoding="utf-8").rstrip()
        marker = "## 2026-06-13 AF3 static analysis"
        if marker in text:
            text = text.split(marker)[0].rstrip()
        PROGRESS.write_text(text + "\n\n" + progress_block, encoding="utf-8")
    else:
        PROGRESS.write_text(progress_block, encoding="utf-8")


def run_analysis(workers: int) -> None:
    manifest = pd.read_csv(OUT_DIR / "af3_1000seed_input_manifest.csv", low_memory=False)
    panel = load_panel_with_parents()
    print("Collecting AF3 confidence summaries...", flush=True)
    raw = collect_confidence_rows(manifest)
    if raw.empty:
        raise RuntimeError("No AF3 summary_confidences.json rows were collected.")
    raw.to_csv(OUT_DIR / "af3_1000seed_raw_model_summary.csv", index=False)
    top_cohort(raw).to_csv(OUT_DIR / "af3_1000seed_top_cohort_manifest.csv", index=False)

    seed_best = one_best_sample_per_seed(raw)
    print(f"Collected raw models={len(raw)}; seed representatives={len(seed_best)}", flush=True)
    seed_status = (
        seed_best.groupby(["target", "variant_id"])["seed"]
        .nunique()
        .reset_index(name="successful_seed_count")
    )
    seed_status["seed_status"] = seed_status["successful_seed_count"].map(
        lambda x: "success" if x >= 1000 else ("valid_partial" if x >= 950 else "incomplete")
    )
    seed_status.to_csv(OUT_DIR / "af3_1000seed_seed_status.csv", index=False)

    print(f"Parsing seed-representative interface geometry with workers={workers}...", flush=True)
    seed_best = add_geometry(seed_best, panel, workers)
    seed_best = compute_scores(seed_best)
    seed_best.to_csv(OUT_DIR / "af3_1000seed_seed_best_model_summary.csv", index=False)

    variant_summary, cluster_summary, parent_compare, classification = build_variant_outputs(seed_best)
    variant_summary.to_csv(OUT_DIR / "af3_1000seed_variant_summary.csv", index=False)
    variant_summary.to_csv(OUT_DIR / "best_model_selection_audit.csv", index=False)
    cluster_summary.to_csv(OUT_DIR / "af3_1000seed_cluster_summary.csv", index=False)
    cluster_summary.to_csv(OUT_DIR / "af3_1000seed_interface_summary.csv", index=False)
    parent_compare.to_csv(OUT_DIR / "af3_1000seed_parent_baseline_comparison.csv", index=False)
    classification.to_csv(OUT_DIR / "af3_1000seed_candidate_classification.csv", index=False)
    write_reports(raw, seed_best, variant_summary, cluster_summary, parent_compare, classification)
    print("AF3 static analysis complete.", flush=True)


def write_reports_from_existing() -> None:
    raw = pd.read_csv(OUT_DIR / "af3_1000seed_raw_model_summary.csv", low_memory=False)
    seed_best = pd.read_csv(OUT_DIR / "af3_1000seed_seed_best_model_summary.csv", low_memory=False)
    variant_summary, cluster_summary, parent_compare, classification = build_variant_outputs(seed_best)
    variant_summary.to_csv(OUT_DIR / "af3_1000seed_variant_summary.csv", index=False)
    variant_summary.to_csv(OUT_DIR / "best_model_selection_audit.csv", index=False)
    cluster_summary.to_csv(OUT_DIR / "af3_1000seed_cluster_summary.csv", index=False)
    cluster_summary.to_csv(OUT_DIR / "af3_1000seed_interface_summary.csv", index=False)
    parent_compare.to_csv(OUT_DIR / "af3_1000seed_parent_baseline_comparison.csv", index=False)
    classification.to_csv(OUT_DIR / "af3_1000seed_candidate_classification.csv", index=False)
    write_reports(raw, seed_best, variant_summary, cluster_summary, parent_compare, classification)
    print("AF3 static reports regenerated from existing CSV artifacts.", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--reports-only", action="store_true")
    args = parser.parse_args()
    if args.reports_only:
        write_reports_from_existing()
    else:
        run_analysis(workers=max(1, int(args.workers)))


if __name__ == "__main__":
    main()
