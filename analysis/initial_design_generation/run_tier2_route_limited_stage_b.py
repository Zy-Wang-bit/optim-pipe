#!/usr/bin/env python3
"""Execute and collect Tier2-heavy-lite Stage B route-limited review.

This module keeps the route-limited boundary explicit:
- electrostatics / pKa integrated review is computed from existing local
  heavy-lite PDB and pKa evidence;
- SimpleFold and AF3Score inputs are prepared from the audited panel;
- heavy tool result collectors summarize only fixed audited rows.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser


ROOT = Path(__file__).resolve().parents[2]
ROUTE_DIR = ROOT / "results/initial_design_generation/tier2_route_limited_stage_b"
HEAVY_LITE_DIR = ROOT / "results/initial_design_generation/tier2_heavy_lite"
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"

PANEL = ROUTE_DIR / "tier2_route_limited_panel.csv"
ELECTROSTATICS_LIST = ROUTE_DIR / "electrostatics_pka_review_list.csv"
SIMPLEFOLD_LIST = ROUTE_DIR / "simplefold_candidate_list.csv"
AF3_LIST = ROUTE_DIR / "af3_complex_candidate_list.csv"

OUT_DIR = ROUTE_DIR / "compute"
SIMPLEFOLD_DIR = OUT_DIR / "simplefold"
AF3SCORE_DIR = OUT_DIR / "af3score"

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
}

TARGET_CHAINS = {
    "Ab_1E62": {"antibody": ["H", "L"], "antigen": "C"},
    "Ab_sdAb": {"antibody": ["A"], "antigen": "B"},
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    view = df.copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
    lines = [
        "| " + " | ".join(map(str, view.columns)) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in view.columns) + " |")
    return "\n".join(lines)


def chain_sequence(structure, chain_id: str) -> str:
    if chain_id not in structure[0]:
        raise ValueError(f"chain {chain_id} missing")
    seq = []
    for residue in structure[0][chain_id]:
        if residue.id[0] == " ":
            seq.append(AA3_TO_1.get(residue.get_resname().upper(), "X"))
    return "".join(seq)


def parse_mutation_positions(mutation_list: str) -> list[tuple[str, int, str]]:
    if not isinstance(mutation_list, str) or not mutation_list:
        return []
    out = []
    for mut in mutation_list.split(";"):
        mut = mut.strip()
        match = re.match(r"^([A-Z])([A-Z])?(\d+)([A-Z])$", mut)
        if not match:
            continue
        chain = match.group(1)
        pos = int(match.group(3))
        new_aa = match.group(4)
        out.append((chain, pos, new_aa))
    return out


def atom_coords_for_residue(structure, chain_id: str, pos: int) -> np.ndarray:
    if chain_id not in structure[0]:
        return np.empty((0, 3))
    for residue in structure[0][chain_id]:
        if residue.id[0] == " " and residue.id[1] == pos:
            return np.array([atom.coord for atom in residue.get_atoms()], dtype=float)
    return np.empty((0, 3))


def chain_atom_coords(structure, chain_id: str) -> np.ndarray:
    if chain_id not in structure[0]:
        return np.empty((0, 3))
    coords = []
    for residue in structure[0][chain_id]:
        if residue.id[0] == " ":
            for atom in residue.get_atoms():
                coords.append(atom.coord)
    if not coords:
        return np.empty((0, 3))
    return np.array(coords, dtype=float)


def min_distance(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0 or len(b) == 0:
        return math.nan
    diff = a[:, None, :] - b[None, :, :]
    return float(np.sqrt(np.sum(diff * diff, axis=2)).min())


def integrated_class(row: pd.Series) -> tuple[str, str]:
    pka_score = float(row.get("his_pka_support_t2_score", math.nan))
    validity = float(row.get("pyrosetta_local_structure_validity_t2_score", math.nan))
    delta = float(row.get("pyrosetta_rosetta_delta_score", math.nan))
    his_dist = float(row.get("his_min_antigen_distance", math.nan))
    if pd.isna(pka_score):
        return "not_applicable", "no_His_or_missing_pKa"
    if validity < 0.65 or delta > 120:
        return "pH_mechanism_boundary", "local_structure_limits_interpretation"
    if pka_score >= 0.65 and (pd.isna(his_dist) or his_dist <= 8.0):
        return "pH_mechanism_strong", "high_pKa_support_and_interface_or_near_interface_His"
    if pka_score >= 0.45 and (pd.isna(his_dist) or his_dist <= 12.0):
        return "pH_mechanism_plausible", "moderate_pKa_support_and_interpretable_His_position"
    if pka_score >= 0.25:
        return "pH_mechanism_boundary", "weak_to_moderate_pKa_support"
    return "pH_mechanism_weak", "low_pKa_support"


def run_electrostatics() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(ELECTROSTATICS_LIST)
    parser = PDBParser(QUIET=True)
    rows = []
    for _, row in df.iterrows():
        pdb_path = ROOT / str(row["pyrosetta_pdb_path"])
        target = row["target"]
        chain_cfg = TARGET_CHAINS[target]
        status = "success"
        notes = []
        his_positions = []
        his_coords = []
        his_min_antigen_distance = math.nan
        try:
            structure = parser.get_structure(row["variant_id"], pdb_path)
            antigen_coords = chain_atom_coords(structure, chain_cfg["antigen"])
            for chain_id, pos, new_aa in parse_mutation_positions(str(row.get("mutation_list", ""))):
                if new_aa != "H":
                    continue
                his_positions.append(f"{chain_id}{pos}H")
                coords = atom_coords_for_residue(structure, chain_id, pos)
                if len(coords):
                    his_coords.append(coords)
            if his_coords:
                his_all = np.concatenate(his_coords, axis=0)
                his_min_antigen_distance = min_distance(his_all, antigen_coords)
        except Exception as exc:  # pragma: no cover - report path
            status = "parse_error"
            notes.append(str(exc))

        result = row.to_dict()
        result["electrostatics_status"] = status
        result["his_positions"] = ";".join(his_positions)
        result["his_min_antigen_distance"] = his_min_antigen_distance
        result["his_pka_values"] = row.get("pka_pka_combined_mean", math.nan)
        result["his_pka_support_score"] = row.get("his_pka_support_t2_score", math.nan)
        result["local_charge_environment"] = "not_explicitly_modeled"
        result["interface_charge_context"] = (
            "near_interface" if pd.notna(his_min_antigen_distance) and his_min_antigen_distance <= 8.0
            else "peripheral_or_unresolved"
        )
        cls, cls_note = integrated_class(pd.Series(result))
        result["pka_review_class"] = cls
        result["electrostatic_mechanism_class"] = cls
        result["mechanism_notes"] = ";".join([cls_note] + notes)
        rows.append(result)

    out = pd.DataFrame(rows)
    out_path = OUT_DIR / "electrostatics_pka_integrated_review.csv"
    out.to_csv(out_path, index=False)

    summary = (
        out.groupby(["target", "electrostatic_mechanism_class"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
    )
    seed_summary = (
        out.groupby(["target", "his_seed_set", "electrostatic_mechanism_class"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "his_seed_set", "count"], ascending=[True, True, False])
    )
    report = [
        "# Electrostatics / pKa Integrated Review",
        "",
        "This review integrates the existing local heavy-lite pKa support with local structure validity and approximate designed-His antigen proximity from the PyRosetta complex PDBs.",
        "",
        "## Class Summary",
        "",
        md_table(summary),
        "",
        "## Seed Summary",
        "",
        md_table(seed_summary.groupby(["target", "his_seed_set"], group_keys=False).head(5)),
        "",
        "## Boundary",
        "",
        "This is still a local/integrated interpretation layer. It is not a complex-level AF3 result and does not prove final pH-sensitive behavior.",
        "",
    ]
    (OUT_DIR / "electrostatics_pka_integrated_review.md").write_text("\n".join(report), encoding="utf-8")
    print(f"wrote {out_path} rows={len(out)}")


def prepare_simplefold() -> None:
    parser = PDBParser(QUIET=True)
    df = pd.read_csv(SIMPLEFOLD_LIST)
    fasta_root = SIMPLEFOLD_DIR / "fasta_input"
    if fasta_root.exists():
        shutil.rmtree(fasta_root)
    fasta_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, row in df.iterrows():
        pdb_path = ROOT / str(row["pyrosetta_pdb_path"])
        structure = parser.get_structure(row["variant_id"], pdb_path)
        chains = TARGET_CHAINS[row["target"]]["antibody"]
        seq = "/".join(chain_sequence(structure, chain_id) for chain_id in chains)
        fasta_path = fasta_root / f"{row['variant_id']}.fasta"
        fasta_path.write_text(f">{row['variant_id']}\n{seq}\n", encoding="utf-8")
        rows.append(
            {
                "variant_id": row["variant_id"],
                "target": row["target"],
                "sequence": seq,
                "fasta_path": fasta_path.as_posix(),
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(SIMPLEFOLD_DIR / "simplefold_input_manifest.csv", index=False)
    command = [
        "# SimpleFold route-limited command",
        "",
        "```bash",
        "CUDA_VISIBLE_DEVICES=1 /data/ziyang/mamba/envs/simplefold/bin/simplefold \\",
        "  --simplefold_model simplefold_3B \\",
        "  --ckpt_dir third_party/ml-simplefold/weights \\",
        "  --num_steps 500 \\",
        "  --tau 0.01 \\",
        "  --nsample_per_protein 3 \\",
        "  --plddt \\",
        f"  --fasta_path {fasta_root.as_posix()} \\",
        f"  --output_dir {(SIMPLEFOLD_DIR / 'outputs').as_posix()} \\",
        "  --output_format pdb \\",
        "  --backend torch",
        "```",
        "",
    ]
    (SIMPLEFOLD_DIR / "simplefold_command_plan.md").write_text("\n".join(command), encoding="utf-8")
    print(f"prepared SimpleFold FASTA rows={len(manifest)}")


def sanitize_variant_id(variant_id: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", variant_id).lower()
    return safe


def prepare_af3score() -> None:
    df = pd.read_csv(AF3_LIST)
    input_dir = AF3SCORE_DIR / "input_pdbs"
    if input_dir.exists():
        shutil.rmtree(input_dir)
    input_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, row in df.iterrows():
        safe = sanitize_variant_id(str(row["variant_id"]))
        src = ROOT / str(row["pyrosetta_pdb_path"])
        dst = input_dir / f"{safe}.pdb"
        shutil.copy2(src, dst)
        rows.append(
            {
                "variant_id": row["variant_id"],
                "target": row["target"],
                "af3score_id": safe,
                "input_pdb": dst.as_posix(),
                "mutation_list": row.get("mutation_list", ""),
                "his_seed_set": row.get("his_seed_set", ""),
                "near_duplicate_cluster_id": row.get("near_duplicate_cluster_id", ""),
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(AF3SCORE_DIR / "af3score_input_manifest.csv", index=False)
    command = [
        "# AF3Score route-limited command plan",
        "",
        "This route uses the AF3Score local runner on the fixed representative PDBs.",
        "",
        "Prepare JSON/batch inputs:",
        "",
        "```bash",
        "/data/ziyang/code/af3-local/env/bin/python /data/ziyang/code/af3-local/src/AF3Score/01_prepare_get_json.py \\",
        f"  --input_dir {(AF3SCORE_DIR / 'input_pdbs').as_posix()} \\",
        f"  --output_dir_cif {(AF3SCORE_DIR / 'single_chain_cif').as_posix()} \\",
        f"  --save_csv {(AF3SCORE_DIR / 'single_seq.csv').as_posix()} \\",
        f"  --output_dir_json {(AF3SCORE_DIR / 'json').as_posix()} \\",
        f"  --batch_dir {(AF3SCORE_DIR / 'af3_input_batch').as_posix()} \\",
        "  --num_jobs 7",
        "```",
        "",
        "Then create H5 files for each generated batch and run AF3Score inference with one GPU per batch.",
        "",
    ]
    (AF3SCORE_DIR / "af3score_command_plan.md").write_text("\n".join(command), encoding="utf-8")
    print(f"prepared AF3Score input PDBs rows={len(manifest)}")


def collect_simplefold() -> None:
    output_roots = [
        SIMPLEFOLD_DIR / "outputs",
        SIMPLEFOLD_DIR / "outputs_sharded",
        SIMPLEFOLD_DIR / "outputs_fast",
        SIMPLEFOLD_DIR / "outputs_tail",
        SIMPLEFOLD_DIR / "outputs_tail_single",
    ]
    manifest_path = SIMPLEFOLD_DIR / "simplefold_input_manifest.csv"
    manifest = pd.read_csv(manifest_path)
    rows = []
    for _, row in manifest.iterrows():
        variant_id = row["variant_id"]
        # SimpleFold output naming can differ by version; collect all files
        # containing the variant id so the status remains auditable.
        matches = []
        for outputs in output_roots:
            if outputs.exists():
                matches.extend(outputs.rglob(f"*{variant_id}*"))
        matches = sorted(set(matches))
        pdbs = [p for p in matches if p.suffix.lower() == ".pdb"]
        pdbs = pdbs[:3]
        plddt_values = []
        for pdb in pdbs:
            with pdb.open() as handle:
                for line in handle:
                    if not line.startswith("ATOM"):
                        continue
                    try:
                        plddt_values.append(float(line[60:66]))
                    except ValueError:
                        pass
        median_plddt = float(np.median(plddt_values)) if plddt_values else math.nan
        if not pdbs:
            fold_class = "fold_failed"
            notes = "no_output_pdb_found"
        elif len(pdbs) >= 3 and median_plddt >= 70:
            fold_class = "fold_supported"
            notes = "three_samples_collected_with_good_median_plddt"
        elif median_plddt >= 50:
            fold_class = "fold_boundary"
            notes = "outputs_collected_with_moderate_median_plddt"
        else:
            fold_class = "fold_low_confidence"
            notes = "outputs_collected_with_low_median_plddt"
        rows.append(
            {
                "variant_id": variant_id,
                "target": row["target"],
                "simplefold_status": "success" if pdbs else "missing",
                "num_samples": len(pdbs),
                "median_plddt": median_plddt,
                "pdb_files": ";".join(p.as_posix() for p in pdbs),
                "fold_review_class": fold_class,
                "fold_review_notes": notes,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "simplefold_antibody_ensemble_review.csv", index=False)
    print(f"collected SimpleFold rows={len(out)} success={out.simplefold_status.eq('success').sum()}")


def collect_af3score() -> None:
    manifest = pd.read_csv(AF3SCORE_DIR / "af3score_input_manifest.csv")
    metrics_path = AF3SCORE_DIR / "af3score_metrics.csv"
    if metrics_path.exists():
        metrics = pd.read_csv(metrics_path)
    else:
        metrics = pd.DataFrame()
    if not metrics.empty:
        id_col = "description" if "description" in metrics.columns else metrics.columns[0]
        merged = manifest.merge(metrics, left_on="af3score_id", right_on=id_col, how="left")
    else:
        merged = manifest.copy()
    merged["af3_status"] = np.where(
        metrics_path.exists() & merged.notna().any(axis=1),
        "success_or_partial",
        "missing_metrics",
    )
    merged.to_csv(OUT_DIR / "af3_complex_review_results.csv", index=False)
    print(f"collected AF3Score rows={len(merged)} metrics_exists={metrics_path.exists()}")


def write_status_report() -> None:
    files = {
        "electrostatics": OUT_DIR / "electrostatics_pka_integrated_review.csv",
        "simplefold": OUT_DIR / "simplefold_antibody_ensemble_review.csv",
        "af3": OUT_DIR / "af3_complex_review_results.csv",
    }
    lines = [
        "# Route-limited Stage B Compute Status",
        "",
        "This report tracks the fixed route-limited compute scope. Broad Tier2-heavy, MD, and final 10K remain locked.",
        "",
        "## Files",
        "",
    ]
    for name, path in files.items():
        if path.exists():
            lines.append(f"- {name}: `{path.as_posix()}` rows={len(pd.read_csv(path))} sha256={sha256_file(path)}")
        else:
            lines.append(f"- {name}: not yet collected")
    lines.append("")
    (OUT_DIR / "route_limited_stage_b_compute_status.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action",
        choices=[
            "electrostatics",
            "prepare-simplefold",
            "prepare-af3score",
            "collect-simplefold",
            "collect-af3score",
            "status",
            "prepare-all",
        ],
    )
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.action == "electrostatics":
        run_electrostatics()
        write_status_report()
    elif args.action == "prepare-simplefold":
        prepare_simplefold()
    elif args.action == "prepare-af3score":
        prepare_af3score()
    elif args.action == "collect-simplefold":
        collect_simplefold()
        write_status_report()
    elif args.action == "collect-af3score":
        collect_af3score()
        write_status_report()
    elif args.action == "status":
        write_status_report()
    elif args.action == "prepare-all":
        run_electrostatics()
        prepare_simplefold()
        prepare_af3score()
        write_status_report()


if __name__ == "__main__":
    main()
