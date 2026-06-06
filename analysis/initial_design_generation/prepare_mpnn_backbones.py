#!/usr/bin/env python3
"""Prepare reviewed AF3 backbones for ProteinMPNN P0 runner inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import MMCIFParser, PDBIO, PDBParser

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT
CONFIG = dry.CONFIG
P0_DIR = ROOT / "results/initial_design_generation/p0_mpnn"
BACKBONE_DIR = ROOT / "results/initial_design_generation/p0_mpnn_backbones"
AF3_MANIFEST = ROOT / "results/ph_sensitive_40aa_window/tables/af3_model_manifest.csv"
REFERENCE = ROOT / "results/ph_sensitive_40aa_window/tables/reference_sequence_map.csv"
TARGET_SAFE = {"1E62": "Ab_1E62", "sdAb": "Ab_sdAb"}
TARGET_REVERSE = {"Ab_1E62": "1E62", "Ab_sdAb": "sdAb"}


def normalize_target(value: str) -> str:
    return TARGET_REVERSE.get(str(value), str(value))


def safe_target(value: str) -> str:
    value = normalize_target(value)
    return TARGET_SAFE.get(value, value)


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_config() -> dict:
    with CONFIG.open() as fh:
        return yaml.safe_load(fh)


def structure_from_path(path: Path):
    suffix = path.suffix.lower()
    if suffix == ".cif":
        return MMCIFParser(QUIET=True).get_structure(path.stem, path)
    if suffix == ".pdb":
        return PDBParser(QUIET=True).get_structure(path.stem, path)
    raise ValueError(f"Unsupported structure format: {path}")


def chain_sequences(structure) -> dict[str, str]:
    model = next(structure.get_models())
    seqs: dict[str, str] = {}
    for chain in model:
        seq = "".join(
            protein_letters_3to1.get(res.get_resname().upper(), "X")
            for res in chain
            if res.id[0] == " "
        )
        if seq:
            seqs[chain.id] = seq
    return seqs


def reference_sequence(reference: pd.DataFrame, target: str, chain: str) -> str:
    sub = reference[
        (reference["target"].map(normalize_target) == target)
        & (reference["chain"] == chain)
    ].sort_values("local_pos")
    if sub.empty:
        raise ValueError(f"No reference sequence for {target}/{chain}")
    return "".join(sub["aa"].astype(str).tolist())


def convert_to_pdb(source: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    structure = structure_from_path(source)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output))


def pick_af3_models(top_n: int) -> pd.DataFrame:
    manifest = pd.read_csv(AF3_MANIFEST)
    manifest["target_norm"] = manifest["target"].map(normalize_target)
    manifest = manifest[
        (manifest["quality_status"] == "pass")
        & (~manifest["blocker_flag"].astype(str).str.lower().isin(["true", "1"]))
    ].copy()
    rows = []
    for target in ["1E62", "sdAb"]:
        sub = manifest[manifest["target_norm"] == target].sort_values(
            ["ranking_score", "seed", "model_index"],
            ascending=[False, True, True],
        )
        if len(sub) < top_n:
            raise ValueError(f"Not enough pass AF3 models for {target}: requested {top_n}, found {len(sub)}")
        rows.append(sub.head(top_n))
    return pd.concat(rows, ignore_index=True)


def build_manifest(top_n: int, output_manifest: Path, output_csv: Path) -> tuple[dict, list[dict]]:
    config = read_config()
    reference = pd.read_csv(REFERENCE)
    selected = pick_af3_models(top_n)
    entries = []
    audit_rows = []
    for _, row in selected.iterrows():
        target = normalize_target(row["target"])
        target_cfg = config["targets"][target]
        design_chain = target_cfg["chain_id"]
        source = ROOT / row["structure_path"]
        if not source.exists():
            raise FileNotFoundError(source)
        structure = structure_from_path(source)
        seqs = chain_sequences(structure)
        if design_chain not in seqs:
            raise ValueError(f"{source} missing design chain {design_chain} for {target}")
        expected = reference_sequence(reference, target, design_chain)
        seq_match = seqs[design_chain] == expected
        if not seq_match:
            raise ValueError(
                f"{target}/{design_chain} sequence mismatch for {source}: "
                f"structure_len={len(seqs[design_chain])}, reference_len={len(expected)}"
            )
        visible = [chain for chain in sorted(seqs) if chain != design_chain]
        rank_within_target = (
            selected[selected["target_norm"] == target]
            .sort_values("ranking_score", ascending=False)
            .reset_index(drop=True)
        )
        target_rank = int(rank_within_target.index[rank_within_target["structure_path"] == row["structure_path"]][0] + 1)
        backbone_id = f"af3_{target}_seed{int(row['seed'])}_model{int(row['model_index'])}_rank{target_rank}"
        output_pdb = BACKBONE_DIR / f"{backbone_id}.pdb"
        convert_to_pdb(source, output_pdb)
        entries.append(
            {
                "target": target,
                "backbone_id": backbone_id,
                "structure_path": display_path(output_pdb),
                "source_structure_path": row["structure_path"],
                "source_format": "cif",
                "converted_format": "pdb",
                "project_chain_id": design_chain,
                "mpnn_design_chain_id": design_chain,
                "visible_chain_ids": visible,
                "local_to_mpnn_index_mode": "chain_order",
                "ranking_score": float(row["ranking_score"]),
                "quality_status": row["quality_status"],
                "seed": int(row["seed"]),
                "model_index": int(row["model_index"]),
                "sequence_match_reference": True,
                "chain_lengths": {chain: len(seq) for chain, seq in seqs.items()},
            }
        )
        audit_rows.append(
            {
                "target": safe_target(target),
                "backbone_id": backbone_id,
                "source_structure_path": row["structure_path"],
                "converted_pdb_path": display_path(output_pdb),
                "mpnn_design_chain_id": design_chain,
                "visible_chain_ids": ";".join(visible),
                "ranking_score": float(row["ranking_score"]),
                "quality_status": row["quality_status"],
                "sequence_match_reference": True,
                "chain_lengths": ";".join(f"{chain}:{len(seq)}" for chain, seq in seqs.items()),
            }
        )
    manifest = {
        "version": 1,
        "status": "ready_for_runner_input_generation",
        "source": "af3_model_manifest_top_pass_by_ranking_score",
        "backbones": entries,
    }
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    with output_manifest.open("w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, allow_unicode=True)
    dry.write_csv(pd.DataFrame(audit_rows), output_csv)
    return manifest, audit_rows


def write_report(audit_rows: list[dict], report_path: Path) -> None:
    lines = [
        "# MPNN Backbone Preparation Report",
        "",
        "Status: `ready_for_runner_input_generation`",
        "",
        "AF3 top-pass CIF models were converted to PDB and checked against the reference antibody sequence map.",
        "",
        "| target | backbone_id | ranking_score | design_chain | visible_chains | sequence_match |",
        "|---|---|---:|---|---|---|",
    ]
    for row in audit_rows:
        lines.append(
            f"| {row['target']} | {row['backbone_id']} | {row['ranking_score']:.6f} | "
            f"{row['mpnn_design_chain_id']} | {row['visible_chain_ids']} | {row['sequence_match_reference']} |"
        )
    report_path.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare AF3 backbones for P0 ProteinMPNN runner inputs.")
    parser.add_argument("--top-n", type=int, default=1, help="Number of pass AF3 models per target.")
    parser.add_argument(
        "--manifest-out",
        default=str((P0_DIR / "mpnn_backbone_manifest.yaml").relative_to(ROOT)),
    )
    parser.add_argument(
        "--audit-csv-out",
        default=str((P0_DIR / "mpnn_backbone_manifest_audit.csv").relative_to(ROOT)),
    )
    parser.add_argument(
        "--report-out",
        default=str((P0_DIR / "mpnn_backbone_preparation_report.md").relative_to(ROOT)),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest_path = ROOT / args.manifest_out
    audit_path = ROOT / args.audit_csv_out
    report_path = ROOT / args.report_out
    manifest, audit_rows = build_manifest(args.top_n, manifest_path, audit_path)
    write_report(audit_rows, report_path)
    print(
        {
            "status": manifest["status"],
            "manifest": display_path(manifest_path),
            "audit_csv": display_path(audit_path),
            "report": display_path(report_path),
            "backbone_count": len(audit_rows),
        }
    )


if __name__ == "__main__":
    main()
