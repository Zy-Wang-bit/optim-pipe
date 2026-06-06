#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

if __package__ is None or __package__ == "":
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from analysis.window_selection.common import (
    af3_input_dir,
    json_dump,
    load_inputs,
    make_arg_parser,
    read_fasta,
)


def protein_entry(chain_id: str, sequence: str) -> dict:
    return {
        "protein": {
            "id": chain_id,
            "sequence": sequence,
        }
    }


def target_chains(target_cfg: dict) -> list[tuple[str, str]]:
    antibody_type = target_cfg["antibody_type"]
    antigen_chain = target_cfg.get("antigen_chain")
    if not antigen_chain:
        antigen_chain = "C" if antibody_type == "scFv" else "B"
    if antibody_type == "scFv":
        chains = [
            ("H", read_fasta(target_cfg["heavy_fasta"])[2]),
            ("L", read_fasta(target_cfg["light_fasta"])[2]),
            (antigen_chain, read_fasta(target_cfg["antigen_fasta"], target_cfg["antigen_name"])[2]),
        ]
    elif antibody_type == "VHH":
        chains = [
            ("A", read_fasta(target_cfg["sdab_fasta"])[2]),
            (antigen_chain, read_fasta(target_cfg["antigen_fasta"])[2]),
        ]
    else:
        raise ValueError(f"Unsupported antibody_type={antibody_type!r}")
    return chains


def make_job(target_name: str, seed: int, chains: list[tuple[str, str]]) -> dict:
    job_name = f"{target_name.lower()}_aes_seed_{seed}".replace(" ", "_")
    return {
        "dialect": "alphafold3",
        "version": 1,
        "name": job_name,
        "sequences": [protein_entry(chain_id, sequence) for chain_id, sequence in chains],
        "modelSeeds": [int(seed)],
        "bondedAtomPairs": None,
        "userCCD": None,
    }


def build_af3_inputs() -> list[str]:
    cfg = load_inputs()
    seeds = cfg.get("af3", {}).get("seeds", [])
    if not seeds:
        raise ValueError("No AF3 seeds configured in inputs.yaml")

    written: list[str] = []
    for target_cfg in cfg["targets"]:
        target_name = target_cfg["name"]
        chains = target_chains(target_cfg)
        out_dir = af3_input_dir(target_name)
        for seed in seeds:
            job = make_job(target_name, int(seed), chains)
            out_path = out_dir / f"{job['name']}.json"
            json_dump(job, out_path)
            written.append(out_path.as_posix())
    return written


def main() -> None:
    make_arg_parser("Prepare per-target, per-seed AlphaFold 3 JSON inputs.").parse_args()
    written = build_af3_inputs()
    print(f"Wrote {len(written)} AF3 input JSON files.")


if __name__ == "__main__":
    main()
