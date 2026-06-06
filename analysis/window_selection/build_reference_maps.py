#!/usr/bin/env python
from __future__ import annotations

import pandas as pd

from analysis.window_selection.common import (
    cdr_region,
    load_inputs,
    make_arg_parser,
    read_fasta,
    write_csv,
)


def build_reference_maps() -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = load_inputs()
    rows = []
    antigen_rows = []

    for target in cfg["targets"]:
        name = target["name"]
        if target["antibody_type"] == "scFv":
            chains = [
                ("heavy", "H", target["heavy_fasta"]),
                ("light", "L", target["light_fasta"]),
            ]
        else:
            chains = [("sdab", "A", target["sdab_fasta"])]

        for molecule, chain, fasta in chains:
            _, _, seq = read_fasta(fasta)
            for pos, aa in enumerate(seq, start=1):
                region = cdr_region(chain, pos, name)
                rows.append(
                    {
                        "target": name,
                        "molecule": molecule,
                        "chain": chain,
                        "local_pos": pos,
                        "aa": aa,
                        "global_pos": f"{chain}{pos}",
                        "region": region,
                        "is_design_eligible_base": True,
                        "ab_numbering_scheme": "local_fallback",
                        "ab_number": f"{chain}{pos}",
                        "cdr_scheme": "local_fallback",
                        "pdb_chain": "",
                        "pdb_resseq": "",
                        "pdb_icode": "",
                        "source_position_label": f"{chain}{pos}",
                    }
                )

        antigen_id = target.get("antigen_name") if target["antibody_type"] == "scFv" else None
        if target["antibody_type"] == "scFv":
            _, _, antigen_seq = read_fasta(target["antigen_fasta"], antigen_id)
        else:
            _, _, antigen_seq = read_fasta(target["antigen_fasta"])
        for pos, aa in enumerate(antigen_seq, start=1):
            antigen_rows.append(
                {
                    "target": name,
                    "antigen_pos": pos,
                    "aa": aa,
                    "is_mhr_99_169": 99 <= pos <= 169,
                    "is_antigenic_loop_100_164": 100 <= pos <= 164,
                    "is_a_determinant_124_147": 124 <= pos <= 147,
                    "is_n146_motif_146_148": 146 <= pos <= 148,
                }
            )

    return pd.DataFrame(rows), pd.DataFrame(antigen_rows)


def main() -> None:
    make_arg_parser("Build reference sequence and AeS region maps.").parse_args()
    ref, antigen = build_reference_maps()
    write_csv(ref, "reference_sequence_map.csv")
    write_csv(antigen, "hbsag_aes_regions.csv")
    print(f"Wrote {len(ref)} reference rows and {len(antigen)} antigen rows.")


if __name__ == "__main__":
    main()
