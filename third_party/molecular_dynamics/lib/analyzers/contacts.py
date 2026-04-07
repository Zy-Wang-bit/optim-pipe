# third_party/molecular_dynamics/lib/analyzers/contacts.py
"""界面接触分析：重原子距离 < cutoff 的残基对频率。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import pandas as pd

from .base import BaseAnalyzer

logger = logging.getLogger(__name__)


def _chain_selection(chains: list[str], keyword: str = "segid") -> str:
    parts = [f"{keyword} {c}" for c in chains]
    return "(" + " or ".join(parts) + ")"


def _select_heavy(u: mda.Universe, chains: list[str]) -> mda.AtomGroup:
    """Select heavy atoms (not H) for given chains, segid then chainID."""
    sel = f"not name H* and {_chain_selection(chains, 'segid')}"
    grp = u.select_atoms(sel)
    if len(grp) == 0:
        sel = f"not name H* and {_chain_selection(chains, 'chainID')}"
        grp = u.select_atoms(sel)
    return grp


class ContactsAnalyzer(BaseAnalyzer):
    """界面接触数时间序列 + 残基对接触频率。"""

    name = "contacts"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))
        ab_chains = self.analysis_cfg["antibody_chains"]
        ag_chains = self.analysis_cfg["antigen_chains"]
        cutoff = self.analysis_cfg.get("contact_distance", 4.5)

        ab = _select_heavy(u, ab_chains)
        ag = _select_heavy(u, ag_chains)

        if len(ab) == 0 or len(ag) == 0:
            raise ValueError(f"Empty selection: ab={len(ab)}, ag={len(ag)}")

        # Pre-compute residue info for ab and ag atoms
        ab_res_info = []
        for atom in ab:
            chain = atom.segid if atom.segid.strip() else atom.chainID
            ab_res_info.append((chain, atom.resid, atom.resname))
        ag_res_info = []
        for atom in ag:
            chain = atom.segid if atom.segid.strip() else atom.chainID
            ag_res_info.append((chain, atom.resid, atom.resname))

        n_frames = len(u.trajectory)
        pair_frame_counts: dict[tuple, int] = {}
        ts_data = []

        for frame_i, ts in enumerate(u.trajectory):
            dists = distance_array(ab.positions, ag.positions)
            contact_mask = dists < cutoff

            # Count total contacts this frame
            n_contacts = int(np.sum(contact_mask))
            ts_data.append({"time_ns": ts.time / 1000.0, "n_contacts": n_contacts})

            # Track residue-pair contacts (unique pairs per frame)
            ab_idx, ag_idx = np.where(contact_mask)
            seen_pairs: set[tuple] = set()
            for ai, gi in zip(ab_idx, ag_idx):
                ab_key = (ab_res_info[ai][0], ab_res_info[ai][1])
                ag_key = (ag_res_info[gi][0], ag_res_info[gi][1])
                pair = (ab_key, ag_key)
                if pair not in seen_pairs:
                    seen_pairs.add(pair)
                    full_key = (
                        ab_res_info[ai][0], ab_res_info[ai][1], ab_res_info[ai][2],
                        ag_res_info[gi][0], ag_res_info[gi][1], ag_res_info[gi][2],
                    )
                    pair_frame_counts[full_key] = pair_frame_counts.get(full_key, 0) + 1

        timeseries_df = pd.DataFrame(ts_data)

        # Frequency dataframe
        freq_rows = []
        for key, count in pair_frame_counts.items():
            freq_rows.append({
                "ab_chain": key[0],
                "ab_resid": key[1],
                "ab_resname": key[2],
                "ag_chain": key[3],
                "ag_resid": key[4],
                "ag_resname": key[5],
                "frequency": count / n_frames,
            })
        frequency_df = pd.DataFrame(freq_rows)
        if not frequency_df.empty:
            frequency_df.sort_values("frequency", ascending=False, inplace=True)

        n_contacts_arr = timeseries_df["n_contacts"].values
        results = {
            "timeseries_df": timeseries_df,
            "frequency_df": frequency_df,
            "summary": {
                "n_contacts_mean": float(np.mean(n_contacts_arr)),
                "n_contacts_std": float(np.std(n_contacts_arr)),
                "n_unique_pairs": len(pair_frame_counts),
            },
        }

        self.save(results)
        logger.info(
            "Contacts analysis: mean=%.1f, unique_pairs=%d",
            results["summary"]["n_contacts_mean"],
            results["summary"]["n_unique_pairs"],
        )
        return results

    def save(self, results: dict[str, Any]) -> Path:
        ts_path = self.output_dir / "contacts_timeseries.csv"
        results["timeseries_df"].to_csv(ts_path, index=False, float_format="%.4f")

        freq_path = self.output_dir / "contacts_frequency.csv"
        results["frequency_df"].to_csv(freq_path, index=False, float_format="%.4f")

        logger.info("Saved: %s, %s", ts_path, freq_path)
        return ts_path
