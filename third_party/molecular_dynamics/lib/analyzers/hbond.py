# third_party/molecular_dynamics/lib/analyzers/hbond.py
"""界面氢键占有率分析。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
import pandas as pd

from .base import BaseAnalyzer

logger = logging.getLogger(__name__)


def _chain_selection(chains: list[str], keyword: str = "segid") -> str:
    """Build selection string for multiple chains."""
    parts = [f"{keyword} {c}" for c in chains]
    return "(" + " or ".join(parts) + ")"


class HBondAnalyzer(BaseAnalyzer):
    """界面氢键分析：占有率 + 时间序列。"""

    name = "hbond"

    def _select_groups(self, u: mda.Universe):
        """Select antibody and antigen atoms, trying segid then chainID."""
        ab_chains = self.analysis_cfg["antibody_chains"]
        ag_chains = self.analysis_cfg["antigen_chains"]

        ab = u.select_atoms(_chain_selection(ab_chains, "segid"))
        if len(ab) == 0:
            ab = u.select_atoms(_chain_selection(ab_chains, "chainID"))

        ag = u.select_atoms(_chain_selection(ag_chains, "segid"))
        if len(ag) == 0:
            ag = u.select_atoms(_chain_selection(ag_chains, "chainID"))

        return ab, ag

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))
        ab, ag = self._select_groups(u)

        if len(ab) == 0 or len(ag) == 0:
            logger.error("Could not select antibody (%d) or antigen (%d) atoms",
                         len(ab), len(ag))
            raise ValueError("Empty antibody or antigen selection")

        ab_indices = set(ab.indices)
        ag_indices = set(ag.indices)

        # Run H-bond analysis on entire system
        hbond = HydrogenBondAnalysis(u)
        hbond.run()

        # hbond.results.hbonds columns:
        # [frame, donor_ix, hydrogen_ix, acceptor_ix, distance, angle]
        hbonds = hbond.results.hbonds

        if len(hbonds) == 0:
            logger.warning("No hydrogen bonds found in trajectory")
            results = {
                "occupancy_df": pd.DataFrame(
                    columns=["donor_chain", "donor_resid", "donor_resname",
                             "acceptor_chain", "acceptor_resid", "acceptor_resname",
                             "occupancy"]
                ),
                "timeseries_df": pd.DataFrame(columns=["time_ns", "n_interface_hbonds"]),
                "summary": {"n_hbond_mean": 0.0, "n_hbond_std": 0.0, "n_unique_pairs": 0},
            }
            self.save(results)
            return results

        n_frames = len(u.trajectory)

        # Filter interface H-bonds: donor in ab & acceptor in ag, or vice versa
        interface_mask = np.zeros(len(hbonds), dtype=bool)
        for i, row in enumerate(hbonds):
            donor_ix = int(row[1])
            acceptor_ix = int(row[3])
            d_in_ab = donor_ix in ab_indices
            d_in_ag = donor_ix in ag_indices
            a_in_ab = acceptor_ix in ab_indices
            a_in_ag = acceptor_ix in ag_indices
            if (d_in_ab and a_in_ag) or (d_in_ag and a_in_ab):
                interface_mask[i] = True

        iface_hbonds = hbonds[interface_mask]

        # Per-residue-pair occupancy
        pair_counts: dict[tuple, int] = {}
        pair_info: dict[tuple, dict] = {}
        frame_counts: dict[int, int] = {}

        for row in iface_hbonds:
            frame = int(row[0])
            donor_ix = int(row[1])
            acceptor_ix = int(row[3])

            d_atom = u.atoms[donor_ix]
            a_atom = u.atoms[acceptor_ix]

            d_chain = d_atom.segid if d_atom.segid.strip() else d_atom.chainID
            a_chain = a_atom.segid if a_atom.segid.strip() else a_atom.chainID

            pair_key = (d_chain, d_atom.resid, a_chain, a_atom.resid)
            if pair_key not in pair_counts:
                pair_counts[pair_key] = set()
                pair_info[pair_key] = {
                    "donor_chain": d_chain,
                    "donor_resid": d_atom.resid,
                    "donor_resname": d_atom.resname,
                    "acceptor_chain": a_chain,
                    "acceptor_resid": a_atom.resid,
                    "acceptor_resname": a_atom.resname,
                }
            pair_counts[pair_key].add(frame)
            frame_counts[frame] = frame_counts.get(frame, 0) + 1

        # Occupancy dataframe
        occ_rows = []
        for pair_key, frames_set in pair_counts.items():
            info = pair_info[pair_key]
            info["occupancy"] = len(frames_set) / n_frames
            occ_rows.append(info)
        occupancy_df = pd.DataFrame(occ_rows)
        if not occupancy_df.empty:
            occupancy_df.sort_values("occupancy", ascending=False, inplace=True)

        # Timeseries
        ts_data = []
        for i, ts in enumerate(u.trajectory):
            t_ns = ts.time / 1000.0
            n_hb = frame_counts.get(i, 0)
            ts_data.append({"time_ns": t_ns, "n_interface_hbonds": n_hb})
        timeseries_df = pd.DataFrame(ts_data)

        n_hbond_series = timeseries_df["n_interface_hbonds"].values
        results = {
            "occupancy_df": occupancy_df,
            "timeseries_df": timeseries_df,
            "summary": {
                "n_hbond_mean": float(np.mean(n_hbond_series)),
                "n_hbond_std": float(np.std(n_hbond_series)),
                "n_unique_pairs": len(pair_counts),
            },
        }

        self.save(results)
        logger.info(
            "H-bond analysis: mean=%.1f, unique_pairs=%d",
            results["summary"]["n_hbond_mean"],
            results["summary"]["n_unique_pairs"],
        )
        return results

    def save(self, results: dict[str, Any]) -> Path:
        occ_path = self.output_dir / "hbond_occupancy.csv"
        results["occupancy_df"].to_csv(occ_path, index=False, float_format="%.4f")

        ts_path = self.output_dir / "hbond_timeseries.csv"
        results["timeseries_df"].to_csv(ts_path, index=False, float_format="%.4f")

        logger.info("Saved: %s, %s", occ_path, ts_path)
        return occ_path
