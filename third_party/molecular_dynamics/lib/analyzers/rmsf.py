# third_party/molecular_dynamics/lib/analyzers/rmsf.py
"""Per-residue Cα RMSF 分析 + CDR 区域统计。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import numpy as np
import pandas as pd

from .base import BaseAnalyzer, detect_convergence

logger = logging.getLogger(__name__)


def _normalize_chain(segid: str) -> str:
    """将 GROMACS segid (如 seg_0_Protein_chain_A) 归一化为简单链名 (A)。"""
    if "_chain_" in segid:
        return segid.rsplit("_chain_", 1)[-1]
    return segid


class RMSFAnalyzer(BaseAnalyzer):
    """Per-residue Cα RMSF，自动跳过非收敛区间。"""

    name = "rmsf"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))
        ca = u.select_atoms("name CA")

        # --- Quick RMSD for convergence detection ---
        rmsd_analysis = rms.RMSD(ca, ca, ref_frame=0)
        rmsd_analysis.run()
        time_ns = rmsd_analysis.results.rmsd[:, 1] / 1000.0
        global_rmsd = rmsd_analysis.results.rmsd[:, 2]

        conv_cfg = self.analysis_cfg["convergence"]
        converge_ns = detect_convergence(
            time_ns, global_rmsd,
            window_ns=conv_cfg["window_ns"],
            variance_threshold=conv_cfg["variance_threshold"],
        )
        start_frame = 0
        if converge_ns > 0:
            start_frame = int(np.searchsorted(time_ns, converge_ns))
        logger.info("RMSF: using frames from %d (converge=%.1f ns)", start_frame, converge_ns)

        # --- Align to average structure (post-convergence) ---
        avg = align.AverageStructure(u, u, select="name CA", ref_frame=0)
        avg.run(start=start_frame)
        ref = avg.results.universe

        aligner = align.AlignTraj(u, ref, select="name CA", in_memory=True)
        aligner.run(start=start_frame)

        # --- RMSF ---
        rmsf_analysis = rms.RMSF(ca)
        rmsf_analysis.run(start=start_frame)
        rmsf_vals = rmsf_analysis.results.rmsf

        # Build per-residue dataframe
        chains = []
        resids = []
        resnames = []
        for atom in ca:
            raw = atom.segid if atom.segid.strip() else atom.chainID
            chains.append(_normalize_chain(raw))
            resids.append(atom.resid)
            resnames.append(atom.resname)

        df = pd.DataFrame({
            "chain": chains,
            "resid": resids,
            "resname": resnames,
            "rmsf": rmsf_vals,
        })

        # --- CDR-specific RMSF ---
        cdr_regions = self.analysis_cfg["cdr_regions"]
        cdr_means = {}
        for cdr_name, cdr_def in cdr_regions.items():
            mask = (
                (df["chain"] == cdr_def["chain"])
                & (df["resid"] >= cdr_def["start"])
                & (df["resid"] <= cdr_def["end"])
            )
            if mask.any():
                cdr_means[cdr_name] = float(df.loc[mask, "rmsf"].mean())
            else:
                cdr_means[cdr_name] = float("nan")
                logger.warning("CDR %s: no residues found for RMSF", cdr_name)

        results = {
            "df": df,
            "converge_ns": converge_ns,
            "summary": {
                "rmsf_global_mean": float(np.mean(rmsf_vals)),
                "converge_ns": converge_ns,
            },
        }
        for cdr_name, val in cdr_means.items():
            results["summary"][f"{cdr_name.lower()}_rmsf_mean"] = val

        self.save(results)
        logger.info("RMSF analysis: global_mean=%.3f Å", results["summary"]["rmsf_global_mean"])
        return results

    def save(self, results: dict[str, Any]) -> Path:
        out = self.output_dir / "rmsf_per_residue.csv"
        results["df"].to_csv(out, index=False, float_format="%.4f")
        logger.info("Saved: %s", out)
        return out
