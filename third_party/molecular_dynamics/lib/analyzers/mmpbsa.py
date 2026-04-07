# third_party/molecular_dynamics/lib/analyzers/mmpbsa.py
"""MM-GBSA / MM-PBSA 结合能估算（包装 gmx_MMPBSA）。"""

import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .base import BaseAnalyzer

logger = logging.getLogger(__name__)


class MMPBSAAnalyzer(BaseAnalyzer):
    """Optional MM-GBSA/MM-PBSA analyzer wrapping gmx_MMPBSA."""

    name = "mmpbsa"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        mmpbsa_cfg = self.analysis_cfg.get("mmpbsa", {})
        if not mmpbsa_cfg.get("enabled", False):
            logger.info("MM-PBSA is disabled in config, skipping")
            return {"summary": {"skipped": True}}

        # Check for gmx_MMPBSA binary
        gmx_mmpbsa = shutil.which("gmx_MMPBSA")
        if gmx_mmpbsa is None:
            logger.error("gmx_MMPBSA not found in PATH")
            raise FileNotFoundError("gmx_MMPBSA binary not found")

        method = mmpbsa_cfg.get("method", "gb")  # 'gb' or 'pb'
        frame_interval_ps = mmpbsa_cfg.get("frame_interval_ps", 100)

        # Generate input file
        input_content = self._generate_input(method, frame_interval_ps)
        work_dir = Path(tempfile.mkdtemp(prefix="mmpbsa_", dir=str(self.output_dir)))
        input_file = work_dir / "mmpbsa.in"
        input_file.write_text(input_content)

        # Run gmx_MMPBSA
        output_file = work_dir / "FINAL_RESULTS_MMPBSA.dat"
        cmd = [
            gmx_mmpbsa,
            "-i", str(input_file),
            "-cs", str(self.tpr),
            "-ct", str(self.xtc),
            "-ci", str(self.tpr),   # complex topology
            "-eo", str(work_dir / "energy.csv"),
        ]

        logger.info("Running gmx_MMPBSA: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, capture_output=True, text=True, cwd=str(work_dir), timeout=7200,
        )

        if proc.returncode != 0:
            logger.error("gmx_MMPBSA failed:\nstdout: %s\nstderr: %s",
                         proc.stdout[-2000:], proc.stderr[-2000:])
            raise RuntimeError(f"gmx_MMPBSA exited with code {proc.returncode}")

        # Parse results
        results = self._parse_results(work_dir)
        self.save(results)

        dg = results["summary"].get("dG_bind_mean", float("nan"))
        logger.info("MM-PBSA analysis: dG_bind = %.2f ± %.2f kcal/mol",
                     dg, results["summary"].get("dG_bind_std", float("nan")))
        return results

    def _generate_input(self, method: str, frame_interval_ps: int) -> str:
        """Generate gmx_MMPBSA input file."""
        if method == "pb":
            return (
                "&general\n"
                f"  interval={frame_interval_ps},\n"
                "  verbose=2,\n"
                "/\n"
                "&pb\n"
                "  istrng=0.150,\n"
                "  fillratio=4.0,\n"
                "/\n"
            )
        else:
            # Default: GB (igb=5 is OBC2, commonly used)
            return (
                "&general\n"
                f"  interval={frame_interval_ps},\n"
                "  verbose=2,\n"
                "/\n"
                "&gb\n"
                "  igb=5,\n"
                "  saltcon=0.150,\n"
                "/\n"
            )

    def _parse_results(self, work_dir: Path) -> dict[str, Any]:
        """Parse gmx_MMPBSA output files."""
        # Try to parse the energy CSV if it exists
        energy_csv = work_dir / "energy.csv"
        final_dat = work_dir / "FINAL_RESULTS_MMPBSA.dat"

        dg_values = []

        # Try energy.csv first
        if energy_csv.exists():
            try:
                df = pd.read_csv(energy_csv)
                if "TOTAL" in df.columns:
                    dg_values = df["TOTAL"].dropna().tolist()
                elif "DeltaG" in df.columns:
                    dg_values = df["DeltaG"].dropna().tolist()
            except Exception as e:
                logger.warning("Could not parse energy.csv: %s", e)

        # Fallback: parse FINAL_RESULTS_MMPBSA.dat
        if not dg_values and final_dat.exists():
            text = final_dat.read_text()
            # Look for "DELTA TOTAL" line with mean and std
            match = re.search(r"DELTA\s+TOTAL\s+([-\d.]+)\s+([-\d.]+)", text)
            if match:
                dg_values = [float(match.group(1))]

        if dg_values:
            arr = np.array(dg_values)
            summary = {
                "dG_bind_mean": float(np.mean(arr)),
                "dG_bind_std": float(np.std(arr)),
            }
        else:
            logger.warning("Could not parse any dG values from gmx_MMPBSA output")
            summary = {
                "dG_bind_mean": float("nan"),
                "dG_bind_std": float("nan"),
            }

        return {
            "dg_values": dg_values,
            "work_dir": str(work_dir),
            "summary": summary,
        }

    def save(self, results: dict[str, Any]) -> Path:
        out = self.output_dir / "mmpbsa_summary.csv"
        summary = results.get("summary", {})

        if "skipped" in summary:
            pd.DataFrame([{"status": "skipped"}]).to_csv(out, index=False)
        else:
            dg_values = results.get("dg_values", [])
            if dg_values:
                df = pd.DataFrame({"frame": range(len(dg_values)), "dG_bind": dg_values})
            else:
                df = pd.DataFrame({
                    "metric": ["dG_bind_mean", "dG_bind_std"],
                    "value": [summary.get("dG_bind_mean"), summary.get("dG_bind_std")],
                })
            df.to_csv(out, index=False, float_format="%.4f")

        logger.info("Saved: %s", out)
        return out
