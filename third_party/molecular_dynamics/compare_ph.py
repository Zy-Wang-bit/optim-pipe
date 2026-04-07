#!/usr/bin/env python3
# third_party/molecular_dynamics/compare_ph.py
"""双 pH 差异比较。

用法：
  python compare_ph.py --variant-dir experiments/1E62/R3/md/HE1H/ --base-ph 7.4 --target-ph 6.0
"""

import argparse
import json
import logging
import sys
from pathlib import Path

_MODULE_ROOT = Path(__file__).resolve().parent
if str(_MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(_MODULE_ROOT))

from lib.config import load_config

logger = logging.getLogger(__name__)


def compare_ph(variant_dir: Path, base_ph: float, target_ph: float) -> dict:
    """比较同一变体在两个 pH 下的分析结果差异。"""
    base_dir = variant_dir / f"pH_{base_ph}" / "analysis"
    target_dir = variant_dir / f"pH_{target_ph}" / "analysis"

    if not base_dir.exists():
        raise FileNotFoundError(f"Base pH analysis not found: {base_dir}")
    if not target_dir.exists():
        raise FileNotFoundError(f"Target pH analysis not found: {target_dir}")

    base_summary = _load_summary(base_dir / "summary.json")
    target_summary = _load_summary(target_dir / "summary.json")

    # Δ = target - base
    deltas = {}
    compare_keys = [
        ("rmsd_rmsd_mean", "delta_rmsd"),
        ("rmsf_rmsf_global_mean", "delta_rmsf_global"),
        ("hbond_n_hbond_mean", "delta_n_hbond"),
        ("sasa_buried_sasa_mean", "delta_buried_sasa"),
        ("contacts_n_contacts_mean", "delta_n_contacts"),
        ("mmpbsa_dG_bind_mean", "delta_dG_bind"),
    ]
    for cdr in ("h1", "h2", "h3", "l1", "l2", "l3"):
        compare_keys.append((f"rmsf_rmsf_{cdr}_mean", f"delta_rmsf_{cdr}"))

    for src_key, delta_key in compare_keys:
        base_val = base_summary.get(src_key)
        target_val = target_summary.get(src_key)
        if base_val is not None and target_val is not None:
            try:
                deltas[delta_key] = float(target_val) - float(base_val)
            except (ValueError, TypeError):
                deltas[delta_key] = None

    result = {
        "variant": variant_dir.name,
        "base_ph": base_ph,
        "target_ph": target_ph,
        "base_summary": base_summary,
        "target_summary": target_summary,
        "deltas": deltas,
    }

    out_path = variant_dir / "ph_comparison.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    logger.info("pH comparison saved: %s", out_path)

    logger.info("=== pH %.1f vs %.1f comparison for %s ===", target_ph, base_ph, variant_dir.name)
    for key, val in deltas.items():
        if val is not None:
            direction = "↑" if val > 0 else "↓" if val < 0 else "="
            logger.info("  %s: %+.4f %s", key, val, direction)

    return result


def _load_summary(path: Path) -> dict:
    if not path.exists():
        logger.warning("Summary not found: %s", path)
        return {}
    with open(path) as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description="双 pH 分析比较")
    parser.add_argument("--variant-dir", type=Path, required=True)
    parser.add_argument("--base-ph", type=float, default=7.4)
    parser.add_argument("--target-ph", type=float, default=6.0)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    compare_ph(args.variant_dir, args.base_ph, args.target_ph)


if __name__ == "__main__":
    main()
