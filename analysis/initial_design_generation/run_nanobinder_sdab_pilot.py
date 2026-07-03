#!/usr/bin/env python3
"""Run a local NanoBinder-style Rosetta feature pilot for sdAb AF3 models.

The upstream NanoBinder repository currently exposes training scripts and the
published feature list, but no pretrained model or prediction CLI. This script
therefore runs a conservative local pilot: it prepares sdAb AF3 complex models
in the chain convention described by NanoBinder (nanobody H, antigen A), then
extracts Rosetta interface features that overlap the published NanoBinder
feature set.

Run with the PyRosetta environment:

    /data/ziyang/mamba/envs/pyrosetta/bin/python \
      analysis/initial_design_generation/run_nanobinder_sdab_pilot.py
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
AF3_DIR = ROOT / "results/initial_design_generation/af3_1000seed_validation"
DEFAULT_OUT = ROOT / "results/initial_design_generation/nanobinder_sdab_pilot"

VARIANT_SUMMARY = AF3_DIR / "af3_1000seed_variant_summary.csv"
PARENT_COMPARISON = AF3_DIR / "af3_1000seed_parent_baseline_comparison.csv"
NANOBINDER_DIR = ROOT / "third_party/NanoBinder"

NANOBINDER_FINAL_FEATURES = [
    "complex_normalized",
    "dG_cross",
    "dG_cross/dSASAx100",
    "dSASA_hphobic",
    "dSASA_int",
    "dSASA_polar",
    "delta_unsatHbonds",
    "dslf_fa13",
    "fa_atr",
    "hbond_E_fraction",
    "hbond_bb_sc",
    "hbond_lr_bb",
    "hbond_sc",
    "hbond_sr_bb",
    "hbonds_int",
    "nres_int",
    "omega",
    "per_residue_energy_int",
    "pro_close",
    "rama_prepro",
    "ref",
    "side1_normalized",
    "side1_score",
    "side2_normalized",
    "side2_score",
    "yhh_planarity",
]


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def safe_float(value: object) -> float:
    try:
        if pd.isna(value):
            return math.nan
        return float(value)
    except Exception:
        return math.nan


def run_git_head(path: Path) -> str:
    if not (path / ".git").exists():
        return "not_a_git_checkout"
    try:
        out = subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True)
        return out.strip()
    except Exception as exc:
        return f"unavailable:{type(exc).__name__}"


def select_panel(limit_unsupported: int) -> pd.DataFrame:
    df = pd.read_csv(VARIANT_SUMMARY)
    parent = pd.read_csv(PARENT_COMPARISON)
    parent = parent[["target", "variant_id", "mean_parent_contact_retention", "delta_interpretation"]]
    df = df.merge(parent, on=["target", "variant_id"], how="left")
    sdab = df[df["target"].eq("Ab_sdAb")].copy()

    rows: list[pd.Series] = []
    parent_rows = sdab[sdab["variant_id"].eq("sdAb_parent_WT")]
    if not parent_rows.empty:
        row = parent_rows.iloc[0].copy()
        row["pilot_role"] = "parent_baseline"
        rows.append(row)

    plausible = sdab[sdab["static_classification"].eq("best-plausible")].copy()
    plausible = plausible.sort_values(
        ["mean_parent_contact_retention", "best_model_score"],
        ascending=[False, False],
        na_position="last",
    )
    for _, r in plausible.iterrows():
        row = r.copy()
        row["pilot_role"] = "sdAb_best_plausible"
        rows.append(row)

    known_negative_ids = [
        "sdAb_sdAb_VHH_072_111_04e3eaf858",
        "sdAb_sdAb_VHH_072_111_09d883bc8d",
        "sdAb_sdAb_VHH_072_111_4e063dd35e",
    ]
    selected_ids = {str(r["variant_id"]) for r in rows}
    negatives = []
    for vid in known_negative_ids:
        sub = sdab[sdab["variant_id"].eq(vid)]
        if not sub.empty and vid not in selected_ids:
            row = sub.iloc[0].copy()
            row["pilot_role"] = "known_shifted_or_unsupported_control"
            negatives.append(row)
            selected_ids.add(vid)

    remaining = sdab[
        sdab["static_classification"].eq("unsupported")
        & ~sdab["variant_id"].isin(selected_ids)
    ].copy()
    remaining = remaining.sort_values(
        ["best_model_score", "interface_plausible_fraction"],
        ascending=[False, False],
        na_position="last",
    )
    for _, r in remaining.head(max(0, limit_unsupported - len(negatives))).iterrows():
        row = r.copy()
        row["pilot_role"] = "high_score_unsupported_control"
        negatives.append(row)

    rows.extend(negatives[:limit_unsupported])
    out = pd.DataFrame(rows).drop_duplicates("variant_id", keep="first")
    keep = [
        "target",
        "variant_id",
        "pilot_role",
        "static_classification",
        "parent_comparison_class",
        "delta_interpretation",
        "mean_parent_contact_retention",
        "best_model_score",
        "best_cluster_fraction",
        "interface_plausible_fraction",
        "best_model_path",
    ]
    return out[[c for c in keep if c in out.columns]]


def parse_atom_site(cif_path: Path) -> tuple[list[str], list[list[str]]]:
    columns: list[str] = []
    rows: list[list[str]] = []
    in_atom_site = False
    with cif_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped == "loop_":
                if in_atom_site:
                    break
                columns = []
                continue
            if stripped.startswith("_atom_site."):
                in_atom_site = True
                columns.append(stripped)
                continue
            if in_atom_site:
                if stripped.startswith("_") or stripped.startswith("#"):
                    break
                parts = stripped.split()
                if len(parts) >= len(columns):
                    rows.append(parts[: len(columns)])
    return columns, rows


def pdb_atom_name(name: str, element: str) -> str:
    # PDB atom names are right-aligned for one-letter elements, with common
    # four-character names left intact.
    if len(name) >= 4:
        return name[:4]
    if len(element) == 1:
        return f" {name:<3}"[:4]
    return f"{name:<4}"[:4]


def cif_to_nanobinder_pdb(cif_path: Path, out_pdb: Path) -> dict[str, object]:
    cols, rows = parse_atom_site(cif_path)
    idx = {c.replace("_atom_site.", ""): i for i, c in enumerate(cols)}
    required = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_comp_id",
        "label_asym_id",
        "auth_seq_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
    ]
    missing = [c for c in required if c not in idx]
    if missing:
        raise ValueError(f"missing atom_site columns in {cif_path}: {missing}")

    chain_counts: dict[str, int] = {}
    serial = 1
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    with out_pdb.open("w", encoding="utf-8") as out:
        for parts in rows:
            group = parts[idx["group_PDB"]]
            if group not in {"ATOM", "HETATM"}:
                continue
            src_chain = parts[idx["label_asym_id"]]
            # AF3 sdAb outputs use chain A for nanobody and chain B for antigen.
            if src_chain == "A":
                chain = "H"
            elif src_chain == "B":
                chain = "A"
            else:
                chain = src_chain[:1] or "X"
            chain_counts[chain] = chain_counts.get(chain, 0) + 1
            atom = parts[idx["label_atom_id"]]
            resn = parts[idx["label_comp_id"]]
            elem = parts[idx["type_symbol"]]
            resi = int(float(parts[idx["auth_seq_id"]]))
            x = float(parts[idx["Cartn_x"]])
            y = float(parts[idx["Cartn_y"]])
            z = float(parts[idx["Cartn_z"]])
            occ = float(parts[idx["occupancy"]])
            bfac = float(parts[idx["B_iso_or_equiv"]])
            record = "HETATM" if group == "HETATM" else "ATOM  "
            out.write(
                f"{record}{serial:5d} {pdb_atom_name(atom, elem)} {resn:>3} {chain:1}"
                f"{resi:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}"
                f"          {elem:>2}\n"
            )
            serial += 1
        out.write("END\n")
    return {"atom_count": serial - 1, "chain_atom_counts": chain_counts}


def score_pose_with_pyrosetta(pdb_path: Path) -> dict[str, object]:
    import pyrosetta
    from pyrosetta.rosetta.core.scoring import ScoreType
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

    if not pyrosetta.rosetta.basic.was_init_called():
        pyrosetta.init("-mute all -ignore_unrecognized_res true")

    pose = pyrosetta.pose_from_pdb(str(pdb_path))
    scorefxn = pyrosetta.create_score_function("ref2015")
    total_score = float(scorefxn(pose))
    energies = pose.energies().total_energies()

    iam = InterfaceAnalyzerMover("H_A")
    iam.set_pack_separated(False)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_compute_interface_delta_hbond_unsat(True)
    iam.set_compute_packstat(True)
    iam.apply(pose)

    dsasa = float(iam.get_interface_delta_sasa())
    dG_cross = float(iam.get_crossterm_interface_energy())
    dG_separated = float(iam.get_separated_interface_energy())
    dG = float(iam.get_interface_dG())
    hbonds_int = float(iam.get_total_Hbond_E())
    nres_int = int(iam.get_num_interface_residues())
    side1_score = float(iam.get_side1_score())
    side2_score = float(iam.get_side2_score())
    side1_nres = int(iam.get_side1_nres())
    side2_nres = int(iam.get_side2_nres())
    side1_normalized = side1_score / side1_nres if side1_nres else math.nan
    side2_normalized = side2_score / side2_nres if side2_nres else math.nan
    complex_normalized = total_score / pose.total_residue() if pose.total_residue() else math.nan
    per_residue_energy_int = dG / nres_int if nres_int else math.nan

    def energy(score_type: object) -> float:
        try:
            return float(energies[score_type])
        except Exception:
            return math.nan

    out = {
        "rosetta_status": "success",
        "prepared_pdb_path": rel(pdb_path),
        "pose_total_residue": pose.total_residue(),
        "total_score": total_score,
        "complex_normalized": complex_normalized,
        "dG_cross": dG_cross,
        "dG_cross/dSASAx100": (dG_cross / dsasa * 100.0) if dsasa else math.nan,
        "dG_separated": dG_separated,
        "dG_separated/dSASAx100": (dG_separated / dsasa * 100.0) if dsasa else math.nan,
        "dSASA_int": dsasa,
        "delta_unsatHbonds": float(iam.get_interface_delta_hbond_unsat()),
        "hbonds_int": hbonds_int,
        "hbond_E_fraction": (hbonds_int / dG) if dG else math.nan,
        "nres_int": nres_int,
        "per_residue_energy_int": per_residue_energy_int,
        "side1_score": side1_score,
        "side2_score": side2_score,
        "side1_normalized": side1_normalized,
        "side2_normalized": side2_normalized,
        "packstat": float(iam.get_interface_packstat()),
        "fa_atr": energy(ScoreType.fa_atr),
        "fa_rep": energy(ScoreType.fa_rep),
        "fa_sol": energy(ScoreType.fa_sol),
        "fa_elec": energy(ScoreType.fa_elec),
        "hbond_bb_sc": energy(ScoreType.hbond_bb_sc),
        "hbond_lr_bb": energy(ScoreType.hbond_lr_bb),
        "hbond_sc": energy(ScoreType.hbond_sc),
        "hbond_sr_bb": energy(ScoreType.hbond_sr_bb),
        "omega": energy(ScoreType.omega),
        "pro_close": energy(ScoreType.pro_close),
        "rama_prepro": energy(ScoreType.rama_prepro),
        "ref": energy(ScoreType.ref),
        "dslf_fa13": energy(ScoreType.dslf_fa13),
        "yhh_planarity": energy(ScoreType.yhh_planarity),
    }
    # These are published NanoBinder fields but not available from this simple
    # InterfaceAnalyzer pass without a full RosettaAntibody scorefile workflow.
    out["dSASA_hphobic"] = math.nan
    out["dSASA_polar"] = math.nan
    return out


def add_relative_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    parent_rows = out[out["pilot_role"].eq("parent_baseline")]
    if parent_rows.empty:
        out["relative_to_parent_dG_cross"] = math.nan
        out["relative_to_parent_dG_cross_dSASA"] = math.nan
        out["nanobinder_pilot_interpretation"] = "no_parent_baseline"
        return out
    parent = parent_rows.iloc[0]
    parent_ratio = safe_float(parent.get("dG_cross/dSASAx100"))
    parent_dg = safe_float(parent.get("dG_cross"))
    flags = []
    rel_ratio = []
    rel_dg = []
    for _, row in out.iterrows():
        ratio = safe_float(row.get("dG_cross/dSASAx100"))
        dg = safe_float(row.get("dG_cross"))
        rr = ratio - parent_ratio if not math.isnan(ratio) and not math.isnan(parent_ratio) else math.nan
        rd = dg - parent_dg if not math.isnan(dg) and not math.isnan(parent_dg) else math.nan
        rel_ratio.append(rr)
        rel_dg.append(rd)
        cls = str(row.get("static_classification", ""))
        retention = safe_float(row.get("mean_parent_contact_retention"))
        if row.get("pilot_role") == "parent_baseline":
            flags.append("parent_reference")
        elif cls == "best-plausible" and retention >= 0.30 and (math.isnan(rr) or rr <= 5.0):
            flags.append("rosetta_interface_parent_like")
        elif cls == "best-plausible" and retention >= 0.15:
            flags.append("rosetta_interface_partial_support")
        elif not math.isnan(rr) and rr <= 0 and retention < 0.15:
            flags.append("binder_like_energy_but_shifted_pose")
        else:
            flags.append("rosetta_interface_weak_or_shifted")
    out["relative_to_parent_dG_cross"] = rel_dg
    out["relative_to_parent_dG_cross_dSASA"] = rel_ratio
    out["nanobinder_pilot_interpretation"] = flags
    return out


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_empty_"
    text = df.copy()
    for col in text.columns:
        text[col] = text[col].map(lambda v: "" if pd.isna(v) else str(v))
    widths = {
        col: max(len(str(col)), *(len(str(v)) for v in text[col].tolist()))
        for col in text.columns
    }
    header = "| " + " | ".join(str(col).ljust(widths[col]) for col in text.columns) + " |"
    sep = "| " + " | ".join("-" * widths[col] for col in text.columns) + " |"
    rows = [
        "| " + " | ".join(str(row[col]).ljust(widths[col]) for col in text.columns) + " |"
        for _, row in text.iterrows()
    ]
    return "\n".join([header, sep, *rows])


def write_report(
    panel: pd.DataFrame,
    features: pd.DataFrame,
    out_dir: Path,
    smoke_status: str,
    nanobinder_head: str,
) -> None:
    counts = features.groupby(["pilot_role", "static_classification"], dropna=False).size().reset_index(name="count")
    interp = features.groupby(["nanobinder_pilot_interpretation"], dropna=False).size().reset_index(name="count")
    key_cols = [
        "variant_id",
        "pilot_role",
        "static_classification",
        "mean_parent_contact_retention",
        "dG_cross",
        "dG_cross/dSASAx100",
        "relative_to_parent_dG_cross_dSASA",
        "dSASA_int",
        "hbonds_int",
        "nres_int",
        "nanobinder_pilot_interpretation",
    ]
    lines = [
        "# NanoBinder sdAb Pilot Report",
        "",
        "## Summary",
        "",
        "- Upstream NanoBinder source was installed under `third_party/NanoBinder`.",
        f"- Upstream git HEAD: `{nanobinder_head}`.",
        f"- Upstream Training.py smoke test: `{smoke_status}`.",
        "- No upstream pretrained model or prediction CLI was present in the repository.",
        "- This pilot therefore extracted NanoBinder-style Rosetta interface features, not NanoBinder binding probabilities.",
        "",
        "## Panel",
        "",
        markdown_table(counts),
        "",
        "## Pilot Interpretation Counts",
        "",
        markdown_table(interp),
        "",
        "## Key Rows",
        "",
        markdown_table(features[[c for c in key_cols if c in features.columns]]),
        "",
        "## Practical Reading",
        "",
        "- Treat this as a feasibility check for adding Rosetta/NanoBinder-style evidence to the sdAb branch.",
        "- Rows marked `binder_like_energy_but_shifted_pose` are exactly the risk case we care about: Rosetta-like energies may look acceptable even when AF3 contact geometry is shifted.",
        "- A true NanoBinder probability requires the missing released model or retraining data; this repository currently only exposes the training scripts and feature definitions.",
        "",
        "## Output Files",
        "",
        f"- `{rel(out_dir / 'nanobinder_sdab_pilot_panel.csv')}`",
        f"- `{rel(out_dir / 'nanobinder_sdab_prepared_structures_manifest.csv')}`",
        f"- `{rel(out_dir / 'nanobinder_sdab_rosetta_feature_pilot.csv')}`",
        f"- `{rel(out_dir / 'nanobinder_sdab_pilot_report.md')}`",
    ]
    (out_dir / "nanobinder_sdab_pilot_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    parser.add_argument("--unsupported-controls", type=int, default=6)
    parser.add_argument("--keep-prepared-pdbs", action="store_true", default=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    prepared_dir = out_dir / "prepared_pdb"
    out_dir.mkdir(parents=True, exist_ok=True)

    panel = select_panel(args.unsupported_controls)
    panel.to_csv(out_dir / "nanobinder_sdab_pilot_panel.csv", index=False)

    nanobinder_head = run_git_head(NANOBINDER_DIR)
    smoke_status = "not_run"
    if (NANOBINDER_DIR / "Training.py").exists():
        proc = subprocess.run(
            [sys.executable, str(NANOBINDER_DIR / "Training.py")],
            cwd=str(NANOBINDER_DIR),
            text=True,
            capture_output=True,
            timeout=20,
        )
        if proc.returncode == 0:
            smoke_status = "success"
        else:
            last = (proc.stderr or proc.stdout).strip().splitlines()[-1] if (proc.stderr or proc.stdout) else ""
            smoke_status = f"failed_as_expected:{last[:160]}"

    prepared_rows = []
    feature_rows = []
    for _, row in panel.iterrows():
        variant_id = str(row["variant_id"])
        cif_path = ROOT / str(row["best_model_path"])
        safe_id = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in variant_id)
        pdb_path = prepared_dir / f"{safe_id}_HA.pdb"
        prep = {
            "target": row.get("target"),
            "variant_id": variant_id,
            "pilot_role": row.get("pilot_role"),
            "source_cif_path": rel(cif_path),
            "prepared_pdb_path": rel(pdb_path),
        }
        try:
            prep.update(cif_to_nanobinder_pdb(cif_path, pdb_path))
            prep["prepare_status"] = "success"
            score = score_pose_with_pyrosetta(pdb_path)
        except Exception as exc:
            prep["prepare_status"] = f"failed:{type(exc).__name__}:{exc}"
            score = {"rosetta_status": f"failed:{type(exc).__name__}:{exc}"}
        prepared_rows.append(prep)
        merged = row.to_dict()
        merged.update(score)
        feature_rows.append(merged)

    prepared = pd.DataFrame(prepared_rows)
    features = pd.DataFrame(feature_rows)
    features = add_relative_flags(features)

    prepared.to_csv(out_dir / "nanobinder_sdab_prepared_structures_manifest.csv", index=False)
    features.to_csv(out_dir / "nanobinder_sdab_rosetta_feature_pilot.csv", index=False)

    coverage_rows = []
    for feat in NANOBINDER_FINAL_FEATURES:
        coverage_rows.append(
            {
                "feature": feat,
                "available_in_pilot": feat in features.columns and features[feat].notna().any(),
                "non_null_count": int(features[feat].notna().sum()) if feat in features.columns else 0,
            }
        )
    pd.DataFrame(coverage_rows).to_csv(out_dir / "nanobinder_published_feature_coverage.csv", index=False)

    write_report(panel, features, out_dir, smoke_status, nanobinder_head)
    print(f"Wrote {rel(out_dir / 'nanobinder_sdab_pilot_report.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
