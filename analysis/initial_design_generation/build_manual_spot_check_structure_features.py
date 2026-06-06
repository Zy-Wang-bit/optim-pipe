#!/usr/bin/env python3
"""Build structure-feature table for manual spot-check variants.

The analysis is geometry-only and read-only. It compares each mutant PDB with
the corresponding parent anchor PDB already present in Stage-1 outputs.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


DEFAULT_TRIAGE = Path(".tasks/active/initial-design-generation/manual_spot_check_panel_assisted_triage.csv")
DEFAULT_STAGE1 = Path("results/initial_design_generation/tier2_staged/full_stage1/tier2b_full_stage1_results.csv")
DEFAULT_PKA = Path("results/initial_design_generation/tier2_staged/full_stage1/tier2b_full_pka_detail.csv")
DEFAULT_OUT_DIR = Path("results/initial_design_generation/tier2_stage1_diagnostics")

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
}

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90}
HBOND_ATOMS = {"N", "O", "S"}
ACIDIC_OXYGENS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}

ANTIBODY_CHAINS = {"Ab_1E62": {"H", "L"}, "Ab_sdAb": {"A"}}
ANTIGEN_CHAINS = {"Ab_1E62": {"C"}, "Ab_sdAb": {"B"}}
WINDOW_CHAIN = {"VH": "H", "VL": "L", "VHH": "A"}
CDR_RANGES = {
    "Ab_1E62": [("L", 24, 38)],
    "Ab_sdAb": [("A", 95, 111)],
}
GLYCAN_GUARDRAIL_RESID = 146


@dataclass(frozen=True)
class ResKey:
    chain: str
    resid: int
    icode: str = ""

    def label(self) -> str:
        return f"{self.chain}{self.resid}{self.icode}".strip()


@dataclass(frozen=True)
class Atom:
    index: int
    name: str
    resname: str
    key: ResKey
    element: str
    coord: np.ndarray


@dataclass(frozen=True)
class Mutation:
    chain: str
    wt: str
    resid: int
    mut: str

    @property
    def key(self) -> ResKey:
        return ResKey(self.chain, self.resid)

    def label(self) -> str:
        return f"{self.chain}{self.wt}{self.resid}{self.mut}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--triage", default=str(DEFAULT_TRIAGE))
    parser.add_argument("--stage1", default=str(DEFAULT_STAGE1))
    parser.add_argument("--pka-detail", default=str(DEFAULT_PKA))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    return parser.parse_args()


def clean_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    text = str(value).strip()
    if not text:
        return None
    try:
        val = float(text)
    except ValueError:
        return None
    if math.isnan(val):
        return None
    return val


def infer_element(atom_name: str, element_field: str) -> str:
    element = element_field.strip().upper()
    if element:
        return element
    stripped = re.sub(r"^[0-9]+", "", atom_name.strip()).upper()
    if stripped.startswith("CL"):
        return "CL"
    if stripped.startswith("SE"):
        return "SE"
    return stripped[:1]


def parse_pdb_atoms(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip()
            resid_text = line[22:26].strip()
            if not chain or not resid_text:
                continue
            element = infer_element(atom_name, line[76:78] if len(line) >= 78 else "")
            if element == "H" or atom_name.upper().startswith(("H", "1H", "2H", "3H")):
                continue
            try:
                resid = int(resid_text)
                coord = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                )
            except ValueError:
                continue
            icode = line[26].strip()
            atoms.append(
                Atom(
                    index=len(atoms),
                    name=atom_name,
                    resname=resname,
                    key=ResKey(chain, resid, icode),
                    element=element,
                    coord=coord,
                )
            )
    return atoms


def atoms_by_residue(atoms: Iterable[Atom]) -> dict[ResKey, list[Atom]]:
    by_res: dict[ResKey, list[Atom]] = defaultdict(list)
    for atom in atoms:
        by_res[atom.key].append(atom)
    return dict(by_res)


def ca_coords(atoms: Iterable[Atom]) -> dict[ResKey, np.ndarray]:
    return {atom.key: atom.coord for atom in atoms if atom.name == "CA"}


def parse_mutations(text: object) -> list[Mutation]:
    if text is None or pd.isna(text):
        return []
    muts: list[Mutation] = []
    for part in re.split(r"[;,| ]+", str(text).strip()):
        if not part:
            continue
        match = re.fullmatch(r"([A-Za-z])([A-Z])(\d+)([A-Z])", part)
        if not match:
            continue
        chain, wt, resid, mut = match.groups()
        muts.append(Mutation(chain=chain, wt=wt, resid=int(resid), mut=mut))
    return muts


def parse_window(window: object, target: str) -> tuple[str, int, int] | None:
    if window is None or pd.isna(window) or not str(window).strip():
        return ("L", 1, 40) if target == "Ab_1E62" else ("A", 72, 111)
    match = re.search(r"_(VH|VL|VHH)_(\d+)_(\d+)", str(window))
    if not match:
        return ("L", 1, 40) if target == "Ab_1E62" else ("A", 72, 111)
    region, start, end = match.groups()
    chain = WINDOW_CHAIN.get(region)
    if not chain:
        return ("L", 1, 40) if target == "Ab_1E62" else ("A", 72, 111)
    if target == "Ab_sdAb":
        chain = "A"
    return chain, int(start), int(end)


def kabsch_transform(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    mob = mobile - mobile_center
    tar = target - target_center
    cov = mob.T @ tar
    u, _, vt = np.linalg.svd(cov)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    return rot, target_center - mobile_center @ rot


def apply_transform(coord: np.ndarray, rot: np.ndarray, trans: np.ndarray) -> np.ndarray:
    return coord @ rot + trans


def rmsd_for_keys(
    keys: Iterable[ResKey],
    parent_ca: dict[ResKey, np.ndarray],
    variant_ca: dict[ResKey, np.ndarray],
    rot: np.ndarray,
    trans: np.ndarray,
) -> float | None:
    sq: list[float] = []
    for key in keys:
        if key not in parent_ca or key not in variant_ca:
            continue
        diff = apply_transform(variant_ca[key], rot, trans) - parent_ca[key]
        sq.append(float(diff @ diff))
    if not sq:
        return None
    return math.sqrt(sum(sq) / len(sq))


def distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def vdw_radius(atom: Atom) -> float:
    return VDW.get(atom.element.upper(), 1.70)


def is_covalent_like(a: Atom, b: Atom) -> bool:
    if a.key == b.key:
        return True
    if a.key.chain == b.key.chain and abs(a.key.resid - b.key.resid) <= 1:
        return True
    if a.resname == "CYS" and b.resname == "CYS" and a.name == "SG" and b.name == "SG":
        d = distance(a.coord, b.coord)
        if 1.7 <= d <= 2.7:
            return True
    return False


def clash_pair_id(a: Atom, b: Atom) -> tuple[tuple[str, int, str, str], tuple[str, int, str, str]]:
    left = (a.key.chain, a.key.resid, a.key.icode, a.name)
    right = (b.key.chain, b.key.resid, b.key.icode, b.name)
    return tuple(sorted([left, right]))  # type: ignore[return-value]


def clash_map(atoms: list[Atom], overlap_threshold: float = 0.4) -> dict[tuple[tuple[str, int, str, str], tuple[str, int, str, str]], float]:
    if len(atoms) < 2:
        return {}
    coords = np.array([a.coord for a in atoms])
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=3.6)
    out: dict[tuple[tuple[str, int, str, str], tuple[str, int, str, str]], float] = {}
    for i, j in pairs:
        a = atoms[i]
        b = atoms[j]
        if is_covalent_like(a, b):
            continue
        d = distance(a.coord, b.coord)
        overlap = vdw_radius(a) + vdw_radius(b) - d
        if overlap <= overlap_threshold:
            continue
        out[clash_pair_id(a, b)] = overlap
    return out


def clash_stats(
    atoms: list[Atom],
    shell_keys: set[ResKey],
    parent_clashes: dict[tuple[tuple[str, int, str, str], tuple[str, int, str, str]], float] | None = None,
    overlap_threshold: float = 0.4,
) -> tuple[int, int, float, int, int, float]:
    clashes = clash_map(atoms, overlap_threshold=overlap_threshold)
    total = len(clashes)
    shell = 0
    max_overlap = max(clashes.values()) if clashes else 0.0
    atom_shell: dict[tuple[str, int, str, str], bool] = {}
    for atom in atoms:
        atom_shell[(atom.key.chain, atom.key.resid, atom.key.icode, atom.name)] = atom.key in shell_keys
    for left, right in clashes:
        if atom_shell.get(left, False) or atom_shell.get(right, False):
            shell += 1

    parent_clashes = parent_clashes or {}
    new_items = {k: v for k, v in clashes.items() if k not in parent_clashes}
    new_total = len(new_items)
    new_shell = 0
    for left, right in new_items:
        if atom_shell.get(left, False) or atom_shell.get(right, False):
            new_shell += 1
    max_new_overlap = max(new_items.values()) if new_items else 0.0
    return total, shell, round(max_overlap, 4), new_total, new_shell, round(max_new_overlap, 4)


def min_distance_between(atom_group_a: list[Atom], atom_group_b: list[Atom]) -> float | None:
    if not atom_group_a or not atom_group_b:
        return None
    coords_b = np.array([a.coord for a in atom_group_b])
    tree = cKDTree(coords_b)
    dists, _ = tree.query(np.array([a.coord for a in atom_group_a]), k=1)
    return round(float(np.min(dists)), 4)


def residue_shell(atoms: list[Atom], mutation_keys: set[ResKey], radius: float = 6.0) -> set[ResKey]:
    if not mutation_keys:
        return set()
    mut_atoms = [a for a in atoms if a.key in mutation_keys]
    if not mut_atoms:
        return set()
    coords = np.array([a.coord for a in atoms])
    tree = cKDTree(coords)
    shell: set[ResKey] = set(mutation_keys)
    for atom in mut_atoms:
        for idx in tree.query_ball_point(atom.coord, radius):
            shell.add(atoms[idx].key)
    return shell


def contact_map(
    atoms: list[Atom],
    antibody_chains: set[str],
    antigen_chains: set[str],
    cutoff: float = 4.5,
) -> dict[tuple[ResKey, ResKey], tuple[float, float]]:
    ab_atoms = [a for a in atoms if a.key.chain in antibody_chains]
    ag_atoms = [a for a in atoms if a.key.chain in antigen_chains]
    if not ab_atoms or not ag_atoms:
        return {}
    ag_coords = np.array([a.coord for a in ag_atoms])
    tree = cKDTree(ag_coords)
    contacts: dict[tuple[ResKey, ResKey], tuple[float, float]] = {}
    for a in ab_atoms:
        for idx in tree.query_ball_point(a.coord, cutoff):
            b = ag_atoms[idx]
            d = distance(a.coord, b.coord)
            overlap = vdw_radius(a) + vdw_radius(b) - d
            key = (a.key, b.key)
            if key not in contacts:
                contacts[key] = (d, overlap)
            else:
                min_d, max_overlap = contacts[key]
                contacts[key] = (min(min_d, d), max(max_overlap, overlap))
    return contacts


def count_hbond_and_saltbridge(his_atoms: list[Atom], atoms: list[Atom]) -> tuple[int, int]:
    if not his_atoms:
        return 0, 0
    his_nd_ne = [a for a in his_atoms if a.name in {"ND1", "NE2"}]
    if not his_nd_ne:
        return 0, 0
    other_atoms = [
        a
        for a in atoms
        if a.element in HBOND_ATOMS and all(a.key != h.key for h in his_nd_ne)
    ]
    if not other_atoms:
        return 0, 0
    coords = np.array([a.coord for a in other_atoms])
    tree = cKDTree(coords)
    hbond_pairs: set[tuple[int, int]] = set()
    salt_pairs: set[tuple[int, int]] = set()
    for h in his_nd_ne:
        for idx in tree.query_ball_point(h.coord, 3.5):
            other = other_atoms[idx]
            hbond_pairs.add((h.index, other.index))
        for idx in tree.query_ball_point(h.coord, 4.0):
            other = other_atoms[idx]
            if (other.resname, other.name) in ACIDIC_OXYGENS:
                salt_pairs.add((h.index, other.index))
    return len(hbond_pairs), len(salt_pairs)


def sphere_points(n_points: int = 96) -> np.ndarray:
    points = []
    increment = math.pi * (3.0 - math.sqrt(5.0))
    offset = 2.0 / n_points
    for i in range(n_points):
        y = (i * offset - 1) + (offset / 2)
        r = math.sqrt(max(0.0, 1 - y * y))
        phi = i * increment
        points.append([math.cos(phi) * r, y, math.sin(phi) * r])
    return np.array(points, dtype=float)


SPHERE_POINTS = sphere_points()


def residue_sasa_for_keys(
    atoms: list[Atom],
    residue_keys: set[ResKey],
    probe_radius: float = 1.4,
) -> dict[ResKey, float]:
    """Approximate residue SASA for selected residues only.

    This avoids whole-structure SASA on hydrogen-rich Rosetta PDBs. The values
    are intended for manual triage bins, not for thermodynamic interpretation.
    """

    if not residue_keys:
        return {}
    coords = np.array([a.coord for a in atoms])
    tree = cKDTree(coords)
    max_query = max(VDW.values()) + probe_radius + max(VDW.values()) + probe_radius
    sasa: dict[ResKey, float] = defaultdict(float)
    for atom in atoms:
        if atom.key not in residue_keys:
            continue
        radius = vdw_radius(atom) + probe_radius
        surface = atom.coord + SPHERE_POINTS * radius
        accessible = 0
        for point in surface:
            occluded = False
            for idx in tree.query_ball_point(point, max_query):
                other = atoms[idx]
                if other.index == atom.index:
                    continue
                other_radius = vdw_radius(other) + probe_radius
                if distance(point, other.coord) < other_radius:
                    occluded = True
                    break
            if not occluded:
                accessible += 1
        sasa[atom.key] += 4.0 * math.pi * radius * radius * accessible / len(SPHERE_POINTS)
    return dict(sasa)


def format_float(value: float | None, digits: int = 4) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def load_pka_detail(path: Path) -> dict[str, dict[ResKey, str]]:
    if not path.exists():
        return {}
    detail: dict[str, dict[ResKey, str]] = defaultdict(dict)
    with path.open() as handle:
        for row in csv.DictReader(handle):
            variant_id = row.get("variant_id", "")
            chain = row.get("chain", "")
            resid = row.get("resid", "")
            if not variant_id or not chain or not resid:
                continue
            propka = clean_float(row.get("pKa_propka"))
            pkai = clean_float(row.get("pKa_pkai"))
            vals = [v for v in [propka, pkai] if v is not None]
            combined = sum(vals) / len(vals) if vals else None
            label = f"{chain}{resid}:propka={format_float(propka, 2)},pkai={format_float(pkai, 2)},mean={format_float(combined, 2)}"
            detail[variant_id][ResKey(chain, int(float(resid)))] = label
    return detail


def parent_rows(stage1: pd.DataFrame) -> dict[str, pd.Series]:
    parents: dict[str, pd.Series] = {}
    for _, row in stage1.iterrows():
        if str(row.get("control_type", "")) == "parent_wt_anchor":
            parents[str(row["target"])] = row
    return parents


def structure_feature_row(
    triage_row: pd.Series,
    stage_row: pd.Series,
    parent_row: pd.Series,
    parent_cache: dict[str, dict[str, object]],
    pka_detail: dict[str, dict[ResKey, str]],
) -> dict[str, object]:
    variant_id = str(triage_row["variant_id"])
    target = str(stage_row["target"])
    status: list[str] = []

    pdb_path = Path(str(stage_row.get("pdb_path", "")))
    parent_pdb = Path(str(parent_row.get("pdb_path", "")))
    mutation_list = str(stage_row.get("mutation_list", "") or "")
    mutations = parse_mutations(mutation_list)
    mutation_keys = {m.key for m in mutations}
    his_muts = [m for m in mutations if m.mut == "H"]
    his_keys = {m.key for m in his_muts}

    base = {
        "variant_id": variant_id,
        "target": target,
        "mutation_list": mutation_list,
        "t2_class_current": str(triage_row.get("tier2_class", stage_row.get("tier2_class", ""))),
        "rosetta_delta": stage_row.get("rosetta_delta_score", ""),
        "local_validity_score": stage_row.get("local_structure_validity_t2_score", ""),
        "his_count": len(his_muts),
        "his_positions": ";".join(m.key.label() for m in his_muts),
        "his_pka_if_available": ";".join(pka_detail.get(variant_id, {}).get(m.key, "") for m in his_muts).strip(";"),
    }

    if not pdb_path.exists():
        base.update({"structure_parse_status": "missing_variant_pdb"})
        return base
    if not parent_pdb.exists():
        base.update({"structure_parse_status": "missing_parent_pdb"})
        return base

    cache_key = target
    if cache_key not in parent_cache:
        parent_atoms = parse_pdb_atoms(parent_pdb)
        parent_ca = ca_coords(parent_atoms)
        parent_contacts = contact_map(parent_atoms, ANTIBODY_CHAINS[target], ANTIGEN_CHAINS[target])
        parent_clashes = clash_map(parent_atoms)
        parent_cache[cache_key] = {
            "atoms": parent_atoms,
            "ca": parent_ca,
            "contacts": parent_contacts,
            "clashes": parent_clashes,
            "pdb": parent_pdb,
        }

    parent_atoms = parent_cache[cache_key]["atoms"]  # type: ignore[assignment]
    parent_ca = parent_cache[cache_key]["ca"]  # type: ignore[assignment]
    parent_contacts = parent_cache[cache_key]["contacts"]  # type: ignore[assignment]
    parent_clashes = parent_cache[cache_key]["clashes"]  # type: ignore[assignment]

    try:
        variant_atoms = parse_pdb_atoms(pdb_path)
    except Exception as exc:  # pragma: no cover - defensive status capture
        base.update({"structure_parse_status": f"variant_parse_failed:{type(exc).__name__}"})
        return base

    variant_ca = ca_coords(variant_atoms)
    shared_ca = sorted(set(parent_ca) & set(variant_ca), key=lambda k: (k.chain, k.resid, k.icode))
    if len(shared_ca) >= 3:
        mobile = np.array([variant_ca[k] for k in shared_ca])
        target_coords = np.array([parent_ca[k] for k in shared_ca])
        rot, trans = kabsch_transform(mobile, target_coords)
    else:
        rot, trans = np.eye(3), np.zeros(3)
        status.append("insufficient_ca_for_alignment")

    by_res = atoms_by_residue(variant_atoms)
    parent_by_res = atoms_by_residue(parent_atoms)
    missing_mut_residues = [m.key.label() for m in mutations if m.key not in by_res]
    if missing_mut_residues:
        status.append("missing_mutation_residue:" + ";".join(missing_mut_residues))
    if not mutations:
        status.append("no_mutations")

    shell_keys = residue_shell(variant_atoms, mutation_keys)
    total_clashes, shell_clashes, max_overlap, new_clashes, new_shell_clashes, max_new_overlap = clash_stats(
        variant_atoms,
        shell_keys,
        parent_clashes=parent_clashes,
    )

    antigen_atoms = [a for a in variant_atoms if a.key.chain in ANTIGEN_CHAINS[target]]
    mutation_atoms = [a for a in variant_atoms if a.key in mutation_keys]
    min_mut_ag = min_distance_between(mutation_atoms, antigen_atoms)

    variant_contacts = contact_map(variant_atoms, ANTIBODY_CHAINS[target], ANTIGEN_CHAINS[target])
    parent_contact_keys = set(parent_contacts.keys())
    variant_contact_keys = set(variant_contacts.keys())
    new_contacts = variant_contact_keys - parent_contact_keys
    lost_contacts = parent_contact_keys - variant_contact_keys
    new_bad = 0
    for key in new_contacts:
        min_d, max_contact_overlap = variant_contacts[key]
        if min_d < 2.8 or max_contact_overlap > 0.4:
            new_bad += 1

    cdr_keys = {
        ResKey(chain, resid)
        for chain, start, end in CDR_RANGES.get(target, [])
        for resid in range(start, end + 1)
    }
    window = parse_window(stage_row.get("window", ""), target)
    window_keys: set[ResKey] = set()
    if window:
        chain, start, end = window
        window_keys = {ResKey(chain, resid) for resid in range(start, end + 1)}
    else:
        status.append("window_parse_failed")

    cdr_rmsd = rmsd_for_keys(cdr_keys, parent_ca, variant_ca, rot, trans)
    window_rmsd = rmsd_for_keys(window_keys, parent_ca, variant_ca, rot, trans)
    shell_rmsd = rmsd_for_keys(shell_keys, parent_ca, variant_ca, rot, trans)

    his_atoms = [a for a in variant_atoms if a.key in his_keys]
    his_min_ag = min_distance_between(his_atoms, antigen_atoms)
    hbond_count, saltbridge_count = count_hbond_and_saltbridge(his_atoms, variant_atoms)

    his_sasa_total: float | None = None
    his_sasa_by_pos = ""
    try:
        sasa = residue_sasa_for_keys(variant_atoms, his_keys)
        parts = []
        total = 0.0
        for key in sorted(his_keys, key=lambda k: (k.chain, k.resid, k.icode)):
            val = sasa.get(key)
            if val is not None:
                total += val
                parts.append(f"{key.label()}:{val:.2f}")
        his_sasa_total = total if parts else None
        his_sasa_by_pos = ";".join(parts)
    except Exception as exc:  # pragma: no cover - defensive status capture
        status.append(f"sasa_failed:{type(exc).__name__}")

    glycan_distance = None
    glycan_atoms = [
        a
        for a in antigen_atoms
        if a.key.resid == GLYCAN_GUARDRAIL_RESID
    ]
    if glycan_atoms and mutation_atoms:
        glycan_distance = min_distance_between(mutation_atoms, glycan_atoms)
    elif not glycan_atoms:
        status.append("glycan_guardrail_residue_missing")

    base.update(
        {
            "clash_count_total": total_clashes,
            "clash_count_mutation_shell": shell_clashes,
            "max_clash_overlap": max_overlap,
            "new_clash_count_total": new_clashes,
            "new_clash_count_mutation_shell": new_shell_clashes,
            "max_new_clash_overlap": max_new_overlap,
            "min_mutation_to_antigen_distance": format_float(min_mut_ag),
            "interface_contact_change_count": len(new_contacts) + len(lost_contacts),
            "cdr_rmsd": format_float(cdr_rmsd),
            "window_rmsd": format_float(window_rmsd),
            "mutation_shell_rmsd": format_float(shell_rmsd),
            "his_min_antigen_distance": format_float(his_min_ag),
            "his_sasa": format_float(his_sasa_total, 2),
            "his_sasa_by_position": his_sasa_by_pos,
            "his_hbond_count": hbond_count,
            "his_saltbridge_count": saltbridge_count,
            "new_bad_contact_count": new_bad,
            "lost_parent_contact_count": len(lost_contacts),
            "glycan_risk_distance_if_available": format_float(glycan_distance),
            "structure_parse_status": "ok" if not status else "ok;" + "|".join(status),
        }
    )
    return base


def summarize(features: pd.DataFrame) -> list[str]:
    lines: list[str] = []
    lines.append("# Manual Spot-Check Structure Feature Review")
    lines.append("")
    lines.append("This report summarizes geometry-derived features from the 134-row manual spot-check panel.")
    lines.append("All metrics are computed from existing Stage-1 PyRosetta PDBs and parent-anchor PDBs; no new structure modeling was run.")
    lines.append("His SASA is a selected-residue heavy-atom approximation for triage, not a full thermodynamic SASA calculation.")
    lines.append("")

    def num(col: str) -> pd.Series:
        return pd.to_numeric(features[col], errors="coerce")

    lines.append("## Coverage")
    lines.append("")
    lines.append(f"- Rows analyzed: {len(features)}")
    lines.append(f"- Parse status ok-like rows: {int(features['structure_parse_status'].astype(str).str.startswith('ok').sum())}")
    lines.append(f"- Rows with designed His: {int(num('his_count').gt(0).sum())}")
    lines.append("")

    lines.append("## Question 1: Are clashes concentrated near mutation neighborhoods?")
    lines.append("")
    total_clash = num("clash_count_total")
    shell_clash = num("clash_count_mutation_shell")
    new_clash = num("new_clash_count_total")
    new_shell_clash = num("new_clash_count_mutation_shell")
    with_clash = features[total_clash.gt(0)].copy()
    if not with_clash.empty:
        ratio = pd.to_numeric(with_clash["clash_count_mutation_shell"], errors="coerce") / pd.to_numeric(
            with_clash["clash_count_total"], errors="coerce"
        )
        lines.append(f"- Rows with any geometry clash: {len(with_clash)} / {len(features)}")
        lines.append(f"- Median mutation-shell clash fraction among clash rows: {ratio.median():.3f}")
        lines.append(f"- Rows with >=50% clashes in mutation shell: {int(ratio.ge(0.5).sum())} / {len(with_clash)}")
        lines.append(f"- Rows with >=80% clashes in mutation shell: {int(ratio.ge(0.8).sum())} / {len(with_clash)}")
    else:
        lines.append("- No rows had geometry clashes under the heavy-atom overlap rule.")
    with_new_clash = features[new_clash.gt(0)].copy()
    if not with_new_clash.empty:
        new_ratio = pd.to_numeric(with_new_clash["new_clash_count_mutation_shell"], errors="coerce") / pd.to_numeric(
            with_new_clash["new_clash_count_total"], errors="coerce"
        )
        lines.append(f"- Rows with parent-relative new clashes: {len(with_new_clash)} / {len(features)}")
        lines.append(f"- Median new-clash mutation-shell fraction: {new_ratio.median():.3f}")
        lines.append(f"- Rows with >=50% new clashes in mutation shell: {int(new_ratio.ge(0.5).sum())} / {len(with_new_clash)}")
        lines.append(f"- Rows with >=80% new clashes in mutation shell: {int(new_ratio.ge(0.8).sum())} / {len(with_new_clash)}")
    else:
        lines.append("- No parent-relative new clashes were detected.")
    lines.append("")

    lines.append("## Question 2: Does high Rosetta delta coincide with low local validity?")
    lines.append("")
    delta = num("rosetta_delta")
    local = num("local_validity_score")
    high_delta = delta.gt(600)
    if high_delta.any():
        corr = delta.corr(local)
        lines.append(f"- Rosetta delta > 600 rows: {int(high_delta.sum())}")
        lines.append(f"- Among high-delta rows, local validity <0.35: {int((high_delta & local.lt(0.35)).sum())}")
        lines.append(f"- Among high-delta rows, local validity <0.55: {int((high_delta & local.lt(0.55)).sum())}")
        lines.append(f"- Pearson correlation between Rosetta delta and local validity: {corr:.3f}")
    else:
        lines.append("- No rows had Rosetta delta > 600.")
    lines.append("")

    lines.append("## Question 3: Are designed His residues near interface/regulatory positions?")
    lines.append("")
    his_rows = features[num("his_count").gt(0)].copy()
    his_dist = pd.to_numeric(his_rows["his_min_antigen_distance"], errors="coerce")
    his_sasa = pd.to_numeric(his_rows["his_sasa"], errors="coerce")
    if not his_rows.empty:
        lines.append(f"- His-containing rows: {len(his_rows)}")
        lines.append(f"- His within 5 A of antigen: {int(his_dist.le(5.0).sum())}")
        lines.append(f"- His within 5-10 A of antigen: {int(his_dist.gt(5.0).mul(his_dist.le(10.0)).sum())}")
        lines.append(f"- His farther than 10 A from antigen: {int(his_dist.gt(10.0).sum())}")
        lines.append(f"- His SASA median: {his_sasa.median():.2f} A^2")
        lines.append(f"- His with very low SASA (<10 A^2 total): {int(his_sasa.lt(10).sum())}")
        lines.append(f"- His with high SASA (>80 A^2 total): {int(his_sasa.gt(80).sum())}")
        lines.append(f"- Rows with at least one His hbond-like contact: {int(pd.to_numeric(his_rows['his_hbond_count'], errors='coerce').gt(0).sum())}")
        lines.append(f"- Rows with at least one His saltbridge-like contact: {int(pd.to_numeric(his_rows['his_saltbridge_count'], errors='coerce').gt(0).sum())}")
    else:
        lines.append("- No designed His rows were present.")
    lines.append("")

    lines.append("## Question 4: Are CDR/window RMSDs large?")
    lines.append("")
    for col in ["cdr_rmsd", "window_rmsd", "mutation_shell_rmsd"]:
        vals = num(col)
        lines.append(
            f"- {col}: median={vals.median():.4f} A, max={vals.max():.4f} A, >0.5 A={int(vals.gt(0.5).sum())}, >1.0 A={int(vals.gt(1.0).sum())}"
        )
    lines.append("")

    lines.append("## Question 5: Are parent interface contacts disrupted?")
    lines.append("")
    lost = num("lost_parent_contact_count")
    changed = num("interface_contact_change_count")
    new_bad = num("new_bad_contact_count")
    lines.append(f"- Rows losing at least one parent antibody-antigen contact: {int(lost.gt(0).sum())} / {len(features)}")
    lines.append(f"- Median lost parent contact count: {lost.median():.1f}; max: {int(lost.max()) if lost.notna().any() else 0}")
    lines.append(f"- Median total interface contact change count: {changed.median():.1f}; max: {int(changed.max()) if changed.notna().any() else 0}")
    lines.append(f"- Rows with new bad interface contacts: {int(new_bad.gt(0).sum())} / {len(features)}")
    lines.append("")

    lines.append("## By-target quick view")
    lines.append("")
    lines.append("| target | rows | high_delta_gt600 | low_local_lt0.35 | any_lost_parent_contact | his_within_5A | cdr_rmsd_gt0.5 |")
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for target, sub in features.groupby("target", dropna=False):
        sub_delta = pd.to_numeric(sub["rosetta_delta"], errors="coerce")
        sub_local = pd.to_numeric(sub["local_validity_score"], errors="coerce")
        sub_lost = pd.to_numeric(sub["lost_parent_contact_count"], errors="coerce")
        sub_his_dist = pd.to_numeric(sub["his_min_antigen_distance"], errors="coerce")
        sub_cdr = pd.to_numeric(sub["cdr_rmsd"], errors="coerce")
        lines.append(
            f"| {target} | {len(sub)} | {int(sub_delta.gt(600).sum())} | {int(sub_local.lt(0.35).sum())} | "
            f"{int(sub_lost.gt(0).sum())} | {int(sub_his_dist.le(5.0).sum())} | {int(sub_cdr.gt(0.5).sum())} |"
        )
    lines.append("")
    return lines


def main() -> None:
    args = parse_args()
    triage_path = Path(args.triage)
    stage1_path = Path(args.stage1)
    pka_path = Path(args.pka_detail)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    triage = pd.read_csv(triage_path, low_memory=False)
    stage1 = pd.read_csv(stage1_path, low_memory=False)
    stage1_by_variant = stage1.set_index("variant_id", drop=False)
    parents = parent_rows(stage1)
    pka_detail = load_pka_detail(pka_path)

    rows: list[dict[str, object]] = []
    parent_cache: dict[str, dict[str, object]] = {}
    for _, triage_row in triage.iterrows():
        variant_id = str(triage_row["variant_id"])
        if variant_id not in stage1_by_variant.index:
            rows.append(
                {
                    "variant_id": variant_id,
                    "target": triage_row.get("target", ""),
                    "structure_parse_status": "missing_stage1_row",
                }
            )
            continue
        stage_row = stage1_by_variant.loc[variant_id]
        target = str(stage_row["target"])
        parent_row = parents.get(target)
        if parent_row is None:
            rows.append(
                {
                    "variant_id": variant_id,
                    "target": target,
                    "structure_parse_status": "missing_parent_anchor_row",
                }
            )
            continue
        rows.append(structure_feature_row(triage_row, stage_row, parent_row, parent_cache, pka_detail))

    features = pd.DataFrame(rows)
    required_order = [
        "variant_id",
        "target",
        "mutation_list",
        "t2_class_current",
        "rosetta_delta",
        "local_validity_score",
        "clash_count_total",
        "clash_count_mutation_shell",
        "max_clash_overlap",
        "new_clash_count_total",
        "new_clash_count_mutation_shell",
        "max_new_clash_overlap",
        "min_mutation_to_antigen_distance",
        "interface_contact_change_count",
        "cdr_rmsd",
        "window_rmsd",
        "mutation_shell_rmsd",
        "his_count",
        "his_positions",
        "his_min_antigen_distance",
        "his_sasa",
        "his_hbond_count",
        "his_saltbridge_count",
        "his_pka_if_available",
        "new_bad_contact_count",
        "lost_parent_contact_count",
        "glycan_risk_distance_if_available",
        "structure_parse_status",
    ]
    extra = [c for c in features.columns if c not in required_order]
    features = features[[c for c in required_order if c in features.columns] + extra]

    feature_path = out_dir / "manual_spot_check_structure_features.csv"
    report_path = out_dir / "manual_spot_check_structure_feature_review.md"
    features.to_csv(feature_path, index=False)
    report_path.write_text("\n".join(summarize(features)) + "\n", encoding="utf-8")

    print(f"Wrote {feature_path} ({len(features)} rows)")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
