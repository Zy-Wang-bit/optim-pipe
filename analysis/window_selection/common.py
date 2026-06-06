from __future__ import annotations

import argparse
import json
import math
import re
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from Bio import SeqIO


REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIG_DIR = PACKAGE_ROOT / "config"
SCHEMA_DIR = PACKAGE_ROOT / "schemas"


def repo_path(path: str | Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else REPO_ROOT / path


def load_yaml(path: str | Path) -> dict[str, Any]:
    with repo_path(path).open() as handle:
        data = yaml.safe_load(handle) or {}
    return data


def load_inputs() -> dict[str, Any]:
    return load_yaml(CONFIG_DIR / "inputs.yaml")


def load_decision_rules() -> dict[str, Any]:
    return load_yaml(CONFIG_DIR / "decision_rules.yaml")


def load_capacity_config() -> dict[str, Any]:
    return load_yaml(CONFIG_DIR / "library_capacity_config.yaml")


def load_source_policy() -> dict[str, Any]:
    return load_yaml(CONFIG_DIR / "source_policy.yaml")


def output_root() -> Path:
    return repo_path(load_inputs()["output_root"])


def tables_dir() -> Path:
    path = output_root() / "tables"
    path.mkdir(parents=True, exist_ok=True)
    return path


def af3_input_dir(target: str) -> Path:
    path = output_root() / "structures" / "af3" / target / "input"
    path.mkdir(parents=True, exist_ok=True)
    return path


def af3_output_dir(target: str) -> Path:
    path = output_root() / "structures" / "af3" / target / "output"
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_fasta(path: str | Path, record_id: str | None = None) -> tuple[str, str, str]:
    records = list(SeqIO.parse(str(repo_path(path)), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    chosen = None
    if record_id is None:
        chosen = records[0]
    else:
        for record in records:
            if record.id == record_id or record.description.split()[0] == record_id:
                chosen = record
                break
    if chosen is None:
        raise ValueError(f"Record {record_id!r} not found in {path}")
    seq = str(chosen.seq).strip()
    if seq.endswith("*"):
        seq = seq[:-1]
    return chosen.id, chosen.description, seq


def cdr_region(chain: str, pos: int, target: str) -> str:
    if target == "1E62" and chain == "L":
        if 24 <= pos <= 34:
            return "CDR1"
        if 50 <= pos <= 56:
            return "CDR2"
        if 89 <= pos <= 97:
            return "CDR3"
        return "FR"
    if chain in {"H", "A"}:
        if 26 <= pos <= 35:
            return "CDR1"
        if 50 <= pos <= 65:
            return "CDR2"
        if 95 <= pos <= 113:
            return "CDR3"
        return "FR"
    return "unknown"


def diff_mutations(wt: str, seq: str, chain: str) -> list[str]:
    if len(wt) != len(seq):
        return [f"{chain}:length_{len(wt)}->{len(seq)}"]
    mutations = []
    for i, (a, b) in enumerate(zip(wt, seq, strict=True), start=1):
        if a != b:
            mutations.append(f"{chain}{a}{i}{b}")
    return mutations


def mutation_positions(mutations: list[str]) -> list[int]:
    out = []
    for mut in mutations:
        match = re.search(r"(\d+)", mut)
        if match:
            out.append(int(match.group(1)))
    return out


def split_list(value: str | float | int | None) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in re.split(r"[;,]", text) if item.strip()]


def parse_constraint_positions(row: pd.Series) -> list[int]:
    values: list[int] = []
    pos = row.get("position")
    if pd.notna(pos) and str(pos).strip():
        try:
            values.append(int(float(pos)))
        except ValueError:
            pass
    pair = row.get("position_pair")
    if pd.notna(pair):
        values.extend(int(x) for x in re.findall(r"(\d+)", str(pair)))
    return sorted(set(values))


def write_csv(df: pd.DataFrame, filename: str) -> Path:
    path = tables_dir() / filename
    df.to_csv(path, index=False)
    return path


def read_table(filename: str) -> pd.DataFrame:
    return pd.read_csv(tables_dir() / filename, dtype=str, keep_default_na=False)


def json_dump(data: Any, path: str | Path) -> None:
    path = repo_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(data, handle, indent=2)


def make_arg_parser(description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--output-root",
        default=str(output_root()),
        help="Output root; config value is used by default.",
    )
    return parser


def site_option_count(classes: str, wt_aa: str, capacity_config: dict[str, Any]) -> int:
    class_names = [x for x in str(classes).split("|") if x]
    residues = set()
    for cls in class_names:
        residues.update(capacity_config["mutation_classes"].get(cls, []))
    residues.discard(wt_aa)
    return len(residues)


def count_variant_space(option_counts: list[int], max_order: int) -> int:
    total = 0
    indexed = [c for c in option_counts if c > 0]
    for order in range(1, min(max_order, len(indexed)) + 1):
        for combo in combinations(indexed, order):
            product = 1
            for count in combo:
                product *= count
            total += product
    return int(total)


def ensure_required_paths(paths: list[str]) -> list[str]:
    missing = []
    for path in paths:
        if not repo_path(path).exists():
            missing.append(path)
    return missing
