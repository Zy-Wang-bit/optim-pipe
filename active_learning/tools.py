from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from collections.abc import Iterable
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from analysis.naming.convert import parse_unified


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CORE_TRAINING_COLUMNS = {
    "sequence_id",
    "sequence_sha256",
    "A7",
    "A6",
    "noise",
    "utility_recomputed",
    "experiment_status",
    "label_is_direct_measurement",
    "usable_for_training",
}
DIRECT_TRAINING_SOURCES = {"raw_elisa", "retrospective_aggregate", "dose_curve_aggregate", "wet_observation"}
DEFAULT_CONFIG: dict[str, Any] = {
    "naming": {
        "fasta_record_id": "A",
        "mutation_chain_label": "H",
        "ignore_chain_for_single_sequence": True,
    },
    "objective": {
        "name": "ph_switch",
        "formula": "A7 * (1 - A6)",
        "input_clip": True,
        "noise_lambda": 0.0,
        "normalization": "reference_clip01",
    },
    "candidate_pool": {
        "mode": "full_single",
        "excluded_positions": [],
        "excluded_mutations": [],
    },
    "embedding": {
        "embedding_model_name": "ESMC-600M",
        "embedding_model_version": "unknown",
        "embedding_layer": "final",
        "pooling": "mean",
        "embedding_dim": None,
    },
    "model": {
        "type": "random_forest_regressor",
        "n_estimators": 100,
        "criterion": "friedman_mse",
        "min_samples_leaf": 1,
        "bootstrap": True,
        "random_state": 42,
    },
    "acquisition": {
        "policy": "topk",
        "kappa": 0.0,
        "top_k": 24,
    },
}


# config / IO
def _merge_defaults(config: dict[str, Any], defaults: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(defaults)
    for key, value in config.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_defaults(value, merged[key])
        else:
            merged[key] = value
    return merged


def load_round_config(path: Path) -> dict[str, Any]:
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    raw = json.loads(text) if path.suffix.lower() == ".json" else yaml.safe_load(text)
    return _merge_defaults(raw or {}, DEFAULT_CONFIG)


def read_optional_csv(path: str | None) -> pd.DataFrame | None:
    return pd.read_csv(path) if path else None


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reject_output_overwrite(output: Path, inputs: list[Path | None]) -> None:
    output = Path(output)
    input_paths = {Path(path).resolve() for path in inputs if path is not None}
    if output.resolve() in input_paths:
        raise ValueError(f"Refusing to overwrite input file: {output}")


def reject_output_paths_overwrite_inputs(paths: dict[str, Path], inputs: list[Path]) -> None:
    input_paths = {path.resolve() for path in inputs}
    for label, output in paths.items():
        if output.resolve() in input_paths:
            raise ValueError(f"Refusing to overwrite input {label} CSV: {output}")


def _stable_hash(value: Any) -> str:
    return sha256_text(json.dumps(value, sort_keys=True, separators=(",", ":")))


def _git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return "unknown"
    return result.stdout.strip() or "unknown"


def build_round_provenance(config: dict[str, Any], paths: dict[str, Path]) -> dict[str, Any]:
    path_map = {key: Path(value) for key, value in paths.items()}
    for path in path_map.values():
        if not path.exists():
            raise FileNotFoundError(path)

    embedding = config.get("embedding", {})
    model = config.get("model", {})
    objective = config.get("objective", {})
    acquisition = config.get("acquisition", {})
    wt_path = path_map.get("wt_fasta")
    record_id = config.get("naming", {}).get("fasta_record_id")
    wt_sequence = parse_fasta_first_sequence(wt_path, record_id) if wt_path else ""

    return {
        "system_id": config.get("system_id"),
        "round_id": config.get("round_id"),
        "wt_fasta_path": str(config.get("wt_fasta_path", wt_path or "")),
        "wt_sequence_sha256": sha256_text(wt_sequence) if wt_sequence else None,
        "training_csv_sha256": sha256_file(path_map["training_csv"]) if "training_csv" in path_map else None,
        "candidate_pool_sha256": sha256_file(path_map["candidate_pool"]) if "candidate_pool" in path_map else None,
        "embedding_csv_sha256": sha256_file(path_map["embedding_csv"]) if "embedding_csv" in path_map else None,
        "embedding_model_name": embedding.get("embedding_model_name"),
        "embedding_model_version": embedding.get("embedding_model_version"),
        "embedding_layer": embedding.get("embedding_layer"),
        "embedding_pooling": embedding.get("pooling"),
        "embedding_dim": embedding.get("embedding_dim"),
        "objective_config": objective,
        "objective_config_hash": _stable_hash(objective),
        "model_config_hash": _stable_hash(model),
        "acquisition_config": acquisition,
        "random_state": model.get("random_state"),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "code_version_or_git_commit": _git_commit(),
    }


# mutation / sequence
def parse_fasta_first_sequence(path: Path, record_id: str | None = None) -> str:
    current_id = None
    current_lines: list[str] = []
    first_sequence: str | None = None

    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                sequence = "".join(current_lines)
                if record_id is None and first_sequence is None:
                    first_sequence = sequence
                if current_id == record_id:
                    return sequence
            current_id = line[1:].split()[0]
            current_lines = []
            continue
        current_lines.append(line)

    if current_id is not None:
        sequence = "".join(current_lines)
        if current_id == record_id:
            return sequence
        if record_id is None and first_sequence is None:
            return sequence
    if record_id is None and first_sequence is not None:
        return first_sequence
    raise ValueError(f"FASTA record not found: {record_id}")


def apply_mutation_string(
    sequence: str,
    mutation_string: str,
    expected_chain: str | None = None,
    ignore_chain_for_single_sequence: bool = True,
) -> str:
    tokens = [token.strip() for token in re.split(r"[;,]", str(mutation_string or "")) if token.strip()]
    if not tokens or tokens == ["WT"]:
        return sequence

    chars = list(sequence)
    for token in tokens:
        try:
            chain, orig, pos, new = parse_unified(token)
        except ValueError as exc:
            raise ValueError(f"invalid mutation token {token}") from exc
        if expected_chain and chain != expected_chain and not ignore_chain_for_single_sequence:
            raise ValueError(f"unexpected chain in mutation token {token}")
        index = pos - 1
        if index < 0 or index >= len(chars):
            raise ValueError(f"mutation position out of range in {token}")
        if chars[index] != orig:
            raise ValueError(f"original amino acid mismatch in {token}: expected {chars[index]}")
        chars[index] = new
    return "".join(chars)


def sequence_id_from_mutations(prefix: str, mutation_string: str) -> str:
    mutation_string = str(mutation_string or "").strip() or "WT"
    return f"{prefix}:{mutation_string}"


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


# labels
def clip01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def compute_switch_utility(a7: float, a6: float, noise: float = 0.0, lambda_noise: float = 0.0) -> float:
    return float(a7) * (1.0 - float(a6)) - float(lambda_noise) * float(noise)


def normalize_by_reference(value: float, blank: float, reference: float) -> float:
    denominator = float(reference) - float(blank)
    if denominator == 0:
        raise ValueError("reference must differ from blank")
    return clip01((float(value) - float(blank)) / denominator)


def _require_core_columns(rows: pd.DataFrame) -> None:
    missing = sorted(CORE_TRAINING_COLUMNS - set(rows.columns))
    if missing:
        raise ValueError(f"missing core training columns: {', '.join(missing)}")


def aggregate_training_replicates(rows: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["sequence_id", "sequence_sha256"]
    required = set(group_cols) | {"A7", "A6", "noise"}
    missing = sorted(required - set(rows.columns))
    if missing:
        raise ValueError(f"missing replicate aggregation columns: {', '.join(missing)}")

    records = []
    for _, group in rows.groupby(group_cols, sort=False, dropna=False):
        record = group.iloc[0].to_dict()
        a7 = float(group["A7"].mean())
        a6 = float(group["A6"].mean())
        noise = float(group["noise"].mean())
        record.update(
            {
                "A7": a7,
                "A6": a6,
                "A7_mean": a7,
                "A6_mean": a6,
                "noise": noise,
                "n_replicates": int(len(group)),
                "utility_recomputed": compute_switch_utility(a7, a6, noise),
            }
        )
        records.append(record)
    return pd.DataFrame(records)


def filter_training_rows(rows: pd.DataFrame) -> pd.DataFrame:
    _require_core_columns(rows)
    mask = (
        (rows["experiment_status"] == "measured")
        & (rows["label_is_direct_measurement"] == True)
        & (rows["usable_for_training"] == True)
    )
    if "label_source_type" in rows.columns:
        source = rows["label_source_type"].fillna("").astype(str)
        mask &= source.isin(DIRECT_TRAINING_SOURCES)
    return rows.loc[mask].copy()


# embeddings
def _embedding_columns(columns: list[str] | pd.Index) -> list[str]:
    emb_cols = [col for col in columns if re.match(r"^emb_\d+$", str(col))]
    if not emb_cols:
        raise ValueError("embedding table requires emb_<n> columns")
    return sorted(emb_cols, key=lambda col: int(str(col).split("_")[1]))


def load_embedding_table(path: Path) -> pd.DataFrame:
    embeddings = pd.read_csv(path)
    missing = {"sequence_id", "sequence_sha256"} - set(embeddings.columns)
    if missing:
        raise ValueError(f"missing embedding columns: {', '.join(sorted(missing))}")
    _embedding_columns(embeddings.columns)
    return embeddings


def join_embeddings(rows: pd.DataFrame, embeddings: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    missing = {"sequence_id", "sequence_sha256"} - set(rows.columns)
    if missing:
        raise ValueError(f"missing row columns: {', '.join(sorted(missing))}")
    emb_cols = _embedding_columns(embeddings.columns)
    if embeddings["sequence_id"].duplicated().any():
        raise ValueError("duplicate sequence_id in embedding table")

    ids = rows["sequence_id"].tolist()
    embedding_index = embeddings.set_index("sequence_id", drop=False)
    missing_ids = [seq_id for seq_id in ids if seq_id not in embedding_index.index]
    if missing_ids:
        raise ValueError(f"Missing embeddings for sequence_id: {', '.join(map(str, missing_ids))}")

    selected = embedding_index.loc[ids].reset_index(drop=True)
    expected_hashes = rows["sequence_sha256"].reset_index(drop=True)
    actual_hashes = selected["sequence_sha256"].reset_index(drop=True)
    if not expected_hashes.eq(actual_hashes).all():
        bad = rows.loc[~expected_hashes.eq(actual_hashes), "sequence_id"].astype(str).tolist()
        raise ValueError(f"sequence_sha256 mismatch for sequence_id: {', '.join(bad)}")

    joined = pd.concat([rows.reset_index(drop=True), selected[emb_cols].reset_index(drop=True)], axis=1)
    return joined, joined[emb_cols].to_numpy(dtype=float)


# candidates
def _bool_series(values: pd.Series, default: bool = True) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(default)
    return values.fillna(default).map(lambda value: str(value).lower() in {"true", "1", "yes"})


def generate_single_substitution_pool(
    sequence: str,
    allowed_positions: Iterable[int] | None = None,
    mutation_chain_label: str = "H",
    sequence_id_prefix: str = "candidate",
) -> pd.DataFrame:
    positions = list(range(1, len(sequence) + 1)) if allowed_positions is None else list(allowed_positions)
    rows = []
    order = 1
    for pos in positions:
        if pos < 1 or pos > len(sequence):
            raise ValueError(f"invalid allowed position: {pos}")
        orig = sequence[pos - 1]
        for new in AMINO_ACIDS:
            if new == orig:
                continue
            mutated = sequence[: pos - 1] + new + sequence[pos:]
            mutation_string = f"{mutation_chain_label}{orig}{pos}{new}"
            rows.append(
                {
                    "sequence_id": sequence_id_from_mutations(sequence_id_prefix, mutation_string),
                    "sequence": mutated,
                    "sequence_sha256": sequence_sha256(mutated),
                    "mutation_string": mutation_string,
                    "candidate_origin": "full_single",
                    "order": order,
                    "valid": True,
                    "invalid_reason": "",
                    "tested": False,
                    "failed_synthesis": False,
                    "previously_selected": False,
                    "available_pre_embedding": True,
                    "unavailable_reason_pre_embedding": "",
                }
            )
            order += 1
    return pd.DataFrame(rows)


def _matches(candidates: pd.DataFrame, rows: pd.DataFrame | None) -> pd.Series:
    mask = pd.Series(False, index=candidates.index)
    if rows is None or rows.empty:
        return mask
    for column in ["sequence_id", "sequence_sha256", "mutation_string"]:
        if column in candidates.columns and column in rows.columns:
            values = set(rows[column].dropna().astype(str))
            mask |= candidates[column].astype(str).isin(values)
    return mask


def _tested_rows(training: pd.DataFrame) -> pd.DataFrame:
    if training is None or training.empty:
        return pd.DataFrame()
    measured = (
        training["experiment_status"].eq("measured")
        if "experiment_status" in training.columns
        else pd.Series(False, index=training.index)
    )
    legacy_summary = (
        training["label_source_type"].eq("legacy_training_summary")
        if "label_source_type" in training.columns
        else pd.Series(False, index=training.index)
    )
    return training.loc[measured | legacy_summary]


def mark_available_candidates(
    candidates: pd.DataFrame,
    training: pd.DataFrame,
    failed_synthesis: pd.DataFrame | None = None,
    invalid_candidates: pd.DataFrame | None = None,
    previously_selected: pd.DataFrame | None = None,
) -> pd.DataFrame:
    marked = candidates.copy()
    existing_available = (
        _bool_series(marked["available_pre_embedding"])
        if "available_pre_embedding" in marked.columns
        else pd.Series(True, index=marked.index)
    )
    existing_reason = (
        marked["unavailable_reason_pre_embedding"].fillna("").astype(str)
        if "unavailable_reason_pre_embedding" in marked.columns
        else pd.Series("", index=marked.index)
    )

    marked["tested"] = _matches(marked, _tested_rows(training))
    marked["failed_synthesis"] = _matches(marked, failed_synthesis)
    marked["invalid_sequence"] = _matches(marked, invalid_candidates)
    marked["previously_selected"] = _matches(marked, previously_selected)
    marked["unavailable_reason_pre_embedding"] = existing_reason
    missing_existing_reason = ~existing_available & marked["unavailable_reason_pre_embedding"].eq("")
    marked.loc[missing_existing_reason, "unavailable_reason_pre_embedding"] = "invalid_sequence"
    for column, reason in [
        ("tested", "tested"),
        ("failed_synthesis", "failed_synthesis"),
        ("invalid_sequence", "invalid_sequence"),
        ("previously_selected", "previously_selected"),
    ]:
        mask = marked["unavailable_reason_pre_embedding"].eq("") & marked[column]
        marked.loc[mask, "unavailable_reason_pre_embedding"] = reason
    marked["available_pre_embedding"] = existing_available & marked["unavailable_reason_pre_embedding"].eq("")
    return marked


def build_candidate_pool_from_config(
    config: dict,
    wt_fasta: Path,
    system_id: str | None = None,
    round_id: str | None = None,
    candidate_mode: str | None = None,
) -> pd.DataFrame:
    mode = candidate_mode or config["candidate_pool"]["mode"]
    if mode != "full_single":
        raise ValueError(f"unsupported candidate mode: {mode}")

    naming = config["naming"]
    sequence = parse_fasta_first_sequence(Path(wt_fasta), naming.get("fasta_record_id"))
    excluded = set(config["candidate_pool"].get("excluded_positions", []))
    allowed_positions = [pos for pos in range(1, len(sequence) + 1) if pos not in excluded]
    pool = generate_single_substitution_pool(
        sequence,
        allowed_positions=allowed_positions,
        mutation_chain_label=naming.get("mutation_chain_label", "H"),
        sequence_id_prefix=system_id or config.get("system_id") or "candidate",
    )
    excluded_mutations = set(config["candidate_pool"].get("excluded_mutations", []))
    if excluded_mutations:
        pool = pool.loc[~pool["mutation_string"].isin(excluded_mutations)].reset_index(drop=True)
        pool["order"] = range(1, len(pool) + 1)
    pool.insert(0, "round_id", round_id or config.get("round_id"))
    pool.insert(0, "system_id", system_id or config.get("system_id"))
    return pool


# sdAb warm-start adapter
def _mutation_value(value) -> str:
    if pd.isna(value):
        return ""
    value = str(value).strip()
    return "" if value.upper() in {"WT", "NONE", "NAN"} else value


def _row_base(config: dict, wt_sequence: str, mutation_string: str) -> dict:
    sequence = apply_mutation_string(
        wt_sequence,
        mutation_string,
        expected_chain=config["naming"].get("mutation_chain_label", "H"),
        ignore_chain_for_single_sequence=config["naming"].get("ignore_chain_for_single_sequence", True),
    )
    return {
        "system_id": config.get("system_id"),
        "sequence_id": sequence_id_from_mutations(config.get("system_id", "sdab"), mutation_string),
        "sequence": sequence,
        "sequence_sha256": sequence_sha256(sequence),
        "mutation_string": mutation_string,
    }


def _source_round_id(row: pd.Series, path: Path, config: dict) -> str | None:
    value = row.get("round_id")
    if value is not None and not pd.isna(value) and str(value).strip():
        return str(value)
    for part in reversed(path.parts):
        if re.fullmatch(r"R\d+", part, flags=re.IGNORECASE):
            return part.upper()
    return config.get("round_id")


def _legacy_rows(path: Path, config: dict, wt_sequence: str) -> list[dict]:
    rows = pd.read_csv(path)
    raw_a7 = rows["log_pH74"].map(math.exp) if "log_pH74" in rows.columns else rows.get("A7", 0.0)
    raw_a6 = rows["log_pH6"].map(math.exp) if "log_pH6" in rows.columns else rows.get("A6", 0.0)
    wt_mask = rows.get("mutations_unified", pd.Series("", index=rows.index)).fillna("").astype(str).str.strip().isin(["", "WT"])
    wt_a7 = float(raw_a7[wt_mask].iloc[0]) if wt_mask.any() else None
    wt_a6 = float(raw_a6[wt_mask].iloc[0]) if wt_mask.any() else None

    records = []
    for idx, row in rows.iterrows():
        mutation_string = _mutation_value(row.get("mutations_unified", row.get("mutation_string", "")))
        a7_raw = float(raw_a7.loc[idx])
        a6_raw = float(raw_a6.loc[idx])
        a7 = normalize_by_reference(a7_raw, 0.0, wt_a7) if wt_a7 else clip01(a7_raw)
        a6 = normalize_by_reference(a6_raw, 0.0, wt_a6) if wt_a6 else clip01(a6_raw)
        record = _row_base(config, wt_sequence, mutation_string)
        record.update(
            {
                "A7": a7,
                "A6": a6,
                "noise": 0.0,
                "utility_recomputed": compute_switch_utility(a7, a6),
                "round_id": _source_round_id(row, path, config),
                "experiment_status": "excluded",
                "label_source_type": "legacy_training_summary",
                "label_source_file": str(path),
                "label_source_row_id": idx,
                "label_is_direct_measurement": False,
                "usable_for_training": False,
                "source": row.get("source", ""),
                "group": row.get("group", ""),
                "batch_id": row.get("batch_id", ""),
                "assay_id": row.get("assay_id", ""),
                "collapsed": row.get("collapsed", ""),
                "leverage_flag": row.get("leverage_flag", ""),
            }
        )
        records.append(record)
    return records


def _direct_rows(path: Path, config: dict, wt_sequence: str) -> list[dict]:
    rows = pd.read_csv(path)
    records = []
    for idx, row in rows.iterrows():
        mutation_string = _mutation_value(row.get("mutation_string", row.get("mutations_unified", "")))
        a7 = clip01(row["A7"])
        a6 = clip01(row["A6"])
        noise = float(row.get("noise", 0.0))
        confirmed = str(row.get("confirmed", False)).lower() in {"true", "1", "yes"}
        record = _row_base(config, wt_sequence, mutation_string)
        record.update(
            {
                "A7": a7,
                "A6": a6,
                "noise": noise,
                "utility_recomputed": compute_switch_utility(a7, a6, noise),
                "round_id": _source_round_id(row, path, config),
                "experiment_status": "measured" if confirmed else "excluded",
                "label_source_type": "raw_elisa",
                "label_source_file": str(path),
                "label_source_row_id": idx,
                "label_is_direct_measurement": confirmed,
                "usable_for_training": confirmed,
                "batch_id": row.get("batch_id", ""),
                "assay_id": row.get("assay_id", ""),
            }
        )
        records.append(record)
    return records


def build_sdab_warm_start_rows(
    config: dict,
    wt_fasta: Path,
    legacy_r4_training: Path,
    direct_wet_labels: Path | None = None,
) -> pd.DataFrame:
    wt_sequence = parse_fasta_first_sequence(Path(wt_fasta), config["naming"].get("fasta_record_id"))
    records = _legacy_rows(Path(legacy_r4_training), config, wt_sequence)
    if direct_wet_labels:
        records.extend(_direct_rows(Path(direct_wet_labels), config, wt_sequence))
    return pd.DataFrame(records)
