#!/usr/bin/env python3
"""Build and audit the 64x1000 AF3 validation workflow.

This module is an execution scaffold for the post-design validation stage. It
does not implicitly launch AF3 or MD. The heavy compute steps are represented as
audited input manifests, shard schedules, and command plans so they can be
unlocked and run explicitly.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

import pandas as pd

try:
    from Bio.Data.IUPACData import protein_letters_3to1_extended
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.Polypeptide import is_aa
except Exception:  # pragma: no cover - import failure is reported at runtime.
    MMCIFParser = None
    PDBParser = None
    is_aa = None
    protein_letters_3to1_extended = {}


ROOT = Path(__file__).resolve().parents[2]
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"
OUT_DIR = ROOT / "results/initial_design_generation/af3_1000seed_validation"

PANEL_1E62 = (
    ROOT
    / "results/initial_design_generation/per_target_15k_candidate_pools/"
    / "final_candidate_pool_1E62_15k_draft.csv"
)
PANEL_SDAB = (
    ROOT
    / "results/initial_design_generation/per_target_15k_candidate_pools/"
    / "final_candidate_pool_sdAb_15k_draft.csv"
)
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"
PROGRESS = TASK_DIR / "progress.md"

AF3_PYTHON = "/data/ziyang/code/af3-local/env/bin/python"
AF3_RUNNER = "/data/ziyang/code/af3-local/src/AF3Score/run_af3score.py"
AF3_MODEL_DIR = "/data/ziyang/code/af3-local/models"
AF3_DB_DIR = "/data/ziyang/code/af3-local/databases/af3_database"
AF3_HMMER_BIN_DIR = "/data/ziyang/code/af3-local/env/bin"
AF3_JACKHMMER = f"{AF3_HMMER_BIN_DIR}/jackhmmer"
AF3_HMMALIGN = f"{AF3_HMMER_BIN_DIR}/hmmalign"
AF3_HMMBUILD = f"{AF3_HMMER_BIN_DIR}/hmmbuild"
AF3_HMMSEARCH = f"{AF3_HMMER_BIN_DIR}/hmmsearch"
AF3_NHMMER = f"{AF3_HMMER_BIN_DIR}/nhmmer"
AF3_CUDA_DATA_DIR = "/usr/local/cuda-12.6"
AF3_INFERENCE_ENV_PREFIX = (
    f"XLA_FLAGS=--xla_gpu_cuda_data_dir={AF3_CUDA_DATA_DIR} "
    f"PATH={AF3_CUDA_DATA_DIR}/bin:$PATH "
)

AF3_ALLOWED_GPUS = list(range(8))
AF3_EXCLUDED_GPUS: List[int] = []
AF3_SEEDS = list(range(1, 1001))
SHARD_SIZE = 50
MSA_MAX_PARALLEL = 2
MSA_CORE_SETS = ["0-15", "16-31"]
MD_CORE_SETS = [
    "0-5",
    "6-11",
    "12-17",
    "18-23",
    "24-29",
    "30-35",
    "36-41",
    "42-47",
    "48-53",
    "54-59",
    "60-65",
    "66-71",
]

CONTACT_CUTOFF = 5.0
CLASH_CUTOFF = 2.0


TARGETS = {
    "Ab_1E62": {
        "short": "1E62",
        "panel_csv": PANEL_1E62,
        "antibody_type": "scFv",
        "heavy_fasta": ROOT / "experiments/1E62/data/heavy.fasta",
        "light_fasta": ROOT / "experiments/1E62/data/light.fasta",
        "antigen_fasta": ROOT / "experiments/1E62/data/antigen_genotypes.fasta",
        "antigen_name": "AeS",
        "chains": {"heavy": "H", "design": "L", "antigen": "C"},
        "design_chain": "L",
        "antigen_chain": "C",
        "panel_size": 32,
        "md_mutant_size": 16,
        "top_seed_cap": 8,
        "cluster_cap": 2,
        "required_seeds": [
            "LK24H;LQ35H",
            "LK24H;LY38H",
            "LY31H;LQ35H",
            "LK24H;LY31H",
            "LQ35H;LY38H",
        ],
    },
    "Ab_sdAb": {
        "short": "sdAb",
        "panel_csv": PANEL_SDAB,
        "antibody_type": "VHH",
        "sdab_fasta": ROOT / "experiments/sdab/R2/data/sdab.fasta",
        "antigen_fasta": ROOT / "experiments/sdab/R2/data/antigen.fasta",
        "antigen_name": "AeS",
        "chains": {"design": "A", "antigen": "B"},
        "design_chain": "A",
        "antigen_chain": "B",
        "panel_size": 32,
        "md_mutant_size": 16,
        "top_seed_cap": 8,
        "cluster_cap": 2,
        "required_seeds": [
            "AD110H",
            "AQ100H",
            "AY111H",
            "AD110H;AY111H",
            "AQ100H;AD110H",
            "AG102H",
            "AV105H",
        ],
    },
}

BEST_MODEL_SELECTION_VERSION = "best_model_score_v1"
BEST_MODEL_WEIGHTS = {
    "normalized_complex_confidence": 0.30,
    "interface_plausibility_score": 0.25,
    "parent_contact_retention_score": 0.20,
    "his_interpretability_score": 0.15,
    "new_clash_or_bad_contact_penalty": -0.10,
}


def ensure_outdir() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for sub in [
        "inputs/json",
        "inputs/msa_json",
        "inputs/processed_json_shards",
        "msa",
        "runs",
        "summaries",
        "top_models",
        "cluster_representatives",
        "md_input_structures",
    ]:
        (OUT_DIR / sub).mkdir(parents=True, exist_ok=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    view = df.copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
    lines = [
        "| " + " | ".join(map(str, view.columns)) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in view.columns) + " |")
    return "\n".join(lines)


def read_fasta(path: Path, preferred_name: Optional[str] = None) -> Tuple[str, str]:
    records: List[Tuple[str, str]] = []
    current_name = ""
    current_seq: List[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    records.append((current_name, "".join(current_seq).replace("*", "")))
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_name:
        records.append((current_name, "".join(current_seq).replace("*", "")))
    for name, seq in records:
        if preferred_name is None or name.split()[0] == preferred_name:
            return name, seq
    if not records:
        raise ValueError(f"Could not read FASTA sequence from {path} preferred={preferred_name}")
    raise ValueError(f"Could not find FASTA sequence from {path} preferred={preferred_name}")


def safe_id(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return text[:180]


def truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def as_int(value: Any, default: int = 0) -> int:
    try:
        if value is None or str(value).strip() == "":
            return default
        return int(float(value))
    except Exception:
        return default


def as_float(value: Any, default: float = math.nan) -> float:
    try:
        if value is None or str(value).strip() == "":
            return default
        return float(value)
    except Exception:
        return default


def parse_mutations(text: Any) -> List[Dict[str, Any]]:
    parsed = []
    for token in str(text or "").split(";"):
        token = token.strip()
        if not token or token.lower() in {"nan", "none"}:
            continue
        match = re.match(r"^([A-Z])([A-Z])(\d+)([A-Z])$", token)
        if not match:
            parsed.append({"token": token, "status": "parse_error"})
            continue
        chain, old, pos, new = match.groups()
        parsed.append(
            {
                "token": token,
                "chain": chain,
                "old": old,
                "pos": int(pos),
                "new": new,
                "status": "ok",
            }
        )
    return parsed


def row_score(row: pd.Series, target: str) -> float:
    score = 0.0
    quality = str(row.get("candidate_quality_tier", ""))
    match = str(row.get("template_match_tier", ""))
    template_status = str(row.get("assigned_template_status", ""))
    source = str(row.get("source_pool", ""))
    stage2a = str(row.get("stage2a_final_class", ""))
    role = str(row.get("candidate_type", ""))

    score += {
        "bank_reviewed": 100,
        "high_support": 80,
        "supported": 60,
        "mpnn_favorable_pH_weak": 45,
        "neutral_boundary_or_high_risk": 20,
    }.get(quality, 0)
    score += {
        "template_seed": 80,
        "A_strong": 60,
        "local_expansion_A": 52,
        "B_medium": 42,
        "local_expansion_B": 35,
        "C_seed_only": 10,
    }.get(match, 0)
    score += {
        "active_primary_template": 35,
        "active_primary_backfill_template": 25,
        "secondary_template_complex_weak": 30,
        "secondary_template_complex_unchecked": 22,
        "low_priority_secondary_template": 5,
        "limited_boundary_template": 3,
        "boundary_representative_only": 3,
    }.get(template_status, 0)
    score += {
        "tier2_candidate_bank": 80,
        "stage2a_candidate_list": 45,
        "tier2_candidate_snapshot": 35,
        "sdab_recovery_passlike_supplement": 35,
        "tier1_ranked_production_pool": 25,
        "constrained_local_expansion_v2": 18,
        "constrained_low_frequency_seed_local_expansion_v2": 15,
        "control_anchor_panel": 8,
    }.get(source, 0)
    score += {
        "T2A_high_confidence": 40,
        "T2A_good": 30,
        "T2A_supported_boundary_confirmed": 15,
    }.get(stage2a, 0)
    if role == "control_anchor":
        score += 5

    score += min(as_float(row.get("expanded_pool_score"), 0.0), 100.0) * 0.05
    score += max(0.0, 100.0 - as_float(row.get("mpnn_score_percentile_within_target"), 100.0)) * 0.03
    score += as_float(row.get("his_pka_support_t2_score"), 0.0) * 5.0

    mutation_count = as_int(row.get("mutation_count"), 0)
    if target == "Ab_1E62" and mutation_count == 4:
        score += 6
    if target == "Ab_sdAb" and mutation_count <= 3:
        score += 5
    return score


def sort_candidates(df: pd.DataFrame, target: str) -> pd.DataFrame:
    out = df.copy()
    out["_representative_score"] = out.apply(lambda row: row_score(row, target), axis=1)
    out["_mutation_count_sort"] = pd.to_numeric(out.get("mutation_count", 0), errors="coerce").fillna(0)
    return out.sort_values(
        [
            "_representative_score",
            "candidate_quality_tier",
            "template_match_tier",
            "_mutation_count_sort",
            "variant_id",
        ],
        ascending=[False, True, True, False, True],
    )


def choose_one(
    ranked: pd.DataFrame,
    selected: List[pd.Series],
    selected_ids: Set[str],
    seed_counts: Counter,
    cluster_counts: Counter,
    reason: str,
    seed_cap: int,
    cluster_cap: int,
    mask: Optional[pd.Series] = None,
    allow_cap_relax: bool = False,
) -> bool:
    sub = ranked if mask is None else ranked[mask]
    for _, row in sub.iterrows():
        variant_id = str(row["variant_id"])
        if variant_id in selected_ids:
            continue
        seed = str(row.get("his_seed_set", ""))
        cluster = str(row.get("near_duplicate_cluster_id", ""))
        if not allow_cap_relax and (seed_counts[seed] >= seed_cap or cluster_counts[cluster] >= cluster_cap):
            continue
        row = row.copy()
        row["panel_selection_reason"] = reason
        selected.append(row)
        selected_ids.add(variant_id)
        seed_counts[seed] += 1
        cluster_counts[cluster] += 1
        return True
    return False


def select_target_panel(target: str, cfg: Dict[str, Any]) -> pd.DataFrame:
    df = pd.read_csv(cfg["panel_csv"], low_memory=False)
    df = df[df["target"].eq(target)].copy()
    if "variant_id" not in df.columns:
        raise ValueError(f"{cfg['panel_csv']} has no variant_id column")
    df["mutation_count"] = pd.to_numeric(df.get("mutation_count", 0), errors="coerce").fillna(0).astype(int)
    ranked = sort_candidates(df, target)
    main_ranked = ranked[ranked.get("candidate_type", "main_candidate").ne("control_anchor")]
    control_ranked = ranked[ranked.get("candidate_type", "").eq("control_anchor")]

    selected: List[pd.Series] = []
    selected_ids: Set[str] = set()
    seed_counts: Counter = Counter()
    cluster_counts: Counter = Counter()
    seed_cap = int(cfg["top_seed_cap"])
    cluster_cap = int(cfg["cluster_cap"])

    # Required seeds, including risk controls for sdAb.
    for seed in cfg["required_seeds"]:
        pool = control_ranked if seed in {"AG102H", "AV105H"} else main_ranked
        choose_one(
            pool,
            selected,
            selected_ids,
            seed_counts,
            cluster_counts,
            reason=f"required_seed:{seed}",
            seed_cap=seed_cap,
            cluster_cap=cluster_cap,
            mask=pool["his_seed_set"].astype(str).eq(seed),
            allow_cap_relax=True,
        )

    # Evidence and source strata.
    strata = [
        ("bank_reviewed_template_seed", main_ranked["template_match_tier"].astype(str).eq("template_seed")),
        ("A_strong_or_local_A", main_ranked["template_match_tier"].astype(str).isin({"A_strong", "local_expansion_A"})),
        ("B_medium_or_local_B", main_ranked["template_match_tier"].astype(str).isin({"B_medium", "local_expansion_B"})),
        ("existing_backfill", main_ranked["source_pool"].astype(str).isin({"tier1_ranked_production_pool", "tier2_candidate_snapshot", "stage2a_candidate_list", "sdab_recovery_passlike_supplement"})),
        ("local_expansion", main_ranked["source_pool"].astype(str).str.contains("local_expansion", na=False)),
        ("boundary_or_high_risk", main_ranked["candidate_quality_tier"].astype(str).eq("neutral_boundary_or_high_risk")),
        ("seed_only", main_ranked["template_match_tier"].astype(str).eq("C_seed_only")),
    ]
    if target == "Ab_1E62":
        strata.extend(
            [
                ("1E62_four_mut", main_ranked["mutation_count"].ge(4)),
                ("1E62_non_four_mut", main_ranked["mutation_count"].le(3)),
                ("1E62_non_LK24H", ~main_ranked["his_seed_set"].astype(str).str.contains("LK24H", na=False)),
            ]
        )
    if target == "Ab_sdAb":
        strata.extend(
            [
                ("sdAb_complex_weak", main_ranked["assigned_template_status"].astype(str).str.contains("complex_weak", na=False)),
                ("sdAb_complex_unchecked", main_ranked["assigned_template_status"].astype(str).str.contains("complex_unchecked", na=False)),
                ("sdAb_low_priority_secondary", main_ranked["assigned_template_status"].astype(str).eq("low_priority_secondary_template")),
            ]
        )

    for reason, mask in strata:
        if len(selected) >= cfg["panel_size"]:
            break
        choose_one(
            main_ranked,
            selected,
            selected_ids,
            seed_counts,
            cluster_counts,
            reason=reason,
            seed_cap=seed_cap,
            cluster_cap=cluster_cap,
            mask=mask,
        )

    def selected_count(predicate) -> int:
        return sum(1 for row in selected if predicate(row))

    def fill_minimum(reason: str, minimum: int, mask: pd.Series, predicate) -> None:
        while len(selected) < cfg["panel_size"] and selected_count(predicate) < minimum:
            if not choose_one(
                main_ranked,
                selected,
                selected_ids,
                seed_counts,
                cluster_counts,
                reason=f"minimum_coverage:{reason}",
                seed_cap=seed_cap,
                cluster_cap=cluster_cap + 1,
                mask=mask,
            ):
                break

    fill_minimum(
        "local_expansion",
        5,
        main_ranked["source_pool"].astype(str).str.contains("local_expansion", na=False),
        lambda row: "local_expansion" in str(row.get("source_pool", "")),
    )
    fill_minimum(
        "existing_backfill",
        6,
        main_ranked["source_pool"].astype(str).isin(
            {
                "tier1_ranked_production_pool",
                "tier2_candidate_snapshot",
                "stage2a_candidate_list",
                "sdab_recovery_passlike_supplement",
            }
        ),
        lambda row: str(row.get("source_pool", ""))
        in {
            "tier1_ranked_production_pool",
            "tier2_candidate_snapshot",
            "stage2a_candidate_list",
            "sdab_recovery_passlike_supplement",
        },
    )
    fill_minimum(
        "boundary_or_high_risk",
        3,
        main_ranked["candidate_quality_tier"].astype(str).eq("neutral_boundary_or_high_risk"),
        lambda row: str(row.get("candidate_quality_tier", "")) == "neutral_boundary_or_high_risk",
    )
    fill_minimum(
        "seed_only",
        2,
        main_ranked["template_match_tier"].astype(str).eq("C_seed_only"),
        lambda row: str(row.get("template_match_tier", "")) == "C_seed_only",
    )
    if target == "Ab_1E62":
        fill_minimum(
            "1E62_non_LK24H",
            4,
            ~main_ranked["his_seed_set"].astype(str).str.contains("LK24H", na=False),
            lambda row: "LK24H" not in str(row.get("his_seed_set", "")),
        )
        fill_minimum(
            "1E62_four_mut",
            6,
            main_ranked["mutation_count"].ge(4),
            lambda row: as_int(row.get("mutation_count"), 0) >= 4,
        )
    if target == "Ab_sdAb":
        fill_minimum(
            "sdAb_complex_weak",
            5,
            main_ranked["assigned_template_status"].astype(str).str.contains("complex_weak", na=False),
            lambda row: "complex_weak" in str(row.get("assigned_template_status", "")),
        )
        fill_minimum(
            "sdAb_complex_unchecked",
            5,
            main_ranked["assigned_template_status"].astype(str).str.contains("complex_unchecked", na=False),
            lambda row: "complex_unchecked" in str(row.get("assigned_template_status", "")),
        )

    # Add a small number of controls/risk controls where available.
    control_target = min(cfg["panel_size"], len(selected) + 2)
    while len(selected) < control_target:
        if not choose_one(
            control_ranked,
            selected,
            selected_ids,
            seed_counts,
            cluster_counts,
            reason="control_or_risk_control",
            seed_cap=seed_cap,
            cluster_cap=max(cluster_cap, 3),
            allow_cap_relax=True,
        ):
            break

    # Main fill with diversity caps, then conservative cluster relaxation.
    non_bank_ranked = main_ranked[~main_ranked["source_pool"].astype(str).eq("tier2_candidate_bank")]
    fill_pools = [non_bank_ranked, main_ranked]
    for fill_pool in fill_pools:
        for cap in [cluster_cap, cluster_cap + 1, cluster_cap + 2]:
            for _, row in fill_pool.iterrows():
                if len(selected) >= cfg["panel_size"]:
                    break
                choose_one(
                    fill_pool,
                    selected,
                    selected_ids,
                    seed_counts,
                    cluster_counts,
                    reason=f"diverse_ranked_fill_cluster_cap_{cap}",
                    seed_cap=seed_cap,
                    cluster_cap=cap,
                    mask=fill_pool["variant_id"].astype(str).eq(str(row["variant_id"])),
                )
            if len(selected) >= cfg["panel_size"]:
                break
        if len(selected) >= cfg["panel_size"]:
            break

    out = pd.DataFrame([row.to_dict() for row in selected[: cfg["panel_size"]]])
    out["panel_role"] = out.get("candidate_type", "main_candidate").map(
        lambda value: "risk_control" if value == "control_anchor" else "representative_mutant"
    )
    out["panel_target_limit"] = cfg["panel_size"]
    out["panel_not_top_score_only"] = True
    out["af3_1000seed_required"] = True
    return out


def parent_baseline_rows() -> pd.DataFrame:
    rows = []
    for target, cfg in TARGETS.items():
        if cfg["antibody_type"] == "scFv":
            _, hseq = read_fasta(cfg["heavy_fasta"])
            _, lseq = read_fasta(cfg["light_fasta"])
            _, agseq = read_fasta(cfg["antigen_fasta"], cfg["antigen_name"])
            design_sequence = lseq
            chains = [
                {"chain_id": "H", "sequence": hseq, "role": "heavy"},
                {"chain_id": "L", "sequence": lseq, "role": "design"},
                {"chain_id": "C", "sequence": agseq, "role": "antigen"},
            ]
        else:
            _, aseq = read_fasta(cfg["sdab_fasta"])
            _, agseq = read_fasta(cfg["antigen_fasta"])
            design_sequence = aseq
            chains = [
                {"chain_id": "A", "sequence": aseq, "role": "design"},
                {"chain_id": "B", "sequence": agseq, "role": "antigen"},
            ]
        rows.append(
            {
                "variant_id": f"{cfg['short']}_parent_WT",
                "target": target,
                "panel_role": "parent_baseline",
                "mutation_list": "",
                "his_seed_set": "parent_WT",
                "near_duplicate_cluster_id": f"{target}_parent_WT",
                "sequence": design_sequence,
                "chain_sequences_json": json.dumps(chains, sort_keys=True),
                "af3_1000seed_required": True,
                "parent_not_counted_in_mutant_limit": True,
            }
        )
    return pd.DataFrame(rows)


def chain_entries_for_row(row: pd.Series) -> List[Dict[str, Any]]:
    target = row["target"]
    cfg = TARGETS[target]
    if str(row.get("panel_role", "")) == "parent_baseline":
        chains = json.loads(row["chain_sequences_json"])
    elif cfg["antibody_type"] == "scFv":
        _, hseq = read_fasta(cfg["heavy_fasta"])
        _, agseq = read_fasta(cfg["antigen_fasta"], cfg["antigen_name"])
        chains = [
            {"chain_id": "H", "sequence": hseq, "role": "heavy"},
            {"chain_id": "L", "sequence": str(row["sequence"]), "role": "design"},
            {"chain_id": "C", "sequence": agseq, "role": "antigen"},
        ]
    else:
        _, agseq = read_fasta(cfg["antigen_fasta"])
        chains = [
            {"chain_id": "A", "sequence": str(row["sequence"]), "role": "design"},
            {"chain_id": "B", "sequence": agseq, "role": "antigen"},
        ]
    return [
        {
            "protein": {
                "id": entry["chain_id"],
                "sequence": entry["sequence"],
            }
        }
        for entry in chains
    ]


def make_af3_json(row: pd.Series, job_name: str, seeds: List[int]) -> Dict[str, Any]:
    variant_id = safe_id(str(row["variant_id"]))
    return {
        "dialect": "alphafold3",
        "version": 1,
        "name": job_name,
        "sequences": chain_entries_for_row(row),
        "modelSeeds": [int(seed) for seed in seeds],
        "bondedAtomPairs": None,
        "userCCD": None,
    }


def build_panel() -> None:
    ensure_outdir()
    panels = []
    for target, cfg in TARGETS.items():
        panels.append(select_target_panel(target, cfg))
    panel = pd.concat(panels, ignore_index=True)
    parent = parent_baseline_rows()

    panel_path = OUT_DIR / "validation_panel_64.csv"
    parent_path = OUT_DIR / "validation_panel_parent_baselines.csv"
    panel.to_csv(panel_path, index=False)
    parent.to_csv(parent_path, index=False)

    lines = [
        "# AF3 1000-Seed Validation Panel Audit",
        "",
        "This is a representative validation panel, not a top-score panel.",
        "",
        "Parent / WT baselines are not counted in the 64 mutant limit.",
        "",
        "## Inputs",
        "",
        f"- 1E62 pool: `{PANEL_1E62.relative_to(ROOT)}` sha256={sha256_file(PANEL_1E62)}",
        f"- sdAb pool: `{PANEL_SDAB.relative_to(ROOT)}` sha256={sha256_file(PANEL_SDAB)}",
        "",
        "## Counts",
        "",
    ]
    counts = panel.groupby(["target", "panel_role"], dropna=False).size().reset_index(name="count")
    lines.append(md_table(counts))
    lines.extend(["", "## Seed Distribution", ""])
    lines.append(md_table(panel.groupby(["target", "his_seed_set"], dropna=False).size().reset_index(name="count")))
    lines.extend(["", "## Cluster Distribution Top 12", ""])
    cluster_rows = []
    for target, sub in panel.groupby("target"):
        top = sub["near_duplicate_cluster_id"].value_counts().head(12)
        for cluster, count in top.items():
            cluster_rows.append({"target": target, "near_duplicate_cluster_id": cluster, "count": count})
    lines.append(md_table(pd.DataFrame(cluster_rows)))
    lines.extend(["", "## Mutation Count Distribution", ""])
    lines.append(md_table(panel.groupby(["target", "mutation_count"], dropna=False).size().reset_index(name="count")))
    lines.extend(["", "## Source / Evidence Distribution", ""])
    lines.append(md_table(panel.groupby(["target", "source_pool", "template_match_tier", "candidate_quality_tier"], dropna=False).size().reset_index(name="count")))
    lines.extend(["", "## Parent Baselines", ""])
    lines.append(md_table(parent[["variant_id", "target", "his_seed_set", "parent_not_counted_in_mutant_limit"]]))
    (OUT_DIR / "validation_panel_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote panel rows={len(panel)} parent_rows={len(parent)} to {OUT_DIR}")


def load_panel_with_parents() -> pd.DataFrame:
    panel_path = OUT_DIR / "validation_panel_64.csv"
    parent_path = OUT_DIR / "validation_panel_parent_baselines.csv"
    if not panel_path.exists() or not parent_path.exists():
        raise FileNotFoundError("Run build-panel before this action.")
    panel = pd.read_csv(panel_path, low_memory=False)
    parent = pd.read_csv(parent_path, low_memory=False)
    common_cols = sorted(set(panel.columns).union(parent.columns))
    for df in [panel, parent]:
        for col in common_cols:
            if col not in df.columns:
                df[col] = ""
    return pd.concat([panel[common_cols], parent[common_cols]], ignore_index=True)


def write_text_yaml(path: Path, data: Dict[str, Any]) -> None:
    def render(value: Any, indent: int = 0) -> List[str]:
        prefix = " " * indent
        if isinstance(value, dict):
            lines = []
            for key, item in value.items():
                if isinstance(item, (dict, list)):
                    lines.append(f"{prefix}{key}:")
                    lines.extend(render(item, indent + 2))
                else:
                    lines.append(f"{prefix}{key}: {json.dumps(item) if isinstance(item, str) else item}")
            return lines
        if isinstance(value, list):
            lines = []
            for item in value:
                if isinstance(item, (dict, list)):
                    lines.append(f"{prefix}-")
                    lines.extend(render(item, indent + 2))
                else:
                    lines.append(f"{prefix}- {json.dumps(item) if isinstance(item, str) else item}")
            return lines
        return [f"{prefix}{value}"]

    path.write_text("\n".join(render(data)) + "\n", encoding="utf-8")


def prepare_af3_inputs() -> None:
    ensure_outdir()
    df = load_panel_with_parents()
    json_dir = OUT_DIR / "inputs/json"
    msa_json_dir = OUT_DIR / "inputs/msa_json"
    processed_shard_dir = OUT_DIR / "inputs/processed_json_shards"
    json_dir.mkdir(parents=True, exist_ok=True)
    msa_json_dir.mkdir(parents=True, exist_ok=True)
    processed_shard_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    shard_rows = []
    resource_rows = []
    msa_rows = []
    for i, (_, row) in enumerate(df.iterrows()):
        target = row["target"]
        variant_id = safe_id(str(row["variant_id"]))
        variant_run_dir = OUT_DIR / "runs" / target / variant_id
        msa_dir = OUT_DIR / "msa" / target / variant_id
        msa_dir.mkdir(parents=True, exist_ok=True)
        msa_job_name = f"{variant_id}_msa_1000seed"
        msa_json_path = msa_json_dir / f"{msa_job_name}.json"
        msa_json = make_af3_json(row, msa_job_name, AF3_SEEDS)
        msa_json_path.write_text(json.dumps(msa_json, indent=2), encoding="utf-8")
        processed_json_path = msa_dir / f"{msa_job_name}_data.json"
        seeds_by_shard = [AF3_SEEDS[j : j + SHARD_SIZE] for j in range(0, len(AF3_SEEDS), SHARD_SIZE)]
        msa_job_id = f"msa_{target}_{variant_id}"
        msa_core_set = MSA_CORE_SETS[i % len(MSA_CORE_SETS)]
        msa_command = (
            f"taskset -c {msa_core_set} {AF3_PYTHON} {AF3_RUNNER} "
            f"--json_path {msa_json_path.as_posix()} "
            f"--output_dir {msa_dir.as_posix()} "
            f"--model_dir {AF3_MODEL_DIR} "
            f"--db_dir {AF3_DB_DIR} "
            f"--jackhmmer_binary_path {AF3_JACKHMMER} "
            f"--hmmalign_binary_path {AF3_HMMALIGN} "
            f"--hmmbuild_binary_path {AF3_HMMBUILD} "
            f"--hmmsearch_binary_path {AF3_HMMSEARCH} "
            f"--nhmmer_binary_path {AF3_NHMMER} "
            "--jackhmmer_n_cpu 8 --nhmmer_n_cpu 8 --run_data_pipeline=true --run_inference=false "
            "--write_fold_input_json_file=true"
        )
        msa_rows.append(
            {
                "msa_job_id": msa_job_id,
                "target": target,
                "variant_id": row["variant_id"],
                "max_parallel_msa_jobs": MSA_MAX_PARALLEL,
                "cpu_core_set": msa_core_set,
                "msa_json_path": msa_json_path.relative_to(ROOT).as_posix(),
                "msa_output_dir": msa_dir.relative_to(ROOT).as_posix(),
                "processed_json_path": processed_json_path.relative_to(ROOT).as_posix(),
                "command": msa_command,
                "notes": "MSA concurrency hard-capped at 2 due database IO throughput.",
            }
        )
        resource_rows.append(
            {
                "task_type": "af3_msa_data_pipeline",
                "target": target,
                "variant_id": row["variant_id"],
                "shard_id": "",
                "gpu_id": "",
                "cpu_core_set": msa_core_set,
                "resource_policy": "MSA/data pipeline runs once per variant; max two concurrent jobs.",
                "start_time": "",
                "end_time": "",
            }
        )
        for shard_id, seeds in enumerate(seeds_by_shard):
            shard_name = f"shard_{shard_id:03d}"
            shard_dir = variant_run_dir / shard_name
            shard_dir.mkdir(parents=True, exist_ok=True)
            shard_job_name = f"{variant_id}_{shard_name}"
            job = make_af3_json(row, shard_job_name, seeds)
            json_path = json_dir / f"{variant_id}_{shard_name}.json"
            json_path.write_text(json.dumps(job, indent=2), encoding="utf-8")
            processed_shard_path = processed_shard_dir / f"{variant_id}_{shard_name}.json"
            gpu_id = AF3_ALLOWED_GPUS[(len(shard_rows)) % len(AF3_ALLOWED_GPUS)]
            output_dir = shard_dir / "output"
            command = (
                f"CUDA_VISIBLE_DEVICES={gpu_id} {AF3_PYTHON} {AF3_RUNNER} "
                f"--json_path {processed_shard_path.as_posix()} "
                f"--output_dir {output_dir.as_posix()} "
                f"--model_dir {AF3_MODEL_DIR} "
                f"--db_dir {AF3_DB_DIR} "
                "--run_data_pipeline=false --run_inference=true --write_fold_input_json_file=false"
            )
            rec = {
                "target": target,
                "variant_id": row["variant_id"],
                "panel_role": row.get("panel_role", ""),
                "shard_id": shard_id,
                "shard_name": shard_name,
                "seed_start": min(seeds),
                "seed_end": max(seeds),
                "seed_count": len(seeds),
                "raw_json_path": json_path.relative_to(ROOT).as_posix(),
                "msa_processed_json_path": processed_json_path.relative_to(ROOT).as_posix(),
                "processed_shard_json_path": processed_shard_path.relative_to(ROOT).as_posix(),
                "output_dir": output_dir.relative_to(ROOT).as_posix(),
                "gpu_id": gpu_id,
                "excluded_gpu_ids": ";".join(map(str, AF3_EXCLUDED_GPUS)),
                "msa_job_id": msa_job_id,
                "command": command,
            }
            shard_rows.append(rec)
            resource_rows.append(
                {
                    "task_type": "af3_complex_prediction",
                    "target": target,
                    "variant_id": row["variant_id"],
                    "shard_id": shard_id,
                    "gpu_id": gpu_id,
                    "cpu_core_set": "",
                    "resource_policy": "AF3 GPU selection is based on real-time GPU availability; no GPU is permanently excluded.",
                    "start_time": "",
                    "end_time": "",
                }
            )
            rows.append(rec)

    pd.DataFrame(rows).to_csv(OUT_DIR / "af3_1000seed_input_manifest.csv", index=False)
    pd.DataFrame(shard_rows).to_csv(OUT_DIR / "af3_shard_schedule.csv", index=False)
    pd.DataFrame(msa_rows).to_csv(OUT_DIR / "msa_job_schedule.csv", index=False)
    pd.DataFrame(resource_rows).to_csv(OUT_DIR / "resource_allocation_manifest.csv", index=False)

    write_text_yaml(
        OUT_DIR / "af3_1000seed_run_manifest.yaml",
        {
            "stage": "af3_1000seed_validation",
            "status": "input_prepared_manual_compute_required",
            "af3_static_prediction_is_not_ph_specific": True,
            "seed_count_per_variant": 1000,
            "shard_size": SHARD_SIZE,
            "shards_per_variant": len(AF3_SEEDS) // SHARD_SIZE,
            "allowed_gpus_for_af3": AF3_ALLOWED_GPUS,
            "excluded_gpus_for_af3": AF3_EXCLUDED_GPUS,
            "msa_max_parallel_jobs": MSA_MAX_PARALLEL,
            "resource_policy": {
                "prefer_idle_gpu": True,
                "af3_msa_md_cpu_sets_separate": True,
                "msa_max_parallel": MSA_MAX_PARALLEL,
                "md_max_tasks_per_idle_gpu": 2,
                "do_not_stack_md_on_gpu_running_af3": True,
            },
            "inputs": {
                "validation_panel": (OUT_DIR / "validation_panel_64.csv").relative_to(ROOT).as_posix(),
                "parent_baselines": (OUT_DIR / "validation_panel_parent_baselines.csv").relative_to(ROOT).as_posix(),
            },
        },
    )
    write_text_yaml(
        OUT_DIR / "best_model_selection_config.yaml",
        {
            "version": BEST_MODEL_SELECTION_VERSION,
            "weights": BEST_MODEL_WEIGHTS,
            "best_model_meaning": "most_favorable_static_structure_hypothesis",
            "sampling_support_meaning": "confidence_weight_for_best_model_hypothesis",
        },
    )
    write_text_yaml(
        OUT_DIR / "interface_cluster_config.yaml",
        {
            "version": "interface_cluster_v1",
            "features": [
                "antibody_antigen_contact_fingerprint",
                "his_to_antigen_distance_vector",
                "cdr_antigen_contact_pattern",
                "interface_residue_pair_contact_map",
                "binding_pose_orientation",
                "interface_rmsd_to_parent_or_best_when_available",
            ],
            "jaccard_threshold": 0.50,
            "contact_cutoff_angstrom": CONTACT_CUTOFF,
        },
    )

    commands = [
        "# AF3 1000-seed validation command plan",
        "",
        "AF3 static prediction is complex conformational sampling only; it does not model pH 7.4 or pH 6.0.",
        "",
        "Allowed AF3 GPUs: `0,1,5,6,7`. Excluded GPUs: `2,3,4`.",
        "",
        "Run MSA/data-pipeline once per variant before inference. The preferred launcher respects the shared MSA/database cap and only starts this task when global slots are available:",
        "",
        "```bash",
        f"{Path('/data/ziyang/mamba/envs/optim-pipe/bin/python').as_posix()} analysis/initial_design_generation/run_af3_1000seed_validation.py run-msa-queue --max-parallel 2 --respect-global-msa-cap",
        "```",
        "",
        "Only use the non-global-cap MSA command if external AF3 MSA jobs have been stopped or a manual policy override has been made:",
        "",
        "```bash",
        f"{Path('/data/ziyang/mamba/envs/optim-pipe/bin/python').as_posix()} analysis/initial_design_generation/run_af3_1000seed_validation.py run-msa-queue --max-parallel 2",
        "```",
        "",
        "After MSA completes, run:",
        "",
        "```bash",
        f"{Path('/data/ziyang/mamba/envs/optim-pipe/bin/python').as_posix()} analysis/initial_design_generation/run_af3_1000seed_validation.py split-processed-json",
        "```",
        "",
        "Then run inference-only shard commands from `af3_shard_schedule.csv`.",
        "",
        "Example MSA command:",
        "",
        "```bash",
    ]
    if msa_rows:
        commands.append(msa_rows[0]["command"])
    commands.extend(["```", "", "Example inference command:", "", "```bash"])
    if rows:
        commands.append(rows[0]["command"])
    commands.extend(["```", ""])
    (OUT_DIR / "af3_1000seed_command_plan.md").write_text("\n".join(commands), encoding="utf-8")
    print(f"Prepared AF3 inputs shards={len(rows)} variants={len(df)} allowed_gpus={AF3_ALLOWED_GPUS}")


def split_processed_json() -> None:
    """Split each processed MSA JSON into inference-only seed-shard JSONs."""
    ensure_outdir()
    sched_path = OUT_DIR / "af3_shard_schedule.csv"
    msa_path = OUT_DIR / "msa_job_schedule.csv"
    if not sched_path.exists() or not msa_path.exists():
        raise FileNotFoundError("Run prepare-af3-inputs before split-processed-json.")
    sched = pd.read_csv(sched_path, low_memory=False)
    msa = pd.read_csv(msa_path, low_memory=False)

    def resolve_processed_json(path_value: str) -> Optional[Path]:
        expected = ROOT / str(path_value)
        if expected.exists():
            return expected
        if expected.parent.exists():
            matches = sorted(expected.parent.rglob("*_data.json"))
            if matches:
                return matches[0]
        return None

    msa_lookup = {
        (str(row["target"]), str(row["variant_id"])): resolve_processed_json(str(row["processed_json_path"]))
        for _, row in msa.iterrows()
    }
    rows = []
    for _, row in sched.iterrows():
        key = (str(row["target"]), str(row["variant_id"]))
        processed = msa_lookup.get(key)
        shard_json = ROOT / str(row["processed_shard_json_path"])
        shard_json.parent.mkdir(parents=True, exist_ok=True)
        seed_start = as_int(row["seed_start"])
        seed_end = as_int(row["seed_end"])
        seeds = list(range(seed_start, seed_end + 1))
        if processed is None or not processed.exists():
            rows.append(
                {
                    "target": row["target"],
                    "variant_id": row["variant_id"],
                    "shard_id": row["shard_id"],
                    "processed_json_path": "" if processed is None else processed.relative_to(ROOT).as_posix(),
                    "processed_shard_json_path": row["processed_shard_json_path"],
                    "status": "missing_msa_processed_json",
                    "seed_count": len(seeds),
                }
            )
            continue
        data = json.loads(processed.read_text(encoding="utf-8"))
        data["modelSeeds"] = seeds
        data["name"] = f"{safe_id(str(row['variant_id']))}_shard_{int(row['shard_id']):03d}"
        shard_json.write_text(json.dumps(data, indent=2), encoding="utf-8")
        rows.append(
            {
                "target": row["target"],
                "variant_id": row["variant_id"],
                "shard_id": row["shard_id"],
                "processed_json_path": processed.relative_to(ROOT).as_posix(),
                "processed_shard_json_path": row["processed_shard_json_path"],
                "status": "written",
                "seed_count": len(seeds),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "processed_json_shard_status.csv", index=False)
    written = int(out["status"].eq("written").sum()) if not out.empty else 0
    missing = int(out["status"].ne("written").sum()) if not out.empty else 0
    print(f"Processed JSON shards written={written} missing={missing}")


def run_command_queue(
    jobs: List[Dict[str, Any]],
    max_parallel: int,
    log_dir: Path,
    status_path: Path,
    skip_done=None,
    max_parallel_getter: Optional[Callable[[], int]] = None,
) -> None:
    """Run shell commands with a fixed concurrency and audit status."""
    log_dir.mkdir(parents=True, exist_ok=True)
    status_rows: List[Dict[str, Any]] = []
    existing_completed_job_ids: Set[str] = set()
    if status_path.exists() and status_path.stat().st_size:
        try:
            existing_status = pd.read_csv(status_path, low_memory=False)
        except Exception:
            existing_status = pd.DataFrame()
        if not existing_status.empty:
            status_rows = existing_status.to_dict("records")
            if {"job_id", "status"}.issubset(existing_status.columns):
                completed = existing_status["status"].isin(["success", "skipped_existing_output"])
                existing_completed_job_ids = set(existing_status.loc[completed, "job_id"].astype(str))
    active: List[Dict[str, Any]] = []
    pending = list(jobs)

    def launch(job: Dict[str, Any]) -> None:
        safe = safe_id(str(job.get("job_id", job.get("variant_id", "job"))))
        log_path = log_dir / f"{safe}.log"
        handle = log_path.open("a", encoding="utf-8")
        handle.write(f"\n===== started {time.strftime('%Y-%m-%d %H:%M:%S')} =====\n")
        handle.flush()
        proc = subprocess.Popen(
            str(job["command"]),
            shell=True,
            cwd=str(ROOT),
            stdout=handle,
            stderr=subprocess.STDOUT,
            executable="/bin/bash",
        )
        job["_proc"] = proc
        job["_log_handle"] = handle
        job["_log_path"] = log_path
        job["_start_time"] = time.strftime("%Y-%m-%d %H:%M:%S")
        active.append(job)

    def finish(job: Dict[str, Any], returncode: int) -> None:
        handle = job.pop("_log_handle")
        handle.write(f"\n===== finished {time.strftime('%Y-%m-%d %H:%M:%S')} returncode={returncode} =====\n")
        handle.close()
        row = {key: value for key, value in job.items() if not key.startswith("_")}
        row.update(
            {
                "start_time": job.get("_start_time", ""),
                "end_time": time.strftime("%Y-%m-%d %H:%M:%S"),
                "returncode": returncode,
                "status": "success" if returncode == 0 else "failed",
                "log_path": job["_log_path"].relative_to(ROOT).as_posix(),
            }
        )
        status_rows.append(row)
        if returncode == 0 and row.get("job_id") is not None:
            existing_completed_job_ids.add(str(row["job_id"]))
        pd.DataFrame(status_rows).to_csv(status_path, index=False)

    while pending or active:
        current_max_parallel = max_parallel
        if max_parallel_getter is not None:
            current_max_parallel = max(0, min(max_parallel, int(max_parallel_getter())))
        while pending and len(active) < current_max_parallel:
            job = pending.pop(0)
            job_id = str(job.get("job_id", ""))
            if job_id in existing_completed_job_ids:
                continue
            if skip_done is not None and skip_done(job):
                row = {key: value for key, value in job.items() if not key.startswith("_")}
                row.update(
                    {
                        "start_time": "",
                        "end_time": time.strftime("%Y-%m-%d %H:%M:%S"),
                        "returncode": 0,
                        "status": "skipped_existing_output",
                        "log_path": "",
                    }
                )
                status_rows.append(row)
                if row.get("job_id") is not None:
                    existing_completed_job_ids.add(str(row["job_id"]))
                pd.DataFrame(status_rows).to_csv(status_path, index=False)
                continue
            launch(job)
        time.sleep(10)
        still_active = []
        for job in active:
            proc = job["_proc"]
            rc = proc.poll()
            if rc is None:
                still_active.append(job)
            else:
                finish(job, rc)
        active = still_active


def count_real_af3_msa_jobs() -> int:
    """Count live AF3 data-pipeline jobs without matching this script text."""
    count = 0
    proc_root = Path("/proc")
    for entry in proc_root.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            parts = (entry / "cmdline").read_bytes().split(b"\0")
        except Exception:
            continue
        args = [part.decode("utf-8", "ignore") for part in parts if part]
        if not args:
            continue
        has_af3_runner = any(arg.endswith("run_alphafold.py") for arg in args)
        has_no_inference = "--norun_inference" in args or (
            "--run_data_pipeline=true" in args and "--run_inference=false" in args
        )
        if has_af3_runner and has_no_inference:
            count += 1
    return count


def run_msa_queue(max_parallel: int = MSA_MAX_PARALLEL, respect_global_cap: bool = False) -> None:
    ensure_outdir()
    msa_path = OUT_DIR / "msa_job_schedule.csv"
    if not msa_path.exists():
        raise FileNotFoundError("Run prepare-af3-inputs before run-msa-queue.")
    msa = pd.read_csv(msa_path, low_memory=False)
    max_parallel = min(int(max_parallel), MSA_MAX_PARALLEL)
    jobs = []
    for _, row in msa.iterrows():
        jobs.append(
            {
                "job_id": row["msa_job_id"],
                "target": row["target"],
                "variant_id": row["variant_id"],
                "cpu_core_set": row["cpu_core_set"],
                "processed_json_path": row["processed_json_path"],
                "command": row["command"],
            }
        )

    def done(job: Dict[str, Any]) -> bool:
        expected = ROOT / str(job["processed_json_path"])
        if expected.exists():
            return True
        return expected.parent.exists() and bool(list(expected.parent.rglob("*_data.json")))

    def available_msa_slots() -> int:
        external_jobs = count_real_af3_msa_jobs()
        return max(0, MSA_MAX_PARALLEL - external_jobs)

    run_command_queue(
        jobs=jobs,
        max_parallel=max_parallel,
        log_dir=OUT_DIR / "logs/msa",
        status_path=OUT_DIR / "msa_run_status.csv",
        skip_done=done,
        max_parallel_getter=available_msa_slots if respect_global_cap else None,
    )
    split_processed_json()


def parse_gpu_ids(value: Optional[str]) -> List[int]:
    if not value:
        return list(AF3_ALLOWED_GPUS)
    gpu_ids = []
    for token in str(value).split(","):
        token = token.strip()
        if not token:
            continue
        gpu_id = int(token)
        if gpu_id in AF3_EXCLUDED_GPUS:
            raise ValueError(f"GPU {gpu_id} is excluded by current AF3 policy.")
        gpu_ids.append(gpu_id)
    if not gpu_ids:
        raise ValueError("No GPU ids provided.")
    return gpu_ids


def rewrite_cuda_visible_devices(command: str, gpu_id: int) -> str:
    if re.search(r"\bCUDA_VISIBLE_DEVICES=\S+", command):
        command = re.sub(r"\bCUDA_VISIBLE_DEVICES=\S+", f"CUDA_VISIBLE_DEVICES={gpu_id}", command, count=1)
    else:
        command = f"CUDA_VISIBLE_DEVICES={gpu_id} {command}"
    if "xla_gpu_cuda_data_dir" not in command:
        command = f"{AF3_INFERENCE_ENV_PREFIX}{command}"
    return command


def run_inference_queue(
    max_parallel_per_gpu: int = 1,
    limit: Optional[int] = None,
    gpu_ids: Optional[List[int]] = None,
    max_total_parallel: Optional[int] = None,
    only_variants: Optional[Set[str]] = None,
    status_file: Optional[str] = None,
) -> None:
    ensure_outdir()
    shard_path = OUT_DIR / "af3_shard_schedule.csv"
    split_path = OUT_DIR / "processed_json_shard_status.csv"
    needs_split = True
    if shard_path.exists() and split_path.exists() and split_path.stat().st_size:
        try:
            shard_count = len(pd.read_csv(shard_path, usecols=["target", "variant_id", "shard_id"]))
            split_status = pd.read_csv(split_path, usecols=["status"])
            needs_split = int(split_status["status"].eq("written").sum()) < shard_count
        except Exception:
            needs_split = True
    if needs_split:
        split_processed_json()
    if not shard_path.exists() or not split_path.exists():
        raise FileNotFoundError("Run prepare-af3-inputs and split-processed-json first.")
    shard = pd.read_csv(shard_path, low_memory=False)
    split = pd.read_csv(split_path, low_memory=False)
    ready = split[split["status"].eq("written")][["target", "variant_id", "shard_id"]]
    shard = shard.merge(ready, on=["target", "variant_id", "shard_id"], how="inner")
    if only_variants:
        shard = shard[shard["variant_id"].astype(str).isin(only_variants)]
    if limit:
        shard = shard.head(int(limit))
    active_gpus = gpu_ids or list(AF3_ALLOWED_GPUS)
    jobs = []
    for _, row in shard.iterrows():
        gpu_id = active_gpus[len(jobs) % len(active_gpus)]
        command = rewrite_cuda_visible_devices(str(row["command"]), gpu_id)
        jobs.append(
            {
                "job_id": f"{row['target']}_{safe_id(str(row['variant_id']))}_shard_{int(row['shard_id']):03d}",
                "target": row["target"],
                "variant_id": row["variant_id"],
                "shard_id": row["shard_id"],
                "gpu_id": gpu_id,
                "output_dir": row["output_dir"],
                "command": command,
            }
        )

    def done(job: Dict[str, Any]) -> bool:
        out = ROOT / str(job["output_dir"])
        return out.exists() and any(out.rglob("summary_confidences.json"))

    status_path = OUT_DIR / (status_file or "af3_inference_run_status.csv")
    # Conservative default: one AF3 inference per allowed GPU.
    queue_parallel = max(1, len(active_gpus) * int(max_parallel_per_gpu))
    if max_total_parallel is not None:
        queue_parallel = max(1, min(queue_parallel, int(max_total_parallel)))
    run_command_queue(
        jobs=jobs,
        max_parallel=queue_parallel,
        log_dir=OUT_DIR / "logs/inference",
        status_path=status_path,
        skip_done=done,
    )


def confidence_files_for(path: Path) -> List[Path]:
    candidates = []
    for base in [path.parent, *path.parents[:4]]:
        if base.exists():
            candidates.extend(base.glob("*summary_confidences.json"))
            candidates.extend(base.glob("*confidences.json"))
            candidates.extend(base.glob("*ranking*.json"))
    return sorted(set(candidates))


def read_confidence(path: Path) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for conf in confidence_files_for(path):
        try:
            data = json.loads(conf.read_text())
        except Exception:
            continue
        if isinstance(data, dict):
            for key in [
                "ranking_score",
                "ranking_confidence",
                "confidence_score",
                "ptm",
                "iptm",
                "fraction_disordered",
            ]:
                if key in data and key not in out:
                    out[key] = data[key]
    return out


def parse_structure(path: Path):
    if path.suffix.lower() in {".cif", ".mmcif"}:
        if MMCIFParser is None:
            raise RuntimeError("Bio.PDB MMCIFParser unavailable")
        return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    if PDBParser is None:
        raise RuntimeError("Bio.PDB PDBParser unavailable")
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))


def residue_pos(residue) -> Optional[int]:
    hetflag, resseq, _icode = residue.id
    if str(hetflag).strip():
        return None
    return int(resseq)


def heavy_atom_coords(residue) -> List[Any]:
    coords = []
    for atom in residue.get_atoms():
        element = (atom.element or "").upper()
        if element == "H" or atom.get_name().upper().startswith("H"):
            continue
        coords.append(atom.coord)
    return coords


def min_distance(coords_a: List[Any], coords_b: List[Any]) -> float:
    best = math.inf
    for a in coords_a:
        for b in coords_b:
            dist = float(((a - b) ** 2).sum() ** 0.5)
            if dist < best:
                best = dist
    return best


def residue_aa(residue) -> str:
    name = residue.get_resname().strip()
    return protein_letters_3to1_extended.get(name.title(), name if name else "UNK")


def compute_structure_metrics(path: Path, target: str, mutation_list: str) -> Dict[str, Any]:
    cfg = TARGETS[target]
    design_chain_id = cfg["design_chain"]
    antigen_chain_id = cfg["antigen_chain"]
    try:
        structure = parse_structure(path)
    except Exception as exc:
        return {
            "structure_parse_status": f"parse_failed:{exc}",
            "interface_contact_count": math.nan,
            "contact_signature": "",
            "his_min_antigen_distance": math.nan,
            "clash_count": math.nan,
            "interface_plausibility_score": 0.0,
            "his_interpretability_score": 0.0,
            "new_clash_or_bad_contact_penalty": 1.0,
        }

    chains = {chain.id: chain for chain in structure.get_chains()}
    if design_chain_id not in chains or antigen_chain_id not in chains:
        return {
            "structure_parse_status": "missing_design_or_antigen_chain",
            "interface_contact_count": math.nan,
            "contact_signature": "",
            "his_min_antigen_distance": math.nan,
            "clash_count": math.nan,
            "interface_plausibility_score": 0.0,
            "his_interpretability_score": 0.0,
            "new_clash_or_bad_contact_penalty": 1.0,
        }

    design_res = [res for res in chains[design_chain_id] if is_aa(res, standard=False)]
    antigen_res = [res for res in chains[antigen_chain_id] if is_aa(res, standard=False)]
    antigen_coords_by_pos = {}
    for res in antigen_res:
        pos = residue_pos(res)
        coords = heavy_atom_coords(res)
        if pos is not None and coords:
            antigen_coords_by_pos[pos] = coords

    contacts = []
    clash_count = 0
    for dres in design_res:
        dpos = residue_pos(dres)
        dcoords = heavy_atom_coords(dres)
        if dpos is None or not dcoords:
            continue
        for apos, acoords in antigen_coords_by_pos.items():
            dist = min_distance(dcoords, acoords)
            if dist <= CONTACT_CUTOFF:
                contacts.append((dpos, apos))
            if dist <= CLASH_CUTOFF:
                clash_count += 1

    his_positions = [m["pos"] for m in parse_mutations(mutation_list) if m.get("status") == "ok" and m.get("new") == "H"]
    his_distances = []
    for dres in design_res:
        dpos = residue_pos(dres)
        if dpos not in his_positions and residue_aa(dres) != "H":
            continue
        dcoords = heavy_atom_coords(dres)
        if not dcoords:
            continue
        best = min((min_distance(dcoords, acoords) for acoords in antigen_coords_by_pos.values()), default=math.inf)
        if math.isfinite(best):
            his_distances.append(best)
    his_min = min(his_distances) if his_distances else math.nan
    if not math.isfinite(his_min):
        his_score = 0.50 if not his_positions else 0.0
    elif his_min <= 6.0:
        his_score = 1.0
    elif his_min <= 10.0:
        his_score = 0.75
    elif his_min <= 14.0:
        his_score = 0.35
    else:
        his_score = 0.05

    signature_pairs = sorted(set(contacts))
    antigen_signature = sorted({apos for _dpos, apos in signature_pairs})
    interface_score = min(len(signature_pairs) / 25.0, 1.0)
    return {
        "structure_parse_status": "ok",
        "interface_contact_count": len(signature_pairs),
        "contact_signature": ";".join(f"{d}:{a}" for d, a in signature_pairs[:200]),
        "antigen_contact_signature": ";".join(map(str, antigen_signature[:200])),
        "his_min_antigen_distance": his_min,
        "clash_count": clash_count,
        "interface_plausibility_score": interface_score,
        "his_interpretability_score": his_score,
        "new_clash_or_bad_contact_penalty": min(clash_count / 10.0, 1.0),
    }


def find_structure_files(shard_output_dir: Path) -> List[Path]:
    if not shard_output_dir.exists():
        return []
    return sorted(
        path
        for path in shard_output_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in {".cif", ".mmcif", ".pdb"}
    )


def seed_from_path(path: Path, fallback: int) -> int:
    text = path.as_posix()
    matches = re.findall(r"(?:seed|model_seed)[_-]?(\d+)", text, flags=re.IGNORECASE)
    if matches:
        return int(matches[-1])
    return fallback


def contact_set(signature: str) -> Set[str]:
    if not isinstance(signature, str) or not signature:
        return set()
    return set(item for item in signature.split(";") if item)


def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / float(len(a | b))


def assign_clusters(df: pd.DataFrame) -> pd.DataFrame:
    out = df.sort_values("best_model_score_raw", ascending=False).copy()
    clusters: List[Dict[str, Any]] = []
    assignments: Dict[int, str] = {}
    for idx, row in out.iterrows():
        sig = contact_set(row.get("contact_signature", ""))
        assigned = None
        for cluster in clusters:
            if jaccard(sig, cluster["signature"]) >= 0.50:
                assigned = cluster["id"]
                break
        if assigned is None:
            assigned = f"ifc_{len(clusters):04d}"
            clusters.append({"id": assigned, "signature": sig})
        assignments[idx] = assigned
    out["cluster_id"] = pd.Series(assignments)
    return out.sort_index()


def collect_af3_results() -> None:
    ensure_outdir()
    manifest_path = OUT_DIR / "af3_1000seed_input_manifest.csv"
    if not manifest_path.exists():
        raise FileNotFoundError("Run prepare-af3-inputs before collect-af3-results.")
    manifest = pd.read_csv(manifest_path, low_memory=False)
    panel = load_panel_with_parents()
    meta = {
        (str(row["target"]), str(row["variant_id"])): row
        for _, row in panel.iterrows()
    }

    seed_status_rows = []
    raw_rows = []
    for _, shard in manifest.iterrows():
        target = str(shard["target"])
        variant_id = str(shard["variant_id"])
        meta_row = meta.get((target, variant_id))
        mutation_list = "" if meta_row is None else str(meta_row.get("mutation_list", ""))
        output_dir = ROOT / str(shard["output_dir"])
        structures = find_structure_files(output_dir)
        seeds = list(range(as_int(shard["seed_start"]), as_int(shard["seed_end"]) + 1))
        seed_seen: Counter[int] = Counter()
        for fallback_i, struct_path in enumerate(structures):
            seed = seed_from_path(struct_path, seeds[min(fallback_i, len(seeds) - 1)])
            seed_seen[seed] += 1
            conf = read_confidence(struct_path)
            metrics = compute_structure_metrics(struct_path, target, mutation_list)
            iptm = as_float(conf.get("iptm", conf.get("ranking_score", math.nan)))
            ptm = as_float(conf.get("ptm", math.nan))
            ranking = as_float(conf.get("ranking_score", conf.get("ranking_confidence", iptm)))
            norm_conf = max(0.0, min(1.0, iptm if math.isfinite(iptm) else ranking))
            parent_retention = 0.5
            best_score_raw = (
                BEST_MODEL_WEIGHTS["normalized_complex_confidence"] * norm_conf
                + BEST_MODEL_WEIGHTS["interface_plausibility_score"] * metrics["interface_plausibility_score"]
                + BEST_MODEL_WEIGHTS["parent_contact_retention_score"] * parent_retention
                + BEST_MODEL_WEIGHTS["his_interpretability_score"] * metrics["his_interpretability_score"]
                + BEST_MODEL_WEIGHTS["new_clash_or_bad_contact_penalty"] * metrics["new_clash_or_bad_contact_penalty"]
            )
            raw_rows.append(
                {
                    "target": target,
                    "variant_id": variant_id,
                    "panel_role": shard.get("panel_role", ""),
                    "shard_id": shard["shard_id"],
                    "seed": seed,
                    "model_path": struct_path.relative_to(ROOT).as_posix(),
                    "ptm": ptm,
                    "iptm": iptm,
                    "ranking_score": ranking,
                    "normalized_complex_confidence": norm_conf,
                    "parent_contact_retention_score": parent_retention,
                    "best_model_score_raw": best_score_raw,
                    **metrics,
                }
            )
        for seed in seeds:
            seed_status_rows.append(
                {
                    "target": target,
                    "variant_id": variant_id,
                    "shard_id": shard["shard_id"],
                    "seed": seed,
                    "seed_status": "success" if seed_seen[seed] else "missing_or_failed",
                    "model_count": seed_seen[seed],
                }
            )

    seed_status = pd.DataFrame(seed_status_rows)
    raw = pd.DataFrame(raw_rows)
    seed_status.to_csv(OUT_DIR / "af3_1000seed_seed_status.csv", index=False)
    if raw.empty:
        pd.DataFrame(
            columns=[
                "target",
                "variant_id",
                "panel_role",
                "shard_id",
                "seed",
                "model_path",
                "ptm",
                "iptm",
                "ranking_score",
                "normalized_complex_confidence",
                "parent_contact_retention_score",
                "best_model_score_raw",
                "structure_parse_status",
                "interface_contact_count",
                "contact_signature",
                "his_min_antigen_distance",
                "clash_count",
                "interface_plausibility_score",
                "his_interpretability_score",
                "new_clash_or_bad_contact_penalty",
            ]
        ).to_csv(OUT_DIR / "af3_1000seed_raw_model_summary.csv", index=False)
        write_static_empty_report(seed_status)
        print("No AF3 model outputs found; wrote seed status only.")
        return

    clustered_parts = []
    for (_target, _variant), sub in raw.groupby(["target", "variant_id"], dropna=False):
        clustered_parts.append(assign_clusters(sub))
    raw = pd.concat(clustered_parts, ignore_index=True)
    raw.to_csv(OUT_DIR / "af3_1000seed_raw_model_summary.csv", index=False)

    parent_signatures = parent_contact_signatures(raw)
    raw["parent_contact_retention_score"] = raw.apply(
        lambda row: contact_retention(row, parent_signatures), axis=1
    )
    raw["best_model_score_raw"] = raw.apply(best_model_score, axis=1)
    raw.to_csv(OUT_DIR / "af3_1000seed_raw_model_summary.csv", index=False)

    best_rows = []
    cluster_rows = []
    top_rows = []
    parent_compare_rows = []
    classifications = []
    for (target, variant_id), sub in raw.groupby(["target", "variant_id"], dropna=False):
        success_count = int(seed_status[seed_status["target"].eq(target) & seed_status["variant_id"].eq(variant_id)]["seed_status"].eq("success").sum())
        sorted_sub = sub.sort_values("best_model_score_raw", ascending=False)
        best = sorted_sub.iloc[0].copy()
        best_cluster = best["cluster_id"]
        cluster_sub = sub[sub["cluster_id"].eq(best_cluster)]
        interface_plausible = sub["interface_plausibility_score"].fillna(0).ge(0.50)
        best_cluster_fraction = len(cluster_sub) / float(max(len(sub), 1))
        interface_fraction = float(interface_plausible.mean()) if len(sub) else 0.0
        parent_delta = parent_delta_class(target, variant_id, sub, parent_signatures)
        static_class = classify_static(best, best_cluster_fraction, interface_fraction, parent_delta)
        best_rows.append(
            {
                "target": target,
                "variant_id": variant_id,
                "num_successful_models": success_count,
                "variant_completion_status": "variant_full"
                if success_count >= 1000
                else ("variant_valid" if success_count >= 950 else "variant_incomplete"),
                "best_model_path": best["model_path"],
                "best_model_seed": best["seed"],
                "best_model_rank": 1,
                "best_model_score": best["best_model_score_raw"],
                "best_model_reason": "max_best_model_score_v1",
                "best_model_selection_version": BEST_MODEL_SELECTION_VERSION,
                "best_model_for_static_review": best["model_path"],
                "best_model_for_MD_start": best["model_path"],
                "best_cluster_id": best_cluster,
                "best_cluster_size": len(cluster_sub),
                "best_cluster_fraction": best_cluster_fraction,
                "interface_plausible_fraction": interface_fraction,
                "iptm_max": sub["iptm"].max(),
                "iptm_top10_mean": sorted_sub["iptm"].head(10).mean(),
                "iptm_top50_mean": sorted_sub["iptm"].head(50).mean(),
                "iptm_median": sub["iptm"].median(),
                "ptm_top10_mean": sorted_sub["ptm"].head(10).mean(),
                "his_min_antigen_distance_best": best["his_min_antigen_distance"],
                "static_classification": static_class,
                "parent_comparison_class": parent_delta,
            }
        )
        for cluster_id, csub in sub.groupby("cluster_id", dropna=False):
            rep = csub.sort_values("best_model_score_raw", ascending=False).iloc[0]
            cluster_rows.append(
                {
                    "target": target,
                    "variant_id": variant_id,
                    "cluster_id": cluster_id,
                    "cluster_size": len(csub),
                    "cluster_fraction": len(csub) / float(max(len(sub), 1)),
                    "cluster_median_iptm": csub["iptm"].median(),
                    "cluster_interface_plausible_fraction": csub["interface_plausibility_score"].fillna(0).ge(0.50).mean(),
                    "cluster_contact_signature": rep.get("contact_signature", ""),
                    "cluster_representative_model_path": rep["model_path"],
                }
            )
        for rank, (_, row) in enumerate(sorted_sub.head(50).iterrows(), start=1):
            top_rows.append(
                {
                    "target": target,
                    "variant_id": variant_id,
                    "top_rank": rank,
                    "model_path": row["model_path"],
                    "seed": row["seed"],
                    "best_model_score": row["best_model_score_raw"],
                    "iptm": row["iptm"],
                    "ptm": row["ptm"],
                    "cluster_id": row["cluster_id"],
                }
            )
        parent_compare_rows.append(parent_compare_summary(target, variant_id, sub, parent_signatures))
        classifications.append(
            {
                "target": target,
                "variant_id": variant_id,
                "static_classification": static_class,
                "parent_comparison_class": parent_delta,
                "best_cluster_fraction": best_cluster_fraction,
                "interface_plausible_fraction": interface_fraction,
                "md_eligibility": md_eligibility(static_class),
            }
        )

    variant_summary = pd.DataFrame(best_rows)
    cluster_summary = pd.DataFrame(cluster_rows)
    top_manifest = pd.DataFrame(top_rows)
    parent_compare = pd.DataFrame(parent_compare_rows)
    classification = pd.DataFrame(classifications)
    variant_summary.to_csv(OUT_DIR / "af3_1000seed_variant_summary.csv", index=False)
    variant_summary.to_csv(OUT_DIR / "best_model_selection_audit.csv", index=False)
    cluster_summary.to_csv(OUT_DIR / "af3_1000seed_cluster_summary.csv", index=False)
    cluster_summary.to_csv(OUT_DIR / "af3_1000seed_interface_summary.csv", index=False)
    top_manifest.to_csv(OUT_DIR / "af3_1000seed_top_cohort_manifest.csv", index=False)
    parent_compare.to_csv(OUT_DIR / "af3_1000seed_parent_baseline_comparison.csv", index=False)
    classification.to_csv(OUT_DIR / "af3_1000seed_candidate_classification.csv", index=False)
    write_static_review_reports(variant_summary, cluster_summary, parent_compare, classification)
    print(f"Collected AF3 models={len(raw)} variants={len(variant_summary)}")


def best_model_score(row: pd.Series) -> float:
    return (
        BEST_MODEL_WEIGHTS["normalized_complex_confidence"] * as_float(row.get("normalized_complex_confidence"), 0.0)
        + BEST_MODEL_WEIGHTS["interface_plausibility_score"] * as_float(row.get("interface_plausibility_score"), 0.0)
        + BEST_MODEL_WEIGHTS["parent_contact_retention_score"] * as_float(row.get("parent_contact_retention_score"), 0.5)
        + BEST_MODEL_WEIGHTS["his_interpretability_score"] * as_float(row.get("his_interpretability_score"), 0.0)
        + BEST_MODEL_WEIGHTS["new_clash_or_bad_contact_penalty"] * as_float(row.get("new_clash_or_bad_contact_penalty"), 1.0)
    )


def parent_contact_signatures(raw: pd.DataFrame) -> Dict[str, Set[str]]:
    out = {}
    parents = raw[raw["variant_id"].astype(str).str.contains("parent_WT")]
    for target, sub in parents.groupby("target"):
        if sub.empty:
            continue
        best = sub.sort_values("best_model_score_raw", ascending=False).iloc[0]
        out[target] = contact_set(best.get("contact_signature", ""))
    return out


def contact_retention(row: pd.Series, parent_signatures: Dict[str, Set[str]]) -> float:
    parent = parent_signatures.get(row["target"], set())
    current = contact_set(row.get("contact_signature", ""))
    if not parent:
        return 0.5
    return len(parent & current) / float(max(len(parent), 1))


def parent_delta_class(
    target: str,
    variant_id: str,
    sub: pd.DataFrame,
    parent_signatures: Dict[str, Set[str]],
) -> str:
    if "parent_WT" in str(variant_id):
        frac = float(sub["interface_plausibility_score"].fillna(0).ge(0.50).mean()) if len(sub) else 0.0
        return "parent_supported" if frac >= 0.20 else "parent_weak"
    parent = parent_signatures.get(target)
    if not parent:
        return "parent_baseline_missing"
    mutant_frac = float(sub["interface_plausibility_score"].fillna(0).ge(0.50).mean()) if len(sub) else 0.0
    retention = float(sub["parent_contact_retention_score"].mean()) if "parent_contact_retention_score" in sub else 0.5
    parent_supported = bool(parent)
    if parent_supported and mutant_frac >= 0.20 and retention >= 0.30:
        return "parent_supported_mutant_supported"
    if parent_supported and (mutant_frac < 0.05 or retention < 0.15):
        return "parent_supported_mutant_weak"
    if mutant_frac >= 0.20:
        return "parent_weak_mutant_supported_manual_review"
    return "parent_weak_mutant_weak"


def classify_static(best: pd.Series, cluster_fraction: float, interface_fraction: float, parent_delta: str) -> str:
    score = as_float(best.get("best_model_score_raw"), 0.0)
    parse_ok = str(best.get("structure_parse_status", "")) == "ok"
    severe = as_float(best.get("new_clash_or_bad_contact_penalty"), 1.0) >= 0.8
    if not parse_ok or severe or parent_delta == "parent_supported_mutant_weak":
        return "unsupported"
    if score >= 0.45 and cluster_fraction >= 0.20 and interface_fraction >= 0.20:
        return "best-supported"
    if score >= 0.35 and (cluster_fraction >= 0.05 or interface_fraction >= 0.05):
        return "best-plausible"
    if score >= 0.35:
        return "best-rare"
    return "unsupported"


def md_eligibility(static_class: str) -> str:
    if static_class in {"best-supported", "best-plausible"}:
        return "main_md_candidate"
    if static_class == "best-rare":
        return "boundary_or_manual_review"
    return "control_only_or_reject"


def parent_compare_summary(target: str, variant_id: str, sub: pd.DataFrame, parent_signatures: Dict[str, Set[str]]) -> Dict[str, Any]:
    return {
        "target": target,
        "variant_id": variant_id,
        "parent_signature_available": target in parent_signatures,
        "mean_parent_contact_retention": sub.get("parent_contact_retention_score", pd.Series(dtype=float)).mean(),
        "delta_interpretation": parent_delta_class(target, variant_id, sub, parent_signatures),
    }


def write_static_empty_report(seed_status: pd.DataFrame) -> None:
    lines = [
        "# AF3 1000-Seed Static Review",
        "",
        "No AF3 outputs were collected yet. Input and seed-status manifests are present.",
        "",
        "## Seed Status Counts",
        "",
        md_table(seed_status.groupby("seed_status").size().reset_index(name="count")) if not seed_status.empty else "_No rows._",
    ]
    (OUT_DIR / "af3_1000seed_static_review.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_static_review_reports(
    variant_summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    parent_compare: pd.DataFrame,
    classification: pd.DataFrame,
) -> None:
    lines = [
        "# AF3 1000-Seed Static Review",
        "",
        "AF3 1000 seeds are complex conformational sampling only; they are not pH 7.4 / pH 6.0 predictions.",
        "",
        "Best model is the most favorable static structure hypothesis; sampling support is its confidence weight.",
        "",
        "## Classification Counts",
        "",
        md_table(classification.groupby(["target", "static_classification"], dropna=False).size().reset_index(name="count")),
        "",
        "## Best Model Summary",
        "",
        md_table(
            variant_summary[
                [
                    "target",
                    "variant_id",
                    "num_successful_models",
                    "variant_completion_status",
                    "best_model_seed",
                    "best_model_score",
                    "best_cluster_fraction",
                    "interface_plausible_fraction",
                    "static_classification",
                    "parent_comparison_class",
                ]
            ]
        ),
        "",
        "## Parent Baseline Comparison",
        "",
        md_table(parent_compare),
    ]
    (OUT_DIR / "af3_1000seed_static_review.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    sdab = classification[classification["target"].eq("Ab_sdAb")].copy()
    sdab_lines = [
        "# sdAb Specific Interpretation Report",
        "",
        "sdAb remains a secondary / exploratory branch unless parent baseline and mutant comparison clearly support upgrade.",
        "",
        "## sdAb Classification Counts",
        "",
        md_table(sdab.groupby(["static_classification", "parent_comparison_class"], dropna=False).size().reset_index(name="count")) if not sdab.empty else "_No sdAb rows._",
    ]
    (OUT_DIR / "sdAb_specific_interpretation_report.md").write_text("\n".join(sdab_lines) + "\n", encoding="utf-8")

    parent_lines = [
        "# Parent Baseline Comparison Report",
        "",
        "Parent / WT 1000-seed baselines are required before strong mutant conclusions.",
        "",
        md_table(parent_compare),
    ]
    (OUT_DIR / "parent_baseline_comparison_report.md").write_text("\n".join(parent_lines) + "\n", encoding="utf-8")


def build_md_panel() -> None:
    ensure_outdir()
    class_path = OUT_DIR / "af3_1000seed_candidate_classification.csv"
    summary_path = OUT_DIR / "af3_1000seed_variant_summary.csv"
    panel_path = OUT_DIR / "validation_panel_64.csv"
    parent_path = OUT_DIR / "validation_panel_parent_baselines.csv"
    if not class_path.exists() or not summary_path.exists():
        raise FileNotFoundError("Run collect-af3-results with completed AF3 outputs before build-md-panel.")
    classification = pd.read_csv(class_path, low_memory=False)
    summary = pd.read_csv(summary_path, low_memory=False)
    panel = pd.read_csv(panel_path, low_memory=False)
    parents = pd.read_csv(parent_path, low_memory=False)
    merged = panel.merge(classification, on=["target", "variant_id"], how="left")
    merged = merged.merge(summary[["target", "variant_id", "best_model_for_MD_start", "best_model_score", "best_cluster_fraction", "interface_plausible_fraction"]], on=["target", "variant_id"], how="left")
    selected = []
    for target, cfg in TARGETS.items():
        sub = merged[merged["target"].eq(target)].copy()
        sub["_priority"] = sub["static_classification"].map({"best-supported": 0, "best-plausible": 1, "best-rare": 2, "unsupported": 9}).fillna(9)
        sub = sub.sort_values(["_priority", "best_model_score", "variant_id"], ascending=[True, False, True])
        main = sub[sub["static_classification"].isin({"best-supported", "best-plausible"})].head(10)
        boundary = sub[sub["static_classification"].eq("best-rare")].head(4)
        controls = sub[sub["panel_role"].astype(str).str.contains("control", na=False)].head(2)
        take = pd.concat([main, boundary, controls]).drop_duplicates("variant_id").head(cfg["md_mutant_size"])
        selected.append(take)
    md_panel = pd.concat(selected, ignore_index=True)
    md_panel["md_panel_role"] = md_panel["static_classification"].map(
        {
            "best-supported": "main_md_candidate",
            "best-plausible": "main_md_candidate",
            "best-rare": "boundary_md_candidate",
            "unsupported": "negative_or_risk_control_only",
        }
    ).fillna("manual_review")
    md_panel["parent_not_counted_in_mutant_limit"] = False
    parent_md = parents.copy()
    parent_md["static_classification"] = "parent_baseline"
    parent_md["md_panel_role"] = "parent_baseline"
    parent_md["best_model_for_MD_start"] = ""
    parent_md["parent_not_counted_in_mutant_limit"] = True
    for col in sorted(set(md_panel.columns).union(parent_md.columns)):
        if col not in md_panel.columns:
            md_panel[col] = ""
        if col not in parent_md.columns:
            parent_md[col] = ""
    out = pd.concat([md_panel, parent_md[md_panel.columns]], ignore_index=True)
    out.to_csv(OUT_DIR / "md_panel_32.csv", index=False)

    lines = [
        "# MD Panel Audit",
        "",
        "Parent / WT baselines are added separately and do not count toward the 32 mutant limit.",
        "",
        "## Counts",
        "",
        md_table(out.groupby(["target", "md_panel_role"], dropna=False).size().reset_index(name="count")),
        "",
        "## Static Class Composition",
        "",
        md_table(out.groupby(["target", "static_classification"], dropna=False).size().reset_index(name="count")),
    ]
    (OUT_DIR / "md_panel_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote MD panel rows={len(out)} mutant_rows={len(md_panel)}")


def prepare_md_inputs() -> None:
    ensure_outdir()
    md_panel_path = OUT_DIR / "md_panel_32.csv"
    if not md_panel_path.exists():
        raise FileNotFoundError("Run build-md-panel before prepare-md-inputs.")
    panel = pd.read_csv(md_panel_path, low_memory=False)
    protonation_rows = []
    prep_rows = []
    schedule_rows = []
    resource_rows = []
    for i, (_, row) in enumerate(panel.iterrows()):
        target = row["target"]
        variant_id = safe_id(str(row["variant_id"]))
        mutations = parse_mutations(row.get("mutation_list", ""))
        his_positions = sorted({m["pos"] for m in mutations if m.get("status") == "ok" and m.get("new") == "H"})
        if not his_positions:
            his_positions = []
        for pos in his_positions:
            protonation_rows.append(
                {
                    "variant_id": row["variant_id"],
                    "target": target,
                    "his_position": pos,
                    "pH7p4_state": "neutral_HIE_or_HID",
                    "pH6p0_state": "protonated_HIP",
                    "reason": "default_design_His_fixed_protonation_policy",
                    "source": "AF3_1000seed_validation_plan_v3",
                }
            )
        if not his_positions:
            protonation_rows.append(
                {
                    "variant_id": row["variant_id"],
                    "target": target,
                    "his_position": "",
                    "pH7p4_state": "no_design_His_detected",
                    "pH6p0_state": "no_design_His_detected",
                    "reason": "parent_or_non_His_control",
                    "source": "AF3_1000seed_validation_plan_v3",
                }
            )
        pdb_path = str(row.get("best_model_for_MD_start", ""))
        prep_rows.append(
            {
                "variant_id": row["variant_id"],
                "target": target,
                "input_model_path": pdb_path,
                "severe_clash_check": "pending",
                "missing_residue_check": "pending",
                "chain_naming_check": "pending",
                "antibody_antigen_completeness_check": "pending",
                "his_protonation_local_conflict_check": "pending",
                "energy_minimization_status": "pending",
                "nvt_equilibration_status": "pending",
                "npt_equilibration_status": "pending",
            }
        )
        gpu_id = AF3_ALLOWED_GPUS[(i // 2) % len(AF3_ALLOWED_GPUS)]
        core_set = MD_CORE_SETS[i % len(MD_CORE_SETS)]
        for ph in [7.4, 6.0]:
            schedule_rows.append(
                {
                    "variant_id": row["variant_id"],
                    "target": target,
                    "ph": ph,
                    "duration_ns": 50,
                    "gpu_id": gpu_id,
                    "cpu_core_set": core_set,
                    "max_md_tasks_per_idle_gpu": 2,
                    "protonation": "fixed",
                    "notes": "Do not share CPU core sets between MD tasks; do not stack MD on a GPU running AF3.",
                }
            )
            resource_rows.append(
                {
                    "task_type": "md_fixed_protonation",
                    "target": target,
                    "variant_id": row["variant_id"],
                    "shard_id": "",
                    "gpu_id": gpu_id,
                    "cpu_core_set": core_set,
                    "resource_policy": "max two MD tasks per fully idle GPU; CPU core sets must not overlap",
                    "start_time": "",
                    "end_time": "",
                }
            )
    pd.DataFrame(protonation_rows).to_csv(OUT_DIR / "md_protonation_state_table.csv", index=False)
    pd.DataFrame(schedule_rows).to_csv(OUT_DIR / "md_resource_schedule.csv", index=False)
    pd.DataFrame(resource_rows).to_csv(OUT_DIR / "md_resource_allocation_manifest.csv", index=False)

    prep_df = pd.DataFrame(prep_rows)
    lines = [
        "# MD Input Structure Preparation Report",
        "",
        "This report is generated before MD. All checks must be completed before production MD is unlocked.",
        "",
        md_table(prep_df),
    ]
    (OUT_DIR / "md_input_structure_preparation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Prepared MD protonation/schedule rows={len(schedule_rows)}")


def update_task_files() -> None:
    report = OUT_DIR / "validation_panel_audit.md"
    lines = [
        "# AF3 1000-Seed Validation Workflow Implementation",
        "",
        "Verdict: `WORKFLOW_SCAFFOLD_IMPLEMENTED__COMPUTE_LOCKED_UNTIL_MANUAL_RUN`.",
        "",
        "This implementation adds the validation-stage workflow for representative panel construction, AF3 1000-seed input/shard scheduling, static result collection, best-model classification, and MD panel/protonation preparation.",
        "",
        "No AF3 or MD production compute is run by the implementation step.",
        "",
        "## Current Resource Policy",
        "",
        "- AF3 stage allowed GPUs: `0,1,5,6,7`.",
        "- AF3 stage excluded GPUs: `2,3,4`.",
        "- MSA/data-pipeline concurrency hard cap: 2 jobs.",
        "- MD tasks must use non-overlapping CPU core sets.",
        "- A fully idle GPU can host at most two MD tasks.",
        "- Do not stack MD on a GPU currently running AF3 or another heavy GPU job.",
        "",
        "## Main Outputs",
        "",
        f"- `{(OUT_DIR / 'validation_panel_64.csv').relative_to(ROOT)}`",
        f"- `{(OUT_DIR / 'validation_panel_parent_baselines.csv').relative_to(ROOT)}`",
        f"- `{(OUT_DIR / 'validation_panel_audit.md').relative_to(ROOT)}`",
        f"- `{(OUT_DIR / 'af3_1000seed_run_manifest.yaml').relative_to(ROOT)}`",
        f"- `{(OUT_DIR / 'af3_shard_schedule.csv').relative_to(ROOT)}`",
        f"- `{(OUT_DIR / 'msa_job_schedule.csv').relative_to(ROOT)}`",
        f"- `{(OUT_DIR / 'resource_allocation_manifest.csv').relative_to(ROOT)}`",
        "",
        "## Compute Status",
        "",
        "- AF3 1000-seed compute: locked until manual run.",
        "- MD: locked until AF3 static classification and MD panel audit.",
        "- Final synthesis/order list: locked.",
    ]
    if report.exists():
        lines.extend(["", "## Panel Audit Snapshot", "", report.read_text(encoding="utf-8")])
    CURRENT_STAGE_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    progress_block = (
        "## 2026-06-05 AF3 1000-seed validation scaffold\n\n"
        "- Implemented representative validation panel and AF3/MD validation workflow scaffold.\n"
        "- AF3 compute remains locked; inputs and resource schedules are prepared for manual execution.\n"
        "- AF3 GPU policy excludes GPUs 2,3,4 and uses GPUs 0,1,5,6,7.\n"
    )
    if PROGRESS.exists():
        text = PROGRESS.read_text(encoding="utf-8")
        marker = "## 2026-06-05 AF3 1000-seed validation scaffold"
        if marker in text:
            text = text.split(marker)[0].rstrip() + "\n\n" + progress_block
        else:
            text = text.rstrip() + "\n\n" + progress_block
    else:
        text = progress_block
    PROGRESS.write_text(text.rstrip() + "\n", encoding="utf-8")


def build_all_inputs() -> None:
    build_panel()
    prepare_af3_inputs()
    update_task_files()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=[
            "build-panel",
            "prepare-af3-inputs",
            "split-processed-json",
            "run-msa-queue",
            "run-inference-queue",
            "collect-af3-results",
            "build-md-panel",
            "prepare-md-inputs",
            "build-all-inputs",
            "update-task-report",
        ],
    )
    parser.add_argument("--max-parallel", type=int, default=None)
    parser.add_argument("--max-total-parallel", type=int, default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--gpus", type=str, default=None)
    parser.add_argument("--only-variant", action="append", default=None)
    parser.add_argument("--status-file", type=str, default=None)
    parser.add_argument("--respect-global-msa-cap", action="store_true")
    args = parser.parse_args()
    if args.action == "build-panel":
        build_panel()
    elif args.action == "prepare-af3-inputs":
        prepare_af3_inputs()
    elif args.action == "split-processed-json":
        split_processed_json()
    elif args.action == "run-msa-queue":
        run_msa_queue(
            max_parallel=args.max_parallel or MSA_MAX_PARALLEL,
            respect_global_cap=args.respect_global_msa_cap,
        )
    elif args.action == "run-inference-queue":
        run_inference_queue(
            max_parallel_per_gpu=args.max_parallel or 1,
            limit=args.limit,
            gpu_ids=parse_gpu_ids(args.gpus),
            max_total_parallel=args.max_total_parallel,
            only_variants=set(args.only_variant or []),
            status_file=args.status_file,
        )
    elif args.action == "collect-af3-results":
        collect_af3_results()
    elif args.action == "build-md-panel":
        build_md_panel()
    elif args.action == "prepare-md-inputs":
        prepare_md_inputs()
    elif args.action == "build-all-inputs":
        build_all_inputs()
    elif args.action == "update-task-report":
        update_task_files()


if __name__ == "__main__":
    main()
