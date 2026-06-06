#!/usr/bin/env python3
"""Build executable ProteinMPNN runner inputs for P0 integration.

The P0 integration scaffold intentionally does not run ProteinMPNN. This
command bridges the next step: given an explicit backbone manifest, it creates
ProteinMPNN JSONL inputs with the correct constraint semantics.

Without a manifest, it writes a manifest template and a report explaining why
runner input generation remains blocked.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml

from analysis.initial_design_generation import run_dry_run as dry
from analysis.initial_design_generation import run_p0_mpnn_integration as p0


ROOT = dry.ROOT
CONFIG = dry.CONFIG
P0_DIR = ROOT / "results/initial_design_generation/p0_mpnn"
DEFAULT_OUT = ROOT / "results/initial_design_generation/p0_mpnn_runner_inputs"
REFERENCE = ROOT / "results/ph_sensitive_40aa_window/tables/reference_sequence_map.csv"
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
TARGET_REVERSE = {"Ab_1E62": "1E62", "Ab_sdAb": "sdAb"}
TARGET_SAFE = {"1E62": "Ab_1E62", "sdAb": "Ab_sdAb"}


def load_yaml(path: Path) -> dict:
    with path.open() as fh:
        return yaml.safe_load(fh)


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_config() -> dict:
    return load_yaml(CONFIG)


def normalize_target(value: str) -> str:
    return TARGET_REVERSE.get(str(value), str(value))


def safe_target(value: str) -> str:
    value = normalize_target(value)
    return TARGET_SAFE.get(value, value)


def split_positions(text: str | float | int | None) -> list[int]:
    if text is None or pd.isna(text):
        return []
    out: list[int] = []
    for item in str(text).split(";"):
        item = item.strip()
        if item:
            out.append(int(item))
    return out


def split_alphabet(text: str | float | None) -> set[str]:
    if text is None or pd.isna(text):
        return set()
    values = [x.strip() for x in str(text).replace("|", ";").split(";") if x.strip()]
    return {x for x in values if len(x) == 1 and x in ALPHABET}


def truthy(value: object) -> bool:
    if value is None or pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_seed_positions(seed_set: str | float | None) -> list[int]:
    if seed_set is None or pd.isna(seed_set):
        return []
    positions: list[int] = []
    for token in str(seed_set).split(";"):
        match = re.match(r"[A-Z][A-Z](\d+)H$", token.strip())
        if match:
            positions.append(int(match.group(1)))
    return positions


def parse_pdb_chains(path: Path) -> dict[str, list[dict]]:
    chains: dict[str, OrderedDict] = OrderedDict()
    with path.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")) or len(line) < 27:
                continue
            chain = line[21].strip()
            if not chain:
                continue
            key = (line[22:26].strip(), line[26].strip())
            chains.setdefault(chain, OrderedDict()).setdefault(
                key,
                {
                    "pdb_resseq": key[0],
                    "pdb_icode": key[1],
                    "resname": line[17:20].strip(),
                },
            )
    return {chain: list(items.values()) for chain, items in chains.items()}


def write_manifest_template(config: dict, path: Path) -> None:
    entries = []
    for target, target_cfg in config["targets"].items():
        entries.append(
            {
                "target": target,
                "backbone_id": f"TODO_{target}_backbone_001",
                "structure_path": "TODO/path/to/backbone.pdb",
                "project_chain_id": target_cfg["chain_id"],
                "mpnn_design_chain_id": "TODO",
                "visible_chain_ids": ["TODO"],
                "local_to_mpnn_index_mode": "chain_order",
                "notes": "Use a reviewed parent complex backbone; do not use old structures unless explicitly approved.",
            }
        )
    data = {
        "version": 1,
        "status": "template_not_runner_ready",
        "backbones": entries,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False, allow_unicode=True)


def validate_manifest(manifest: dict, config: dict) -> tuple[list[dict], list[str]]:
    errors: list[str] = []
    entries = manifest.get("backbones", []) if isinstance(manifest, dict) else []
    if not entries:
        return [], ["manifest has no backbones"]
    valid: list[dict] = []
    for i, entry in enumerate(entries):
        prefix = f"backbones[{i}]"
        target = normalize_target(entry.get("target", ""))
        if target not in config["targets"]:
            errors.append(f"{prefix}: unknown target {entry.get('target')!r}")
            continue
        structure = ROOT / str(entry.get("structure_path", ""))
        if not structure.exists():
            errors.append(f"{prefix}: missing structure_path {structure}")
            continue
        if structure.suffix.lower() != ".pdb":
            errors.append(f"{prefix}: only PDB is currently supported for runner JSONL generation")
            continue
        chains = parse_pdb_chains(structure)
        design_chain = str(entry.get("mpnn_design_chain_id", "")).strip()
        if design_chain not in chains:
            errors.append(f"{prefix}: design chain {design_chain!r} not found in {structure}")
            continue
        visible = [str(x).strip() for x in entry.get("visible_chain_ids", []) if str(x).strip()]
        missing_visible = [x for x in visible if x not in chains]
        if missing_visible:
            errors.append(f"{prefix}: visible chains not found: {missing_visible}")
            continue
        target_cfg = config["targets"][target]
        max_needed = max(
            int(target_cfg["window"].split(":")[1].split("-")[1]),
            max(int(x) for x in target_cfg["hard_protect_positions"]),
            max(int(x) for x in target_cfg["his_seed_positions"]),
        )
        if len(chains[design_chain]) < max_needed:
            errors.append(
                f"{prefix}: design chain length {len(chains[design_chain])} < required local position {max_needed}"
            )
            continue
        valid.append({**entry, "target": target, "_chains": chains})
    return valid, errors


def group_omit_rules(position_to_omit: dict[int, str]) -> list[list[object]]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for pos, omitted in position_to_omit.items():
        grouped[omitted].append(int(pos))
    return [[sorted(positions), omitted] for omitted, positions in sorted(grouped.items())]


def allowed_by_position(evidence: pd.DataFrame, target: str) -> dict[int, set[str]]:
    sub = evidence[evidence["target"].map(normalize_target) == target]
    out: dict[int, set[str]] = {}
    for _, row in sub.iterrows():
        pos = int(row["position"])
        allowed = split_alphabet(row.get("restricted_alphabet"))
        if not allowed:
            allowed = STANDARD_AA.copy()
        out[pos] = allowed
    return out


def reference_sequences() -> dict[tuple[str, str], str]:
    ref = pd.read_csv(REFERENCE)
    out: dict[tuple[str, str], str] = {}
    for (target, chain), sub in ref.groupby([ref["target"].map(normalize_target), "chain"]):
        out[(target, chain)] = "".join(sub.sort_values("local_pos")["aa"].astype(str).tolist())
    return out


def nxs_t_sites(seq: str) -> set[int]:
    sites: set[int] = set()
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in {"S", "T"}:
            sites.add(i + 1)
    return sites


def creates_new_nxs_t(parent: str, position: int, aa: str) -> bool:
    if position < 1 or position > len(parent):
        return False
    mutated = parent[: position - 1] + aa + parent[position:]
    return bool(nxs_t_sites(mutated) - nxs_t_sites(parent))


def rescue_allowed_set(allowed_map: dict[int, set[str]], parent: str, position: int) -> set[str]:
    allowed = set(allowed_map.get(position, STANDARD_AA.copy()))
    allowed.discard("C")
    allowed.discard("H")
    allowed.discard("X")
    allowed = {aa for aa in allowed if not creates_new_nxs_t(parent, position, aa)}
    return allowed


def ranked_sparse_rescue_positions(
    evidence: pd.DataFrame,
    target: str,
    rescue_positions: list[int],
    seed_positions: list[int],
    hard_positions: set[int],
    allowed_map: dict[int, set[str]],
    parent: str,
) -> list[int]:
    sub = evidence[evidence["target"].map(normalize_target) == target].copy()
    by_pos = {int(row.position): row for row in sub.itertuples(index=False)}
    ranked: list[tuple[float, int]] = []
    for pos in sorted(set(rescue_positions)):
        if pos in hard_positions or pos in seed_positions:
            continue
        row = by_pos.get(pos)
        if row is None:
            continue
        if truthy(getattr(row, "protected_from_mutation", False)):
            continue
        if truthy(getattr(row, "known_failure_flag", False)):
            continue
        if not truthy(getattr(row, "rescue_eligible", False)):
            continue
        if not truthy(getattr(row, "applies_to_generation", False)):
            continue
        if not rescue_allowed_set(allowed_map, parent, pos):
            continue
        distance = min(abs(pos - seed) for seed in seed_positions) if seed_positions else 99
        region = str(getattr(row, "region", ""))
        score = -float(distance)
        if "CDR" in region:
            score += 2.0
        if str(getattr(row, "wetlab_support", "")).strip().lower() not in {"", "nan"}:
            score += 2.5
        if str(getattr(row, "structural_support", "")).strip().lower() not in {"", "nan"}:
            score += 1.0
        if str(getattr(row, "interface_support", "")).strip().lower() not in {"", "nan"}:
            score += 1.0
        ranked.append((score, pos))
    ranked.sort(key=lambda x: (-x[0], x[1]))
    return [pos for _, pos in ranked]


def sparse_rescue_subsets(target: str, ranked_positions: list[int]) -> list[tuple[int, ...]]:
    if target == "sdAb":
        return [(pos,) for pos in ranked_positions[:4]]
    subsets: list[tuple[int, ...]] = [(pos,) for pos in ranked_positions[:3]]
    if len(ranked_positions) >= 2:
        subsets.append((ranked_positions[0], ranked_positions[1]))
    return subsets[:4]


def make_symlink(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    os.symlink(src.resolve(), dst)


def build_constrained_inputs(
    manifest_entries: list[dict],
    output_dir: Path,
    *,
    sparse_smoke: bool = False,
    sparse_full: bool = False,
    smoke_samples_per_job: int = 32,
    smoke_seed_limit_per_target: int = 4,
) -> dict:
    if sparse_smoke and sparse_full:
        raise ValueError("sparse_smoke and sparse_full are mutually exclusive")
    config = read_config()
    evidence = pd.read_csv(ROOT / "results/initial_design_generation/tables/evidence_ledger.csv")
    plan = pd.read_csv(P0_DIR / "constrained_mpnn_rescue_plan.csv")
    refs = reference_sequences()

    pdb_dir = output_dir / "pdbs"
    jsonl_dir = output_dir / "jsonl"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    jsonl_dir.mkdir(parents=True, exist_ok=True)

    chain_id: dict[str, list[list[str]]] = {}
    fixed_positions: dict[str, dict[str, list[int]]] = {}
    omit_aa: dict[str, dict[str, list[list[object]]]] = {}
    job_rows: list[dict] = []

    for entry in manifest_entries:
        target = normalize_target(entry["target"])
        target_plan = plan[plan["target"].map(normalize_target) == target].copy()
        if target_plan.empty:
            continue
        if sparse_smoke:
            seed_order = target_plan["his_seed_set"].astype(str).drop_duplicates().head(smoke_seed_limit_per_target)
            target_plan = target_plan[
                target_plan["his_seed_set"].astype(str).isin(set(seed_order))
                & target_plan["mpnn_temperature"].astype(float).isin({0.1, 0.2})
            ].copy()
        design_chain = str(entry["mpnn_design_chain_id"])
        chains = entry["_chains"]
        all_chains = list(chains)
        visible = [x for x in entry.get("visible_chain_ids", []) if x in chains and x != design_chain]
        for chain in all_chains:
            if chain != design_chain and chain not in visible:
                visible.append(chain)
        design_len = len(chains[design_chain])
        allowed_map = allowed_by_position(evidence, target)
        structure = ROOT / str(entry["structure_path"])
        parent = refs[(target, design_chain)]
        hard_positions = set(int(x) for x in config["targets"][target]["hard_protect_positions"])

        for row_index, row in target_plan.iterrows():
            seed_positions = split_positions(row["fixed_his_positions"])
            full_rescue_positions = split_positions(row["rescue_design_positions"])
            if sparse_smoke or sparse_full:
                ranked_rescue = ranked_sparse_rescue_positions(
                    evidence,
                    target,
                    full_rescue_positions,
                    seed_positions,
                    hard_positions,
                    allowed_map,
                    parent,
                )
                rescue_subsets = sparse_rescue_subsets(target, ranked_rescue)
            else:
                rescue_subsets = [tuple(full_rescue_positions)]
            for subset_index, rescue_subset in enumerate(rescue_subsets):
                rescue_positions = list(rescue_subset)
                design_positions = sorted(set(seed_positions) | set(rescue_positions))
                fixed_design = [pos for pos in range(1, design_len + 1) if pos not in set(design_positions)]
                position_to_omit: dict[int, str] = {}
                for pos in seed_positions:
                    position_to_omit[pos] = "".join(aa for aa in ALPHABET if aa != "H")
                for pos in rescue_positions:
                    allowed = rescue_allowed_set(allowed_map, parent, pos)
                    if not allowed:
                        continue
                    position_to_omit[pos] = "".join(aa for aa in ALPHABET if aa not in allowed)

                suffix = f"row{int(row_index):04d}"
                if sparse_smoke or sparse_full:
                    suffix += f"__subset{subset_index:02d}"
                job_name = f"{entry['backbone_id']}__constrained__{suffix}"
                make_symlink(structure, pdb_dir / f"{job_name}.pdb")
                chain_id[job_name] = [[design_chain], visible]
                fixed_positions[job_name] = {design_chain: fixed_design}
                for chain in visible:
                    fixed_positions[job_name][chain] = []
                omit_aa[job_name] = {design_chain: group_omit_rules(position_to_omit)}
                for chain in visible:
                    omit_aa[job_name][chain] = []
                if sparse_smoke:
                    planned_samples = smoke_samples_per_job
                    input_policy = "sparse_smoke"
                elif sparse_full:
                    planned_samples = max(1, int(round(int(row["planned_raw_samples"]) / len(rescue_subsets))))
                    if target == "sdAb" and planned_samples > 144:
                        planned_samples = 144
                    input_policy = "sparse_full"
                else:
                    planned_samples = int(row["planned_raw_samples"])
                    input_policy = "full_rescue_region"
                job_rows.append(
                    {
                        "job_name": job_name,
                        "target": safe_target(target),
                        "backbone_id": entry["backbone_id"],
                        "structure_path": str(structure),
                        "mpnn_design_chain_id": design_chain,
                        "visible_chain_ids": ";".join(visible),
                        "his_seed_set": row["his_seed_set"],
                        "seed_positions_forced_to_H": ";".join(str(x) for x in seed_positions),
                        "rescue_design_positions": ";".join(str(x) for x in rescue_positions),
                        "fixed_design_chain_position_count": len(fixed_design),
                        "designable_position_count": len(design_positions),
                        "mpnn_temperature": row["mpnn_temperature"],
                        "planned_raw_samples": planned_samples,
                        "runner_ready": True,
                        "constrained_input_policy": input_policy,
                    }
                )

    with (jsonl_dir / "chain_id.jsonl").open("w") as fh:
        fh.write(json.dumps(chain_id) + "\n")
    with (jsonl_dir / "fixed_positions.jsonl").open("w") as fh:
        fh.write(json.dumps(fixed_positions) + "\n")
    with (jsonl_dir / "omit_AA.jsonl").open("w") as fh:
        fh.write(json.dumps(omit_aa) + "\n")
    parsed_jsonl = jsonl_dir / "parsed_pdbs.jsonl"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "third_party/ProteinMPNN/helper_scripts/parse_multiple_chains.py"),
            "--input_path",
            str(pdb_dir),
            "--output_path",
            str(parsed_jsonl),
        ],
        check=True,
    )
    parsed_by_name = {
        json.loads(line)["name"]: json.loads(line)
        for line in parsed_jsonl.read_text().splitlines()
        if line.strip()
    }
    jobs = pd.DataFrame(job_rows)
    dry.write_csv(jobs, output_dir / "mpnn_runner_jobs_constrained.csv")
    constrained_commands: list[str] = [
        "# ProteinMPNN constrained-generation commands grouped by planned sampling temperature.",
        "# Run on an appropriate GPU node.",
        "# Each command uses jobs with the same planned temperature and num_seq_per_target.",
        "",
    ]
    command_rows: list[dict] = []
    temp_group_count = 0
    command_count = 0
    for (temp, planned_samples), sub in jobs.groupby(["mpnn_temperature", "planned_raw_samples"], sort=True):
        temp_tag = str(temp).replace(".", "p")
        sample_count = int(planned_samples)
        group_tag = f"T_{temp_tag}_N{sample_count}"
        temp_dir = jsonl_dir / "by_temperature" / group_tag
        temp_dir.mkdir(parents=True, exist_ok=True)
        names = sub["job_name"].astype(str).tolist()
        with (temp_dir / "chain_id.jsonl").open("w") as fh:
            fh.write(json.dumps({name: chain_id[name] for name in names}) + "\n")
        with (temp_dir / "fixed_positions.jsonl").open("w") as fh:
            fh.write(json.dumps({name: fixed_positions[name] for name in names}) + "\n")
        with (temp_dir / "omit_AA.jsonl").open("w") as fh:
            fh.write(json.dumps({name: omit_aa[name] for name in names}) + "\n")
        with (temp_dir / "parsed_pdbs.jsonl").open("w") as fh:
            for name in names:
                fh.write(json.dumps(parsed_by_name[name]) + "\n")
        temp_out = output_dir / "mpnn_out" / group_tag
        constrained_commands.append(
            " ".join(
                [
                    "/data/ziyang/mamba/envs/proteinmpnn/bin/python",
                    "third_party/ProteinMPNN/protein_mpnn_run.py",
                    f"--jsonl_path {temp_dir / 'parsed_pdbs.jsonl'}",
                    f"--chain_id_jsonl {temp_dir / 'chain_id.jsonl'}",
                    f"--fixed_positions_jsonl {temp_dir / 'fixed_positions.jsonl'}",
                    f"--omit_AA_jsonl {temp_dir / 'omit_AA.jsonl'}",
                    f"--out_folder {temp_out}",
                    f"--num_seq_per_target {sample_count}",
                    f"--sampling_temp {temp}",
                    f"--seed {config['random_seed']}",
                    "--batch_size 1",
                ]
            )
        )
        command_rows.append(
            {
                "command_index": command_count,
                "mpnn_temperature": temp,
                "num_seq_per_target": sample_count,
                "job_count": int(len(sub)),
                "planned_raw_samples": int(sub["planned_raw_samples"].sum()),
                "parsed_pdbs_jsonl": str(temp_dir / "parsed_pdbs.jsonl"),
                "out_folder": str(temp_out),
            }
        )
        command_count += 1
        temp_group_count += 1
    (output_dir / "constrained_generation_runner_commands.sh").write_text(
        "\n".join(constrained_commands) + "\n"
    )
    dry.write_csv(pd.DataFrame(command_rows), output_dir / "mpnn_constrained_command_plan.csv")
    commands = [
        "# Parsed PDB JSONL is generated by build_mpnn_runner_inputs.py.",
        f"# Parsed JSONL: {parsed_jsonl}",
        "",
        "# Example ProteinMPNN command for debugging only.",
        "# Production constrained generation should use constrained_generation_runner_commands.sh to preserve planned temperatures.",
        f"/data/ziyang/mamba/envs/proteinmpnn/bin/python third_party/ProteinMPNN/protein_mpnn_run.py --jsonl_path {parsed_jsonl} --chain_id_jsonl {jsonl_dir / 'chain_id.jsonl'} --fixed_positions_jsonl {jsonl_dir / 'fixed_positions.jsonl'} --omit_AA_jsonl {jsonl_dir / 'omit_AA.jsonl'} --out_folder {output_dir / 'mpnn_out'} --num_seq_per_target 1 --sampling_temp 0.1 --seed {config['random_seed']} --batch_size 1",
    ]
    (output_dir / "example_runner_commands.sh").write_text("\n".join(commands) + "\n")
    return {
        "status": "runner_inputs_ready",
        "output_dir": display_path(output_dir),
        "job_count": int(len(jobs)),
        "jsonl_dir": display_path(jsonl_dir),
        "parsed_jsonl": display_path(parsed_jsonl),
        "temperature_group_count": int(temp_group_count),
        "constrained_command_count": int(command_count),
        "planned_constrained_raw_samples": int(jobs["planned_raw_samples"].sum()),
    }


def top_single_his_seed_sets(target: str, limit: int) -> list[str]:
    summary_path = ROOT / "results/initial_design_generation/production_initial_pool/his_seed_set_summary.csv"
    summary = pd.read_csv(summary_path)
    safe = safe_target(target)
    sub = summary[(summary["target"] == safe) & (summary["His_count"].astype(int) == 1)].copy()
    sub = sub.sort_values(["count", "his_seed_set"], ascending=[False, True])
    return sub["his_seed_set"].astype(str).head(limit).tolist()


def build_relaxed_inputs(manifest_entries: list[dict], output_dir: Path) -> dict:
    config = read_config()
    relax_cfg = config["proteinmpnn_integration"]["relaxed_counterfactual_generation"]
    evidence = pd.read_csv(ROOT / "results/initial_design_generation/tables/evidence_ledger.csv")
    plan = pd.read_csv(P0_DIR / "relaxed_mpnn_counterfactual_plan.csv")

    pdb_dir = output_dir / "relaxed_pdbs"
    jsonl_dir = output_dir / "jsonl" / "relaxed"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    jsonl_dir.mkdir(parents=True, exist_ok=True)

    chain_id: dict[str, list[list[str]]] = {}
    fixed_positions: dict[str, dict[str, list[int]]] = {}
    omit_aa: dict[str, dict[str, list[list[object]]]] = {}
    job_rows: list[dict] = []

    temperatures = [float(x) for x in relax_cfg.get("temperatures", [0.1, 0.2])]
    top_seed_limit = int(relax_cfg.get("his_seeded_top_single_seed_sets_per_target", 4))

    for entry in manifest_entries:
        target = normalize_target(entry["target"])
        target_plan = plan[plan["target"].map(normalize_target) == target].copy()
        if target_plan.empty:
            continue
        design_chain = str(entry["mpnn_design_chain_id"])
        chains = entry["_chains"]
        all_chains = list(chains)
        visible = [x for x in entry.get("visible_chain_ids", []) if x in chains and x != design_chain]
        for chain in all_chains:
            if chain != design_chain and chain not in visible:
                visible.append(chain)
        design_len = len(chains[design_chain])
        allowed_map = allowed_by_position(evidence, target)
        structure = ROOT / str(entry["structure_path"])
        hard_positions = set(int(x) for x in config["targets"][target]["hard_protect_positions"])

        for plan_index, row in target_plan.iterrows():
            branch = str(row["branch"])
            design_positions = [pos for pos in split_positions(row["mpnn_designed_positions"]) if pos not in hard_positions]
            if branch.endswith("His_optional"):
                seed_sets = [""]
            else:
                seed_sets = top_single_his_seed_sets(target, top_seed_limit)
            if not seed_sets:
                continue
            raw_total = int(relax_cfg["raw_samples_per_target_by_mode"].get(branch, 5000))
            samples_per_job = max(1, raw_total // (len(seed_sets) * len(temperatures)))

            for seed_set in seed_sets:
                seed_positions = parse_seed_positions(seed_set)
                for temp in temperatures:
                    temp_tag = str(temp).replace(".", "p")
                    seed_tag = seed_set.replace(";", "_") if seed_set else "HisOptional"
                    branch_tag = branch.replace("relaxed_MPNN_", "")
                    job_name = f"{entry['backbone_id']}__relaxed__{branch_tag}__{seed_tag}__T{temp_tag}"
                    make_symlink(structure, pdb_dir / f"{job_name}.pdb")
                    design_set = set(design_positions) | set(seed_positions)
                    fixed_design = [pos for pos in range(1, design_len + 1) if pos not in design_set]
                    position_to_omit: dict[int, str] = {}
                    for pos in design_set:
                        if pos in seed_positions:
                            position_to_omit[pos] = "".join(aa for aa in ALPHABET if aa != "H")
                            continue
                        allowed = set(allowed_map.get(pos, STANDARD_AA.copy()))
                        allowed.discard("C")
                        allowed.discard("X")
                        if not allowed:
                            allowed = STANDARD_AA - {"C"}
                        position_to_omit[pos] = "".join(aa for aa in ALPHABET if aa not in allowed)

                    chain_id[job_name] = [[design_chain], visible]
                    fixed_positions[job_name] = {design_chain: fixed_design}
                    for chain in visible:
                        fixed_positions[job_name][chain] = []
                    omit_aa[job_name] = {design_chain: group_omit_rules(position_to_omit)}
                    for chain in visible:
                        omit_aa[job_name][chain] = []
                    job_rows.append(
                        {
                            "job_name": job_name,
                            "target": safe_target(target),
                            "branch": branch,
                            "backbone_id": entry["backbone_id"],
                            "structure_path": str(structure),
                            "mpnn_design_chain_id": design_chain,
                            "visible_chain_ids": ";".join(visible),
                            "his_seed_set": seed_set,
                            "seed_positions_forced_to_H": ";".join(str(x) for x in seed_positions),
                            "relaxed_design_positions": ";".join(str(x) for x in sorted(design_set)),
                            "fixed_design_chain_position_count": len(fixed_design),
                            "designable_position_count": len(design_set),
                            "mpnn_temperature": temp,
                            "planned_raw_samples": samples_per_job,
                            "final_allocation_status": "audit_only_no_preallocated_quota",
                            "runner_ready": True,
                        }
                    )

    with (jsonl_dir / "chain_id.jsonl").open("w") as fh:
        fh.write(json.dumps(chain_id) + "\n")
    with (jsonl_dir / "fixed_positions.jsonl").open("w") as fh:
        fh.write(json.dumps(fixed_positions) + "\n")
    with (jsonl_dir / "omit_AA.jsonl").open("w") as fh:
        fh.write(json.dumps(omit_aa) + "\n")
    parsed_jsonl = jsonl_dir / "parsed_pdbs.jsonl"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "third_party/ProteinMPNN/helper_scripts/parse_multiple_chains.py"),
            "--input_path",
            str(pdb_dir),
            "--output_path",
            str(parsed_jsonl),
        ],
        check=True,
    )
    parsed_by_name = {
        json.loads(line)["name"]: json.loads(line)
        for line in parsed_jsonl.read_text().splitlines()
        if line.strip()
    }
    jobs = pd.DataFrame(job_rows)
    dry.write_csv(jobs, output_dir / "mpnn_runner_jobs_relaxed.csv")
    relaxed_commands: list[str] = [
        "# ProteinMPNN relaxed/counterfactual commands grouped by branch, temperature, and sample count.",
        "# Run on an appropriate GPU node. Outputs are audit-only until manual unlock.",
        "",
    ]
    command_rows: list[dict] = []
    command_count = 0
    if not jobs.empty:
        for (branch, temp, planned_samples), sub in jobs.groupby(
            ["branch", "mpnn_temperature", "planned_raw_samples"], sort=True
        ):
            temp_tag = str(temp).replace(".", "p")
            branch_tag = str(branch).replace("relaxed_MPNN_", "")
            sample_count = int(planned_samples)
            group_tag = f"{branch_tag}_T_{temp_tag}_N{sample_count}"
            group_dir = jsonl_dir / "by_mode" / group_tag
            group_dir.mkdir(parents=True, exist_ok=True)
            names = sub["job_name"].astype(str).tolist()
            with (group_dir / "chain_id.jsonl").open("w") as fh:
                fh.write(json.dumps({name: chain_id[name] for name in names}) + "\n")
            with (group_dir / "fixed_positions.jsonl").open("w") as fh:
                fh.write(json.dumps({name: fixed_positions[name] for name in names}) + "\n")
            with (group_dir / "omit_AA.jsonl").open("w") as fh:
                fh.write(json.dumps({name: omit_aa[name] for name in names}) + "\n")
            with (group_dir / "parsed_pdbs.jsonl").open("w") as fh:
                for name in names:
                    fh.write(json.dumps(parsed_by_name[name]) + "\n")
            out_dir = output_dir / "relaxed_mpnn_out" / group_tag
            relaxed_commands.append(
                " ".join(
                    [
                        "/data/ziyang/mamba/envs/proteinmpnn/bin/python",
                        "third_party/ProteinMPNN/protein_mpnn_run.py",
                        f"--jsonl_path {group_dir / 'parsed_pdbs.jsonl'}",
                        f"--chain_id_jsonl {group_dir / 'chain_id.jsonl'}",
                        f"--fixed_positions_jsonl {group_dir / 'fixed_positions.jsonl'}",
                        f"--omit_AA_jsonl {group_dir / 'omit_AA.jsonl'}",
                        f"--out_folder {out_dir}",
                        f"--num_seq_per_target {sample_count}",
                        f"--sampling_temp {temp}",
                        f"--seed {config['random_seed']}",
                        "--batch_size 1",
                    ]
                )
            )
            command_rows.append(
                {
                    "command_index": command_count,
                    "branch": branch,
                    "mpnn_temperature": temp,
                    "num_seq_per_target": sample_count,
                    "job_count": int(len(sub)),
                    "planned_raw_samples": int(sub["planned_raw_samples"].sum()),
                    "parsed_pdbs_jsonl": str(group_dir / "parsed_pdbs.jsonl"),
                    "out_folder": str(out_dir),
                    "final_allocation_status": "audit_only_no_preallocated_quota",
                }
            )
            command_count += 1
    (output_dir / "relaxed_generation_runner_commands.sh").write_text(
        "\n".join(relaxed_commands) + "\n"
    )
    dry.write_csv(pd.DataFrame(command_rows), output_dir / "mpnn_relaxed_command_plan.csv")
    planned = int(jobs["planned_raw_samples"].sum()) if not jobs.empty else 0
    return {
        "relaxed_job_count": int(len(jobs)),
        "relaxed_command_count": int(command_count),
        "planned_relaxed_raw_samples": planned,
    }


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def build_score_inputs(manifest_entries: list[dict], output_dir: Path, shard_size: int = 5000) -> dict:
    pool = pd.read_csv(ROOT / "results/initial_design_generation/production_initial_pool/production_initial_pool_candidates_all.csv")
    fasta_dir = output_dir / "score_fastas"
    score_out_dir = output_dir / "mpnn_score_only_out"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    score_out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    commands: list[str] = [
        "# ProteinMPNN score-only commands for current production pool.",
        "# These commands are compute-heavy; run on an appropriate GPU node.",
        "",
    ]
    for entry in manifest_entries:
        target = normalize_target(entry["target"])
        safe = safe_target(target)
        sub = pool[pool["target"].map(normalize_target) == target].copy()
        if sub.empty:
            continue
        design_chain = str(entry["mpnn_design_chain_id"])
        chain_len = len(entry["_chains"][design_chain])
        bad_len = sub[sub["sequence"].astype(str).str.len() != chain_len]
        if not bad_len.empty:
            raise ValueError(
                f"{target}: {len(bad_len)} candidate sequences do not match design chain length {chain_len}"
            )
        backbone = ROOT / str(entry["structure_path"])
        for shard_index, start in enumerate(range(0, len(sub), shard_size)):
            shard = sub.iloc[start : start + shard_size]
            fasta_path = fasta_dir / f"{safe}_score_shard_{shard_index:04d}.fa"
            records = [
                (str(row.variant_id), str(row.sequence))
                for row in shard.itertuples(index=False)
            ]
            write_fasta(records, fasta_path)
            shard_out = score_out_dir / f"{safe}_score_shard_{shard_index:04d}"
            rows.append(
                {
                    "target": safe,
                    "backbone_id": entry["backbone_id"],
                    "structure_path": str(backbone),
                    "mpnn_design_chain_id": design_chain,
                    "fasta_path": str(fasta_path),
                    "score_output_dir": str(shard_out),
                    "candidate_count": len(shard),
                    "shard_index": shard_index,
                    "runner_ready": True,
                    "score_status": "pending_mpnn_score_only_run",
                }
            )
            commands.append(
                " ".join(
                    [
                        "/data/ziyang/mamba/envs/proteinmpnn/bin/python",
                        "third_party/ProteinMPNN/protein_mpnn_run.py",
                        "--score_only 1",
                        f"--path_to_fasta {fasta_path}",
                        f"--pdb_path {backbone}",
                        f"--pdb_path_chains {design_chain}",
                        f"--out_folder {shard_out}",
                        "--num_seq_per_target 1",
                        "--sampling_temp 0.1",
                        "--seed 20260523",
                        "--batch_size 1",
                    ]
                )
            )
    score_manifest = pd.DataFrame(rows)
    dry.write_csv(score_manifest, output_dir / "mpnn_score_only_input_manifest.csv")
    (output_dir / "score_only_runner_commands.sh").write_text("\n".join(commands) + "\n")
    return {
        "score_shard_count": int(len(score_manifest)),
        "score_candidate_count": int(score_manifest["candidate_count"].sum()) if not score_manifest.empty else 0,
    }


def write_blocked_report(output_dir: Path, errors: list[str] | None = None) -> None:
    lines = [
        "# MPNN Runner Input Generation Report",
        "",
        "Status: `blocked_pending_backbone_manifest`",
        "",
        "No executable ProteinMPNN JSONL inputs were generated.",
        "",
        "Provide a reviewed backbone manifest before running constrained MPNN.",
    ]
    if errors:
        lines += ["", "## Manifest Errors", ""]
        lines += [f"- {err}" for err in errors]
    (output_dir / "mpnn_runner_input_generation_report.md").write_text("\n".join(lines) + "\n")


def write_ready_report(output_dir: Path, summary: dict) -> None:
    is_sparse = str(summary.get("runner_mode", "")).startswith("constrained_sparse_")
    lines = [
        "# MPNN Runner Input Generation Report",
        "",
        "Status: `runner_inputs_ready`" if not is_sparse else f"Status: `{summary.get('runner_mode')}_inputs_ready`",
        "",
        f"- Job count: {summary['job_count']}",
        f"- JSONL directory: `{summary['jsonl_dir']}`",
        f"- Parsed PDB JSONL: `{summary['parsed_jsonl']}`",
        f"- Constrained temperature groups: {summary.get('temperature_group_count', 0)}",
        f"- Constrained command count: {summary.get('constrained_command_count', 0)}",
        f"- Planned constrained raw samples: {summary.get('planned_constrained_raw_samples', 0)}",
        f"- Relaxed counterfactual jobs: {summary.get('relaxed_job_count', 0)}",
        f"- Relaxed command count: {summary.get('relaxed_command_count', 0)}",
        f"- Planned relaxed raw samples: {summary.get('planned_relaxed_raw_samples', 0)}",
        f"- Score-only FASTA shards: {summary.get('score_shard_count', 0)}",
        f"- Score-only candidates: {summary.get('score_candidate_count', 0)}",
        "",
        "These inputs encode full fixed-position complements and per-position omit rules for constrained rescue generation.",
    ]
    if is_sparse:
        lines += [
            "This directory contains constrained sparse-rescue inputs only; it does not replace the original P0 score-only or relaxed runner inputs.",
        ]
    else:
        lines += [
            "They also include relaxed/MPNN-only counterfactual runner inputs and score-only FASTA shards for the current production pool.",
        ]
    lines += [
        "Run constrained generation with `constrained_generation_runner_commands.sh` after assigning appropriate compute resources.",
        "Run relaxed counterfactual generation with `relaxed_generation_runner_commands.sh`; outputs remain audit-only until evidence review and manual unlock." if not is_sparse else "No relaxed or score-only runner inputs are included in this sparse constrained directory.",
        "`example_runner_commands.sh` is for debugging only and does not preserve the planned temperature split.",
    ]
    (output_dir / "mpnn_runner_input_generation_report.md").write_text("\n".join(lines) + "\n")


def write_p0_ready_preflight(output_dir: Path, summary: dict) -> None:
    requirements = {
        "status": "runner_inputs_ready",
        "runner_ready": True,
        "safe_to_submit_to_proteinmpnn": True,
        "compute_boundary": "run only on appropriate compute node; builder did not execute ProteinMPNN",
        "runner_input_dir": display_path(output_dir),
        "constrained_generation": {
            "runner_jobs": int(summary["job_count"]),
            "command_count": int(summary.get("constrained_command_count", 0)),
            "planned_raw_samples": int(summary.get("planned_constrained_raw_samples", 0)),
            "command_file": display_path(output_dir / "constrained_generation_runner_commands.sh"),
        },
        "relaxed_counterfactual_generation": {
            "runner_jobs": int(summary.get("relaxed_job_count", 0)),
            "command_count": int(summary.get("relaxed_command_count", 0)),
            "planned_raw_samples": int(summary.get("planned_relaxed_raw_samples", 0)),
            "command_file": display_path(output_dir / "relaxed_generation_runner_commands.sh"),
            "final_allocation_status": "audit_only_no_preallocated_quota",
        },
        "score_only": {
            "fasta_shards": int(summary.get("score_shard_count", 0)),
            "candidate_count": int(summary.get("score_candidate_count", 0)),
            "command_file": display_path(output_dir / "score_only_runner_commands.sh"),
        },
        "runner_jsonl_inputs": {
            "parsed_pdbs_jsonl": display_path(output_dir / "jsonl/parsed_pdbs.jsonl"),
            "chain_id_jsonl": display_path(output_dir / "jsonl/chain_id.jsonl"),
            "fixed_positions_jsonl": display_path(output_dir / "jsonl/fixed_positions.jsonl"),
            "omit_AA_jsonl": display_path(output_dir / "jsonl/omit_AA.jsonl"),
        },
        "execution_notes": [
            "Use constrained_generation_runner_commands.sh, not example_runner_commands.sh, for true constrained generation.",
            "Run collect_p0_mpnn_results.py after ProteinMPNN jobs finish.",
            "Run validate_p0_mpnn_integration.py --require-mpnn-results before P0 evidence review.",
        ],
    }
    with (P0_DIR / "mpnn_runner_input_requirements.yaml").open("w") as fh:
        yaml.safe_dump(requirements, fh, sort_keys=False)

    lines = [
        "# ProteinMPNN Runner Input Preflight",
        "",
        "Status: `runner_inputs_ready`",
        "",
        "Reviewed AF3-derived backbones and runner-specific JSONL inputs have been generated.",
        "",
        "## Runner Inputs",
        "",
        f"- Runner input directory: `{display_path(output_dir)}`",
        f"- Constrained runner jobs: {summary['job_count']}",
        f"- Constrained command count: {summary.get('constrained_command_count', 0)}",
        f"- Planned constrained raw samples: {summary.get('planned_constrained_raw_samples', 0)}",
        f"- Relaxed counterfactual jobs: {summary.get('relaxed_job_count', 0)}",
        f"- Relaxed command count: {summary.get('relaxed_command_count', 0)}",
        f"- Planned relaxed raw samples: {summary.get('planned_relaxed_raw_samples', 0)}",
        f"- Score-only FASTA shards: {summary.get('score_shard_count', 0)}",
        f"- Score-only candidates: {summary.get('score_candidate_count', 0)}",
        "",
        "## Execution Boundary",
        "",
        "- These files are ready for ProteinMPNN execution on an appropriate compute node.",
        "- True ProteinMPNN score-only inference and constrained generation have not been run by this builder.",
        "- Relaxed counterfactual outputs remain audit-only and have no preallocated final-library quota.",
        "- Use `constrained_generation_runner_commands.sh`, not the debug-only example command, for constrained generation.",
    ]
    (P0_DIR / "mpnn_runner_preflight_report.md").write_text("\n".join(lines) + "\n")


def update_p0_comparison_runner_status(output_dir: Path, summary: dict) -> None:
    report = P0_DIR / "p0_mpnn_comparison_report.md"
    if not report.exists():
        return
    text = report.read_text()
    old = "- Scaffold plan rows are not direct ProteinMPNN runner inputs; runner-specific JSONL must be generated separately with reviewed backbones before compute."
    new = "\n".join(
        [
            f"- Scaffold plan rows are not direct ProteinMPNN runner inputs; runner-specific JSONL has now been generated under `{display_path(output_dir)}`.",
            f"- Current runner command plan encodes {summary.get('planned_constrained_raw_samples', 0)} planned constrained raw samples across {summary.get('constrained_command_count', 0)} commands.",
            f"- Relaxed counterfactual command plan encodes {summary.get('planned_relaxed_raw_samples', 0)} audit-only raw samples across {summary.get('relaxed_command_count', 0)} commands.",
        ]
    )
    changed = False
    if old in text:
        text = text.replace(old, new)
        changed = True
    elif "## Runner Input Status" not in text:
        text += "\n\n## Runner Input Status\n\n" + new + "\n"
        changed = True
    elif "runner-specific JSONL has now been generated" in text:
        changed = True
    text = text.replace(
        "- His-seeded relaxed mode requires explicit seed-set enumeration in runner inputs before execution.",
        "- His-seeded relaxed runner inputs enumerate top single-His seed sets from the production pool; outputs remain audit-only.",
    )
    text = text.replace(
        "- His-seeded relaxed mode still requires explicit seed-set enumeration before execution.",
        "- His-seeded relaxed runner inputs enumerate top single-His seed sets from the production pool; outputs remain audit-only.",
    )
    if changed:
        report.write_text(text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build ProteinMPNN runner JSONL inputs.")
    parser.add_argument("--manifest", default=None, help="Backbone input manifest YAML.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUT.relative_to(ROOT)))
    parser.add_argument("--write-template", action="store_true", help="Write manifest template.")
    parser.add_argument(
        "--constrained-sparse-smoke",
        action="store_true",
        help="Build only sparse constrained-MPNN smoke inputs with small rescue subsets.",
    )
    parser.add_argument(
        "--constrained-sparse",
        action="store_true",
        help="Build only repaired full sparse constrained-MPNN inputs with all constrained plan rows.",
    )
    parser.add_argument("--smoke-samples-per-job", type=int, default=32)
    parser.add_argument("--smoke-seed-limit-per-target", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = read_config()
    output_dir = ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    template_path = P0_DIR / "mpnn_backbone_manifest_template.yaml"
    if args.write_template or not args.manifest:
        write_manifest_template(config, template_path)
    if not args.manifest:
        write_blocked_report(output_dir)
        print(
            {
                "status": "blocked_pending_backbone_manifest",
                "manifest_template": display_path(template_path),
                "output_dir": display_path(output_dir),
            }
        )
        return

    manifest = load_yaml(ROOT / args.manifest)
    entries, errors = validate_manifest(manifest, config)
    if errors:
        write_blocked_report(output_dir, errors)
        raise RuntimeError("; ".join(errors))
    if args.constrained_sparse_smoke and args.constrained_sparse:
        raise SystemExit("--constrained-sparse-smoke and --constrained-sparse are mutually exclusive")
    summary = build_constrained_inputs(
        entries,
        output_dir,
        sparse_smoke=args.constrained_sparse_smoke,
        sparse_full=args.constrained_sparse,
        smoke_samples_per_job=args.smoke_samples_per_job,
        smoke_seed_limit_per_target=args.smoke_seed_limit_per_target,
    )
    if args.constrained_sparse_smoke or args.constrained_sparse:
        summary.update(
            {
                "relaxed_job_count": 0,
                "relaxed_command_count": 0,
                "planned_relaxed_raw_samples": 0,
                "score_shard_count": 0,
                "score_candidate_count": 0,
                "runner_mode": "constrained_sparse_smoke"
                if args.constrained_sparse_smoke
                else "constrained_sparse_full",
            }
        )
    else:
        summary.update(build_relaxed_inputs(entries, output_dir))
        summary.update(build_score_inputs(entries, output_dir))
    write_ready_report(output_dir, summary)
    if not args.constrained_sparse_smoke and not args.constrained_sparse:
        write_p0_ready_preflight(output_dir, summary)
        update_p0_comparison_runner_status(output_dir, summary)
    print(summary)


if __name__ == "__main__":
    main()
