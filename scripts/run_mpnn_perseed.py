#!/usr/bin/env python3
"""Per-seed MPNN design: 每个 His seed 组合作为固定模版, MPNN 设计其余位置.

用法:
  python scripts/run_mpnn_perseed.py configs/config_sdab_r2_v2.yaml \
    --num-per-seed 100 --shards 5
"""
import os, sys, json, glob, argparse, shutil, yaml
import numpy as np
import pandas as pd
from itertools import combinations
import random

# ProteinMPNN alphabet: H is at index 6
MPNN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'
HIS_IDX = MPNN_ALPHABET.index('H')  # 6


def load_wt_seq(cfg, chain):
    wt_dir = cfg["paths"]["wt_dir"]
    mapping = cfg["paths"].get("wt_files", {})
    p = os.path.join(wt_dir, mapping[chain]) if chain in mapping else None
    if not p or not os.path.exists(p):
        raise FileNotFoundError(f"WT fasta for chain {chain} not found")
    seq = ""
    with open(p) as f:
        for line in f:
            if line.startswith(">"): continue
            seq += line.strip()
    return seq


def generate_seeds(cfg):
    """Generate 142 His seed combinations based on user-specified rules."""
    high = {51, 52, 55, 105, 110}       # pH ratio >= 2.0
    low  = {27, 28, 31, 53, 56}          # pH ratio < 1.5
    exp_all = {27, 28, 29, 30, 31, 51, 52, 53, 55, 56, 105, 108, 109, 110}
    exp_keep = exp_all - low              # 9 positions
    top3 = {52, 105, 51}

    seeds = []

    # S1: exp_keep 单点 (9)
    for p in sorted(exp_keep):
        seeds.append({"positions": [p], "source": "his_seed_S1"})

    # S2: exp_keep 全配对 (36) + (61,99) (1)
    for combo in combinations(sorted(exp_keep), 2):
        seeds.append({"positions": list(combo), "source": "his_seed_S2"})
    seeds.append({"positions": [61, 99], "source": "his_seed_S2"})

    # S3: exp_keep 三组合, >=2 高分 + (61,82,99)
    for combo in combinations(sorted(exp_keep), 3):
        if sum(1 for p in combo if p in high) >= 2:
            seeds.append({"positions": list(combo), "source": "his_seed_S3"})
    seeds.append({"positions": [61, 82, 99], "source": "his_seed_S3"})

    # S4: exp_keep 四组合, >=3 高分 + 含 top3
    for combo in combinations(sorted(exp_keep), 4):
        n_high = sum(1 for p in combo if p in high)
        has_top3 = any(p in top3 for p in combo)
        if n_high >= 3 and has_top3:
            seeds.append({"positions": list(combo), "source": "his_seed_S4"})

    # Count
    counts = {}
    for s in seeds:
        counts[s["source"]] = counts.get(s["source"], 0) + 1
    for src, n in sorted(counts.items()):
        print(f"  {src}: {n}")
    print(f"  Total seeds: {len(seeds)}")

    return seeds


def setup_perseed_pdb_dir(pdb_src, seeds, out_dir):
    """Create symlinks: sdab_seed_NNN.pdb -> original PDB, one per seed."""
    pdb_dir = os.path.join(out_dir, "_perseed_pdbs")
    os.makedirs(pdb_dir, exist_ok=True)
    src_abs = os.path.abspath(pdb_src)
    seed_pdbs = []
    for i, seed in enumerate(seeds):
        name = f"seed_{i:04d}"
        dst = os.path.join(pdb_dir, f"{name}.pdb")
        if os.path.exists(dst):
            os.remove(dst)
        os.symlink(src_abs, dst)
        seed_pdbs.append((name, dst, seed))
    return pdb_dir, seed_pdbs


def build_perseed_jsonl(seed_pdbs, design_chain, antigen_chain, chain_lengths, out_dir):
    """Build JSONL files for per-seed MPNN: chain_id, fixed_positions, bias_by_res."""
    jsonl_dir = os.path.join(out_dir, "_mpnn_config")
    os.makedirs(jsonl_dir, exist_ok=True)

    chain_jsonl = os.path.join(jsonl_dir, "chain_id.jsonl")
    fixed_jsonl = os.path.join(jsonl_dir, "fixed_positions.jsonl")
    bias_jsonl  = os.path.join(jsonl_dir, "bias_by_res.jsonl")

    CJ = {}  # {name: [masked_chains, visible_chains]}
    FP = {}
    BR = {}

    L_design = chain_lengths[design_chain]
    L_antigen = chain_lengths[antigen_chain]

    for name, pdb_path, seed in seed_pdbs:
        his_positions = seed["positions"]  # 1-based

        # chain_id: masked=design chain, visible=antigen
        CJ[name] = [[design_chain], [antigen_chain]]

        # fixed_positions: fix antigen entirely, chain A is all designable
        FP[name] = {
            design_chain: [],  # nothing fixed in design chain
            antigen_chain: list(range(1, L_antigen + 1)),  # fix all antigen
        }

        # bias_by_res: +100 at His positions for amino acid H
        bias_A = []
        for pos in range(1, L_design + 1):
            vec = [0.0] * 21
            if pos in his_positions:
                vec[HIS_IDX] = 100.0
            bias_A.append(vec)

        bias_B = [[0.0] * 21 for _ in range(L_antigen)]

        BR[name] = {design_chain: bias_A, antigen_chain: bias_B}

    # Write all as single JSON dicts (MPNN expected format)
    with open(chain_jsonl, "w") as f:
        f.write(json.dumps(CJ) + "\n")

    with open(fixed_jsonl, "w") as f:
        f.write(json.dumps(FP) + "\n")

    with open(bias_jsonl, "w") as f:
        f.write(json.dumps(BR) + "\n")

    print(f"[OK] JSONL -> {jsonl_dir}")
    return os.path.abspath(chain_jsonl), os.path.abspath(fixed_jsonl), os.path.abspath(bias_jsonl)


def get_chain_lengths(pdb_path):
    """Get chain lengths from PDB file."""
    from collections import OrderedDict
    chains = OrderedDict()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                ch = line[21]
                resid = line[22:27].strip()
                chains.setdefault(ch, set()).add(resid)
    return {ch: len(resids) for ch, resids in chains.items()}


def run_mpnn(pdb_dir, chain_jsonl, fixed_jsonl, bias_jsonl, out_dir,
             runner, num_per_seed, shards, temps, base_seed):
    """Run ProteinMPNN across all seed PDBs.

    Step 1: parse PDBs into JSONL using helper script.
    Step 2: run MPNN with --jsonl_path.
    """
    import subprocess

    # Step 1: Parse PDB directory into JSONL
    parser_script = os.path.join(os.path.dirname(runner), "helper_scripts", "parse_multiple_chains.py")
    parsed_jsonl = os.path.join(out_dir, "_parsed_pdbs.jsonl")
    print(f"[MPNN] Parsing {pdb_dir} -> {parsed_jsonl}")
    subprocess.run([
        sys.executable, parser_script,
        "--input_path", pdb_dir,
        "--output_path", parsed_jsonl,
    ], check=True)

    # Step 2: Run MPNN per shard
    per_shard = max(1, num_per_seed // shards)
    for i in range(shards):
        temp = temps[i % len(temps)]
        seed = base_seed + i
        shard_dir = os.path.join(out_dir, f"shard_{i:03d}")
        os.makedirs(shard_dir, exist_ok=True)

        json.dump({"sampling_temp": temp, "seed": seed, "per_shard": per_shard},
                  open(os.path.join(shard_dir, "marker.json"), "w"))

        cmd = [
            sys.executable, runner,
            "--jsonl_path", parsed_jsonl,
            "--out_folder", shard_dir,
            "--chain_id_jsonl", chain_jsonl,
            "--fixed_positions_jsonl", fixed_jsonl,
            "--bias_by_res_jsonl", bias_jsonl,
            "--num_seq_per_target", str(per_shard),
            "--sampling_temp", str(temp),
            "--seed", str(seed),
            "--batch_size", "1",
        ]
        print(f"[MPNN] shard {i}: temp={temp}, seed={seed}, per_seed={per_shard}")
        subprocess.run(cmd, check=True)


def collect_results(out_dir, seeds, wt_seq, design_chain):
    """Collect MPNN FASTA outputs, verify His positions, build CSV."""
    rows = []
    for shard_dir in sorted(glob.glob(os.path.join(out_dir, "shard_*"))):
        shard_name = os.path.basename(shard_dir)
        meta = {}
        mpath = os.path.join(shard_dir, "marker.json")
        if os.path.exists(mpath):
            meta = json.load(open(mpath))

        for fa in glob.glob(os.path.join(shard_dir, "**", "*.fa"), recursive=True):
            # Parse seed index from filename
            fname = os.path.basename(fa)
            # Find which seed this belongs to
            parent = os.path.basename(os.path.dirname(fa))
            if parent == "seqs":
                parent = os.path.basename(os.path.dirname(os.path.dirname(fa)))

            # Read sequences
            seqs = []
            with open(fa) as f:
                s = ""
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if s: seqs.append(s)
                        s = ""
                    else:
                        s += line
                if s: seqs.append(s)

            # Match to seed by filename
            seed_idx = None
            for i, seed in enumerate(seeds):
                seed_name = f"seed_{i:04d}"
                if seed_name in fa:
                    seed_idx = i
                    break

            if seed_idx is None:
                continue

            seed = seeds[seed_idx]
            his_positions = seed["positions"]

            for seq in seqs[1:]:  # skip first (WT)
                # Extract chain A sequence (first chain)
                chain_a_seq = seq.split("/")[0] if "/" in seq else seq

                # Verify His at seed positions
                his_ok = all(chain_a_seq[p-1] == "H" for p in his_positions
                             if p-1 < len(chain_a_seq))
                if not his_ok:
                    continue  # skip sequences where His wasn't placed

                # Generate mutation codes
                muts = []
                for j, (w, m) in enumerate(zip(wt_seq, chain_a_seq)):
                    if w != m:
                        muts.append(f"{w}{design_chain}{j+1}{m}")

                rows.append({
                    "sequence": chain_a_seq,
                    "mutations": ",".join(muts),
                    "source": seed["source"],
                    "seed_id": seed_idx,
                    "his_positions": ";".join(str(p) for p in his_positions),
                    "shard": shard_name,
                    "sampling_temp": meta.get("sampling_temp"),
                })

    df = pd.DataFrame(rows)
    n_raw = len(df)
    df = df.drop_duplicates(subset=["sequence"], keep="first")
    print(f"[Collect] raw={n_raw}, dedup={len(df)}, dup_rate={100*(1-len(df)/max(n_raw,1)):.1f}%")

    # Stats per source
    for src, grp in df.groupby("source"):
        print(f"  {src}: {len(grp)} seqs")

    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--num-per-seed", type=int, default=100)
    ap.add_argument("--shards", type=int, default=5)
    ap.add_argument("--temps", default="0.10,0.15,0.20,0.25,0.30")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    design_chain = cfg["design"]["chain"]
    ag_chains = cfg["interface"]["ag_chains"]
    antigen_chain = ag_chains[0] if isinstance(ag_chains, list) else ag_chains
    out_root = cfg["paths"]["mpnn_out_dir"]
    pdb_dir = cfg["paths"]["pdb_dir"]
    os.makedirs(out_root, exist_ok=True)

    # Load WT
    wt_seq = load_wt_seq(cfg, design_chain)
    print(f"[WT] chain={design_chain}, length={len(wt_seq)}")

    # Find template PDB
    pdbs = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    if not pdbs:
        raise FileNotFoundError(f"No PDB files in {pdb_dir}")
    template_pdb = pdbs[0]
    print(f"[Template] {template_pdb}")

    # Get chain lengths
    chain_lengths = get_chain_lengths(template_pdb)
    print(f"[Chain lengths] {chain_lengths}")

    # Generate seeds
    print("\n[Generating seeds]")
    seeds = generate_seeds(cfg)

    # Save seed definitions
    seed_csv = os.path.join(out_root, "seed_definitions.csv")
    seed_rows = []
    for i, s in enumerate(seeds):
        seed_rows.append({
            "seed_id": i,
            "positions": ";".join(str(p) for p in s["positions"]),
            "source": s["source"],
            "mutations": ",".join(f"{wt_seq[p-1]}{design_chain}{p}H" for p in s["positions"]
                                  if p-1 < len(wt_seq) and wt_seq[p-1] != "H"),
        })
    pd.DataFrame(seed_rows).to_csv(seed_csv, index=False)
    print(f"[OK] Seed definitions -> {seed_csv}")

    # Setup PDB symlinks
    print("\n[Setting up per-seed PDB directory]")
    perseed_pdb_dir, seed_pdbs = setup_perseed_pdb_dir(template_pdb, seeds, out_root)
    print(f"  Created {len(seed_pdbs)} PDB symlinks in {perseed_pdb_dir}")

    # Build JSONL files
    print("\n[Building JSONL files]")
    chain_jsonl, fixed_jsonl, bias_jsonl = build_perseed_jsonl(
        seed_pdbs, design_chain, antigen_chain, chain_lengths, out_root)

    # Run MPNN
    runner = cfg["paths"].get("mpnn_runner", "third_party/ProteinMPNN/protein_mpnn_run.py")
    temps = [float(x) for x in args.temps.split(",")]

    print(f"\n[Running MPNN] {len(seeds)} seeds × {args.num_per_seed} seqs = ~{len(seeds)*args.num_per_seed} total")
    run_mpnn(perseed_pdb_dir, chain_jsonl, fixed_jsonl, bias_jsonl, out_root,
             runner, args.num_per_seed, args.shards, temps, args.seed)

    # Collect results
    print("\n[Collecting results]")
    df = collect_results(out_root, seeds, wt_seq, design_chain)

    # Save
    out_csv = os.path.join(out_root, "mpnn_perseed_all.csv")
    df.to_csv(out_csv, index=False)
    print(f"\n[DONE] {len(df)} sequences -> {out_csv}")

    # Also save to screening dir for pipeline compatibility
    screen_dir = os.path.join(cfg["paths"]["results_dir"], "screening")
    os.makedirs(screen_dir, exist_ok=True)
    shutil.copy2(out_csv, os.path.join(screen_dir, "mpnn_all.csv"))


if __name__ == "__main__":
    main()
