#!/usr/bin/env python3
"""Phase B+C only: AnalyseComplex + Summary for completed BM batches.

Usage:
  python scripts/run_ac_only.py configs/config_sdab_r2_v2_node0.yaml
"""
import os, sys, glob, yaml, csv
from concurrent.futures import ProcessPoolExecutor

# Reuse functions from run_foldx_batch.py
sys.path.insert(0, os.path.dirname(__file__))
from run_foldx_batch import (
    load_cfg, find_batches, copy_template_if_needed, list_mutant_pdbs,
    find_interaction_file, parse_ac_tsv, load_wt_cache, run,
)


def _ac_one(args):
    """Module-level AC worker for ProcessPoolExecutor."""
    bdir, mpdb, ph, foldx, groups, env, left_set, right_set = args
    mpdb_base = os.path.basename(mpdb)
    out_tag = f"AC_{os.path.splitext(mpdb_base)[0]}_{ph}"

    cmd = [
        foldx,
        "--command=AnalyseComplex",
        f"--pdb={mpdb_base}",
        f"--analyseComplexChains={groups}",
        f"--pH={ph}",
        f"--output-file={out_tag}",
    ]
    try:
        run(cmd, cwd=bdir, env=env, log=os.path.join(bdir, out_tag))
    except Exception:
        return (bdir, mpdb_base, ph, None)

    fx = find_interaction_file(bdir, os.path.splitext(mpdb_base)[0], ph)
    if not fx:
        return (bdir, mpdb_base, ph, None)
    val = parse_ac_tsv(fx, mpdb_base, left_set, right_set)
    return (bdir, mpdb_base, ph, val)


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "configs/config_sdab_r2_v2.yaml"
    cfg = load_cfg(cfg_path)

    foldx_dir = cfg["paths"]["foldx_dir"]
    foldx_bin = os.path.abspath(cfg["paths"]["foldx_bin"])
    groups = cfg["foldx"].get("analyse_groups", "A;B")
    ph_points = cfg["foldx"].get("ph_points", [7.4, 6.0])
    omp = str(cfg.get("resources", {}).get("omp_threads", 1))
    max_procs = int(cfg.get("resources", {}).get("foldx_max_procs", 128))
    env = {"OMP_NUM_THREADS": omp}

    left_set = set(groups.split(";")[0].split(","))
    right_set = set(groups.split(";")[1].split(","))

    # Find batches with optional range filter
    batches_root = os.path.join(foldx_dir, "batches")
    batches = find_batches(batches_root)

    batch_range = cfg.get("foldx_batch_range")
    if batch_range:
        lo, hi = int(batch_range[0]), int(batch_range[1])
        batches = [b for b in batches
                   if lo <= int(os.path.basename(b).replace("batch_", "")) <= hi]
        print(f"[AC] batch 范围过滤: [{lo}, {hi}] → {len(batches)} batches")

    # Filter to completed BM batches only
    complete_batches = []
    skipped = 0
    for b in batches:
        pid = os.path.basename(os.path.dirname(b))
        ml = os.path.join(b, "individual_list.txt")
        if not os.path.exists(ml):
            skipped += 1
            continue
        muts = sum(1 for _ in open(ml))
        pdbs = len(glob.glob(os.path.join(b, "*.pdb")))
        expected = muts * 2 + 1
        if pdbs >= expected and expected > 1:
            complete_batches.append(b)
        else:
            skipped += 1

    print(f"[AC] 完成 batch: {len(complete_batches)}, 跳过: {skipped}")

    # Load WT cache
    wt_cache, wt_cache_csv, repaired_dir = load_wt_cache(cfg)
    print(f"[AC] WT 基线: {len(wt_cache)} 项")

    # Build AC task list
    ac_tasks = []
    batch_mpdbs = {}  # batch -> [mpdb_base, ...]
    for b in complete_batches:
        pid = os.path.basename(os.path.dirname(b))
        copy_template_if_needed(b, repaired_dir, pid)
        mpdbs = list_mutant_pdbs(b, pid)
        batch_mpdbs[b] = mpdbs
        for ph in ph_points:
            for mpdb in mpdbs:
                ac_tasks.append((b, mpdb, ph, foldx_bin, groups, env, left_set, right_set))

    print(f"[AC] 总 AC 任务: {len(ac_tasks)} ({len(complete_batches)} batch × {len(ph_points)} pH)")
    print(f"[AC] 并发: {max_procs} 进程")

    # Phase B: Run all AC in parallel
    ac_results = {}
    done_count = 0
    with ProcessPoolExecutor(max_workers=max_procs) as pool:
        for bdir, mpdb_base, ph, val in pool.map(_ac_one, ac_tasks):
            if val is not None:
                ac_results[(bdir, mpdb_base, ph)] = val
            done_count += 1
            if done_count % 500 == 0 or done_count == len(ac_tasks):
                print(f"[AC] 进度: {done_count}/{len(ac_tasks)}")

    print(f"[AC] Phase B 完成: {len(ac_results)} 有效结果")

    # Phase C: Generate per-batch summary CSV
    for b in complete_batches:
        pid = os.path.basename(os.path.dirname(b))
        mpdbs = batch_mpdbs.get(b, [])
        out_csv = os.path.join(b, "foldx_summary.csv")
        with open(out_csv, "w") as w:
            w.write("mpdb,dG_pH7_4,dG_pH6_0,WT_dG_pH7_4,WT_dG_pH6_0,ddG_pH7_4,ddG_pH6_0,delta,delta_wt,delta_delta\n")
            for mpdb in mpdbs:
                v74 = ac_results.get((b, mpdb, 7.4), "")
                v60 = ac_results.get((b, mpdb, 6.0), "")
                w74 = wt_cache.get((pid, 7.4), "")
                w60 = wt_cache.get((pid, 6.0), "")
                dd74 = (v74 - w74) if v74 != "" and w74 != "" else ""
                dd60 = (v60 - w60) if v60 != "" and w60 != "" else ""
                delt = (v60 - v74) if v60 != "" and v74 != "" else ""
                delw = (w60 - w74) if w60 != "" and w74 != "" else ""
                ddd = (delt - delw) if delt != "" and delw != "" else ""
                w.write(f"{mpdb},{v74},{v60},{w74},{w60},{dd74},{dd60},{delt},{delw},{ddd}\n")

    print(f"[AC] Phase C 完成: {len(complete_batches)} summary CSVs")
    print("[AC] ALL DONE")


if __name__ == "__main__":
    main()
