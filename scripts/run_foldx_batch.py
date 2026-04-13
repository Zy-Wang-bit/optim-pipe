#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, glob, yaml, subprocess, csv
from concurrent.futures import ProcessPoolExecutor, as_completed

def load_cfg(p): 
    with open(p, "r") as f: 
        return yaml.safe_load(f)

def ensure_dir(p): 
    os.makedirs(p, exist_ok=True); return p

def run(cmd, cwd=None, env=None, log=None):
    env2 = dict(os.environ, **(env or {}))
    p = subprocess.run(cmd, cwd=cwd, env=env2,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if log:
        open(log+".out","w").write(p.stdout or "")
        open(log+".err","w").write(p.stderr or "")
    p.check_returncode()
    return p

def find_batches(root):
    return [p for p in sorted(glob.glob(os.path.join(root, "*", "batch_*"))) if os.path.isdir(p)]

def copy_template_if_needed(bdir, repaired_root, pid):
    src = os.path.join(repaired_root, f"{pid}_Repair.pdb")
    dst = os.path.join(bdir, f"{pid}_Repair.pdb")
    if not os.path.exists(src):
        raise FileNotFoundError(f"缺少模板: {src}，请先运行 repair_pdbs.py")
    if not os.path.exists(dst):
        import shutil; shutil.copy2(src, dst)
    return os.path.basename(dst)

def list_mutant_pdbs(bdir, pid):
    muts = []
    for x in glob.glob(os.path.join(bdir, f"{pid}_Repair_*.pdb")):
        base = os.path.basename(x)
        if base.startswith("WT_"): 
            continue
        if base.endswith("_Repair.pdb"):
            continue
        muts.append(base)
    return sorted(set(muts))

def find_interaction_file(bdir, mpdb_noext, ph):
    import re, glob
    labels=[]; s=str(ph)
    labels += [s, f"{ph:.1f}", s.replace(".","_"), str(int(round(ph))), s.replace(".","")]
    labels = list(dict.fromkeys(labels))
    pats=[]
    for lab in labels:
        pats += [
            f"Interaction_AC_{mpdb_noext}_{lab}_AC.fxout",
            f"Interaction_AC_{mpdb_noext}_{lab}.fxout",
            f"Interaction_*{mpdb_noext}*_{lab}*_AC.fxout",
            f"Interaction_*{mpdb_noext}*_{lab}*.fxout",
        ]
    num = re.sub(r"\D","",str(ph))
    pats.append(f"Interaction_*{mpdb_noext}*{num}*.fxout")
    for pat in pats:
        hits = sorted(glob.glob(os.path.join(bdir, pat)))
        if hits: return hits[0]
    return None

def parse_ac_tsv(fp, mpdb_base, left_set, right_set):
    with open(fp, encoding="utf-8", errors="ignore") as f:
        lines = [l.rstrip("\n") for l in f]
    hdr_i=None
    for i, line in enumerate(lines):
        if "\t" in line and "Pdb" in line and "Group1" in line and "Group2" in line and "Interaction Energy" in line:
            hdr_i=i; break
    if hdr_i is None: return None
    header = lines[hdr_i].split("\t")
    c_pdb = header.index("Pdb"); c_g1 = header.index("Group1"); c_g2 = header.index("Group2")
    c_ie  = header.index("Interaction Energy")
    total=0.0; found=False
    for line in lines[hdr_i+1:]:
        if not line.strip(): continue
        parts = line.split("\t")
        if len(parts)<=max(c_pdb,c_g1,c_g2,c_ie): continue
        pdb_field = os.path.basename(parts[c_pdb].strip()).lstrip("./")
        if pdb_field != mpdb_base: continue
        g1 = parts[c_g1].strip(); g2 = parts[c_g2].strip()
        if not ((g1 in left_set and g2 in right_set) or (g2 in left_set and g1 in right_set)):
            continue
        try: val = float(parts[c_ie])
        except: continue
        total += val; found=True
    return total if found else None

def load_wt_cache(cfg):
    repaired_dir = os.path.join(cfg["paths"]["foldx_dir"], "repaired")
    cache_csv = os.path.join(repaired_dir, "WT_ac.csv")
    cache = {}
    if os.path.exists(cache_csv):
        with open(cache_csv) as f:
            r = csv.DictReader(f)
            for row in r:
                try:
                    cache[(row["pdb_id"], float(row["pH"]))] = float(row["dG"])
                except Exception:
                    continue
    return cache, cache_csv, repaired_dir

def append_wt_cache(cache_csv, pid, ph, dG):
    header_needed = not os.path.exists(cache_csv)
    with open(cache_csv, "a") as f:
        if header_needed:
            f.write("pdb_id,pH,dG\n")
        f.write(f"{pid},{ph},{dG}\n")

def process_one_batch(bdir, cfg, wt_cache, wt_cache_csv, repaired_dir):
    foldx = os.path.abspath(cfg["paths"]["foldx_bin"])
    omp   = str(cfg.get("resources", {}).get("omp_threads", 1))
    env   = {"OMP_NUM_THREADS": omp}

    pid = os.path.basename(os.path.dirname(bdir))  # foldx/batches/<pid>/batch_xxx
    template_pdb  = copy_template_if_needed(bdir, repaired_dir, pid)
    template_noext = os.path.splitext(template_pdb)[0]          # e.g. n1-0_Repair
    template_base  = os.path.basename(template_pdb)             # e.g. n1-0_Repair.pdb

    mutfile = os.path.join(bdir, "individual_list.txt")
    if not os.path.exists(mutfile) or os.path.getsize(mutfile) == 0:
        print(f"[SKIP] {bdir}: individual_list.txt 缺失或为空"); 
        return

    ph_points = cfg["foldx"].get("ph_points", [7.4, 6.0])
    nruns     = int(cfg["foldx"].get("number_of_runs", 1))
    groups    = cfg["foldx"].get("analyse_groups", "A,B;C")
    left_set  = set(groups.split(";")[0].split(","))
    right_set = set(groups.split(";")[1].split(","))

    # 1) BuildModel（按 pH）——生成变体
    for ph in ph_points:
        tag = f"build_pH{ph}"
        cmd = [
            foldx,
            f"--command=BuildModel",
            f"--pdb={template_pdb}",
            f"--mutant-file={os.path.basename(mutfile)}",
            f"--pH={ph}",
            f"--numberOfRuns={nruns}",
            f"--output-file={tag}",
        ]
        try:
            run(cmd, cwd=bdir, env=env, log=os.path.join(bdir, tag))
        except subprocess.CalledProcessError as e:
            print(f"[WARN] BuildModel 失败: {bdir} pH={ph}: {e}")

    # 2) 列出突变体
    mpdbs = list_mutant_pdbs(bdir, pid)
    if not mpdbs:
        print(f"[WARN] 未发现突变体PDB：{bdir}（检查 BuildModel 输出）")

    # 3) WT（模板）——优先读缓存；缺失才计算并写回缓存（只算一次/每PDB）
    wt_ac = {}
    for ph in ph_points:
        key = (pid, float(ph))
        if key in wt_cache:
            wt_ac[ph] = wt_cache[key]
            continue
        # 缓存没有：在 repaired 目录计算一次并缓存
        out_tag = f"AC_{pid}_WT_{ph}"
        cmd = [
            foldx,
            f"--command=AnalyseComplex",
            f"--pdb={template_base}",
            f"--analyseComplexChains={groups}",
            f"--pH={ph}",
            f"--output-file={out_tag}",
        ]
        try:
            run(cmd, cwd=repaired_dir, env=env, log=os.path.join(repaired_dir, out_tag))
        except subprocess.CalledProcessError as e:
            print(f"[WARN] WT AnalyseComplex 失败: {pid} pH={ph}: {e}")
            continue
        fx = find_interaction_file(repaired_dir, template_noext, ph)
        if not fx:
            print(f"[WARN] 未找到 WT Interaction 文件: {pid} pH={ph}")
            continue
        val = parse_ac_tsv(fx, template_base, left_set, right_set)
        if val is None:
            print(f"[WARN] 未能解析 WT Interaction TSV: {fx}")
            continue
        wt_ac[ph] = val
        wt_cache[key] = val
        append_wt_cache(wt_cache_csv, pid, ph, val)

    # 4) 各变体在两种 pH 下的相互作用能 —— 并行化 AnalyseComplex
    ac = {}  # (mpdb_base, ph) -> energy
    _ac_max_procs = int(cfg.get("resources", {}).get("foldx_max_procs", 24))

    def _run_one_ac(args_tuple):
        """单个 AnalyseComplex 任务（子进程入口）"""
        mpdb, ph, bdir_, foldx_, groups_, env_, left_set_, right_set_ = args_tuple
        mpdb_base = os.path.basename(mpdb)
        out_tag   = f"AC_{os.path.splitext(mpdb_base)[0]}_{ph}"
        cmd = [
            foldx_,
            f"--command=AnalyseComplex",
            f"--pdb={mpdb}",
            f"--analyseComplexChains={groups_}",
            f"--pH={ph}",
            f"--output-file={out_tag}",
        ]
        try:
            run(cmd, cwd=bdir_, env=env_, log=os.path.join(bdir_, out_tag))
        except subprocess.CalledProcessError:
            return (mpdb_base, ph, None)
        fx = find_interaction_file(bdir_, os.path.splitext(mpdb_base)[0], ph)
        if not fx:
            return (mpdb_base, ph, None)
        val = parse_ac_tsv(fx, mpdb_base, left_set_, right_set_)
        return (mpdb_base, ph, val)

    ac_tasks = []
    for ph in ph_points:
        for mpdb in mpdbs:
            ac_tasks.append((mpdb, ph, bdir, foldx, groups, env, left_set, right_set))

    if ac_tasks:
        n_workers = min(_ac_max_procs, len(ac_tasks))
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            for mpdb_base, ph, val in pool.map(_run_one_ac, ac_tasks):
                if val is not None:
                    ac[(mpdb_base, ph)] = val

    # 5) 汇总（含 WT/ΔΔ 项）
    out_csv = os.path.join(bdir, "foldx_summary.csv")
    with open(out_csv, "w") as w:
        w.write("mpdb,dG_pH7_4,dG_pH6_0,WT_dG_pH7_4,WT_dG_pH6_0,ddG_pH7_4,ddG_pH6_0,delta,delta_wt,delta_delta\n")
        for mpdb in mpdbs:
            mp = os.path.basename(mpdb)
            v74 = ac.get((mp, 7.4), "")
            v60 = ac.get((mp, 6.0), "")
            w74 = wt_ac.get(7.4, "")
            w60 = wt_ac.get(6.0, "")
            dd74 = (v74 - w74) if v74 != "" and w74 != "" else ""
            dd60 = (v60 - w60) if v60 != "" and w60 != "" else ""
            delt = (v60 - v74) if v60 != "" and v74 != "" else ""
            delw = (w60 - w74) if w60 != "" and w74 != "" else ""
            ddd  = (delt - delw) if delt != "" and delw != "" else ""
            w.write(f"{mp},{v74},{v60},{w74},{w60},{dd74},{dd60},{delt},{delw},{ddd}\n")
    print(f"[OK] summary -> {out_csv}")

def _build_one_toplevel(args_tuple):
    """模块级：单个 BuildModel 任务。"""
    bdir_, ph_, cfg_ = args_tuple
    foldx_ = os.path.abspath(cfg_["paths"]["foldx_bin"])
    omp_   = str(cfg_.get("resources", {}).get("omp_threads", 1))
    env_   = {"OMP_NUM_THREADS": omp_}
    repaired_ = os.path.join(cfg_["paths"]["foldx_dir"], "repaired")
    pid_ = os.path.basename(os.path.dirname(bdir_))
    template_pdb_ = copy_template_if_needed(bdir_, repaired_, pid_)
    mutfile_ = os.path.join(bdir_, "individual_list.txt")
    if not os.path.exists(mutfile_) or os.path.getsize(mutfile_) == 0:
        return (bdir_, ph_, "skip")
    nruns_ = int(cfg_["foldx"].get("number_of_runs", 1))
    tag_ = f"build_pH{ph_}"
    cmd_ = [foldx_, f"--command=BuildModel", f"--pdb={template_pdb_}",
            f"--mutant-file={os.path.basename(mutfile_)}",
            f"--pH={ph_}", f"--numberOfRuns={nruns_}", f"--output-file={tag_}"]
    try:
        run(cmd_, cwd=bdir_, env=env_, log=os.path.join(bdir_, tag_))
        return (bdir_, ph_, "ok")
    except Exception as e:
        return (bdir_, ph_, str(e))


def _ac_one_toplevel(args_tuple):
    """模块级：单个 AnalyseComplex 任务。"""
    mpdb_, ph_, bdir_, foldx_, groups_, env_, left_set_, right_set_ = args_tuple
    mpdb_base_ = os.path.basename(mpdb_)
    out_tag_ = f"AC_{os.path.splitext(mpdb_base_)[0]}_{ph_}"
    cmd_ = [foldx_, f"--command=AnalyseComplex", f"--pdb={mpdb_}",
            f"--analyseComplexChains={groups_}", f"--pH={ph_}",
            f"--output-file={out_tag_}"]
    try:
        run(cmd_, cwd=bdir_, env=env_, log=os.path.join(bdir_, out_tag_))
    except subprocess.CalledProcessError:
        return (bdir_, mpdb_base_, ph_, None)
    fx_ = find_interaction_file(bdir_, os.path.splitext(mpdb_base_)[0], ph_)
    if not fx_:
        return (bdir_, mpdb_base_, ph_, None)
    val_ = parse_ac_tsv(fx_, mpdb_base_, left_set_, right_set_)
    return (bdir_, mpdb_base_, ph_, val_)


def _worker(args):
    """子进程入口：处理单个 batch。"""
    bdir, cfg, repaired_dir = args
    # 每个子进程独立加载 WT 缓存（避免共享状态竞争）
    wt_cache, wt_cache_csv, _ = load_wt_cache(cfg)
    try:
        process_one_batch(bdir, cfg, wt_cache, wt_cache_csv, repaired_dir)
        return (bdir, None)
    except Exception as e:
        return (bdir, str(e))


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    cfg = load_cfg(cfg_path)
    batches_root = os.path.join(cfg["paths"]["foldx_dir"], "batches")
    repaired_root = os.path.join(cfg["paths"]["foldx_dir"], "repaired")
    ensure_dir(batches_root); ensure_dir(repaired_root)

    batches = find_batches(batches_root)
    if not batches:
        print(f"[ERR] 未发现批次目录：{batches_root}/<pid>/batch_*，请先运行 make_mutlist_chunk.py")
        sys.exit(1)

    # 支持 batch 范围过滤（多节点并行时，每个节点只处理自己的 batch）
    batch_range = cfg.get("foldx_batch_range")
    if batch_range:
        lo, hi = int(batch_range[0]), int(batch_range[1])
        batches = [b for b in batches
                   if lo <= int(os.path.basename(b).replace("batch_", "")) <= hi]
        print(f"[FoldX] batch 范围过滤: [{lo}, {hi}] → {len(batches)} batches")

    max_procs = int(cfg.get("resources", {}).get("foldx_max_procs", 1))

    # 先预计算所有 PDB 的 WT 基线（串行，只需做一次）
    wt_cache, wt_cache_csv, repaired_dir = load_wt_cache(cfg)
    pdb_ids = set()
    for b in batches:
        pid = os.path.basename(os.path.dirname(b))
        pdb_ids.add(pid)

    foldx = os.path.abspath(cfg["paths"]["foldx_bin"])
    omp   = str(cfg.get("resources", {}).get("omp_threads", 1))
    env   = {"OMP_NUM_THREADS": omp}
    groups = cfg["foldx"].get("analyse_groups", "A,B;C")
    left_set  = set(groups.split(";")[0].split(","))
    right_set = set(groups.split(";")[1].split(","))

    for pid in sorted(pdb_ids):
        template_base = f"{pid}_Repair.pdb"
        template_noext = f"{pid}_Repair"
        for ph in cfg["foldx"].get("ph_points", [7.4, 6.0]):
            key = (pid, float(ph))
            if key in wt_cache:
                continue
            out_tag = f"AC_{pid}_WT_{ph}"
            cmd = [foldx, f"--command=AnalyseComplex",
                   f"--pdb={template_base}", f"--analyseComplexChains={groups}",
                   f"--pH={ph}", f"--output-file={out_tag}"]
            try:
                run(cmd, cwd=repaired_root, env=env, log=os.path.join(repaired_root, out_tag))
            except subprocess.CalledProcessError as e:
                print(f"[WARN] WT AnalyseComplex 失败: {pid} pH={ph}: {e}")
                continue
            fx = find_interaction_file(repaired_root, template_noext, ph)
            if fx:
                val = parse_ac_tsv(fx, template_base, left_set, right_set)
                if val is not None:
                    wt_cache[key] = val
                    append_wt_cache(wt_cache_csv, pid, ph, val)

    print(f"[FoldX] WT 基线就绪: {len(wt_cache)} 项")
    print(f"[FoldX] 批次: {len(batches)} | 并发: {max_procs} 进程")

    # ── Phase A: BuildModel —— 并行（每个 batch×pH 一个 FoldX 进程）───────────
    print(f"[FoldX] Phase A: BuildModel ({len(batches)} batches × {len(cfg['foldx'].get('ph_points',[7.4,6.0]))} pH)")

    pass  # _build_one_toplevel defined at module level

    ph_points = cfg["foldx"].get("ph_points", [7.4, 6.0])
    build_tasks = []
    for b in batches:
        for ph in ph_points:
            build_tasks.append((b, ph, cfg))

    n_build_workers = min(max_procs, len(build_tasks))
    with ProcessPoolExecutor(max_workers=n_build_workers) as pool:
        for i, (bdir_, ph_, status) in enumerate(pool.map(_build_one_toplevel, build_tasks)):
            if status not in ("ok", "skip"):
                print(f"[WARN] BuildModel {os.path.basename(bdir_)} pH={ph_}: {status}")
            if (i + 1) % 10 == 0 or i + 1 == len(build_tasks):
                print(f"[FoldX] BuildModel 进度: {i+1}/{len(build_tasks)}")

    print("[FoldX] Phase A 完成: BuildModel 全部就绪")

    # ── Phase B: AnalyseComplex —— 并行（所有 batch 的所有变体×pH 汇总后统一并行）
    print("[FoldX] Phase B: AnalyseComplex")
    groups = cfg["foldx"].get("analyse_groups", "A,B;C")
    left_set  = set(groups.split(";")[0].split(","))
    right_set = set(groups.split(";")[1].split(","))

    foldx_abs = os.path.abspath(cfg["paths"]["foldx_bin"])
    ac_tasks = []
    batch_mpdbs = {}  # bdir -> [mpdb list]
    for b in batches:
        pid = os.path.basename(os.path.dirname(b))
        mpdbs = list_mutant_pdbs(b, pid)
        batch_mpdbs[b] = mpdbs
        for ph in ph_points:
            for mpdb in mpdbs:
                ac_tasks.append((mpdb, ph, b, foldx_abs, groups, {"OMP_NUM_THREADS": "1"}, left_set, right_set))

    print(f"[FoldX] AC 任务: {len(ac_tasks)} | 并发: {max_procs}")

    # {(bdir, mpdb_base, ph) -> energy}
    ac_results = {}
    if ac_tasks:
        n_ac_workers = min(max_procs, len(ac_tasks))
        done_count = 0
        with ProcessPoolExecutor(max_workers=n_ac_workers) as pool:
            for bdir_, mpdb_base_, ph_, val_ in pool.map(_ac_one_toplevel, ac_tasks):
                if val_ is not None:
                    ac_results[(bdir_, mpdb_base_, ph_)] = val_
                done_count += 1
                if done_count % 500 == 0 or done_count == len(ac_tasks):
                    print(f"[FoldX] AC 进度: {done_count}/{len(ac_tasks)}")

    print("[FoldX] Phase B 完成: AnalyseComplex 全部就绪")

    # ── Phase C: 汇总每个 batch 的 summary CSV ─────────────────────────────────
    for b in batches:
        pid = os.path.basename(os.path.dirname(b))
        mpdbs = batch_mpdbs.get(b, [])
        out_csv = os.path.join(b, "foldx_summary.csv")
        with open(out_csv, "w") as w:
            w.write("mpdb,dG_pH7_4,dG_pH6_0,WT_dG_pH7_4,WT_dG_pH6_0,ddG_pH7_4,ddG_pH6_0,delta,delta_wt,delta_delta\n")
            for mpdb in mpdbs:
                mp = os.path.basename(mpdb)
                v74 = ac_results.get((b, mp, 7.4), "")
                v60 = ac_results.get((b, mp, 6.0), "")
                w74 = wt_cache.get((pid, 7.4), "")
                w60 = wt_cache.get((pid, 6.0), "")
                dd74 = (v74 - w74) if v74 != "" and w74 != "" else ""
                dd60 = (v60 - w60) if v60 != "" and w60 != "" else ""
                delt = (v60 - v74) if v60 != "" and v74 != "" else ""
                delw = (w60 - w74) if w60 != "" and w74 != "" else ""
                ddd  = (delt - delw) if delt != "" and delw != "" else ""
                w.write(f"{mp},{v74},{v60},{w74},{w60},{dd74},{dd60},{delt},{delw},{ddd}\n")

    print(f"[FoldX] 全部完成: {len(batches)} 批次")


if __name__ == "__main__":
    main()