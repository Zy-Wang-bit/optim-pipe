#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, glob, yaml, subprocess, csv

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
    foldx = cfg["paths"]["foldx_bin"]          # 例如 "foldx"
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

    # 4) 各变体在两种 pH 下的相互作用能
    ac = {}  # (mpdb_base, ph) -> energy
    for ph in ph_points:
        for mpdb in mpdbs:
            mpdb_base = os.path.basename(mpdb)
            out_tag   = f"AC_{os.path.splitext(mpdb_base)[0]}_{ph}"
            cmd = [
                foldx,
                f"--command=AnalyseComplex",
                f"--pdb={mpdb}",
                f"--analyseComplexChains={groups}",
                f"--pH={ph}",
                f"--output-file={out_tag}",
            ]
            try:
                run(cmd, cwd=bdir, env=env, log=os.path.join(bdir, out_tag))
            except subprocess.CalledProcessError as e:
                print(f"[WARN] AnalyseComplex 失败: {bdir} {mpdb_base} pH={ph}: {e}")
                continue

            fx = find_interaction_file(bdir, os.path.splitext(mpdb_base)[0], ph)
            if not fx:
                print(f"[WARN] 未找到 Interaction 文件: {bdir} {mpdb_base} pH={ph}")
                continue
            val = parse_ac_tsv(fx, mpdb_base, left_set, right_set)
            if val is None:
                print(f"[WARN] 未能解析 Interaction TSV: {fx}")
                continue
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

def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    cfg = load_cfg(cfg_path)
    batches_root = os.path.join(cfg["paths"]["foldx_dir"], "batches")
    repaired_root = os.path.join(cfg["paths"]["foldx_dir"], "repaired")
    ensure_dir(batches_root); ensure_dir(repaired_root)

    wt_cache, wt_cache_csv, repaired_dir = load_wt_cache(cfg)

    batches = find_batches(batches_root)
    if not batches:
        print(f"[ERR] 未发现批次目录：{batches_root}/<pid>/batch_*，请先运行 make_mutlist_chunk.py")
        sys.exit(1)

    for b in batches:
        try:
            process_one_batch(b, cfg, wt_cache, wt_cache_csv, repaired_dir)
        except Exception as e:
            print("[ERROR]", b, e)

if __name__ == "__main__":
    main()