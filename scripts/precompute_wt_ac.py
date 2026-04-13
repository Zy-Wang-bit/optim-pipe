#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, glob, yaml, subprocess, csv

def run(cmd, cwd=None, env=None, log=None):
    env2 = dict(os.environ)            # 继承当前环境（含 PATH）
    if env:
        env2.update(env)               # 只覆盖给定键，例如 OMP_NUM_THREADS
    p = subprocess.run(cmd, cwd=cwd, env=env2, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if log:
        open(log+".out","w").write(p.stdout or "")
        open(log+".err","w").write(p.stderr or "")
    p.check_returncode()
    return p

def find_interaction_file_wt(dirpath, pid, ph):
    import glob, os, re, time

    # 生成 pH 备选标签
    s = str(ph)
    labels = [s, f"{ph:.1f}", s.replace(".", "_"), str(int(round(ph))), s.replace(".", "")]
    labels = [x for x in dict.fromkeys(labels) if x]  # 去重&去空

    patterns = []
    # 最常见：Interaction_AC_<pid>_WT_<lab>_AC.fxout
    for lab in labels:
        patterns.append(os.path.join(dirpath, f"Interaction_AC_{pid}_WT_{lab}_AC.fxout"))
    # 更宽松的备选（不同 FoldX 组合顺序）
    for lab in labels:
        patterns.append(os.path.join(dirpath, f"Interaction_*{pid}*WT*{lab}*_AC.fxout"))
        patterns.append(os.path.join(dirpath, f"Interaction_*{pid}*WT*{lab}*.fxout"))

    hits = []
    for pat in patterns:
        for fp in glob.glob(pat):
            try:
                hits.append((os.path.getmtime(fp), fp))
            except FileNotFoundError:
                pass
    if not hits:
        return None
    # 取最新
    hits.sort(reverse=True)
    return hits[0][1]

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

def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    foldx = os.path.abspath(cfg["paths"]["foldx_bin"])
    repaired_dir = os.path.join(cfg["paths"]["foldx_dir"], "repaired")
    os.makedirs(repaired_dir, exist_ok=True)

    groups = cfg["foldx"].get("analyse_groups","A,B;C")
    left_set  = set(groups.split(";")[0].split(","))
    right_set = set(groups.split(";")[1].split(","))
    ph_points = cfg["foldx"].get("ph_points",[7.4,6.0])

    cache_csv = os.path.join(repaired_dir, "WT_ac.csv")
    cache = {}  # (pid, ph) -> dG
    if os.path.exists(cache_csv):
        with open(cache_csv) as f:
            r = csv.DictReader(f)
            for row in r:
                cache[(row["pdb_id"], float(row["pH"]))] = float(row["dG"])

    # 遍历 repaired 下的模板
    pdbs = sorted(glob.glob(os.path.join(repaired_dir, "*_Repair.pdb")))
    for pdb_path in pdbs:
        pid = os.path.basename(pdb_path).replace("_Repair.pdb","")
        for ph in ph_points:
            key = (pid, float(ph))
            if key in cache:  # 已有缓存
                continue
            out_tag = f"AC_{pid}_WT_{ph}"
            cmd = [
                foldx,
                f"--command=AnalyseComplex",
                f"--pdb={os.path.basename(pdb_path)}",
                f"--analyseComplexChains={groups}",
                f"--pH={ph}",
                f"--output-file={out_tag}",
            ]
            try:
                run(cmd, cwd=repaired_dir, env={"OMP_NUM_THREADS":str(cfg.get("resources",{}).get("omp_threads",1))}, log=os.path.join(repaired_dir,out_tag))
            except subprocess.CalledProcessError as e:
                print(f"[WARN] WT AnalyseComplex 失败: {pid} pH={ph}: {e}")
                continue
            fx = find_interaction_file_wt(repaired_dir, pid, ph)
            if not fx:
                print(f"[WARN] 未找到 WT Interaction 文件: {pid} pH={ph}")
                continue
            val = parse_ac_tsv(fx, f"{pid}_Repair.pdb", left_set, right_set)
            if val is None:
                print(f"[WARN] 解析失败: {fx}")
                continue
            cache[key]=val
            # 追加写入
            mode = "a" if os.path.exists(cache_csv) else "w"
            with open(cache_csv, mode) as f:
                if mode=="w":
                    f.write("pdb_id,pH,dG\n")
                f.write(f"{pid},{ph},{val}\n")
            print(f"[WT] {pid} pH={ph} dG={val:.3f}")

    print(f"[OK] WT 缓存 -> {cache_csv}")

if __name__=="__main__":
    cfg = sys.argv[1] if len(sys.argv)>1 else "configs/config.yaml"
    main(cfg)