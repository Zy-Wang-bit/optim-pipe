#!/usr/bin/env python3
import os, sys, glob, yaml, subprocess, shutil

def run(bin_path, args, env=None, cwd=None):
    env = dict(os.environ, **(env or {}))
    env.setdefault("OMP_NUM_THREADS", "1")
    try:
        return subprocess.run(
            [bin_path] + args,
            env=env,
            cwd=cwd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"[FoldX ERROR] cmd: {bin_path} {' '.join(args)} (cwd={cwd})\n")
        if e.stdout:
            sys.stderr.write(e.stdout.decode(errors="ignore") + "\n")
        if e.stderr:
            sys.stderr.write(e.stderr.decode(errors="ignore") + "\n")
        raise

def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    foldx    = os.path.abspath(cfg["paths"]["foldx_bin"])
    pdb_dir  = cfg["paths"]["pdb_dir"]
    out_dir  = os.path.join(cfg["paths"]["foldx_dir"], "repaired")
    os.makedirs(out_dir, exist_ok=True)

    # 可选 rotabase
    rotabase = None
    if "foldx" in cfg and isinstance(cfg["foldx"], dict):
        rotabase = cfg["foldx"].get("rotabase")
    env = {}
    if rotabase:
        env["FOLDX_ROTABASE"] = os.path.abspath(rotabase)

    for src in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
        base = os.path.basename(src)                # 仅文件名（FoldX需要这样）
        tag  = os.path.splitext(base)[0]
        dst  = os.path.join(out_dir, base)

        # 复制到工作目录（FoldX只在cwd查找PDB）
        if os.path.abspath(src) != os.path.abspath(dst):
            shutil.copy2(src, dst)

        print(f"[Repair] {tag}")
        args = [
            "--command=RepairPDB",
            f"--pdb={base}",                        # 关键：只给文件名
            f"--output-file={tag}",
        ]
        run(foldx, args, env=env, cwd=out_dir)

        # （可选）验证输出是否存在
        repaired = os.path.join(out_dir, f"{tag}_Repair.pdb")
        if not os.path.exists(repaired):
            sys.stderr.write(f"[WARN] 未发现输出 {repaired}，请检查 FoldX 输出信息。\n")

    print("[OK] repaired ->", out_dir)

if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    main(cfg)