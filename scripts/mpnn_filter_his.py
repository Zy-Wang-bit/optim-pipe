#!/usr/bin/env python3
import os, sys, csv, yaml, glob

def keep_his(seq, region, min_his=1):
    a,b = region
    seg = seq[a-1:b]
    return seg.count("H") >= min_his

def process(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    in_dir = cfg["paths"]["mpnn_out_dir"]
    out_dir = os.path.join(in_dir, "filtered")
    os.makedirs(out_dir, exist_ok=True)
    region = cfg["design"]["region"]
    min_his = cfg["his_bias"]["min_his_in_region"]

    for fn in sorted(glob.glob(os.path.join(in_dir, "*.csv"))):
        out_fn = os.path.join(out_dir, os.path.basename(fn).replace(".csv","_his.csv"))
        kept = 0; total = 0
        with open(fn) as f, open(out_fn,"w",newline="") as g:
            r = csv.DictReader(f)
            w = csv.DictWriter(g, fieldnames=r.fieldnames)
            w.writeheader()
            for row in r:
                total += 1
                if keep_his(row["sequence"], region, min_his):
                    w.writerow(row); kept += 1
        print(f"[{os.path.basename(fn)}] kept {kept}/{total} -> {out_fn}")

if __name__=="__main__":
    cfg = sys.argv[1] if len(sys.argv)>1 else "configs/config.yaml"
    process(cfg)