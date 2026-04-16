#!/usr/bin/env python3
"""关键突变重要性评分 + 最小化变体构建

基于 MD 数据 (contacts, hbond, RMSF, pairwise_dist) 对 69 个突变位点评分，
分层分类，并构建只含关键突变的最小化变体。
"""

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
CAND_CSV = REPO / "workspace/sdab_r2_v2/results/final_candidates.csv"
MD_DIR = REPO / "experiments/sdab_variants/md"
OUTPUT_DIR = MD_DIR / "analysis_plots"

CDR = {
    "FR1": (1, 25), "CDR1": (26, 35), "FR2": (36, 49), "CDR2": (50, 65),
    "FR3": (66, 94), "CDR3": (95, 112), "FR4": (113, 122),
}

# WT sequence
WT_AA = {}


def get_region(pos):
    for name, (s, e) in CDR.items():
        if s <= pos <= e:
            return name
    return "?"


def load_wt():
    global WT_AA
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("wt", str(REPO / "experiments/sdab/repaired/WT_Repair.pdb"))
    aa3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
              'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
              'TYR': 'Y', 'VAL': 'V'}
    for res in s[0]["A"].get_residues():
        if res.id[0] != " ":
            continue
        WT_AA[res.id[1]] = aa3to1.get(res.resname, "X")


def parse_mutations(mut_str):
    muts = {}
    for m in str(mut_str).split(","):
        match = re.match(r"([A-Z])A(\d+)([A-Z])", m.strip())
        if match:
            muts[int(match.group(2))] = (match.group(1), match.group(3))
    return muts


def load_contacts(vid, ph):
    path = MD_DIR / vid / f"pH_{ph}" / "analysis" / "contacts_frequency.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    ab = df[df["ab_chain"].str.contains("chain_A", na=False)]
    return ab.groupby("ab_resid")["frequency"].sum().to_dict()


def load_hbond(vid, ph):
    """Return {resid: total_occupancy} for antibody residues in interface hbonds."""
    path = MD_DIR / vid / f"pH_{ph}" / "analysis" / "hbond_occupancy.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    result = {}
    for col_chain, col_resid in [("donor_chain", "donor_resid"), ("acceptor_chain", "acceptor_resid")]:
        if col_chain in df.columns:
            ab = df[df[col_chain].str.contains("chain_A", na=False)]
            for resid, occ in ab.groupby(col_resid)["occupancy"].sum().items():
                result[int(resid)] = result.get(int(resid), 0) + occ
    return result


def load_rmsf(vid, ph):
    path = MD_DIR / vid / f"pH_{ph}" / "analysis" / "rmsf_per_residue.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    df_a = df[df["chain"].str.contains("chain_A", na=False)]
    return dict(zip(df_a["resid"].astype(int), df_a["rmsf"]))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    load_wt()

    cand = pd.read_csv(CAND_CSV)
    variants = [f"sdab_v2_{i:04d}" for i in range(15)]

    # Parse mutations
    all_positions = Counter()
    variant_muts = {}
    for i, (_, r) in enumerate(cand.iterrows()):
        vid = variants[i]
        muts = parse_mutations(r["mutations"])
        variant_muts[vid] = muts
        for pos in muts:
            all_positions[pos] += 1

    # ---- Score each position ----
    rows = []
    for pos in sorted(all_positions.keys()):
        freq = all_positions[pos]
        region = get_region(pos)
        is_cdr = "CDR" in region
        wt_aa = WT_AA.get(pos, "?")

        # Most common mutation at this position
        mut_counter = Counter()
        for vid in variants:
            if pos in variant_muts[vid]:
                mut_counter[variant_muts[vid][pos][1]] += 1
        common_mut = mut_counter.most_common(1)[0][0] if mut_counter else "?"
        is_his = common_mut == "H"

        # MD metrics across variants that have this mutation
        contact_74_vals = []
        contact_60_vals = []
        d_contact_vals = []
        hbond_74_vals = []
        d_hbond_vals = []
        rmsf_74_vals = []
        d_rmsf_vals = []

        for vid in variants:
            if pos not in variant_muts[vid]:
                continue
            c74 = load_contacts(vid, "7.4").get(pos, 0)
            c60 = load_contacts(vid, "6.0").get(pos, 0)
            h74 = load_hbond(vid, "7.4").get(pos, 0)
            h60 = load_hbond(vid, "6.0").get(pos, 0)
            r74 = load_rmsf(vid, "7.4").get(pos, None)
            r60 = load_rmsf(vid, "6.0").get(pos, None)

            contact_74_vals.append(c74)
            contact_60_vals.append(c60)
            d_contact_vals.append(c60 - c74)
            hbond_74_vals.append(h74)
            d_hbond_vals.append(h60 - h74)
            if r74 is not None and r60 is not None:
                rmsf_74_vals.append(r74)
                d_rmsf_vals.append(r60 - r74)

        # Aggregate
        contact_74_med = np.median(contact_74_vals) if contact_74_vals else 0
        d_contact_med = np.median(d_contact_vals) if d_contact_vals else 0
        hbond_74_med = np.median(hbond_74_vals) if hbond_74_vals else 0
        d_hbond_med = np.median(d_hbond_vals) if d_hbond_vals else 0
        rmsf_74_med = np.median(rmsf_74_vals) if rmsf_74_vals else None
        d_rmsf_med = np.median(d_rmsf_vals) if d_rmsf_vals else None

        is_interface = contact_74_med > 0.05 or (np.median(contact_60_vals) > 0.05 if contact_60_vals else False)

        # Importance tier
        if is_his:
            tier = "必须保留"
        elif is_interface and abs(d_contact_med) > 0.1:
            tier = "必须保留"
        elif is_interface and hbond_74_med > 0.05:
            tier = "必须保留"
        elif is_interface:
            tier = "推荐保留"
        elif rmsf_74_med is not None and rmsf_74_med < 1.0:
            tier = "推荐保留"  # buried core
        elif freq >= 10 and is_cdr:
            tier = "推荐保留"
        elif freq < 5 and not is_interface:
            tier = "可去除"
        else:
            tier = "可选"

        rows.append({
            "position": pos,
            "wt_aa": wt_aa,
            "common_mut": common_mut,
            "mutation": f"{wt_aa}{pos}{common_mut}",
            "region": region,
            "is_cdr": is_cdr,
            "is_his": is_his,
            "frequency": f"{freq}/15",
            "freq_num": freq,
            "contact_74_median": round(contact_74_med, 3),
            "d_contact_median": round(d_contact_med, 3),
            "hbond_74_median": round(hbond_74_med, 3),
            "d_hbond_median": round(d_hbond_med, 3),
            "rmsf_74_median": round(rmsf_74_med, 2) if rmsf_74_med is not None else None,
            "d_rmsf_median": round(d_rmsf_med, 2) if d_rmsf_med is not None else None,
            "is_interface": is_interface,
            "tier": tier,
        })

    score_df = pd.DataFrame(rows).sort_values(["tier", "position"])
    score_df.to_csv(OUTPUT_DIR / "critical_mutations_report.csv", index=False)

    # ---- Summary by tier ----
    print("=" * 80)
    print("突变位点重要性分层")
    print("=" * 80)
    for tier in ["必须保留", "推荐保留", "可选", "可去除"]:
        sub = score_df[score_df["tier"] == tier]
        positions = sub["position"].tolist()
        cdr_count = sub["is_cdr"].sum()
        fr_count = len(sub) - cdr_count
        print(f"\n[{tier}] {len(sub)} 位点 (CDR:{cdr_count}, FR:{fr_count})")
        for _, r in sub.iterrows():
            interface_tag = "界面" if r["is_interface"] else ""
            his_tag = "His" if r["is_his"] else ""
            tags = " ".join(filter(None, [his_tag, interface_tag]))
            print(f"  {r['mutation']:>6s}  {r['region']:<5s}  freq={r['frequency']}  "
                  f"contact={r['contact_74_median']:.2f}  Δcontact={r['d_contact_median']:+.2f}  "
                  f"hbond={r['hbond_74_median']:.2f}  {tags}")

    # ---- CDR/FR distribution ----
    print("\n" + "=" * 80)
    print("突变分布 (CDR vs FR)")
    print("=" * 80)
    for region in CDR:
        sub = score_df[score_df["region"] == region]
        if len(sub) == 0:
            continue
        tiers = sub["tier"].value_counts().to_dict()
        tier_str = ", ".join(f"{t}:{n}" for t, n in tiers.items())
        print(f"  {region:5s}: {len(sub):2d} 位点  [{tier_str}]")

    # ---- Build minimal variants ----
    # FR2 保守位点: 实验确认不可突变 (F37V 人源化失败)
    FORBIDDEN_POSITIONS = {37, 47, 49, 50}

    essential = set(score_df[score_df["tier"] == "必须保留"]["position"])
    recommended = set(score_df[score_df["tier"] == "推荐保留"]["position"])
    minimal_positions = (essential | recommended) - FORBIDDEN_POSITIONS

    print(f"\n最小化变体突变数: {len(minimal_positions)} (必须:{len(essential)} + 推荐:{len(recommended)})")

    # For each of top-5 ranked variants, build minimal version
    ranking = pd.read_csv(OUTPUT_DIR / "final_ranking.csv")
    top5 = ranking.head(5)

    print("\n" + "=" * 80)
    print("最小化变体 (Top 5 设计的精简版)")
    print("=" * 80)

    minimal_rows = []
    for _, r in top5.iterrows():
        vid = r["variant_id"]
        muts = variant_muts[vid]
        full_n = len(muts)
        # Keep only essential + recommended mutations
        kept = {pos: aa for pos, aa in muts.items() if pos in minimal_positions}
        removed = {pos: aa for pos, aa in muts.items() if pos not in minimal_positions}
        # Build sequence
        seq = list(WT_AA.get(i, "X") for i in range(1, 123))
        for pos, (wt, mut) in kept.items():
            seq[pos - 1] = mut
        seq_str = "".join(seq)
        kept_muts = ",".join(f"{wt}{pos}{mut}" for pos, (wt, mut) in sorted(kept.items()))

        minimal_rows.append({
            "variant_id": f"{vid}_minimal",
            "base_variant": vid,
            "sequence": seq_str,
            "mutations": kept_muts,
            "n_mutations": len(kept),
            "n_removed": len(removed),
            "identity_to_wt": round((122 - len(kept)) / 122 * 100, 1),
        })
        print(f"  {vid}_minimal: {len(kept)} muts (去除 {len(removed)}, WT identity {(122-len(kept))/122*100:.1f}%)")

    minimal_df = pd.DataFrame(minimal_rows)
    minimal_df.to_csv(OUTPUT_DIR / "minimal_variants.csv", index=False)

    print(f"\n输出: {OUTPUT_DIR / 'critical_mutations_report.csv'}")
    print(f"输出: {OUTPUT_DIR / 'minimal_variants.csv'}")


if __name__ == "__main__":
    main()
