#!/usr/bin/env python3
"""ESM 评分 batch2 (90k 双链序列)"""
import os, torch
import pandas as pd
import esm

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Device: {device}")
if device == "cuda":
    print(f"GPU: {torch.cuda.get_device_name(0)}")

model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
model = model.eval().to(device)
bc = alphabet.get_batch_converter()

df = pd.read_csv("workspace/1e62_r4_dual/mpnn_outputs_batch2/mpnn_batch2_dedup.csv")
print(f"Loaded {len(df)} sequences")

def extract_ab(seq):
    parts = str(seq).split("/")
    return "".join(parts[:2]) if len(parts) >= 2 else seq

seqs = [extract_ab(s) for s in df["sequence"]]

bs = 64
scores = []
for i in range(0, len(seqs), bs):
    batch = seqs[i:i+bs]
    data = [(f"seq_{j}", s) for j, s in enumerate(batch)]
    _, _, tokens = bc(data)
    tokens = tokens.to(device)
    with torch.no_grad():
        results = model(tokens, repr_layers=[33])
    logits = results["logits"]
    for k in range(len(batch)):
        L = len(batch[k])
        logp = torch.log_softmax(logits[k, 1:L+1], dim=-1)
        token_ids = tokens[k, 1:L+1]
        avg = logp.gather(1, token_ids.unsqueeze(1)).mean().item()
        scores.append(avg)
    if (i // bs) % 100 == 0:
        print(f"  ESM progress: {min(i+bs, len(seqs))}/{len(seqs)}")

df["esm_avg_logprob"] = scores
out = "workspace/1e62_r4_dual/esm_scores/mpnn_batch2_dedup_esm.csv"
df.to_csv(out, index=False)
print(f"Done: {out} ({len(df)} rows)")
