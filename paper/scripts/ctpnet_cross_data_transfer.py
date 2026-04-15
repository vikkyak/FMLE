#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from scipy.stats import gmean
from pathlib import Path
from torch.utils.data import DataLoader, TensorDataset
import re
from sklearn.model_selection import train_test_split

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
if hasattr(torch, "set_float32_matmul_precision"):
    torch.set_float32_matmul_precision("high")

# -------------------------
# Args
# -------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--train_dir", required=True)
ap.add_argument("--test_dir",  required=True)
ap.add_argument("--out_dir",   required=True)
ap.add_argument("--panel",     required=True, help="gene_panel_2000.csv from R")
ap.add_argument("--seed", type=int, default=4905)
ap.add_argument("--source_tag", default="A", help="e.g., A or B")
ap.add_argument("--target_tag", default="B", help="e.g., B or A or C")
args = ap.parse_args()

np.random.seed(args.seed)
torch.manual_seed(args.seed)

train_dir = Path(args.train_dir).expanduser()
test_dir  = Path(args.test_dir).expanduser()
out_dir   = Path(args.out_dir).expanduser()
out_dir.mkdir(parents=True, exist_ok=True)

# -------------------------
# Load data
# -------------------------
# =============================================================================
# rna_train = pd.read_csv(train_dir / "rna_train_full.csv", index_col=0)
# adt_train = pd.read_csv(train_dir / "adt_train.csv", index_col=0)
# rna_test  = pd.read_csv(test_dir  / "rna_test_full.csv",  index_col=0)
# 
# =============================================================================
rna_train = pd.read_csv(train_dir / "rna_train.csv", index_col=0)
adt_train = pd.read_csv(train_dir / "adt_train.csv", index_col=0)
rna_test  = pd.read_csv(test_dir  / "rna_test.csv",  index_col=0)


assert (rna_train.index == adt_train.index).all()

genes_use = [g for g in rna_train.columns if g in set(rna_test.columns)]
print("genes_use:", len(genes_use))  # will be 1051 here
assert list(rna_train[genes_use].columns) == list(rna_test[genes_use].columns)
assert len(set(rna_train.columns)) == rna_train.shape[1]
assert len(set(rna_test.columns))  == rna_test.shape[1]
assert (rna_train.index == adt_train.index).all()
print("Train genes:", rna_train[genes_use].shape)
print("Test genes :", rna_test[genes_use].shape)

print("Train mean/std:",
      rna_train[genes_use].values.mean(),
      rna_train[genes_use].values.std())

print("Test mean/std:",
      rna_test[genes_use].values.mean(),
      rna_test[genes_use].values.std())
corr = np.corrcoef(
    rna_train[genes_use].mean(0),
    rna_test[genes_use].mean(0)
)[0,1]

print("Gene-wise mean correlation:", corr)

# -------------------------
# Leakage check
# -------------------------
overlap = set(rna_train.index).intersection(set(rna_test.index))
print("Overlapping cell IDs:", len(overlap))
if len(overlap) > 0:
    print(list(overlap)[:10])
    raise RuntimeError("Train/Test cell leakage detected")

# -------------------------
# FIXED GENE PANEL ALIGNMENT (DO NOT ERROR; PAD WITH 0)
# -------------------------
panel = pd.read_csv(args.panel)["gene"].astype(str).tolist()
print("Panel size:", len(panel))

panel = pd.read_csv(args.panel)["gene"].astype(str).tolist()
train_cols = list(rna_train.columns)

missing_train = sorted(set(panel) - set(train_cols))
print("Missing train examples:", missing_train[:30])
test_cols = list(rna_test.columns)
missing_test = sorted(set(panel) - set(test_cols))
print("Missing test examples:", missing_test[:30])


missing_tr = [g for g in panel if g not in rna_train.columns]
missing_te = [g for g in panel if g not in rna_test.columns]
print("Missing in train:", len(missing_tr))
print("Missing in test :", len(missing_te))

# Pad missing genes with zeros (transfer-safe)
rna_train = rna_train.reindex(columns=panel, fill_value=0)
rna_test  = rna_test.reindex(columns=panel, fill_value=0)

assert list(rna_train.columns) == panel
assert list(rna_test.columns)  == panel
print("Final aligned genes:", rna_train.shape[1])

# -------------------------
# Align RNA/ADT on train
# -------------------------
common_cells = rna_train.index.intersection(adt_train.index)
if len(common_cells) < 100:
    raise RuntimeError(f"Too few common cells: {len(common_cells)}")

rna_train = rna_train.loc[common_cells]
adt_train = adt_train.loc[common_cells]
assert (rna_train.index == adt_train.index).all()

protein_names = list(adt_train.columns)
gene_names    = list(rna_train.columns)

train_genes = set(pd.read_csv(train_dir/"rna_train.csv", nrows=1).columns[1:])
test_genes  = set(pd.read_csv(test_dir/"rna_test.csv",  nrows=1).columns[1:])
panel_genes = set(pd.read_csv(args.panel)["gene"].astype(str))

print("raw train∩test:", len(train_genes & test_genes))
print("panel∩train  :", len(panel_genes & train_genes))
print("panel∩test   :", len(panel_genes & test_genes))



print("Proteins:", len(protein_names))
print("Done preprocessing ✓")

# -------------------------
# cTPnet-style ADT transform (CLR per cell) - TRAIN ONLY
# -------------------------
def norm_adt_cellsxprot(Y_cellsxprot: np.ndarray) -> np.ndarray:
    Yp = Y_cellsxprot.astype(np.float64) + 1.0
    gm = gmean(Yp, axis=1, keepdims=True)  # per cell
    return np.log(Yp / gm)

# IMPORTANT: adt_train is cells x proteins already (from your CSVs)
Y_all = norm_adt_cellsxprot(adt_train.values)

# -------------------------
# Train/Val split (TRAIN ONLY)
# -------------------------
idx = np.arange(len(rna_train))
tr_idx, va_idx = train_test_split(
    idx, test_size=0.15, random_state=args.seed, shuffle=True
)

X_tr = torch.tensor(rna_train.values[tr_idx], dtype=torch.float32)
Y_tr = torch.tensor(Y_all[tr_idx], dtype=torch.float32)

X_va = torch.tensor(rna_train.values[va_idx], dtype=torch.float32)
Y_va = torch.tensor(Y_all[va_idx], dtype=torch.float32)

X_train_full = torch.tensor(rna_train.values, dtype=torch.float32)
X_test_full  = torch.tensor(rna_test.values,  dtype=torch.float32)

# -------------------------
# Model (multi-head)
# -------------------------
def safe_key(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", name)

class CTPNetMultiHead(nn.Module):
    def __init__(self, n_genes: int, protein_names):
        super().__init__()
        self.protein_names = list(protein_names)
        self.protein_keys  = {p: safe_key(p) for p in self.protein_names}

        self.fc1 = nn.Linear(n_genes, 1000)
        self.fc2 = nn.Linear(1000, 256)
        self.fc3 = nn.ModuleDict({self.protein_keys[p]: nn.Linear(256, 64) for p in self.protein_names})
        self.fc4 = nn.ModuleDict({self.protein_keys[p]: nn.Linear(64, 1)  for p in self.protein_names})

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        outs = []
        for p in self.protein_names:
            k = self.protein_keys[p]
            h = F.relu(self.fc3[k](x))
            outs.append(self.fc4[k](h))
        return torch.cat(outs, dim=1)

model = CTPNetMultiHead(n_genes=X_tr.shape[1], protein_names=protein_names).to(DEVICE)
optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-3, amsgrad=True)

# -------------------------
# Loaders
# -------------------------
train_loader = DataLoader(TensorDataset(X_tr, Y_tr), batch_size=64, shuffle=True, drop_last=False)
val_loader   = DataLoader(TensorDataset(X_va, Y_va), batch_size=64, shuffle=False, drop_last=False)

# -------------------------
# Train + early stop on mean val MSE
# -------------------------
patience = 30
best_val = np.inf
bad = 0
best_state = None

def eval_mean_mse(loader):
    model.eval()
    mses = []
    with torch.no_grad():
        for xb, yb in loader:
            xb = xb.to(DEVICE)
            yb = yb.to(DEVICE)
            pred = model(xb)
            mses.append(F.mse_loss(pred, yb).item())
    return float(np.mean(mses)) if len(mses) else np.inf

max_epochs = 200
for ep in range(max_epochs):
    model.train()
    tr_losses = []
    for xb, yb in train_loader:
        xb = xb.to(DEVICE)
        yb = yb.to(DEVICE)

        optimizer.zero_grad()
        pred = model(xb)
        loss = F.mse_loss(pred, yb)
        loss.backward()
        optimizer.step()
        tr_losses.append(loss.item())

    val_mse = eval_mean_mse(val_loader)
    if ep % 10 == 0:
        print(f"epoch {ep:3d} | train_mse={np.mean(tr_losses):.4f} | val_mse={val_mse:.4f}")

    if val_mse < best_val - 1e-4:
        best_val = val_mse
        bad = 0
        best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    else:
        bad += 1
        if bad >= patience:
            print(f"Early stop @ epoch {ep}, best_val={best_val:.4f}")
            break

if best_state is not None:
    model.load_state_dict(best_state)

# -------------------------
# Predict full train + test
# -------------------------
def predict_full(X_t: torch.Tensor, index, batch_size=256):
    loader = DataLoader(TensorDataset(X_t, torch.zeros(len(X_t), 1)), batch_size=batch_size, shuffle=False)
    preds = []
    model.eval()
    with torch.no_grad():
        for xb, _ in loader:
            xb = xb.to(DEVICE)
            pred = model(xb).cpu().numpy()
            preds.append(pred)
    pred = np.vstack(preds)  # (N,P)
    return pd.DataFrame(pred, index=index, columns=protein_names)

pred_train_df = predict_full(X_train_full, rna_train.index.astype(str))
pred_test_df  = predict_full(X_test_full,  rna_test.index.astype(str))

# -------------------------
# Save as proteins × cells (R expects proteins as rownames)
# -------------------------
train_out = out_dir / f"ctpnet_pred_train_source{args.source_tag}.csv"
test_out  = out_dir / f"ctpnet_pred_test_target{args.target_tag}.csv"
ckpt_out  = out_dir / f"ctpnet_source{args.source_tag}.pt"

(pred_train_df.T).to_csv(train_out)
(pred_test_df.T ).to_csv(test_out)

ckpt = {
    "state_dict": model.state_dict(),
    "protein_names": protein_names,
    "gene_names": gene_names,   # <- panel-aligned genes (NOT genes_common)
    "adt_norm": "CLR: log((y+1)/gmean(y+1)) per-cell (TRAIN ONLY)",
    "train_dir": str(train_dir),
    "test_dir": str(test_dir),
}
torch.save(ckpt, ckpt_out)

print("Saved:", train_out, test_out, ckpt_out)








# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_trasnfer_1/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_trasnfer_3/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_A_to_C \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag A --target_tag C
# =============================================================================

# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_trasnfer_3/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_trasnfer_1/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_C_to_A \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag C --target_tag A
# =============================================================================



# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_trasnfer_1/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_trasnfer_2/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_A_to_B \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag A --target_tag B
# =============================================================================


# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_trasnfer_2/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_trasnfer_1/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_B_to_A \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag B --target_tag A
# =============================================================================




# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_trasnfer_2/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_trasnfer_3/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_B_to_C \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag B --target_tag C
# =============================================================================



# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_trasnfer_3/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_trasnfer_2/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_C_to_B \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag C --target_tag B
# =============================================================================

































# =============================================================================
#  python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_1/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_3/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_A_to_C \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag A --target_tag C
# =============================================================================


# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_2/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_1/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_B_to_A \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag B --target_tag A
# =============================================================================


# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_2/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_3/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_B_to_C \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag B --target_tag C
# =============================================================================


# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_3/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_1/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_C_to_A \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag C --target_tag A
# =============================================================================

# =============================================================================
# python ctpnet_transfer.py \
#   --train_dir ~/Desktop/FMLE/benchmarks_3/ctp \
#   --test_dir  ~/Desktop/FMLE/benchmarks_2/ctp \
#   --out_dir   ~/Desktop/FMLE/transfer_preds/ctp_C_to_B \
#   --panel     ~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv \
#   --source_tag C --target_tag B
# =============================================================================












