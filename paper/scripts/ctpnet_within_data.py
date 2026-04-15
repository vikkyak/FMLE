#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 12:37:57 2026

@author: vikas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
torch.set_float32_matmul_precision("high") if hasattr(torch, "set_float32_matmul_precision") else None

# -------------------------
# Paths + Load
# -------------------------
# =============================================================================
# base = Path("~/Desktop/FMLE/benchmarks_3/ctp").expanduser()
# =============================================================================
# =============================================================================
# base = Path("~/Desktop/FMLE/benchmarks_c_doror_1/ctp").expanduser()
# =============================================================================
# =============================================================================
# base = Path("~/Desktop/FMLE//benchmarks_spatial_pbmc10_2/ctp").expanduser()
# =============================================================================
# =============================================================================
# base = Path("~/Desktop/FMLE/COVID_PBMC_CITE_seq/ctp").expanduser()
# =============================================================================
# =============================================================================
# base = Path("~/Desktop/FMLE/Bone_marrow_Ab-seq/ctp").expanduser()
# =============================================================================
base = Path("/your/root/path") / "ctp"
base.mkdir(parents=True, exist_ok=True)

rna_train = pd.read_csv(base / "rna_train.csv", index_col=0)
adt_train = pd.read_csv(base / "adt_train.csv", index_col=0)
rna_test  = pd.read_csv(base / "rna_test.csv",  index_col=0)
adt_test  = pd.read_csv(base / "adt_test.csv",  index_col=0)

# =============================================================================
# for spatial trasnfer
# =============================================================================
# =============================================================================
# rna_test  = pd.read_csv(base / "rna_spatial_test.csv",  index_col=0)
# adt_test  = pd.read_csv(base / "adt_spatial_test.csv",  index_col=0)
# =============================================================================
assert rna_test.shape[1] == 2000

# sanity: same cell order between RNA and ADT
assert (rna_train.index == adt_train.index).all()
assert (rna_test.index  == adt_test.index).all()

protein_names = list(adt_train.columns)
gene_names    = list(rna_train.columns)

# -------------------------
# cTPnet-style ADT transform (CLR per cell)
# -------------------------
def norm_adt_cellsxprot(Y_cellsxprot: np.ndarray) -> np.ndarray:
    Yp = Y_cellsxprot.astype(np.float64) + 1.0
    gm = gmean(Yp, axis=1, keepdims=True)  # per cell
    return np.log(Yp / gm)

Y_all = norm_adt_cellsxprot(adt_train.values)  # TRAIN ONLY for training/val split

# -------------------------
# Train/Val split (from TRAIN only)
# -------------------------
idx = np.arange(len(rna_train))
tr_idx, va_idx = train_test_split(idx, test_size=0.15, random_state=4905, shuffle=True)

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
    return re.sub(r'[^A-Za-z0-9_]', '_', name)

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
            outs.append(self.fc4[k](h))          # (B,1)
        return torch.cat(outs, dim=1)            # (B,P)

model = CTPNetMultiHead(n_genes=X_tr.shape[1], protein_names=protein_names).to(DEVICE)
optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-3, amsgrad=True)

# -------------------------
# Loaders (IMPORTANT: Y is real, not zeros)
# -------------------------
batch_size = 64
train_loader = DataLoader(TensorDataset(X_tr, Y_tr), batch_size=batch_size, shuffle=True, drop_last=False)
val_loader   = DataLoader(TensorDataset(X_va, Y_va), batch_size=batch_size, shuffle=False, drop_last=False)

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
# Save in protein × cell (R expects proteins as rownames)
# -------------------------
(pred_train_df.T).to_csv(base / "kaggle_ctpnet_pred_train.csv")
(pred_test_df.T ).to_csv(base / "kaggle_ctpnet_pred_test.csv")
(pd.concat([pred_train_df, pred_test_df], axis=0).T).to_csv(base / "kaggle_ctpnet_pred_all.csv")

# checkpoint (optional)
ckpt = {
    "state_dict": model.state_dict(),
    "protein_names": protein_names,
    "gene_names": gene_names,
    "adt_norm": "CLR: log((y+1)/gmean(y+1)) per-cell",
}
torch.save(ckpt, base / "ctpnet.pt")

print("Saved:",
      base / "kaggle_ctpnet_pred_train.csv",
      base / "kaggle_ctpnet_pred_test.csv",
      base / "kaggle_ctpnet_pred_all.csv")
