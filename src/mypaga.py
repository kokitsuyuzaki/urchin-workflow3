# paga.py
# -*- coding: utf-8 -*-
import sys
import scanpy as sc
import pandas as pd

args = sys.argv
infile = args[1]
outfile = args[2]

# Load
adata = sc.read_h5ad(infile)

# Neighbors（PAGA に必須）
sc.pp.neighbors(adata)

# Cluster ベースの PAGA
sc.tl.paga(adata, groups="celltype")

# 予約語 '_index' を消す
if "_index" in adata.var.columns:
    adata.var = adata.var.loc[:, [c for c in adata.var.columns if c != "_index"]]
if "_index" in adata.obs.columns:
    adata.obs = adata.obs.loc[:, [c for c in adata.obs.columns if c != "_index"]]

# raw ごと消す
adata.raw = None

adata.write(outfile)