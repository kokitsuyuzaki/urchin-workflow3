# -*- coding: utf-8 -*-
import sys
import scvelo as scv
import scanpy
import pandas as pd

args = sys.argv
infile = args[1]
outfile = args[2]

# Parameters
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# Loading
adata = scv.read(infile)

# Embedding
scv.pl.velocity_embedding_stream(
    adata, 
    basis='umap', 
    save=outfile, 
    dpi=500, 
    figsize=(10,10),
    color='clusters',       # ← celltypeで色分け
    linewidth=2.5,
    size=80,
    alpha=0.8,
    density=1.5
)