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

# Heatmap
top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:100]
scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color="clusters", n_convolve=100, save=outfile, figsize=(10,8))
