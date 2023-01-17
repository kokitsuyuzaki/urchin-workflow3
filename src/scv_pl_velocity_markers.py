# -*- coding: utf-8 -*-
import sys
import scvelo as scv
import scanpy
import pandas as pd

args = sys.argv
infile = args[1]
outfile1 = args[2]
outfile2 = args[3]
outfile3 = args[4]

# Parameters
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# Loading
adata = scv.read(infile)

# Markers
top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:100]
scv.pl.velocity(adata, top_genes[range(0, 27)], ncols=3, save=outfile1)
scv.pl.velocity(adata, top_genes[range(27, 54)], ncols=3, save=outfile2)
scv.pl.velocity(adata, top_genes[range(54, 81)], ncols=3, save=outfile3)
