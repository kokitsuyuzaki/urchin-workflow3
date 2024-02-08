# -*- coding: utf-8 -*-
import sys
import shutil
import loompy

args = sys.argv

outfile = args[1]
infiles = ["output/hpbase/" + i + "/velocyto/" + i + ".loom" for i in ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']]

# cp
infiles2 = ["output/hpbase/" + i + "/velocyto/" + i + ".tmp" for i in ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']]
for x in range(len(infiles)):
    shutil.copyfile(infiles[x], infiles2[x])

# Merge
loompy.combine(infiles2, outfile, key="Accession")
