# -*- coding: utf-8 -*-
import sys
import loompy

args = sys.argv

db = args[1]
outfile = args[2]
infiles = ["output/" + db + "/" + i + "/velocyto/" + i + ".loom" for i in ['cont-24h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']]

loompy.combine(infiles, outfile, key="Accession")
