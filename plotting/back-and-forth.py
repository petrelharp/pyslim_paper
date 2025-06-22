#!/usr/bin/env python

import sys
import tskit, pyslim, tszip
# import xml.dom.minidom # to fix tskit#3144
from plotting import plot_tree

assert len(sys.argv) == 3 or len(sys.argv) == 2, f"Usage: {sys.argv[0]} (trees file) [focal position]"

fname = sys.argv[1]

try:
    ts = tskit.load(fname)
except tskit.FileFormatError:
    ts = tszip.decompress(fname)

if len(sys.argv) < 3:
    focal_pos = ts.site(0).position
else:
    focal_pos = int(sys.argv[2])

outfile = ".".join(fname.split(".")[:-1] + ["muts", "svg"])
plot_tree(ts, focal_pos, outfile)

outfile = ".".join(fname.split(".")[:-1] + ["simp", "muts", "svg"])
plot_tree(ts.simplify(ts.samples(time=0), filter_sites=False), focal_pos, outfile)
