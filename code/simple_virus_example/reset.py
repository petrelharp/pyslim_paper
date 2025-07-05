import tskit
import sys

from util import *

infile = sys.argv[1]
new_time = int(sys.argv[2])
outfile = sys.argv[3]

ts = tskit.load(infile)
rts = reset_time(ts, new_time)
rts.dump(outfile)
