import tskit
import sys

from util import *

infile = sys.argv[1]
host_id = sys.argv[2]
outfile = sys.argv[3]

ts = tskit.load(infile)

new_time = ts.metadata['SLiM']['user_metadata']['FOUNDING_TIME'][0][host_id][0]

rts = reset_time(ts, new_time)
rts.dump(outfile)
