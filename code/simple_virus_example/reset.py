import tskit, pyslim
import sys

infile = sys.argv[1]
host_id = sys.argv[2]
outfile = sys.argv[3]

ts = tskit.load(infile)

tick = ts.metadata['SLiM']['tick']
new_tick = ts.metadata['SLiM']['user_metadata']['FOUNDING_TIME'][0][host_id][0]

founders = [ind.id for ind in ts.individuals()
            if ind.metadata['subpopulation'] == int(host_id)]

rts = pyslim.set_slim_state(ts, time=tick - new_tick, individuals=founders)
rts.dump(outfile)
