import tskit, pyslim
from util import *

ots = tskit.load("sim_00.trees")
sts = reset_time(ots, 10)

for ts in (ots, sts):

    fi = ts.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    for x, st in zip(ts.metadata['SLiM']['user_metadata']['OUTPUT_HOST_IDS'],
                    ts.metadata['SLiM']['user_metadata']['TRANSMISSION_DAYS']):
        t = ts.metadata['SLiM']['tick'] - st
        inds = pyslim.individuals_alive_at(ts, t)
        sids = [ts.individual(i).metadata['pedigree_id'] for i in inds]
        print(x, t, st)
        print("   ", fi[x])
        print("   ", [u in sids for u in fi[x]])
        for ind in ts.individuals():
            assert (ind.id in sids) == (ind.id in fi[x])
