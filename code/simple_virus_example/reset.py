import tskit, pyslim
import sys

infile = sys.argv[1]
host_id = sys.argv[2]
outfile = sys.argv[3]

ts = tskit.load(infile)

tick = ts.metadata['SLiM']['tick']
new_tick = ts.metadata['SLiM']['user_metadata']['FOUNDING_TIME'][0][host_id][0]

# do not resurrect founders with conflicting pedigree IDs
founders = []
num_removed = 0
seen_ids = set()
for ind in ts.individuals():
    if ind.metadata['subpopulation'] == int(host_id):
        sid = ind.metadata['pedigree_id']
        if sid not in seen_ids:
            founders.append(ind.id)
        else:
            print(f"Duplicate pedigree id {sid}; removing.")
            num_removed += 1
        seen_ids.add(sid)
if num_removed > 0:
    print(f"Removed {num_removed} founders due to duplicate IDs.")

rts = pyslim.set_slim_state(ts, time=tick - new_tick, individuals=founders)
rts.dump(outfile)
