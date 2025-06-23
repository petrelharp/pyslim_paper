import tskit, pyslim
import numpy as np

def clear_alive(ts):
    tables = ts.dump_tables()
    tables.individuals.clear()
    for ind in ts.individuals():
        tables.individuals.append(ind.replace(flags=ind.flags & ~pyslim.INDIVIDUAL_ALIVE))
    return tables.tree_sequence()

ts1 = clear_alive(tskit.load("init.trees"))
ts2 = clear_alive(tskit.load("rerun.trees"))

node_map = np.full(ts2.num_nodes, tskit.NULL)
smap = {n.metadata['slim_id'] : n.id for n in ts1.nodes()}
for n in ts2.nodes():
    sid = n.metadata['slim_id']
    if sid in smap:
        node_map[n.id] = smap[sid]

tsu = ts1.union(ts2, node_map, add_populations=False)

tsu.dump("union.trees")
