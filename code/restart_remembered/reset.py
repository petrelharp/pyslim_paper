import tskit, pyslim
ts = tskit.load("init.trees")

def shift_times(ts, t):
    tables = ts.dump_tables()
    tables.nodes.clear()
    tables.mutations.clear()
    for node in ts.nodes():
        tables.nodes.append(node.replace(time=node.time + t))
    for mutation in ts.mutations():
        tables.mutations.append(mutation.replace(time=mutation.time + t))
    return tables

def reset_time(ts, t):
    tables = shift_times(ts, -t)
    tables.individuals.clear()
    for ind in ts.individuals():
        if ind.time == t:
            new_flags = ind.flags | pyslim.INDIVIDUAL_ALIVE
        else:
            new_flags = ind.flags & ~pyslim.INDIVIDUAL_ALIVE
        tables.individuals.append(ind.replace(flags=new_flags))
    md = tables.metadata
    md['SLiM']['tick'] = t
    tables.metadata = md
    return tables.tree_sequence()

rts = reset_time(ts, 5)
rts.dump("reset.trees")
