import tskit, pyslim

def shift_times(ts, dt):
    """
    Add `dt` to the times ago of all nodes and mutations.
    """
    tables = ts.dump_tables()
    tables.nodes.clear()
    tables.mutations.clear()
    for node in ts.nodes():
        tables.nodes.append(node.replace(time=node.time + dt))
    for mutation in ts.mutations():
        tables.mutations.append(mutation.replace(time=mutation.time + dt))
    return tables

def reset_time(ts, t):
    """
    Change times in ts so that it matches what would have been written out
    by SLiM in tick `t`.
    """
    dt = ts.metadata['SLiM']['tick'] - t
    tables = shift_times(ts, -dt)
    tables.individuals.clear()
    alive_inds = pyslim.individuals_alive_at(ts, dt, stage="late", remembered_stage="late")
    for ind in ts.individuals():
        if ind.id in alive_inds:
            new_flags = ind.flags | pyslim.INDIVIDUAL_ALIVE
        else:
            new_flags = ind.flags & ~pyslim.INDIVIDUAL_ALIVE
        tables.individuals.append(ind.replace(flags=new_flags))
    md = tables.metadata
    md['SLiM']['tick'] = t
    md['SLiM']['cycle'] = t
    tables.metadata = md
    return tables.tree_sequence()

