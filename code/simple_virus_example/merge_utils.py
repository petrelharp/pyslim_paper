import numpy as np
import pandas as pd
import scipy.sparse

import tskit, pyslim


def clear_alive(ts):
    tables = ts.dump_tables()
    tables.individuals.clear()
    for ind in ts.individuals():
        tables.individuals.append(ind.replace(flags=ind.flags & ~pyslim.INDIVIDUAL_ALIVE))
    return tables.tree_sequence()


def shift_times(ts, t_to_add):
    ts_tables = ts.dump_tables()
    ts_tables.nodes.clear()
    ts_tables.mutations.clear()
    for node in ts.nodes():
        ts_tables.nodes.append(node.replace(time=node.time + t_to_add))
    for mutation in ts.mutations():
        ts_tables.mutations.append(mutation.replace(time=mutation.time + t_to_add))
    new_ts = ts_tables.tree_sequence()
    return new_ts


def merge_ts(ts1, ts2):
    node_map = node_mapping(ts1, ts2)
    merged = ts1.dump_tables()
    merged.union(ts2.tables, node_map, add_populations = True)
    md = merged.metadata
    md['SLiM']['user_metadata']['FOUNDERS'][0].update(
            ts2.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    )
    merged.metadata = md
    return merged.tree_sequence()


def slim_to_nodes(ts):
    return { i.metadata['pedigree_id'] : i.nodes for i in ts.individuals() }


def shared_founder_nodes(ts1, ts2):
    f1 = ts1.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    f2 = ts2.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    shared = set(f1.keys()).intersection(set(f2.keys()))
    map1 = slim_to_nodes(ts1)
    map2 = slim_to_nodes(ts2)
    sn1 = []
    sn2 = []
    for k in f1:
        if k in f2:
            assert f1[k] == f2[k]
            sn1.extend([u for sid in f1[k] for u in map1[sid]])
            sn2.extend([u for sid in f2[k] for u in map2[sid]])
    t1 = ts1.nodes_time[sn1]
    t2 = ts2.nodes_time[sn2]
    assert np.allclose(t1, t2)
    return sn1, sn2


def node_mapping(ts1, ts2):
    """
    Appropriate for doing ts1.union(ts2, node_mapping)
    """
    sn1, sn2 = shared_founder_nodes(ts1, ts2)
    anc1 = ancestors(ts1, sn1)
    anc2 = ancestors(ts2, sn2)
    assert len(anc1) == len(anc2)
    for a1, a2 in zip(anc1, anc2):
        n1 = ts1.node(a1)
        n2 = ts2.node(a2)
        assert n1.metadata['slim_id'] == n2.metadata['slim_id']
        assert n1.time == n2.time
    node_map = np.full(ts2.num_nodes, tskit.NULL)
    node_map[anc2] = anc1
    return node_map


def ancestors(ts, nodes):
    """
    Returns the list of nodes reachable from `nodes` by following
    child->parent relationships in the edge table.
    """
    out = np.zeros(ts.num_nodes, dtype='int')
    out[nodes] = 1
    x = out.copy()
    A = scipy.sparse.coo_array(
            (
                np.ones(ts.num_edges, dtype='int'),
                (ts.edges_parent, ts.edges_child)
            ), dtype='int', shape=(ts.num_nodes, ts.num_nodes)
    )
    while np.sum(x) > 0:
        x = A @ x
        out += x
    return np.where(out > 0)[0]
