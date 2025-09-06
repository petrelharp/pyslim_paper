import numpy as np
import pandas as pd
import scipy.sparse

import tskit, pyslim


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
    merged = merge_pop_tables(ts1, ts2)
    tables2 = ts2.dump_tables()
    # metadata may differ for founder individuals because they lived for longer
    update_founder_metadata(merged, tables2)
    update_founder_metadata(tables2, merged)
    try:
        merged.union(tables2, node_map, add_populations=False)
    except Exception as err:
        # If this fails because of non-equal overlap we really
        # want to know *what* differs between them; this
        # replicates what happens under the hood in union().
        if "TSK_ERR_UNION_DIFF_HISTORIES" in str(err):
            n2 = np.where(node_map >= 0)[0]
            n1 = node_map[n2]
            sts1 = merged.subset(n1)
            sts2 = tables2.subset(n2)
            sts1.canonicalise()
            sts2.canonicalise()
            sts1.assert_equals(sts2, ignore_provenance=True, ignore_ts_metadata=True)
                               
    md = merged.metadata
    md['SLiM']['user_metadata']['FOUNDERS'][0].update(
            ts2.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    )
    md['SLiM']['user_metadata']['FOUNDING_TIME'][0].update(
            ts2.metadata['SLiM']['user_metadata']['FOUNDING_TIME'][0]
    )
    merged.metadata = md
    return merged.tree_sequence()


def slim_to_nodes(ts):
    return { i.metadata['pedigree_id'] : i.nodes for i in ts.individuals() }

def slim_to_inds(ts):
    return { i.metadata['pedigree_id'] : i.id for i in ts.individuals() }


def merge_pop_tables(ts1, ts2):
    """
    Returns the tables corresponding to ts1, but with any additional
    populations described in ts2 (ie with non-Null metadata) included.
    Note: it would be more natural to return a merged population table here,
    but there's currently not a good way to get a new table with byte-for-byte
    identical population metadata schema, which messes up the shared overlap
    checking.
    """
    tables = ts1.dump_tables()
    tables.populations.clear()
    k = 0
    while k < ts1.num_populations:
        p1 = ts1.population(k)
        md = p1.metadata
        if k < ts2.num_populations:
            p2 = ts2.population(k)
            if p2.metadata is not None:
                if md is not None:
                    assert md == p2.metadata
                else:
                    md = p2.metadata
        tables.populations.append(p1.replace(metadata=md))
        k += 1
    while k < ts2.num_populations:
        p2 = ts2.population(k)
        tables.populations.append(p2)
        k += 1
    return tables


def update_founder_metadata(tables1, tables2):
    """
    Update *in place* the individual metadata in tables1 for founder individuals
    with the metadata from the same individuals in tables2, using whichever
    metadata value has a larger `age`.  Use case: they were remembered in
    tables1 in their youth, but then tables2 recorded them later in life.
    """
    map2 = { i.metadata['pedigree_id'] : k for k, i in enumerate(tables2.individuals) }
    f1 = tables1.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    f2 = tables2.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    slim_ids = [i for k in f1 if k in f2 for i in f2[k]]
    individuals = tables1.individuals.copy()
    tables1.individuals.clear()
    for ind in individuals:
        md = ind.metadata
        sid = md['pedigree_id']
        if sid in slim_ids:
            other = tables2.individuals[map2[sid]].metadata
            if other['age'] > md['age']:
                md = other
        tables1.individuals.append(ind.replace(metadata=md))


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
