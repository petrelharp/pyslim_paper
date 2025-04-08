import tskit
import numpy as np
import pandas as pd

def match_nodes(other, ts, split_time):
    """
    Given SLiM tree sequences `other` and `ts`, builds a numpy array with length
    `other.num_nodes` in which the indexes represent the node id in `other` and the
    entries represent the equivalent node id in `ts`. If a node in `other` has no
    equivalent in `ts`, then the entry takes the value `tskit.NULL` (-1). The
    matching is done by comparing the IDs assigned by SLiM which are kept in
    node metadata. This matching of SLiM IDs is *only* done for nodes with time
    older than the specified `split_time`.
    """
    node_mapping = np.full(other.num_nodes, tskit.NULL)
    sids0 = np.array([n.metadata["slim_id"] for n in ts.nodes()])
    sids1 = np.array([n.metadata["slim_id"] for n in other.nodes()])
    alive_before_split1 = (other.tables.nodes.time >= split_time)
    is_1in0 = np.isin(sids1, sids0)
    both = np.logical_and(alive_before_split1, is_1in0)
    sorted_ids0 = np.argsort(sids0)
    matches = np.searchsorted(
        sids0,
        sids1[both],
        side='left',
        sorter=sorted_ids0
    )
    node_mapping[both] = sorted_ids0[matches]
    return node_mapping

def union_children(parent, df, merged):
    print(f"Going in: {parent}")
    child_rows = df[df.parent == parent]
    assert (len(child_rows) == 2) or (len(childs) == 0)
    if len(child_rows) == 2:
        children = [row.child for _, row in child_rows.iterrows()]
        for child in children:
            if child not in merged:
                union_children(child, df, merged)
        split_time = merged[children[0]]["depth"]
        assert split_time == merged[children[1]]["depth"] # ultrametric
        print(f'Unioning: {children}, Split time: {split_time}')
        ts0 = merged[children[0]]["ts"]
        ts1 = merged[children[1]]["ts"]
        node_map = match_nodes(ts1, ts0, split_time)
        tsu = ts0.union(ts1, node_map, check_shared_equality=True)
        # the time from tip to start of simulation is split_time plus the
        # length of the edge
        parent_edgelength = df[df.child==parent].edgelen.item()
        merged[parent] = {
            "ts": tsu,
            "depth": split_time + parent_edgelength,
            "children": merged[children[0]]["children"] + merged[children[1]]["children"]
        }

merged = {
    row.child : {
        "ts": tskit.load(row.outfile),
        "depth": row.edgelen,
        "children": [row.child]
    }
    for i, row in df[df.is_leaf].iterrows()
}

union_children("root", df, merged)
# union of all three species tree sequences is in the root.
tsu = merged["root"]["ts"]
pops = merged["root"]["children"]
