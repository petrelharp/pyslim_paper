import os, sys
import argparse
import tskit
import numpy as np
import pandas as pd

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

def node_mapping(ts1, ts2):
    ts1_slim_ids = [x.metadata['slim_id'] for x in ts1.nodes()]
    ts1_slim_ids_set = set(ts1_slim_ids)
    ts2_slim_ids = np.array([x.metadata['slim_id'] for x in ts2.nodes()])
    node_map = np.full(ts2.num_nodes, tskit.NULL)
    for i, node in enumerate(ts2.nodes()):
        sid = ts2_slim_ids[i]
        # check for slim_id match
        if sid in ts1_slim_ids_set:
            ts1_match_idx = np.where(ts1_slim_ids == sid)[0][0]
            # check for population id match
            if ts1.nodes_population[ts1_match_idx] == ts2.nodes_population[i]:
                node_map[i] = ts1_match_idx
    return node_map

def relabel_slim_ids(ts, id_dict):
    slim_ids = [x.metadata['slim_id'] for x in ts.nodes()]
    if len(set(slim_ids) & set(id_dict.keys())) > 1:
        nodes_to_update = [x for x, id in enumerate(slim_ids) if id in id_dict and id_dict[id]['pop'] == ts.node(x).population]
        tables = ts.dump_tables()
        tables.nodes.clear()
        for node in ts.nodes():
            if node.id in nodes_to_update:
                new_metadata = node.metadata
                new_metadata['slim_id'] = id_dict[slim_ids[node.id]]['new_slim_id']
                tables.nodes.append(node.replace(metadata=new_metadata))
            else:
                tables.nodes.append(node)
        new_ts = tables.tree_sequence()
    else:
        new_ts = ts
    return new_ts

def main(args):
    argp = argparse.ArgumentParser(description = 'Merges a set of tree sequences that are related to each other')
    argp.add_argument('-i', '--input', type = str, nargs='+', help = 'All tree sequences to be merged', required=True)
    argp.add_argument('-o', '--output', type = str, help = 'Path to desired output file', required=True)
    args = argp.parse_args(args)

    # load tree sequences
    tslist = [tskit.load(x) for x in args.input]
    # shift times for all trees to match up
    ticks = [ts.metadata['SLiM']['tick'] for ts in tslist]
    dt = [max(ticks) - t for t in ticks]
    tslist = [shift_times(ts, t) for ts, t in zip(tslist, dt)]

    # relabel slim ids, if necessary
    working_dir=os.getcwd() + "/"
    slim_update_file = working_dir + "slim_id_update.csv"
    if os.path.isfile(slim_update_file) and os.path.getsize(slim_update_file) > 0:
        slim_update_df = pd.read_csv(slim_update_file, header=None)
        slim_update_df.columns = ['old_slim_id', 'new_slim_id', 'pop']
        slim_update_dict = slim_update_df.set_index(slim_update_df.columns[0]).to_dict(orient='index')
        tslist = [relabel_slim_ids(x, slim_update_dict) for x in tslist]

    # merge trees pairwise
    tsu = tslist[0]
    for k in range(1, len(tslist)):
        ts = tslist[k]
        tsu = merge_ts(tsu, ts)

    tsu.dump(args.output)

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
