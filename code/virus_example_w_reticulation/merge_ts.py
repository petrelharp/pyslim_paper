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

def merge_pop_tables(tables1, tables2):
    pop_table1 = tables1.populations.copy()
    tables1.populations.clear()
    for i in range(pop_table1.num_rows):
        pop_data = pop_table1[i]
        if pop_data.metadata is None:
            pop_data = tables2.populations[i]
        else:
            assert tables2.populations.num_rows < i + 1 or tables2.populations[i].metadata is None or tables2.populations[i] == pop_data
        tables1.populations.append(pop_data)
    if tables2.populations.num_rows > pop_table1.num_rows:
        for j in range(i+1, tables2.populations.num_rows):
        	tables1.populations.append(tables2.populations[j])
    return tables1.populations

def merge_ts(ts1, ts2):
    node_map = node_mapping(ts1, ts2)
    combined_pop_table = merge_pop_tables(ts1.tables, ts2.tables)
    merged = ts1.dump_tables()
    merged.populations.replace_with(combined_pop_table)
    ts2_tables = ts2.dump_tables()
    ts2_tables.populations.replace_with(combined_pop_table)
    merged.union(ts2_tables, node_map, add_populations = False)
    # merge metadata
    md = merged.metadata
    ts1_founders = ts1.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    ts2_founders = ts2.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    md['SLiM']['user_metadata']['FOUNDERS'][0].update(ts2_founders)
    all_infs = set(ts1_founders) | set(ts2_founders)
    # combine founders from both tree sequences
    for inf in all_infs:
        if (inf in ts1_founders) & (inf in ts2_founders):
            md['SLiM']['user_metadata']['FOUNDERS'][0][inf] = list(set(ts1_founders[inf] + ts2_founders[inf]))
    # merge slim id update dictionaries, if necessary
    slim_id_dict = None
    if 'slim_id_update_dict' in md['SLiM']['user_metadata'] and 'slim_id_update_dict' in ts2_tables.metadata['SLiM']['user_metadata']:
        ts1_dict = md['SLiM']['user_metadata']['slim_id_update_dict']
        ts2_dict = ts2_tables.metadata['SLiM']['user_metadata']['slim_id_update_dict']
        if isinstance(ts1_dict, list):
            ts1_dict = delist_dict(ts1_dict[0])
        if isinstance(ts2_dict, list):
            ts2_dict = delist_dict(ts2_dict[0])
        slim_id_dict = merge_update_dicts(ts1_dict, ts2_dict)
    elif 'slim_id_update_dict' in md['SLiM']['user_metadata'] and 'slim_id_update_dict' not in ts2_tables.metadata['SLiM']['user_metadata']:
        slim_id_dict = md['SLiM']['user_metadata']['slim_id_update_dict']
    elif 'slim_id_update_dict' not in md['SLiM']['user_metadata'] and 'slim_id_update_dict' in ts2_tables.metadata['SLiM']['user_metadata']:
        slim_id_dict = ts2_tables.metadata['SLiM']['user_metadata']['slim_id_update_dict']
    if slim_id_dict is not None:
        if isinstance(slim_id_dict, list):
            slim_id_dict = delist_dict(slim_id_dict[0])
        md['SLiM']['user_metadata']['slim_id_update_dict'] = slim_id_dict
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
            ts1_match_idx = np.where(ts1_slim_ids == sid)[0]
            # check for population id match
            for j in ts1_match_idx:
                if ts1.nodes_population[j] == ts2.nodes_population[i]:
                    node_map[i] = j
    return node_map

def relabel_slim_ids(ts, id_dict, slim_id_offset):
    slim_ids = [x.metadata['slim_id'] for x in ts.nodes()]
    nodes_to_update = set([x for x, id in enumerate(slim_ids) if str(id) in id_dict and str(ts.node(x).population) in id_dict[str(id)]])
    inds_to_update_dict = {}
    if len(nodes_to_update) > 0:
        tables = ts.dump_tables()
        tables.nodes.clear()
        for node in ts.nodes():
            if node.id in nodes_to_update:
                old_slim_id = node.metadata['slim_id']
                new_metadata = node.metadata
                new_metadata['slim_id'] = id_dict[str(slim_ids[node.id])][str(node.population)]
                tables.nodes.append(node.replace(metadata=new_metadata))
                if node.individual != tskit.NULL:
                    inds_to_update_dict[node.individual] = old_slim_id
            else:
                tables.nodes.append(node)
        # update individual pedigree ids, if necessary
        if len(inds_to_update_dict) > 0:
            tables.individuals.clear()
            ped_relabel_dict = {}
            for ind in ts.individuals():
                if ind.id in inds_to_update_dict:
                    old_ped_id = ind.metadata['pedigree_id']
                    new_ped_id = int(inds_to_update_dict[ind.id] + ind.population * slim_id_offset)//2
                    new_metadata = ind.metadata
                    new_metadata['pedigree_id'] = new_ped_id
                    tables.individuals.append(ind.replace(metadata=new_metadata))
                    ped_relabel_dict[old_ped_id] = new_ped_id
                else:
                    tables.individuals.append(ind)
            # update founder lists in metadata
            md = tables.metadata
            for host in md['SLiM']['user_metadata']['FOUNDERS'][0]:
                founder_inds = md['SLiM']['user_metadata']['FOUNDERS'][0][host]
                if len(set(founder_inds) & set(ped_relabel_dict)) > 0:
                    founder_inds_new = [x if x not in ped_relabel_dict else ped_relabel_dict[x] for x in founder_inds]
                    md['SLiM']['user_metadata']['FOUNDERS'][0][host] = founder_inds_new
            tables.metadata = md
        ts = tables.tree_sequence()
    return ts

# merge different slim id update dictionaries
def merge_update_dicts(original_dict, updates):
    for old_slim_id, inner_dict in updates.items():
        if old_slim_id not in original_dict:
            original_dict[old_slim_id] = inner_dict
        else:
            for inner_key, value in inner_dict.items():
                if inner_key not in original_dict[old_slim_id]:
                    original_dict[old_slim_id][inner_key] = value
                else:
                    if original_dict[old_slim_id][inner_key] != value:
                        raise ValueError("Conflict at [{old_slim_id}][{inner_key}]: {original_dict[old_slim_id][inner_key]} != {value}")
    return original_dict

def delist_dict(dict):
    return {outer_k: {inner_k: inner_v[0] for d in v for inner_k, inner_v in d.items()} for outer_k, v in dict.items()}

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
    slim_id_update_dict = {}
    for ts in tslist:
        if 'slim_id_update_dict' in ts.metadata['SLiM']['user_metadata']:
            new_dict = ts.metadata['SLiM']['user_metadata']['slim_id_update_dict']
            if isinstance(new_dict, list):
                new_dict = delist_dict(new_dict[0])
            if len(slim_id_update_dict) == 0:
                slim_id_update_dict = new_dict
            else:
                slim_id_update_dict = merge_update_dicts(slim_id_update_dict, new_dict)
    # calculate slim_id_offset
    if len(slim_id_update_dict) > 0:
        old_slim_id =  next(iter(slim_id_update_dict))
        pop = next(iter(slim_id_update_dict[old_slim_id]))
        new_slim_id = slim_id_update_dict[old_slim_id][pop]
        slim_id_offset = int((new_slim_id - int(old_slim_id)) / int(pop))
        tslist = [relabel_slim_ids(x, slim_id_update_dict, slim_id_offset) for x in tslist]

    # merge trees pairwise
    tsu = tslist[0]
    for k in range(1, len(tslist)):
        ts = tslist[k]
        tsu = merge_ts(tsu, ts)

    tsu.dump(args.output)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))