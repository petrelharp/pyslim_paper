import sys
import argparse
import tskit
import pyslim
import numpy as np
from collections import Counter

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
            ts1_match_idx = np.where(ts1_slim_ids == sid)[0][0]
            # check for population id match
            if ts1.nodes_population[ts1_match_idx] == ts2.nodes_population[i]:
                node_map[i] = ts1_match_idx
    return node_map

def shift_ped_ids(ts, pedigree_ids=None, ped_offset_value=None):
    # shift pedigree ids in individual table and top-level metadata, if applicable
    new_ts = update_individual_metadata(ts = ts, ped_offset_value = ped_offset_value,
                                        ids_to_update = set([x.id for x in ts.individuals() if x.metadata['pedigree_id'] in pedigree_ids]))
    # update founder metadata for any hosts, as required
    founder_dict = new_ts.metadata['SLiM']['user_metadata']['FOUNDERS'][0]
    tables = new_ts.dump_tables()
    md = tables.metadata
    for host in founder_dict:
        founder_ped_ids = founder_dict[host]
        if len(set(pedigree_ids) & set(founder_ped_ids)) > 0:
            founder_ped_ids = [int(x + [y.population for y in ts.individuals() if y.metadata['pedigree_id'] == x][0] * ped_offset_value)
                               if x in pedigree_ids else int(x)
                               for x in founder_ped_ids]
            md['SLiM']['user_metadata']['FOUNDERS'][0][host] = founder_ped_ids
    tables.metadata = md
    new_ts = tables.tree_sequence()
    return new_ts

# adjust pedigree id based on the offset value and population id
def update_individual_metadata(ts, ids_to_update, ped_offset_value):
    tables = ts.dump_tables()
    tables.individuals.clear()
    for ind in ts.individuals():
        if ind.id in ids_to_update:
            new_ped_id = int(ind.metadata['pedigree_id'] + ind.population * ped_offset_value)
            new_metadata = ind.metadata
            new_metadata['pedigree_id'] = new_ped_id               
            tables.individuals.append(ind.replace(metadata=new_metadata))
        else:
            tables.individuals.append(ind)
    new_ts = tables.tree_sequence()
    return new_ts

# ensure slim ids match pedigree ids
def fix_slim_ids(ts):
    tables = ts.dump_tables()
    slim_id_dict = {}
    tables.nodes.clear()
    for node in ts.nodes():
        new_metadata = node.metadata
        if node.individual != tskit.NULL:
            ind = ts.individual(node.individual)
            ped_id = ind.metadata['pedigree_id']
            old_slim_id = node.metadata['slim_id']
            if old_slim_id//2 != ped_id:
                if old_slim_id % 2 == 0:
                    new_slim_id = int(ped_id * 2)
                elif old_slim_id % 2 == 1:
                    new_slim_id = int(ped_id * 2 + 1)
                if old_slim_id in slim_id_dict:
                    assert node.population not in slim_id_dict[old_slim_id]
                slim_id_dict.update({str(old_slim_id): {str(node.population): new_slim_id}})
                new_metadata['slim_id'] = new_slim_id
        tables.nodes.append(node.replace(metadata=new_metadata))
    # add slim id update dictionary to top-level metadata, merging entries as required
    md = tables.metadata
    if 'slim_id_update_dict' in md['SLiM']['user_metadata']:
        old_dict = md['SLiM']['user_metadata']['slim_id_update_dict']
        slim_id_dict = merge_update_dicts(old_dict, slim_id_dict)
    md['SLiM']['user_metadata'].update({'slim_id_update_dict': slim_id_dict})
    tables.metadata = md
    new_ts = tables.tree_sequence()
    return new_ts

# shift slim ids if there are duplicates
def shift_dup_slim_ids(ts, slim_offset_value, host_id):
    slim_ids = [x.metadata['slim_id'] for x in ts.nodes()]
    if len(slim_ids) > len(set(slim_ids)):
        tables = ts.dump_tables()
        counts = Counter(slim_ids)
        duplicates = [x for x, count in counts.items() if count > 1]
        nodes_to_update = []
        for dup in duplicates:
            nodes_to_update = nodes_to_update + [int(x) for x in np.where(np.array(slim_ids) == dup)[0]][1:]
        nodes_to_update = set(nodes_to_update)
        inds_to_update = []
        slim_id_dict = {}
        tables.nodes.clear()
        for node in ts.nodes():
            if node.id in nodes_to_update:
                if node.individual is not tskit.NULL:
                    inds_to_update = inds_to_update + [node.individual]
                    assert node.id - 1 in nodes_to_update or node.id + 1 in nodes_to_update
                new_metadata = node.metadata
                old_slim_id = new_metadata['slim_id']
                new_slim_id = int(old_slim_id + node.population * slim_offset_value)
                if old_slim_id in slim_id_dict:
                    assert str(node.population) not in slim_id_dict[str(old_slim_id)]
                slim_id_dict.update({str(old_slim_id): {str(node.population): new_slim_id}})
                new_metadata['slim_id'] = new_slim_id
                tables.nodes.append(node.replace(metadata=new_metadata))
            else:
                tables.nodes.append(node)
        # Adjust individual ids, if necessary
        inds_to_update = set(inds_to_update)
        founders_to_update = {}
        if len(inds_to_update) > 0:
            ped_offset_value = slim_offset_value/2
            all_founders = set(ts.metadata['SLiM']['user_metadata']['FOUNDERS'][0][host_id])
            tables.individuals.clear()
            for ind in ts.individuals():
                if ind.id in inds_to_update:
                    new_metadata = ind.metadata
                    old_ped_id = new_metadata['pedigree_id']
                    new_ped_id = int(old_ped_id + ind.population * ped_offset_value)
                    new_metadata['pedigree_id'] = new_ped_id
                    tables.individuals.append(ind.replace(metadata=new_metadata))
                    if (old_ped_id in all_founders) and (ind.flags & pyslim.INDIVIDUAL_ALIVE):
                        founders_to_update.update({old_ped_id: ind.population})
                else:
                    tables.individuals.append(ind)
        # add slim id update dictionary to top-level metadata, merging entries as required
        md = tables.metadata
        if len(founders_to_update) > 0:
            md['SLiM']['user_metadata']['FOUNDERS'][0][host_id] = [x if x not in founders_to_update
                                                                   else int(x + founders_to_update[x] * ped_offset_value)
                                                                   for x in all_founders]
        if 'slim_id_update_dict' in md['SLiM']['user_metadata']:
            old_dict = md['SLiM']['user_metadata']['slim_id_update_dict']
            slim_id_dict = merge_update_dicts(old_dict, slim_id_dict)
        md['SLiM']['user_metadata']['slim_id_update_dict'] = slim_id_dict
        tables.metadata = md
        ts = tables.tree_sequence()
        # check again for discrepancies just in case only one node was shifted for a particular individual
        ts = fix_slim_ids(ts)
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
    argp = argparse.ArgumentParser(description = 'Merge multiple tree sequences together to generate input for a coinfected host')
    argp.add_argument('-i', '--input', type = str, nargs='+', help = 'All tree sequences to be merged', required=True)
    argp.add_argument('-H', '--host_id', type = str, help = 'ID of host to be transmitted to', required=True)
    argp.add_argument('--slim_id_offset', type = int, default=1e10, help = 'How much to add to slim_id values if duplicates arise')
    argp.add_argument('-o', '--output', type = str, help = 'Desired output file name', required=True)
    args = argp.parse_args(args)

    # load tree sequences
    tslist = [tskit.load(x) for x in args.input]
    # check that all ticks are identical
    ticks = [ts.metadata['SLiM']['tick'] for ts in tslist]
    if len(set(ticks)) > 1:
        raise RuntimeError("Infections were transmitted at different times")

    # merge trees pairwise
    tsu = tslist[0]
    founders0 = set(tsu.metadata['SLiM']['user_metadata']['FOUNDERS'][0][args.host_id])
    for ts in tslist[1:]:
        # check for duplicate founder individuals across tree sequences and resolve
        new_founders = set(ts.metadata['SLiM']['user_metadata']['FOUNDERS'][0][args.host_id])
        dup_ped_ids = founders0 & new_founders
        if len(dup_ped_ids) > 0:
            # shift pedigree and slim ids accordingly
            ts = shift_ped_ids(ts, pedigree_ids=dup_ped_ids, ped_offset_value=args.slim_id_offset/2)
        # merge
        tsu = merge_ts(tsu, ts)

    # update slim ids to conform with shifted pedigree ids, if necessary
    if len(dup_ped_ids) > 0:
        tsu = fix_slim_ids(tsu)

    # check for duplicate slim_ids and shift if necessary (also shifts pedigree ids for associated individuals)
    tsu = shift_dup_slim_ids(tsu, args.slim_id_offset, args.host_id)

    # write out result
    tsu.dump(args.output)

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))