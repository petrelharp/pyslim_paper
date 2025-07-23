import sys
import os
import argparse
import tskit
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
    merged = ts1.dump_tables()
    merged.populations.replace_with(merge_pop_tables(ts1.tables, ts2.tables))
    ts2_tables = ts2.dump_tables()
    ts2_tables.populations.replace_with(merge_pop_tables(ts1.tables, ts2.tables))
    merged.union(ts2_tables, node_map, add_populations = False)
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

def main(args):
    argp = argparse.ArgumentParser(description = 'Merge multiple tree sequences together to generate input for a coinfected host')
    argp.add_argument('-i', '--input', type = str, nargs='+', help = 'All tree sequences to be merged', required=True)
    argp.add_argument('-H', '--host_id', type = str, help = 'ID of host to be transmitted to', required=True)
    argp.add_argument('-o', '--output', type = str, help = 'Desired output file name', required=True)
    args = argp.parse_args(args)

    # load tree sequences
    tslist = [tskit.load(x) for x in args.input]
    # check that all ticks are identical
    ticks = [ts.metadata['SLiM']['tick'] for ts in tslist]
    if len(set(ticks)) > 1:
         stop("Infections were transmitted at different times")

    # merge trees pairwise
    tsu = tslist[0]
    founders_present = tsu.metadata['SLiM']['user_metadata']['FOUNDERS'][0][args.host_id]
    for k in range(1, len(tslist)):
        ts = tslist[k]
        founders_to_add = ts.metadata['SLiM']['user_metadata']['FOUNDERS'][0][args.host_id]
        # make sure none of the founding individuals from different sources share a pedigree id
        common_ids = set(founders_present) & set(founders_to_add)
        if len(common_ids) > 1:
            tables = ts.dump_tables()
            tables.individuals.clear()
            # if they do, kill that individual in ts and remove from founder list
            inds_to_kill = set([i.id for i in ts.individuals() if i.metadata['pedigree_id'] in common_ids])
            for ind in ts.individuals():
                if ind.id in inds_to_kill:
                    new_flags = ind.flags & ~pyslim.INDIVIDUAL_ALIVE
                else:
                    new_flags = ind.flags
                tables.individuals.append(ind.replace(flags=new_flags))
            md = tables.metadata
            new_founder_list = [x for x in founders_to_add if x not in common_ids]
            md['SLiM']['user_metadata']['FOUNDERS'][0][args.host_id] = new_founder_list
            tables.metadata = md
            ts = tables.tree_sequence()
        # merge trees
        tsu = merge_ts(tsu, ts)
    # check for duplicate slim_ids and update if necessary
    slim_ids = [x.metadata['slim_id'] for x in tsu.nodes()]
    if len(slim_ids) > len(set(slim_ids)):
        working_dir=os.getcwd() + "/"
        tables = tsu.dump_tables()
        counts = Counter(slim_ids)
        duplicates = [x for x, count in counts.items() if count > 1]
        nodes_to_update = []
        for dup in duplicates:
            nodes_to_update = nodes_to_update + [int(x) for x in np.where(np.array(slim_ids) == dup)[0]][1:]
        nodes_to_update = set(nodes_to_update)
        tables.nodes.clear()
        for node in tsu.nodes():
            if node.id in nodes_to_update:
                new_metadata = node.metadata
                old_slim_id = node.metadata['slim_id']
                new_slim_id = int(old_slim_id + node.population * 1e12)
                output = [str(old_slim_id), str(new_slim_id), str(node.population)]
                new_metadata['slim_id'] = new_slim_id               
                tables.nodes.append(node.replace(metadata=new_metadata))
                with open(working_dir + "slim_id_update.csv", "a") as file:
                    file.write(",".join(output) + "\n")
            else:
                tables.nodes.append(node)
        tsu = tables.tree_sequence()
    # write out result
    tsu.dump(args.output)

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
