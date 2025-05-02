## SCRIPT: merge_ts.py
## AUTHOR: Shyamalika Gopalan
## OVERVIEW: Merges a set of tree sequences that are phylogenetically related to each other. Conceptually similar to: https://tskit.dev/pyslim/docs/stable/vignette_parallel_phylo.html
## REQUIRED ARGUMENTS:
##		   -i, --input_dir (path to directory containing tree sequences to be merged)
##		   -I, --inf_sequence (path to infection sequence file in .csv format)
##		   -o, --output (desired file path for output tree sequence)
##
## USAGE EXAMPLE:
## $ python merge_ts.py --input_dir ~/test/run/results/ --inf_sequence inf_sequence.csv --output ~/test/run/output/output

# Import necessary libraries
import argparse as ap
import sys
import pandas as pd
import numpy as np
import tskit
import tszip
import os

## SUB-FUNCTIONS ----
## ------------------
# Create SLiM id to node id dictionary
def SLiM_to_node_dict(tree_seq):
	id_map = {}
	for node in tree_seq.nodes():
		slim_id = node.metadata['slim_id']
		assert slim_id not in id_map
		id_map[slim_id] = node.id
	return id_map

# Create node id to SLiM id dictionary
def node_to_SLiM_dict(tree_seq, as_string=False):
	id_map = {}
	for node in tree_seq.nodes():
		slim_id = node.metadata['slim_id']
		if as_string:
			slim_id = str(slim_id)
		id_map[node.id] = slim_id
	return id_map

# Update tree ages
def update_tree_ages(tree_seq, t_to_add):
	ts_tables = tree_seq.dump_tables()
	ts_tables.nodes.clear()
	ts_tables.mutations.clear()
	for node in tree_seq.nodes():
		node.id
		ts_tables.nodes.add_row(flags = node.flags, \
								time = node.time + t_to_add, \
								population = node.population, \
								individual = node.individual, \
								metadata = node.metadata)
	for mutation in tree_seq.mutations():
		ts_tables.mutations.add_row(site = mutation.site, \
			      					node = mutation.node, \
									time = mutation.time + t_to_add, \
									derived_state = mutation.derived_state, \
									parent = mutation.parent, \
									metadata = mutation.metadata)
	# update top level metadata
	top_level_meta = tree_seq.metadata
	top_level_meta['SLiM']['tick'] = top_level_meta['SLiM']['tick'] + t_to_add
	ts_tables.metadata_schema = tskit.MetadataSchema(tree_seq.metadata_schema.schema.copy())
	ts_tables.metadata = top_level_meta
	new_tree_seq = ts_tables.tree_sequence()
	return new_tree_seq

def list_ancestors(inf, inf_sequence):
	get_inf_source = inf_sequence.set_index('inf_id')['inf_source'].to_dict()
	if inf not in [0, '0']:
		current_inf = inf
		source = [get_inf_source[current_inf]]
		while source[len(source) - 1] not in [0, '0']:
			current_inf = source[len(source) - 1]
			source.append(get_inf_source[current_inf])
		return source
	else:
		return []

# Determine most recent common ancestors of all pairs of sample infections
def calc_pairwise_relatedness(samples, ancestor_list):
	pairwise_relatedness_info = []
	for i in range(len(samples)):
		sample_1 = samples[i]
		ancestors_1 = set(ancestor_list[sample_1])
		for j in range(i+1, len(samples)):
			sample_2 = samples[j]
			ancestors_2 = set(ancestor_list[sample_2])
			if sample_1 in ancestor_list[sample_2] or sample_2 in ancestor_list[sample_1]:
				continue
			if sample_1 != sample_2:
				common_ancestors = ancestors_1.intersection(ancestors_2)
				mrca = max(int(x) for x in common_ancestors)
				# append data
				pairwise_relatedness_info.append([sample_1, sample_2, mrca])
	pairwise_relatedness_info = pd.DataFrame(pairwise_relatedness_info, columns = ['sample_1', 'sample_2', 'mrca'])
	pairwise_relatedness_info = pairwise_relatedness_info.sort_values(by='mrca', key=lambda x: x.astype(int), ascending=False).reset_index(drop = True)
	return(pairwise_relatedness_info)

# Add missing edges to merged tree
def add_missing_edges(merged_tree_seq, tree_seq, verbose = False):
	# make slim <> node id maps
	merged_tree_slim_to_node_map = SLiM_to_node_dict(merged_tree_seq)
	merged_tree_node_to_slim_map = node_to_SLiM_dict(merged_tree_seq)
	tree_slim_to_node_map = SLiM_to_node_dict(tree_seq)
	tree_node_to_slim_map = node_to_SLiM_dict(tree_seq)

	# determine what the ultimate roots are
	ultimate_roots = set([tree_node_to_slim_map[x] for x in np.where(tree_seq.nodes_time == tree_seq.max_root_time)[0]] + [merged_tree_node_to_slim_map[x] for x in np.where(merged_tree_seq.nodes_time == merged_tree_seq.max_root_time)[0]])
	# determine which 'false' roots exist in the merged tree
	non_roots = set([merged_tree_node_to_slim_map[i] for t in merged_tree_seq.trees() for i in t.roots]) - ultimate_roots

	# remove false roots that cannot be rooted with the information in this tree_seq
	unrootable = non_roots - set(tree_slim_to_node_map.keys())
	non_roots = non_roots - unrootable

	# add missing edges and iterate until no more non_roots are left
	merged_tree_updated = merged_tree_seq
	OLD_non_root_intervals = None
	OLD_n_non_roots = None
	if len(non_roots) != OLD_n_non_roots and verbose:
		if len(non_roots) == 1:
			print("1 node to root", flush = True)
		else:
			print(str(len(non_roots)), "nodes to root", flush = True)
	while len(non_roots) > 0:
		focal_slim = next(iter(non_roots))
		# for each non_root, determine intervals where merged_tree_seq needs to be rooted
		non_root_intervals = []
		node_id = merged_tree_slim_to_node_map[focal_slim]
		merged_simp = merged_tree_seq.simplify(samples=[node_id], keep_unary = True, keep_input_roots = True, filter_nodes = False)
		for tree in merged_simp.trees():
			if node_id in tree.roots:
				non_root_intervals.append(tree.interval)
		# simplify tree_seq on this node and trim to relevant intervals
		if non_root_intervals != OLD_non_root_intervals:
			tree_trimmed = tree_seq.keep_intervals([[i.left, i.right] for i in non_root_intervals], simplify=False)
		tree_simp = tree_trimmed.simplify(samples=[tree_slim_to_node_map[focal_slim]], keep_unary = True, keep_input_roots = True, filter_nodes = False)
		# add edges where the false root is a child in tree_seq
		edges_to_add = [i for i in tree_simp.edges() if i.child == tree_slim_to_node_map[focal_slim]]
		# if there are no relevant edges, count this node as 'unrootable' and move on
		if len(edges_to_add) == 0:
			unrootable.add(focal_slim)
			non_roots.remove(focal_slim)
			continue
		# add edges from tree_seq to merged_tree_seq
		merged_tables = merged_tree_updated.dump_tables()
		for i in edges_to_add:
			parent_slim = tree_node_to_slim_map[i.parent]
			merged_tables.edges.add_row(left = i.left, \
										right = i.right, \
										parent = merged_tree_slim_to_node_map[parent_slim], \
										child = merged_tree_slim_to_node_map[focal_slim], \
										metadata = i.metadata)
		# is there some where to keep going without doing this? ##
		merged_tables.sort()
		try:
			merged_tree_updated = merged_tables.tree_sequence()
		except:
			unrootable.add(focal_slim)
		# re-compute non_roots
		OLD_n_non_roots = len(non_roots)
		non_roots = set([merged_tree_node_to_slim_map[i] for t in merged_tree_updated.trees() for i in t.roots]) - ultimate_roots - unrootable
		n_non_roots = len(non_roots)
		OLD_non_root_intervals = non_root_intervals
		if n_non_roots != OLD_n_non_roots and verbose:
			if len(non_roots) == 1:
				print("1 node to root", flush = True)
			else:
				print(str(len(non_roots)), "nodes to root", flush = True)
	return merged_tree_updated

# Merge two trees
def merge_ts(tree_seq1, tree_seq2, split_time, reroot_nodes = True, verbose = False):
	# create dataframes to store tree_seq1 and tree_seq2 information
	ts1_dat = pd.DataFrame({'ts1_idx': range(tree_seq1.num_nodes), 'slim_id': [x.metadata['slim_id'] for x in tree_seq1.nodes()], 'time': tree_seq1.nodes_time})
	ts2_dat = pd.DataFrame({'ts2_idx': range(tree_seq2.num_nodes), 'slim_id': [x.metadata['slim_id'] for x in tree_seq2.nodes()], 'time': tree_seq2.nodes_time})
	common_node_dat = pd.merge(ts1_dat, ts2_dat, on=['slim_id', 'time'])
	common_node_dat = common_node_dat[common_node_dat['time'] >= split_time]
	# if there are shared nodes, fill in node map
	if len(common_node_dat) > 0:
		# create null node map
		node_map = np.full(tree_seq1.num_nodes, tskit.NULL)
		# fill in null node map
		node_map[common_node_dat.ts1_idx] = common_node_dat.ts2_idx
	# merge trees
	merged_ts = tree_seq2.union(tree_seq1, node_map, check_shared_equality = False, add_populations = False, record_provenance = False)
	if reroot_nodes:
		if verbose:
			print("Rerooting nodes", flush = True)
		merged_ts = add_missing_edges(merged_ts, tree_seq1, verbose=verbose)
	return merged_ts

## MAIN FUNCTION ----
## ------------------
def main(args):
	argp = ap.ArgumentParser(description = 'Merges a set of tree sequences that are phylogenetically related to each other')
	argp.add_argument('-i', '--input_dir', type = str, help = 'Path to directory containing tree sequences to be merged')
	argp.add_argument('-I', '--inf_sequence', type = str, help = 'Path to infection sequence file in .csv format')
	argp.add_argument('-o', '--output', type = str, help = 'Path to desired output file')
	args = argp.parse_args(args)

	## Load input data ---
	# Load infection sequence file
	try:
		inf_sequence = pd.read_csv(args.inf_sequence)
		sampled_infs = list(inf_sequence.loc[inf_sequence['inf_step'] == max(inf_sequence['inf_step'])]['inf_id'])
	except FileNotFoundError:
		print('Infection sequence file not found', flush = True)

	# Load all sample trees
	input_dir_files = os.listdir(args.input_dir)
	n_failed_samples = 0
	sample_ts = {}
	print('Loading sample trees...', flush = True)
	for i in sampled_infs:
		dat = inf_sequence[inf_sequence.inf_id==i].iloc[0]
		file_name = f"inf{i}_{dat.host_id}_on_inf_day_{dat.transmission_day}_on_overall_day_{dat.overall_day}.ts"
		tree_exists = False
		for file in input_dir_files:
			if file_name in file:
				try:
					tree = tskit.load(args.input_dir + '/' + file)
					tree_exists = True
				except:
					print(' * ' + file + ' could not be loaded', flush = True)
		if tree_exists:
			sample_ts[i] = tree
		else:
			print(' * tree sequence file matching ' + file_name + ' could not be found - are you sure it was created?', flush = True)
			n_failed_samples += + 1

	if len(sample_ts) > 1:
		## Adjust tree ages ---
		print('Adjusting tree ages...', flush = True)
		# Find sample tree with oldest ancestry
		tree_ages = {key: tree.metadata['SLiM']['tick'] for key, tree in sample_ts.items()}
		max_age = max(tree_ages.values())
		# Adjust node ages to match up across all trees
		adjusted_sample_ts = {}
		for ts in sample_ts.keys():
			tree_age = tree_ages[ts]
			t_to_add = max_age - tree_age
			if t_to_add > 0:
				adjusted_sample_ts[ts] = update_tree_ages(sample_ts[ts], t_to_add = t_to_add)
			else:
				adjusted_sample_ts[ts] = sample_ts[ts]

		# Identify ancestors of all sample infections ----
		ancestor_list0 = {}
		print('Finding sample ancestors...', flush = True)
		for i in sampled_infs:
			ancestor_list0[i] = list_ancestors(i, inf_sequence)
		ancestor_list = ancestor_list0.copy()

		# Determine most recent common ancestor between all pairs of infections
		print('Calculating pairwise relatedness metrics...', flush = True)
		pairwise_sample_relatedness0 = calc_pairwise_relatedness(sampled_infs, ancestor_list)
		pairwise_sample_relatedness = pairwise_sample_relatedness0.copy()

		# Determine merge order - resolve groups of nodes from tips to root in order of MRCA age
		print('Determining merge order...', flush = True)
		sampled_infs1 = sampled_infs.copy()
		n_sampled_infs = len(sampled_infs1)
		samples_added = 1
		pairwise_sample_relatedness['sister_pair'] = list(zip(pairwise_sample_relatedness.sample_1, pairwise_sample_relatedness.sample_2))
		while samples_added > 0:
			# determine closest shared ancestors among all pairwise relationships
			while not pairwise_sample_relatedness.empty:
				mrca = pairwise_sample_relatedness.iloc[0]['mrca']
				sister_samples = {pairwise_sample_relatedness.iloc[0]['sample_1'], pairwise_sample_relatedness.iloc[0]['sample_2']}
				# of all other samples that share this mrca, is it the most recent shared ancestor for them too?
				other_samples_df = pairwise_sample_relatedness[pairwise_sample_relatedness['mrca'] == mrca]
				other_samples = set(other_samples_df.sample_1.tolist() + other_samples_df.sample_2.tolist()) - sister_samples
				if len(other_samples) > 0:
					for i in other_samples:
						all_mrcas = set(pairwise_sample_relatedness.mrca[(pairwise_sample_relatedness.sample_1 == i) | (pairwise_sample_relatedness.sample_2 == i)])
						# if yes for any, add those samples to sister_samples group
						if str(max(map(int, all_mrcas))) == mrca:
							sister_samples.add(i)
				# add mrca to sample list if it's not already present
				if mrca not in sampled_infs1:
					sampled_infs1.append(mrca)
					ancestor_list[mrca] = list_ancestors(mrca, inf_sequence)
				# exclude all sister samples from further consideration
				pairwise_sample_relatedness = pairwise_sample_relatedness[~pairwise_sample_relatedness.index.isin(
					pairwise_sample_relatedness.index[pairwise_sample_relatedness.sample_1.isin(sister_samples) | 
					pairwise_sample_relatedness.sample_2.isin(sister_samples)])].reset_index(drop=True)
			samples_added = len(sampled_infs1) - n_sampled_infs
			n_sampled_infs = len(sampled_infs1)
			pairwise_sample_relatedness = calc_pairwise_relatedness(list(sampled_infs1), ancestor_list)
		tmp = pairwise_sample_relatedness.copy()
		# final pass
		merge_order = pd.DataFrame(columns = tmp.columns)
		while not tmp.empty:
			mrca = tmp.iloc[0]['mrca']
			sister_samples = {tmp.iloc[0]['sample_1'], tmp.iloc[0]['sample_2']}
			idx = [0]
			# of all other samples that share this mrca, is it the most recent shared ancestor for them too?
			other_samples_idx = [x for x in tmp.index[tmp.mrca == mrca].tolist() if x not in [0, '0']]
			if other_samples_idx:
				other_samples = set(tmp.loc[other_samples_idx, 'sample_1'].tolist() + tmp.loc[other_samples_idx, 'sample_2'].tolist()) - sister_samples
				for i in other_samples:
					all_mrcas = set(tmp.loc[(tmp['sample_1'] == i) | (tmp['sample_2'] == i), 'mrca'].tolist())
					# if yes for any, add those samples to sister_samples group
					if str(max(map(int, all_mrcas))) == mrca:
						sister_samples.add(i)
						idx.append(other_samples_idx[next(index for index, value in enumerate(other_samples_idx) if tmp.iloc[value]['sample_1'] == i or tmp.iloc[value]['sample_2'] == i)])
			# exclude all sister samples from further consideration
			merge_order = pd.concat([merge_order, tmp.iloc[idx]]).reset_index(drop = True)
			tmp = tmp[~tmp.index.isin(tmp.index[(tmp['sample_1'].isin(sister_samples)) | (tmp['sample_2'].isin(sister_samples))])].reset_index(drop=True)
		del tmp
		inf_day_map = inf_sequence.set_index('inf_id')['overall_day'].to_dict()
		merge_order["split_day"] = merge_order["mrca"].map(inf_day_map)

		# Merge
		unique_mrcas = merge_order.mrca.unique()
		print('Merging a total of ' + str(len(sampled_infs)) + ' sample trees in ' + str(len(unique_mrcas)) +  ' merging steps', flush = True)
		for i in unique_mrcas:
			print('Resolving infection ' + str(i) + ' (' + str(1 + np.where(unique_mrcas == i)[0][0]) + ' of ' + str(len(unique_mrcas)) + ')' , flush = True)
			dat = merge_order[merge_order.mrca == i].reset_index(drop = True)
			focal_samples = set(dat['sample_1']).union(dat['sample_2'])
			samples_added = set()
			# merge samples pairwise
			while len(focal_samples) > 0:
				if len(samples_added) == 0:
					inf1 = list(focal_samples)[0]
					inf2 = list(focal_samples)[1]
				else:
					inf1 = i
					inf2 = list(focal_samples)[0]
				ts1, ts2 = adjusted_sample_ts[inf1], adjusted_sample_ts[inf2]
				adjusted_sample_ts[i] = merge_ts(ts1, ts2, split_time=max_age-dat.split_day[0]-1, reroot_nodes=False, verbose=False)
				for x in [inf1, inf2]:
					if x in focal_samples:
						focal_samples.discard(x)
					if x != i:
						samples_added.add(x)
			merged_ts = adjusted_sample_ts.get(i)
	else:
		print('Only 1 sample tree was found and processed', flush = True)
		merged_ts = sample_ts[list(sample_ts.keys())[0]]

	# Simplify, removing seed nodes, and write out merged tree
	print('Simplifying and writing out final tree...', flush = True)
	# simplify tree
	simple_ts = merged_ts.simplify(keep_input_roots = True)
	# write out
	tszip.compress(output_ts, args.output + '.tsz')

	with open("/Users/shyamag/Desktop/tree_output.txt", "w") as f:
		f.write(simple_ts.draw_text())

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
