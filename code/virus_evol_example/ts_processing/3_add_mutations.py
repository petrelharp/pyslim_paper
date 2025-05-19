## SCRIPT: add_mutations.py
## AUTHOR: Shyamalika Gopalan
## OVERVIEW: Adds mutations to a recapitated tree, ensuring that no new mutations are added after the coalescent portion.
## REQUIRED ARGUMENTS:
##         -i, --input_tree (path to recapitated tree file)
##         -m, --mu (mutation rate)
##         -s, --seed (optional seed)
##         -o, --output (desired file path for output files)
##
## USAGE EXAMPLE:
## $ python add_mutations.py --input_tree ~/results/input_tree.trees.tsz --mu 0.0000000003 --output ~/results/output_tree

# Import necessary libraries
import argparse as ap
import sys
import pyslim
import tszip
import msprime

## MAIN FUNCTION ----
## ------------------
def main(args):
	argp = ap.ArgumentParser(description='Adds mutations to a recapitated tree, ensuring that no new mutations are added after the coalescent portion')
	argp.add_argument('-i', '--input_tree_seq', type=str, help='Name of zipped tree sequence file to add mutations to')
	argp.add_argument('-m', '--mu', type=float, help='Per base pair mutation rate')
	argp.add_argument('-s', '--seed', type=int, help='Optional seed')
	argp.add_argument('-o', '--output', type=str, help='Desired output file path and name')
	args = argp.parse_args(args)

	# Read in tree sequence
	print('Loading tree', flush=True)
	tree_seq = tszip.decompress(args.input_tree_seq)

	# Add mutations to the tree sequence
	print('Adding neutral mutations', flush=True)
	mut_type = max([i.metadata['mutation_list'][0]['mutation_type'] for i in tree_seq.mutations()]) + 1
	if args.seed is None:
		mut_ts = msprime.sim_mutations(tree_seq, rate=float(args.mu), start_time=0, model=msprime.SLiMMutationModel(type = mut_type, next_id=pyslim.next_slim_mutation_id(tree_seq)), keep=True)
	else:
		mut_ts = msprime.sim_mutations(tree_seq, rate=float(args.mu), start_time=0, model=msprime.SLiMMutationModel(type = mut_type, next_id=pyslim.next_slim_mutation_id(tree_seq)), keep=True, random_seed=int(args.seed))

	# Write out tree
	print('Writing out tree', flush=True)
	tszip.compress(ts=mut_ts, destination=args.output + '.tsz')

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))