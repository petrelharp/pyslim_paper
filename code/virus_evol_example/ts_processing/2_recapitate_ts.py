## SCRIPT: recapitate_ts.py
## AUTHOR: Shyamalika Gopalan
## OVERVIEW: Recapitates a tree sequence.
## REQUIRED ARGUMENTS:
##         -i, --input_ts (path to merged tree sequence file)
##         -N, --Ne (effective population size)
##         -s, --seed (optional seed)
##         -o, --output (desired file path for output files)
##
## USAGE EXAMPLE:
## $ python recapitate_ts.py --input_ts ~/results/input_tree.trees.tsz --Ne 30000 --output ~/results/output_tree

# Import necessary libraries
import argparse as ap
import sys
import pyslim
import tszip
import pandas as pd
import msprime
import warnings

## MAIN FUNCTION ----
## ------------------
def main(args):
	argp = ap.ArgumentParser(description='Recapitates a tree sequence with multiple roots')
	argp.add_argument('-i', '--input_ts', type=str, help='Name of zipped tree sequence file to be recapitated')
	argp.add_argument('-N', '--Ne', type=int, help='Effective population size')
	argp.add_argument('-s', '--seed', type=int, help='Optional seed')
	argp.add_argument('-o', '--output', type=str, help='Desired output file path and name')
	args = argp.parse_args(args)

	# Read in tree sequence
	print('Loading tree sequence', flush=True)
	warnings.simplefilter('ignore', UserWarning)
	tree_seq = pyslim.update(tszip.decompress(args.input_ts))

	# Recapitate tree sequence
	print('Recapitating tree sequence', flush=True)
	demography = msprime.Demography()
	demography.add_population(name = "p0", initial_size=int(args.Ne))
	warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)
	if args.seed is None:
		recapped_ts = msprime.sim_ancestry(initial_state=tree_seq, demography = demography)
	else:
		recapped_ts = msprime.sim_ancestry(initial_state=tree_seq, demography = demography, random_seed = int(args.seed))

	# Write out tree sequence
	print('Writing out tree sequence', flush=True)
	tszip.compress(ts=recapped_ts, destination=args.output + '.tsz')

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))