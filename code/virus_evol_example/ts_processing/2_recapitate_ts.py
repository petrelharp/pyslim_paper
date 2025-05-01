## SCRIPT: recapitate_ts.py
## AUTHOR: Shyamalika Gopalan
## OVERVIEW: Recapitates a tree sequence, adds mutations on the branches at the specified rate, and writes out the recapitated tree.
## REQUIRED ARGUMENTS:
##         -i, --input_ts (path to merged tree sequence file)
##         -N, --Ne (effective population size)
##         -r, --recomb_rates (path to recombination rates file)
##         -p, --recomb_positions (path to recombination positions file)
##         -s, --seed (optional seed)
##         -o, --output (desired file path for output files)
##
## USAGE EXAMPLE:
## $ python recapitate_ts.py --input_ts ~/results/input_tree.trees.tsz --Ne 30000 --recomb_rates ~/recombination_rates --recomb_positions ~/recombination_positions --output ~/results/output_tree

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
	argp.add_argument('-r', '--recomb_rates', type=str, help='Path to recombination rates file')	
	argp.add_argument('-p', '--recomb_positions', type=str, help='Path to recombination positions file')
	argp.add_argument('-R', '--recomb_rescale_factor', type=float, default=1., help='Factor by which to rescale recombination rate. Values > 0 and < 1 slow recombination')
	argp.add_argument('-s', '--seed', type=int, help='Optional seed')
	argp.add_argument('-o', '--output', type=str, help='Desired output file path and name')
	args = argp.parse_args(args)

	# Read in tree sequence
	print('Loading tree sequence', flush=True)
	warnings.simplefilter('ignore', UserWarning)
	tree_seq = pyslim.update(tszip.decompress(args.input_ts))

	# Read in recombination information
	rates = pd.read_csv(args.recomb_rates, header=None).iloc[:,0].tolist()
	rates = [float(x) * args.recomb_rescale_factor if x != 0.5 else float(x) for x in rates]
	positions = pd.read_csv(args.recomb_positions, header=None).iloc[:,0].tolist()
	positions.insert(0, 0) # add zero to start of positions
	positions[-1] += 1 # increment all endpoints by 1
	recomb_map = msprime.RateMap(position = positions, rate = rates)

	# Recapitate tree sequence
	print('Recapitating tree sequence', flush=True)
	demography = msprime.Demography()
	demography.add_population(name = "p0", initial_size=int(args.Ne))
	warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)
	if args.seed is None:
		recapped_ts = msprime.sim_ancestry(initial_state=tree_seq, recombination_rate = recomb_map, demography = demography)
	else:
		recapped_ts = msprime.sim_ancestry(initial_state=tree_seq, recombination_rate = recomb_map, demography = demography, random_seed = int(args.seed))

	# Write out tree sequence
	print('Writing out tree sequence', flush=True)
	tszip.compress(ts=recapped_ts, destination=args.output + '.tsz')

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))