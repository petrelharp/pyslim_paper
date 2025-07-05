import sys, glob
import argparse
import tskit
from merge_utils import *

def main(args):
    argp = argparse.ArgumentParser(description = 'Merges a set of tree sequences that are related to each other')
    argp.add_argument('-i', '--input', type = str, nargs='+', help = 'All tree sequences to be merged', required=True)
    argp.add_argument('-o', '--output', type = str, help = 'Path to desired output file', required=True)
    args = argp.parse_args(args)

    input_dir_files = glob.glob(f"{args.input_dir}/*.trees")
    tslist = [clear_alive(tskit.load(x)) for x in input_dir_files]
    ticks = [ts.metadata['SLiM']['tick'] for ts in tslist]
    dt = [max(ticks) - t for t in ticks]
    tslist = [shift_times(ts, t) for ts, t in zip(tslist, dt)]

    tsu = tslist[0]
    for k in range(1, len(tslist)):
        ts = tslist[k]
        tsu = merge_ts(tsu, ts)

    tsu.dump(args.output)


if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
