#!/bin/bash
# first arg: infection sequence table to use (e.g. inf_seq_w_reticulation.csv)
# second arg: desired name for snakefile (e.g. snakefile_w_reticulation)

python make_snakefile.py --input $1 --outfile $2
snakemake -s $2 --rulegraph | dot -Tsvg > DAG.svg
snakemake -s $2 -c 4
