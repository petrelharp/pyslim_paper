#!/bin/bash
python make_snakefile.py --input $1 --outfile $2
snakemake -s $2 --rulegraph | dot -Tsvg > DAG.svg
snakemake -s $2 -c 4
