#!/bin/bash
#python make_snakefile.py --input $1 --outfile $2
#snakemake -s $2 --rulegraph | dot -Tsvg > DAG.svg
rm -f slim_id_update.csv
snakemake -s $2 -c 4
