#!/bin/bash

# SET-UP
param_file="$PWD/sim_parameters.csv"
working_dir="$PWD/results/"

# RUN SLiM SEQUENCES ----
for outpath in ${working_dir}sim*; do
    start=$SECONDS
    snakefile="${outpath}/snakefile"
    echo $snakefile
    if [[ -f ${snakefile} ]];then
        if [[ "$(uname)" == "Darwin" ]]; then sed -i '' 's/\r//g' ${snakefile}; else sed -i 's/\r//g' ${snakefile};fi
        snakemake --snakefile ${snakefile} --unlock --directory ${outpath}
        if [[ ! -f ${outpath}/DAG.png ]];then
            snakemake --snakefile ${snakefile} --rulegraph | dot -Tpng > ${outpath}/DAG.png
        fi
        snakemake --snakefile ${snakefile} --cores 1 --keep-going --rerun-incomplete --directory ${outpath}
    fi
    end=$SECONDS
    echo 'running SLiM for task' ${outpath} 'took' $(($end - $start)) 'seconds to complete'
done
