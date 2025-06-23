#!/bin/bash

results_dir="results/"
plotting_script="../../plotting/back-and-forth.py"

for outpath in results/sim*; do
    if [ -d $outpath ]; then
        start=$SECONDS
        python ${plotting_script} ${outpath}.merged.trees
        end=$SECONDS
        echo 'processing for task' ${outpath} 'took' $(($end - $start)) 'seconds to complete'
    fi
done
