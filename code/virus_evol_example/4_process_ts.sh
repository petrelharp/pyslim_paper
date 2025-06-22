#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --array=1-10

# SET-UP
param_file="sim_parameters.csv"
results_dir="results/"
inf_sequence=${results_dir}"inf_sequence.csv"
merging_script="merge_ts.py"

# MERGE TREE SEQUENCES ----
for outpath in results/sim*; do
    start=$SECONDS
    python ${merging_script} --input_dir ${outpath} --output ${outpath}merged.trees
    end=$SECONDS
    echo 'processing for task' ${outpath} 'took' $(($end - $start)) 'seconds to complete'
done
