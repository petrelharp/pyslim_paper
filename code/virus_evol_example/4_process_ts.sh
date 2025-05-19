#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --array=1-10

# SET-UP
working_dir="/data2/gopalan_lab/shyamag/hybrid_sims_chapter/virus_evol_example/"
param_file=${working_dir}"sim_parameters.csv"
results_dir=${working_dir}"results/"
inf_sequence=${results_dir}"inf_sequence.csv"
merging_script=${working_dir}"ts_processing/1_merge_ts.py"
source ~/.bashrc
conda activate slim-env

# MERGE TREE SEQUENCES ----
task=${SLURM_ARRAY_TASK_ID}
start=$SECONDS
sim_id=$(head -n1 $param_file | cut -f$((1+${task})) -d, | tr -d '\r')
outpath=${results_dir}${sim_id}"/"

python ${merging_script} --input_dir ${outpath} --inf_sequence ${inf_sequence} --output ${outpath}merged.trees

end=$SECONDS
echo 'processing for task' ${task} 'took' $(($end - $start)) 'seconds to complete'
