#!/bin/bash

# SET-UP ----
# set directories
export working_dir="/Users/shyamag/Dropbox/Work/Research/Projects/hybrid_sims_chapter/virus_evol_example/"
export results_dir=${working_dir}"results/"
# point to parameter file
export param_file=${working_dir}"sim_parameters.csv"
# point to infection events file
export SLiM_seq=${results_dir}"inf_sequence.csv"
# point to SLiM file
export SLiM_model=${working_dir}virus_in_host.slim
# specify which version of slim
export which_slim="/Users/shyamag/miniconda3/envs/slim-env/bin/slim"
# specify which runs to launch
export runs_to_launch="$(seq 10 | tr '\n' ',')"

# SIMULATE AN INFECTION SCENARIO
sh 1_model_infection_scenario.sh

# MAKE DAGS AND SNAKEFILES ----
sh 2_make_snakefile.sh

# RUN SLiM SEQUENCES ----
sh 3_run_SLiM.sh

# MERGE FINAL TREE SEQUENCES
sh 4_merge_ts.sh