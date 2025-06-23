#!/bin/bash

# SIMULATE AN INFECTION SCENARIO
sh 1_model_infection_scenario.sh

# MAKE DAGS AND SNAKEFILES ----
sh 2_make_snakefile.sh

# RUN SLiM SEQUENCES ----
sh 3_run_SLiM.sh

# MERGE FINAL TREE SEQUENCES
sh 4_merge_ts.sh

# MAKE PLOTS
sh 5_plot_ts.sh
