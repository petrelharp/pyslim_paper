#!/bin/bash
# SIMULATE AN INFECTION SCENARIO
mkdir -p $results_dir
Rscript ${working_dir}model_infection_scenario.R ${results_dir}