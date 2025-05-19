#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --array=1-10

# SET-UP
param_file="/data2/gopalan_lab/shyamag/hybrid_sims_chapter/virus_evol_example/sim_parameters.csv"
working_dir="/data2/gopalan_lab/shyamag/hybrid_sims_chapter/virus_evol_example/results/"
source ~/.bashrc
conda activate slim-env

# RUN SLiM SEQUENCES ----
task=${SLURM_ARRAY_TASK_ID}
start=$SECONDS
sim_id=$(head -n1 $param_file | cut -f$((1+${task})) -d, | tr -d '\r')
outpath=${working_dir}${sim_id}"/"
snakefile=${outpath}"snakefile"
if [[ -f ${snakefile} ]];then
	if [[ "$(uname)" == "Darwin" ]]; then sed -i '' 's/\r//g' ${snakefile}; else sed -i 's/\r//g' ${snakefile};fi
	snakemake --snakefile ${snakefile} --unlock --directory ${outpath}
	if [[ ! -f ${outpath}/${sim_id}_DAG.png ]];then
		snakemake --snakefile ${snakefile} --rulegraph | dot -Tpng > ${outpath}/${sim_id}_DAG.png
	fi
	snakemake --snakefile ${snakefile} --cores 1 --keep-going --rerun-incomplete --directory ${outpath}
fi
end=$SECONDS
echo 'running SLiM for task' ${task} 'took' $(($end - $start)) 'seconds to complete'
