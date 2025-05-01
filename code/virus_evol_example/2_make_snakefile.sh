#!/bin/bash
#SBATCH --time=1:00:00

# SET-UP
param_file="/data2/gopalan_lab/shyamag/hybrid_sims_chapter/virus_evol_example/sim_parameters.csv"
results_dir="/data2/gopalan_lab/shyamag/hybrid_sims_chapter/virus_evol_example/results/"
which_slim="/home/shyamag/.conda/envs/slim-env/bin/slim"
SLiM_model="/data2/gopalan_lab/shyamag/hybrid_sims_chapter/virus_evol_example/virus_in_host.slim"
SLiM_seq=${results_dir}inf_sequence.csv
runs_to_launch=($(seq 10))

# MAKE DAGS AND SNAKEFILES ----
for task in ${runs_to_launch[@]};do
	start=$SECONDS
	echo Working on task $task
	# specify which simulation to run
	sim_id=$(head -n1 $param_file | cut -f$((1+${task})) -d, | tr -d '\r')
	
	# specify output path
	outpath=${results_dir}"${sim_id}/"
	mkdir -p $outpath
	
	# find input file
	if [[ ! -f ${SLiM_seq} ]];then
		echo "SLiM sequence could not be found"
	fi
	
	# Create snakefile
	if [[ -f ${outpath}snakefile ]];then
		rm ${outpath}snakefile
	fi
	if [[ -f ${outpath}.tmp.all.input ]];then
		rm ${outpath}.tmp.all.input
	fi
	if [[ -f ${outpath}.tmp.snakefile ]];then
		rm ${outpath}.tmp.snakefile
	fi
	if [[ ! -f ${outpath}snakefile ]] && [[ -f ${SLiM_seq} ]];then
		seed=${RANDOM}
		N=$(wc -l ${SLiM_seq} | awk '{print $1}')
		max_inf_step=$(tail -n+2 $SLiM_seq | cut -f1 -d, | sort -n | uniq | tail -n1)
	
		# Create reproducible list of random seeds
		range=($(seq 99999))
		seeds=($(printf "%s\n" "${range[@]}" | awk -v seed="$seed" 'BEGIN {srand(seed)} {print $0, rand()}' | sort -k2,2 | cut -d ' ' -f1 | head -n $((N-1))))
	
		printf "** Generating snakemake file... "
		printf "rule all:\n\tinput:\n" > ${outpath}.tmp.all.input
		for i in `seq 2 $N`;do
			# Extract information for this run
			IFS="," read -a run <<< "$(sed -n ${i}p ${SLiM_seq})"
			inf_step=${run[0]}
			infection_idx=${run[1]}
			infection_source=${run[2]}
			infector_id=${run[3]}
			infected_ids=${run[4]}
			transmission_day=${run[5]}
			overall_day=${run[6]}
	
			# Specify rule name
			target="inf${infection_idx}"
	
			# Add required input
			printf "rule ${target}:\n\tinput:\n\t\t\"${SLiM_model}\",\n" > ${outpath}.tmp.snakefile
			if [[ $infector_id != "seed" ]];then
				n=$(awk '{print $2}' FS=, ${SLiM_seq} | grep -nw $infection_source | awk '{print $1}' FS=:)
				IFS="," read -a source_line <<< "$(sed -n ${n}p ${SLiM_seq})"
				source_transmission_day=${source_line[5]}
				source_overall_day=${source_line[6]}
				infile=$(echo ${outpath}inf${infection_source}_${source_line[3]}_on_inf_day_${source_transmission_day}_on_overall_day_${source_overall_day}.ts)
				printf "\t\t\"${infile}\"\n" >> ${outpath}.tmp.snakefile
			else
				# Remove the trailing comma from the last line of the temporary snakefile - syntax depends on OS
				if [[ "$(uname)" == "Darwin" ]]; then
				    sed -i '' '$s/,$//' ${outpath}.tmp.snakefile # macOS `sed` requires an empty string argument after `-i`
				else
				    sed -i '$s/,$//' ${outpath}.tmp.snakefile # Linux `sed` modifies file in place directly
				fi
				inf_day=0
				infile=NULL
			fi
	
			# Add expected output
			printf "\toutput:\n" >> ${outpath}.tmp.snakefile
			outfile[${x}]=$(echo ${outpath}${target}_${infector_id}_on_inf_day_${transmission_day}_on_overall_day_${overall_day}.ts)
			if [[ ${inf_step} -lt ${max_inf_step} ]];then # if we are NOT in the last time step of the simulation, mark output file as temporary
				printf "\t\ttemporary(\"${outfile}\")\n" >> ${outpath}.tmp.snakefile
			else # if we ARE in the last time step, also add this file to 'rule all'
				printf "\t\t\"${outfile}\"\n" >> ${outpath}.tmp.snakefile
				printf "\t\t\"${outfile}\",\n" >> ${outpath}.tmp.all.input
			fi
			# Add recipe
			# write out each line of recipe...
			printf "\tshell:\n\t\t\"\"\"\n" >> ${outpath}.tmp.snakefile
			recipe=$(printf "%s -l 0 -seed %s -d \"INF_ID='%s'\"" "$which_slim" "${seeds[$((i-2))]}" "$target")
			recipe=$(echo ${recipe} $(printf -- "-d \"INFILE='%s'\"" "${infile}"))
			recipe=$(echo ${recipe} $(printf -- "-d \"INPUT_HOST_ID='%s'\"" "${infector_id}"))
			recipe=$(echo ${recipe} $(printf -- "-d \"OUTPUT_HOST_IDS='%s'\"" "${infected_ids}"))
			recipe=$(echo ${recipe} $(printf -- "-d \"TRANSMISSION_DAY=%s\"" "$transmission_day"))
			recipe=$(echo ${recipe} $(printf -- "-d \"OUTPUT_DIR='%s'\"" "${outpath}"))
			n=$(grep -n '^N_FOUNDERS,' $param_file | awk -F: '{print $1}')
			N_founders=$(sed -n "${n}p" $param_file | cut -f$((1+${task})) -d, | tr -d '\r')
			recipe=$(echo ${recipe} $(printf -- "-d \"N_FOUNDERS=%s\"" "$N_founders"))
			params=($(cut -f1 -d, $param_file | tail -n+2))
			for param in ${params[@]};do
				# if parameter is not 'N_founders', continue
				if [[ $param != "N_FOUNDERS" ]];then
					n=$(grep -n ^${param} $param_file | awk -F: '{print $1}')
					param_val=$(sed -n "${n}p" $param_file | cut -f$((1+${task})) -d, | tr -d '\r')
					# if parameter is not null, add it to the recipe
					if [[ ! -z $param_val ]];then
						recipe=$(echo ${recipe} $(printf -- "-d \"%s=%s\"" "$param" "$param_val"))
					fi
				fi
			done
			recipe=$(echo ${recipe} $(printf "${SLiM_model}\n"))
			echo "$(printf "\t\t\t")${recipe}" >> ${outpath}.tmp.snakefile
			printf "\t\t\"\"\"\n\n" >> ${outpath}.tmp.snakefile
			cat ${outpath}.tmp.snakefile >> ${outpath}snakefile
		done
		if [[ "$(uname)" == "Darwin" ]]; then sed -i '' '$s/,$/\n/' ${outpath}.tmp.all.input; else sed -i '$s/,$/\n/' ${outpath}.tmp.all.input;fi
		cat ${outpath}.tmp.all.input > ${outpath}.tmp.snakefile
		cat ${outpath}snakefile >> ${outpath}.tmp.snakefile
		mv ${outpath}.tmp.snakefile ${outpath}snakefile
		rm ${outpath}.tmp.all.input
		printf "done!\n"
	fi
	end=$SECONDS
	echo 'creating snakemake file for task' ${task} 'took' $(($end - $start)) 'seconds to complete'
done