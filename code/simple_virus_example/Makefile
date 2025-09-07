all.trees : sim_04.trees sim_05.trees sim_03.trees
	python merge_ts.py -o $@ -i $^

sim_00.trees :
	slim -d "HOST_ID='00'" -d "OUTPUT_HOST_IDS=c('01','02')" \
		-d "ABS_TRANSMISSION_DAYS=c(10,20)" -d "ABS_SAMPLING_DAY=20" \
		virus_in_host.slim

sim_00.%.trees : sim_00.trees
	python reset.py $< $* $@

sim_01.%.trees : sim_01.trees
	python reset.py $< $* $@

sim_02.%.trees : sim_02.trees
	python reset.py $< $* $@

sim_01.trees : sim_00.01.trees
	slim -d "INFILE='$<'" -d "HOST_ID='01'" -d "OUTPUT_HOST_IDS=c('03','04')" \
		-d "ABS_TRANSMISSION_DAYS=c(15,35)" -d "ABS_SAMPLING_DAY=36" \
		virus_in_host.slim

sim_02.trees : sim_00.02.trees
	slim -d "INFILE='$<'" -d "HOST_ID='02'" -d "OUTPUT_HOST_IDS='04'" \
		-d "ABS_TRANSMISSION_DAYS=35" -d "ABS_SAMPLING_DAY=60" \
		virus_in_host.slim

sim_03.trees : sim_01.03.trees
	slim -d "INFILE='$<'" -d "HOST_ID='03'" -d "ABS_SAMPLING_DAY=20" virus_in_host.slim

sim_01m02.%.trees : sim_01.%.trees sim_02.%.trees
	python merge_ts.py -o $@ -i $^

sim_01_02.%.trees : sim_01m02.%.trees
	python reset.py $< $* $@

sim_04.trees : sim_01_02.04.trees
	slim -d "INFILE='$<'" -d "HOST_ID='04'" -d "OUTPUT_HOST_IDS='05'" \
		-d "ABS_TRANSMISSION_DAYS=60" -d "ABS_SAMPLING_DAY=60" \
		virus_in_host.slim

sim_05.trees : sim_04.05.trees
	slim -d "INFILE='$<'" -d "HOST_ID='05'" -d "ABS_SAMPLING_DAY=65" virus_in_host.slim

clean : 
	-rm *.trees
