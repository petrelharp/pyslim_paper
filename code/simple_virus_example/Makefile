all.trees : sim_04.trees sim_05.trees sim_03.trees
	python merge_ts.py -o $@ -i $^

sim_00.trees :
	slim -d "HOST_ID='00'" -d "OUTPUT_HOST_IDS=c('01','02')" -d "TRANSMISSION_DAYS=c(10,20)" -d "SAMPLING_DAY=20" virus_in_host.slim

sim_00.%.trees : sim_00.trees
	python reset.py $< $* $@

sim_01.%.trees : sim_01.trees
	python reset.py $< $* $@

sim_02.%.trees : sim_02.trees
	python reset.py $< $* $@

sim_01.trees : sim_00.10.trees
	slim -d "INFILE='$<'" -d "HOST_ID='01'" -d "OUTPUT_HOST_IDS=c('03','04')" -d "TRANSMISSION_DAYS=c(5,25)" -d "SAMPLING_DAY=25" virus_in_host.slim

sim_02.trees : sim_00.20.trees
	slim -d "INFILE='$<'" -d "HOST_ID='02'" -d "OUTPUT_HOST_IDS='05'" -d "TRANSMISSION_DAYS=40" -d "SAMPLING_DAY=40" virus_in_host.slim

sim_03.trees : sim_01.5.trees
	slim -d "INFILE='$<'" -d "HOST_ID='03'" -d "SAMPLING_DAY=25" virus_in_host.slim

sim_04.trees : sim_01.25.trees
	slim -d "INFILE='$<'" -d "HOST_ID='04'" -d "SAMPLING_DAY=35" virus_in_host.slim

sim_05.trees : sim_02.40.trees
	slim -d "INFILE='$<'" -d "HOST_ID='05'" -d "SAMPLING_DAY=15" virus_in_host.slim
