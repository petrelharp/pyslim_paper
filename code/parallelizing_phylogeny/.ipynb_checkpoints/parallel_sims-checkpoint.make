all: root.trees C.trees I.trees B.trees A.trees

root.trees:  ./simulate_branch.slim
	slim -d "infile=''" -d popsize=500 -d "popname='root'" -d num_gens=2000 -d "outfile='root.trees'" {path_to_slimscript}

C.trees: root.trees ./simulate_branch.slim
	slim -d "infile='root.trees'" -d popsize=50 -d "popname='C'" -d num_gens=250 -d "outfile='C.trees'" {path_to_slimscript}

I.trees: root.trees ./simulate_branch.slim
	slim -d "infile='root.trees'" -d popsize=100 -d "popname='I'" -d num_gens=200 -d "outfile='I.trees'" {path_to_slimscript}

B.trees: I.trees ./simulate_branch.slim
	slim -d "infile='I.trees'" -d popsize=70 -d "popname='B'" -d num_gens=50 -d "outfile='B.trees'" {path_to_slimscript}

A.trees: I.trees ./simulate_branch.slim
	slim -d "infile='I.trees'" -d popsize=40 -d "popname='A'" -d num_gens=50 -d "outfile='A.trees'" {path_to_slimscript}

