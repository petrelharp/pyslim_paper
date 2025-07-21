import pandas as pd

# Paths
path_to_tsv = "./phylo.tsv"
path_to_make = "./parallel_sims.make"
path_to_slimscript = "./simulate_branch.slim"

# Reading the phylogeny data frame
df = pd.read_csv(path_to_tsv, sep="\t")
df = df.fillna('')

## Creating intermediate tree sequences filenames
df["infile"] = df.parent + ".trees"
df["outfile"] = df.child + ".trees"
df.loc[df["infile"]==".trees", "infile"] = ""
df["is_leaf"] = ~df.child.isin(df.parent) # setting nodes that are never parents as leaves

# Writing a Makefile 
f = open(path_to_make, "w")
print(f"all: {' '.join(df.outfile.to_list())}\n", file=f)
for i, row in df.iterrows():
    print(f"{row.outfile}: {row.infile} {path_to_slimscript}", file=f)
    print(f"\tslim -d \"infile='{row.infile}'\" -d popsize={row.popsize} "
          f"-d \"popname=\'{row.child}\'\" "
          f"-d num_gens={row.edgelen} " f"-d \"outfile='{row.child}.trees'\" "
          "{path_to_slimscript}\n",
          file=f)
f.close()
