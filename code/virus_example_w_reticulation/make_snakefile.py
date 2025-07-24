import pandas as pd
import numpy as np
import re
import os, sys
import argparse

# define a function to extract and flatten values from a data frame column
def extract_ids(value):
    if pd.isna(value):
        return []
    value = str(value)
    return re.findall(r"['\"]?(\d+)['\"]?", value)

def main(args):
    argp = argparse.ArgumentParser(description = 'Generate snakemake file from a infection sequence table')
    argp.add_argument('-i', '--input', type = str, help = 'Infection sequence table', required=True)
    argp.add_argument('-o', '--outfile', type = str, help = 'Desired output file name', required=True)
    args = argp.parse_args(args)

    working_dir=os.getcwd() + "/"

    # read in infection sequence
    inf_seq = pd.read_csv(working_dir + args.input, dtype=str)

    # collect all host IDs for origin infections
    all_origin_ids = set(sum(inf_seq["origin_id"].apply(extract_ids), []))

    # parse infection sequence
    inf_seq["status"] = "internal"
    inf_seq.loc[inf_seq["origin_id"].isna(), "status"] = "seed"
    inf_seq.loc[~inf_seq["host_id"].isin(all_origin_ids), "status"] = "leaf"

    # write out rule_all
    with open(working_dir + args.outfile, "w") as file:
        file.write("rule all:\n\tinput:\n\t\t\"merged.trees\"\n\n")

    # write out each seed/internal infection rule
    files_to_merge = []
    with open(working_dir + args.outfile, "a") as file:
        for i in inf_seq.index:
            dat_row = inf_seq.iloc[i]
            # write rule name
            file.write("rule H" + dat_row.host_id + ":\n")
            # write input line
            if dat_row.status == "internal" or dat_row.status == "leaf":
                file.write("\tinput:\n")
                inf_sources = dat_row.origin_id.split(",")
                expected_infile = working_dir + '-'.join(inf_sources) + "_to_" + dat_row.host_id + ".trees"
                file.write("\t\t\"" + expected_infile + "\"\n")
            # write output line
            file.write("\toutput:\n")
            expected_outfiles = []
            output_infs = dat_row.output_host_ids.split(",")
            if dat_row.status == "seed" or dat_row.status == "internal":
                for j in output_infs:
                    expected_outfiles = expected_outfiles + [working_dir + dat_row.host_id + "_to_" + j + ".trees"]
                    if j == output_infs[-1]:
                        file.write("\t\t\"" + expected_outfiles[-1] + "\"\n")
                    else:
                        file.write("\t\t\"" + expected_outfiles[-1] + "\",\n")
            elif dat_row.status == "leaf":
                files_to_merge = files_to_merge + [working_dir + dat_row.host_id + ".trees"]
                file.write("\t\t\"" + files_to_merge[-1] + "\"\n")
            # WRITE RECIPE
            file.write("\tshell:\n\t\t\"\"\"\n\t\t\tslim ")
            # figure out what input file name is
            if dat_row.status == "internal" or dat_row.status == "leaf":
                file.write("-d \"INFILE=\'" + expected_infile + "\'\" \\\n")
            elif dat_row.status == "seed":
                file.write("-d \"INFILE=\'\'\" \\\n")
            # print rest of SLiM script - host ID
            file.write("\t\t\t\t-d \"HOST_ID=\'" + dat_row.host_id + "\'\" \\\n")
            # transmission days
            file.write("\t\t\t\t-d \"TRANSMISSION_DAYS=c(" + dat_row.transmission_days + ")\" \\\n")
            # output host IDs
            file.write("\t\t\t\t-d \"OUTPUT_HOST_IDS=c(")
            for j in output_infs:
                if j == output_infs[-1]:
                    file.write("\'" + j + "\'")
                else:
                    file.write("\'" + j + "\',")
            file.write(")\" \\\n\t\t\t\t" + working_dir + "virus_in_host.slim\n")
            # print python script
            if dat_row.status == "seed" or dat_row.status == "internal":
                file.write("\t\t\tpython " + working_dir + "process.py --infile " + "\"" + working_dir + dat_row.host_id + ".trees\"\n")
            file.write("\t\t\"\"\"\n\n")

    # write out rules to generate input tree sequences for coinfections
    with open(working_dir + args.outfile, "a") as file:
        for i in inf_seq.index:
            dat_row = inf_seq.iloc[i]
            if dat_row.origin_id is not np.nan:
                if "," in dat_row.origin_id:
                    file.write("rule make_" + dat_row.host_id + "_infile:\n\tinput:\n")
                    expected_infiles = []
                    inf_sources = dat_row.origin_id.split(",")
                    for j in inf_sources:
                        expected_infiles = expected_infiles + [working_dir + j + "_to_" + dat_row.host_id + ".trees"]
                        if j == inf_sources[-1]:
                            file.write("\t\t\"" + expected_infiles[-1] + "\"\n")
                        else:
                            file.write("\t\t\"" + expected_infiles[-1] + "\",\n")
                    expected_outfile = working_dir + '-'.join(inf_sources) + "_to_" + dat_row.host_id + ".trees"
                    file.write("\toutput:\n\t\t\"" + expected_outfile + "\"\n")
                    file.write("\tshell:\n\t\t\"\"\"\n\t\t\tpython " + working_dir + "generate_input_ts.py \\\n\t\t\t\t--input ")
                    for j in expected_infiles:
                        file.write("\"" + j + "\" \\\n\t\t\t\t")
                    file.write("--host_id \"" + dat_row.host_id + "\" \\\n\t\t\t\t")
                    file.write("--output \"" + expected_outfile + "\"\n\t\t\"\"\"\n\n")

    # write out final merging step
    with open(working_dir + args.outfile, "a") as file:
        file.write("rule merge:\n\tinput:\n")
        for i in files_to_merge:
            if i == files_to_merge[-1]:
                file.write("\t\t\"" + i + "\"\n")
            else:
                file.write("\t\t\"" + i + "\",\n")
        file.write("\toutput:\n\t\t\"merged.trees\"\n")
        file.write("\tshell:\n\t\t\"\"\"\n\t\t\tpython " + working_dir + "merge_ts.py \\\n\t\t\t\t--input ")
        for i in files_to_merge:
            file.write("\"" + i + "\" \\\n\t\t\t\t")
        file.write("\t\t\t\t--output merged.trees\n\t\t\"\"\"")

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))