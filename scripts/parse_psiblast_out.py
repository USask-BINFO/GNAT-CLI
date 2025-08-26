#! /usr/bin/env python3
"""
Script to remove duplicated accessions (keeps accession with lowest e-value) from psi-BLAST output.
"""

import argparse
from decimal import Decimal
from sys import stdout
import os

parser = argparse.ArgumentParser("Remove duplicated psi-BLAST matches and extract relevant information.")

# parser.add_argument("-m", "--master", help="Path to the master match file to update or create. If not specified, "
#                                            "uses current directory.")
parser.add_argument("psi_filename", help="Name of psi-BLAST output file to parse.")
parser.add_argument("-d", "--database", help="Location of the GBFF subdirectory within the database provided to AcrVis.")
# parser.add_argument("-o", "--out", help="Name of output file. If not specified, writes to stdout.")
# Use psi filename + query name + .parsed for file names, output in current directory

args = parser.parse_args()

with open(args.psi_filename, "r") as f:
    lines = [line.rstrip() for line in f]

# If multiple queries are included in the psiblast file, split them into separate lists
separated_lines = {}
for i in range(len(lines)):
    line = lines[i]
    if line.split()[0] in separated_lines.keys():
        separated_lines[line.split()[0]].append(line)
    else:
        separated_lines[line.split()[0]] = [line]

for query_acc in separated_lines.keys():
    lines = separated_lines[query_acc]

    acc_dict = {}
    # Parses for match accession and e values, and keeps the match with the lowest e value upon duplicate accession
    for line in lines:
        acc = line.split()[1].split("|")[2]
        val = acc_dict.get(acc, "fail") # Returns "fail" if accession not in dict, else returns a tuple with (match_accession, e_value)
        if val == "fail":
            acc_dict[acc] = (line, line.split()[5])
        else:
            if Decimal(val[1]) > Decimal(line.split()[5]):
                acc_dict[acc] = (line, line.split()[5])

    outstr = ''
    for a in acc_dict.keys():
        l = acc_dict[a][0].split()
        # GBFF_NAME MATCH_ACCESSION E_VALUE QUERY_ACCESSION PID QUERY_LEN MATCH_LEN MATCH_LEN/QUERY_LEN
        # Deprecated res_line = l[1].split("|")[0] + "\t" + l[1].split("|")[2] + "\t" + l[13] + "\t" + l[0] + "\t" + l[2] + "\t" + \
        #           l[11] + "\t" + l[12] + "\t" + str(round(Decimal(l[12])/Decimal(l[11]), 2))

        res_line = l[1].split("|")[0] + "\t" + l[1].split("|")[2] + "\t" + l[5] + "\t" + l[0] + "\t" + \
                l[2] + "\t" + l[3] + "\t" + l[4] + "\t" + str(round(Decimal(l[4])/Decimal(l[3]), 2))

        outstr += res_line + '\n'

    with open(os.path.basename(args.psi_filename).split('.')[0] + "_" + query_acc + ".parsed", "w") as f:
        f.write(outstr)
    # if args.out:
    #     with open(args.out, "w") as f:
    #         f.write(outstr)
    # else:
    #     stdout.write(outstr)

    # # Master file modification/creation
    # existing_lines = []
    # if args.master:
    #     master = args.master
    #     with open(args.master, "r") as f:
    #         existing_lines = [line.rstrip().split() for line in f]
    # else:
    #     # Check if master file exists in current directory if args.master not specified
    #     try:
    #         with open("all_match_names.txt", "r") as f:
    #             existing_lines = [line.rstrip().split() for line in f]
    #     except FileNotFoundError:
    #         pass
    #     master = "all_match_names.txt"

    # existing_matches = []
    # existing_genes = []
    # existing_evals = []
    # for line in existing_lines:
    #     existing_matches.append(line[0])
    #     existing_genes.append(line[1])
    #     existing_evals.append(line[2])

    # existing_dict = {}
    # for i in range(len(existing_matches)):
    #     existing_dict[existing_matches[i]] = (existing_genes[i], existing_evals[i], i)

    # for key in acc_dict.keys():
    #     line = acc_dict[key][0].split()
    #     e = acc_dict[key][1]
    #     if line[1].split("|")[2] not in existing_matches:
    #         existing_lines.append([line[1].split("|")[2], line[0], line[13]])
    #     else:
    #         # Protein acc. already in the master file from another Acr, keep acc with lowest eval
    #         new_eval = e
    #         old_eval = existing_dict[line[1].split("|")[2]][1]
    #         if Decimal(new_eval) < Decimal(old_eval):
    #             existing_lines[existing_dict[line[1].split("|")[2]][2]] = [line[1].split("|")[2], line[0], line[13]]

    # all_str = ""
    # for line in existing_lines:
    #     all_str += line[0] + "\t" + line[1] + "\t" + line[2] + "\n"

    # with open(master, "w") as f:
    #     f.write(all_str)

