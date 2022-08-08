#! /usr/bin/env python
"""
Script to remove duplicated accessions (keeps accession with lowest e-value) from psi-BLAST output.
"""

import argparse
from decimal import Decimal
from sys import stdout

parser = argparse.ArgumentParser("Remove duplicated psi-BLAST matches and extract relevant information.")

parser.add_argument("-m", "--master", help="Path to the master match file to update or create. If not specified, "
                                           "uses current directory.")
parser.add_argument("psi_filename", help="Name of psi-BLAST output file to parse.")
parser.add_argument("-o", "--out", help="Name of output file. If not specified, writes to stdout.")

args = parser.parse_args()

with open(args.psi_filename, "r") as f:
    lines = [line.rstrip() for line in f]

acc_dict = {}
for line in lines:
    acc = line.split()[1].split("|")[2]
    val = acc_dict.get(acc, "fail")
    if val == "fail":
        acc_dict[acc] = (line, line.split()[13])
    else:
        if Decimal(val[1]) > Decimal(line.split()[13]):
            acc_dict[acc] = (line, line.split()[13])

outstr = ''
for a in acc_dict.keys():
    l = acc_dict[a][0].split()
    # GBFF_NAME MATCH_ACCESSION E_VALUE QUERY_NAME PID QUERY_LEN MATCH_LEN MATCH_LEN/QUERY_LEN
    res_line = l[1].split("|")[0] + "\t" + l[1].split("|")[2] + "\t" + l[13] + "\t" + l[0] + "\t" + l[2] + "\t" + \
               l[11] + "\t" + l[12] + "\t" + str(round(Decimal(l[12])/Decimal(l[11]), 2))
    outstr += res_line + '\n'

if args.out:
    with open(args.out, "w") as f:
        f.write(outstr)
else:
    stdout.write(outstr)

# Master file modification/creation
existing_lines = []
if args.master:
    master = args.master
    with open(args.master, "r") as f:
        existing_lines = [line.rstrip().split() for line in f]
else:
    # Check if master file exists in current directory if args.master not specified
    try:
        with open("all_match_names.txt", "r") as f:
            existing_lines = [line.rstrip().split() for line in f]
    except FileNotFoundError:
        pass
    master = "all_match_names.txt"

existing_matches = []
existing_genes = []
existing_evals = []
for line in existing_lines:
    existing_matches.append(line[0])
    existing_genes.append(line[1])
    existing_evals.append(line[2])

existing_dict = {}
for i in range(len(existing_matches)):
    existing_dict[existing_matches[i]] = (existing_genes[i], existing_evals[i], i)

for key in acc_dict.keys():
    line = acc_dict[key][0].split()
    e = acc_dict[key][1]
    if line[1].split("|")[2] not in existing_matches:
        existing_lines.append([line[1].split("|")[2], line[0], line[13]])
    else:
        # Protein acc. already in the master file from another Acr, keep acc with lowest eval
        new_eval = e
        old_eval = existing_dict[line[1].split("|")[2]][1]
        if Decimal(new_eval) < Decimal(old_eval):
            existing_lines[existing_dict[line[1].split("|")[2]][2]] = [line[1].split("|")[2], line[0], line[13]]

all_str = ""
for line in existing_lines:
    all_str += line[0] + "\t" + line[1] + "\t" + line[2] + "\n"

with open(master, "w") as f:
    f.write(all_str)

