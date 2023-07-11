#! /usr/bin/env python3
"""
Script to filter out redundant (identical) genomic neighbourhoods without relying on time-consuming clinker alignments.
"""

import argparse
from Bio import SeqIO
from os import path, listdir, mkdir
from shutil import move


# Get inputs
parser = argparse.ArgumentParser("Filter out identical genomic neighbourhoods based on amino acid sequence.")
parser.add_argument("gn_dir", help="The name of the directory containing the genomic neighbourhood genbank files.")
parser.add_argument("-m", "--map", help="Name of mapping file to output. Defaults to basename of directory containing"
                                        " GNs.")

args = parser.parse_args()

# Parse out amino acid sequences for all GNs into tuples (acession, seq)
seq_tuples = []
if path.isdir(args.gn_dir):
    genomes = []
    for gn in listdir(args.gn_dir):
        if '.gb' in gn[-3:]:
            genomes.append(gn)
else:
    raise ValueError("Provided argument '%s' not recognized as a directory." % args.gn_dir)

for gn in genomes:
    cds_entries = []
    gb_record = SeqIO.read(open(path.join(args.gn_dir, gn), "r"), "genbank")
    gb_feature = gb_record.features

    gn_name = gn.split("/")[-1]#.split("_")[-2]

    for seq_feature in gb_feature:
        if seq_feature.type == "CDS":
            if 'translation' in seq_feature.qualifiers:
                if seq_feature.qualifiers['translation'][0] != '':
                    cds_entries.append(seq_feature.qualifiers['translation'][0])

    aa_seq = ""
    for seq in cds_entries:
        aa_seq += seq
    rev_seq = ""
    for seq in cds_entries[::-1]:
        rev_seq += seq

    seq_tuples.append((gn_name, aa_seq, rev_seq))


# Create dictionary to find redundant GNs
removed = []
identical_dictionary = {}
for t in seq_tuples:
    check = identical_dictionary.get(t[1], "KeyError")
    if check == "KeyError":
        rev = identical_dictionary.get(t[2], "KeyError")
        if rev == "KeyError":
            identical_dictionary[t[1]] = [t[0]]
        else:
            removed.append(t[0])
            identical_dictionary[t[2]].append(t[0])
    else:
        removed.append(t[0])
        identical_dictionary[t[1]].append(t[0])


# Pick first accession as representative for each sequence in dict
rep_tuples = []

for s in identical_dictionary.keys():
    rep_tuples.append((identical_dictionary[s][0], identical_dictionary[s][1:], len(identical_dictionary[s])))

rem_dir = path.join(args.gn_dir, "redundant")
mkdir(rem_dir)

for rep_tuple in rep_tuples:
    rep = rep_tuple[0]
    rem = rep_tuple[1]
    num_removed = len(rep_tuple[1])
    # Move redundants
    for r in rem:
        move(path.join(args.gn_dir, r), path.join(rem_dir, r))
        pass

    # Rename representative to contain number represented
    new_rep_name = rep.split("_")
    new_rep_name.insert(-2, str(num_removed))
    new_rep_name = "_".join(new_rep_name)

    move(path.join(args.gn_dir, rep), path.join(args.gn_dir, new_rep_name))

# Mapping data:
m = ""
for t in rep_tuples:
    m += ">" + t[0] + "\t" + str(t[2]-1) + "\n"
    for rem in t[1]:
        m += rem + "\n"

if args.map:
    mapname = args.map
else:
    if args.gn_dir[-1] == "/" or args.gn_dir[-1] == "\\":
        args.gn_dir = args.gn_dir[:-1]
    mapname = args.gn_dir.split("/")[-1] + ".map"

with open(path.join(args.gn_dir, mapname), "w") as f:
    f.write(m)

