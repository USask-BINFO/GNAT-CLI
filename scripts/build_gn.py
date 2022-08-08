#! /usr/bin/env python
"""
Program used to isolate the gene clusters of the psiblast matches of Acr/Aca proteins. Can add a command-line arg to
modify the n param in find_cluster()
"""

from Bio import SeqIO
import gzip
import argparse
from os import path

parser = argparse.ArgumentParser("Create genomic neighbourhoods from tab-delimited psi-BLAST output files.")
parser.add_argument("parsed_out", type=str, help="Tab delimited psi-BLAST output file.")
parser.add_argument("-n", help="Number of genes upstream/downstream to be added to create genomic neighbourhoods.",
                    type=int)
parser.add_argument("-d", "--database", help="Path to database folder containing .gbff files (normally called gbff/",
                    required=True)
parser.add_argument("-o", "--out", help="Desired output location (data/query_name/)", required=True)

args = parser.parse_args()

filename = path.basename(args.parsed_out)

with open(args.parsed_out, "r") as f:
    lines = [line for line in f]

if args.n is not None:
    num_genomes = args.n
else:
    num_genomes = 20

# holds the column values from a modified tab delimted psiblast output for each match against a query
genomes = []
accessions = []
e_values = []

# go through each line (match) of the formatted psiblast output and add their values into the corresponding lists
# each element position regards to a specific match across all lists (aka element [0] is 1 match for 3 lists, just
# holding different data)
for line in lines:
    if line.split()[1] in accessions:
        match_already_present_index = accessions.index(line.split()[1])

        if e_values[match_already_present_index] > float(line.split()[2]):
            genomes[match_already_present_index] = line.split()[0]
            accessions[match_already_present_index] = line.split()[1]
            e_values[match_already_present_index] = float(line.split()[2])

    else:
        genomes.append(line.split()[0])
        accessions.append(line.split()[1])
        e_values.append(float(line.split()[2]))


def find_cluster(genomename, accessionpos, n):
    """
    Function to find the clusters of a psiblast match and return the cluster a certain number up and downstream
    from the match found
    Param: genomename = Genome to look through
    Param: accessionpos = accession of the psiblast match
    :param n: number of genes up/downstream to find
    Return: Genbank Record containing a certain number
    """
    # the path for opening the file must point to the directory holding genomes (.gbff.gz files)
    with gzip.open(path.join(args.database, genomename), "rt") as openedgenome:
        recordlist = SeqIO.parse(openedgenome, "gb")
        # use the parser to parse through the genbank record and the features of each record
        for record in recordlist:
            # two lists hold the features n upstream and downstream of the match in record
            n_upstream_plusmatch = []
            n_downstream = []
            matchreached = False

            for feature in record.features:
                if feature.type == "CDS":
                    # if the match is not found in record, and the list is not at n capacity, then add to upstream list
                    if not matchreached:
                        if len(n_upstream_plusmatch) < n + 1:
                            n_upstream_plusmatch.append(feature)

                        # else remove the oldest feature from the list and add the new one (maintain n upstream from
                        # match)
                        else:
                            n_upstream_plusmatch.pop(0)
                            n_upstream_plusmatch.append(feature)

                    # match is found, begin adding to downstream list
                    else:
                        if len(n_downstream) < n:
                            n_downstream.append(feature)
                            if feature == record.features[-1]:
                                print("Reached end of record")  # no more downstream genes in record
                                n_upstream_plusmatch.extend(n_downstream)  # combine the two lists into one
                                gene_cluster = record_match_foundin[(n_upstream_plusmatch[0].location.start): (
                                    n_upstream_plusmatch[-1].location.end)]
                                gene_cluster.annotations["taxonomy"] = record_match_foundin.annotations["taxonomy"]
                                gene_cluster.annotations["organism"] = record_match_foundin.annotations["organism"]
                                return gene_cluster

                        # n in downstream list, create genbank cluster and return
                        else:
                            n_upstream_plusmatch.extend(n_downstream)
                            gene_cluster = record_match_foundin[(n_upstream_plusmatch[0].location.start): (
                                n_upstream_plusmatch[-1].location.end)]
                            gene_cluster.annotations["taxonomy"] = record_match_foundin.annotations["taxonomy"]
                            gene_cluster.annotations["organism"] = record_match_foundin.annotations["organism"]
                            return gene_cluster

                    # if a protein_id exists and it matches the accession number of the match
                    if 'protein_id' in feature.qualifiers and feature.qualifiers['protein_id'][0] == accessionpos:
                        record_match_foundin = record
                        matchreached = True  # match has been found

                        # if the match is the last feature in the record, return the cluster
                        if feature == record.features[-1]:
                            n_upstream_plusmatch.extend(n_downstream)
                            gene_cluster = record_match_foundin[(n_upstream_plusmatch[0].location.start): (
                                n_upstream_plusmatch[-1].location.end)]
                            gene_cluster.annotations["taxonomy"] = record_match_foundin.annotations["taxonomy"]
                            gene_cluster.annotations["organism"] = record_match_foundin.annotations["organism"]
                            print("Match in last Feature")
                            return gene_cluster

                elif feature.type != 'CDS' and feature == record.features[-1] and matchreached:
                    n_upstream_plusmatch.extend(n_downstream)  # combine the two lists into one
                    gene_cluster = record_match_foundin[
                                   (n_upstream_plusmatch[0].location.start): (n_upstream_plusmatch[-1].location.end)]
                    gene_cluster.annotations["taxonomy"] = record_match_foundin.annotations["taxonomy"]
                    gene_cluster.annotations["organism"] = record_match_foundin.annotations["organism"]
                    print("Match in last CDS")
                    return gene_cluster


# Main
print("Creating Clusters for " + args.parsed_out + " as Input file...")

# for each genome in the tab delimited file of matches, find the cluster create genbank file
for each in range(len(genomes)):
    Gene_cluster = find_cluster(genomes[each], accessions[each], num_genomes)

    # Remove special characters from name (add more as needed)
    genbankfile_matchname = Gene_cluster.description.replace(" ", "_").replace(":", "").replace(",", "").replace("|", "_")
    #print(Gene_cluster)
    # print(genbankfile_matchname)
    print("%i GNs generated. Total=%s" % (each+1, str(len(genomes))))

    # Write GN to args.out
    SeqIO.write(Gene_cluster, path.join(args.out, genbankfile_matchname.replace('/', '_') + '_' + str(accessions[each])
                                        + '_' + str(lines[each].split()[2]) + ".gb"), "genbank")

print("------NEXT QUERY-------")

