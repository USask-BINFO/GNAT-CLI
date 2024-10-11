#! /usr/bin/env python3
"""
Program used to isolate the gene clusters of the psiblast matches of Acr/Aca proteins. Can add a command-line arg to
modify the n param in find_cluster()
"""

from Bio import SeqIO
import gzip
import argparse
from os import path

parser = argparse.ArgumentParser("Create genomic neighbourhoods from tab-delimited psi-BLAST output files.")
parser.add_argument("parsed_out", type=str, help="Tab delimited parsed psi-BLAST output file.")
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
    num_genomes = 15

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

def get_num_features(record):
    """
    Get the number of CDS features in a GenBank record.
    :param record: GenBank record to count features in.
    :return: Number of CDS features in the record.
    """
    return len([feature for feature in record.features if feature.type == "CDS"])

def get_feature_index(record, accession):
    """
    Get the index of a given protein accession (feature) in the GenBank record.
    :param record: GenBank record to search for the feature in.
    :param accession: Accession number of the protein in the feature to find.
    :return: Index of the feature in the record.
    """
    for i, feature in enumerate(record.features):
        if 'protein_id' in feature.qualifiers and feature.qualifiers['protein_id'][0] == accession:
            return i
    return -1

def get_flanks(record, split_index):
    """
    Split a record feature list into 3 sublists around the specified index: [before][index][after]
    :param record: GenBank record containing the features to separate
    :param split_index: The index to split flanks off from, usually the match index.
    :return: [list, of, features, before, idx], [idx], [features, after, idx]
    """
    features = [feature for feature in record.features if feature.type == "CDS"]

    before, match, after = features[:split_index], features[split_index], features[split_index+1:]

    return before, match, after

def trim_flank(flank, edge_gene, side, additional_genes_allowed=3):
    """
    Trim the flank list to the edge gene and up to 3 additional overlapping genes.
    :param flank: List of features to trim. (either 'after' or reversed 'before' from get_flanks)
    :param edge_gene: Feature to trim to.
    :param side: +1 or -1 to indicate which side of the match gene the features are on. +1 is after the match gene, -1 is before.
    :param additional_genes_allowed: Number of additional overlapping genes to allow. Default 3.
    :return: Trimmed list of features.
    """
    # Walk through until we hit the edge gene
    edge_position = 0
    for i in range(len(flank)):
        feature = flank[i]
        if feature == edge_gene:
            edge_position = i
            break
    
    # Now that we're at the edge gene, allow up to 3 more genes, only if they overlap with the edge gene
    end_position = edge_position
    for i in range(edge_position + 1, len(flank)):
        if are_features_overlapping(edge_gene, flank[i], side):
            additional_genes_allowed -= 1
            end_position = i
            if additional_genes_allowed == 0:
                break
        else:
            break

    # Trim the list to the edge gene and the additional overlapping genes
    return flank[:end_position + 1]

def are_features_overlapping(feature1, feature2, side):
    """
    Check if two CDS features are overlapping.
    :param feature1: First feature to check (feature closer to the match gene).
    :param feature2: Second feature to check (feature further from the match gene).
    :param side: +1 or -1 to indicate which side of the match gene the features are on. +1 is after the match gene, -1 is before.
    :return: True if the features are overlapping, False otherwise.
    """
    def start_end_to_earliest_latest(start, end):
        """
        Converts the start and end coordinates to earliest and latest to make it orientation-independent.
        """
        return min(start, end), max(start, end)
    
    # Get the start and end coordinates of the features (multiplied by side to make sure conversion to ealiest & latest works properly)
    f1_start, f1_end = feature1.location.start * side, feature1.location.end * side
    f2_start, f2_end = feature2.location.start * side, feature2.location.end * side

    # Convert to earliest, latest (only need f1 early and f2 late)
    _, f1_latest = start_end_to_earliest_latest(f1_start, f1_end)
    f2_earliest, _ = start_end_to_earliest_latest(f2_start, f2_end)

    # If f2 earliest is less than f1 latest, the genes are overlapping
    return f2_earliest < f1_latest

def does_feature_wrap_origin(feature, original_record):
    """
    Checks if a gene wraps around the genome. (i.e. is broken over the start/end of the genome since bacterial genomes are circular).
    We'll eventually remove these genes from the GNs, since Clinker doesn't display them well.
    :param feature: Feature to check for wrapping.
    :param original_record: Original record the feature came from (e.g., full gbff file).
    :return: True if the feature wraps around the genome, False otherwise.
    """
    genome_length = len(original_record.seq)

    # Get the start and end positions of the feature
    start, end = feature.location.start, feature.location.end
    print(start, end)

    if start == 0 and end == genome_length:
        return True
    return False

def build_gn_from_features(original_record, feature_list, match_accession):
    """
    Build a new GenBank record from a list of features extracted by find_cluster().
    Builds the new record from the loci cooridinates of the earliest and latest features,
    so there is a rare edge case where there are more genes than 2n+1 in the GN if either
    of the genes n away from the match have another gene overlapping them, resulting in 2n+2 
    genes.

    :param original_record: Original SeqRecord to use for sequence data and annotations.
    :param feature_list: List of SeqFeatures to include in the new record.
    :param match_accession: Accession of the PSI-BLAST match the GN is being built around.
    :return: New SeqRecord with the specified gene features. 
    """
    if not feature_list:
        raise ValueError("No features to build Genomic Neighbourhood from.")
    
    # Get start and end positions of the new record
    start, end = feature_list[0].location.start, feature_list[-1].location.end

    # No longer using this loop so that the GN cuts off at the start/end of the first/last gene (respectively). It will still include genes fully overlapped by the edge gene, but avoids the rare edge case where a gene within the GN consists of two exons super far apart and results in a GN that is bascially the entire genome.
    # Loop to ensure there aren't overlapping genes that go further than the first/last in the list
    for feature in feature_list:
        if feature.location.start < start:
            start = feature.location.start
        if feature.location.end > end:
            end = feature.location.end

    new_record = original_record[start:end]

    # Check the number of features in the new GN
    if get_num_features(new_record) > (len(feature_list) + 6):
        # If there are more features than expected (plus a little extra for potentially overlapping genes), separate the flanks and start trimming them
        edge_gene_before, edge_gene_after = feature_list[0], feature_list[-1]
        new_record_match_idx = get_feature_index(new_record, match_accession)

        # Get the flanks of the GN
        before, match, after = get_flanks(new_record, new_record_match_idx)

        # Check the lengths of each flank
        if len(before) > num_genomes + 3:
            # Flip the list to get the 'earliest' genes (genes closest to match)
            before = before[::-1]
            before = trim_flank(before, edge_gene_before, -1)
            before = before[::-1]
        
        if len(after) > num_genomes + 3:
            after = trim_flank(after, edge_gene_after, 1)

        # Check the first and last gene features, if either of them overlap to the other end of the genome, remove them
        first_gene = before[0]
        last_gene = after[-1]

        if does_feature_wrap_origin(first_gene, original_record):
            before = before[1:]

        if does_feature_wrap_origin(last_gene, original_record):
            after = after[:-2]

        # Combine the flanks and match gene
        new_features = before + [match] + after

        # Create a new record explicitly with the new features
        start = new_features[0].location.start
        end = new_features[-1].location.end
        new_record = original_record[start:end]
        new_record.features = new_features

    new_record.annotations = original_record.annotations

    return new_record

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
                                # print("Reached end of record")  # no more downstream genes in record
                                n_upstream_plusmatch.extend(n_downstream)  # combine the two lists into one
                                return build_gn_from_features(record_match_foundin, n_upstream_plusmatch, accessionpos)

                        # n in downstream list, create genbank cluster and return
                        else:
                            n_upstream_plusmatch.extend(n_downstream)
                            return build_gn_from_features(record_match_foundin, n_upstream_plusmatch, accessionpos)

                    # if a protein_id exists and it matches the accession number of the match
                    if 'protein_id' in feature.qualifiers and feature.qualifiers['protein_id'][0] == accessionpos:
                        record_match_foundin = record
                        matchreached = True  # match has been found

                        # if the match is the last feature in the record, return the cluster
                        if feature == record.features[-1]:
                            n_upstream_plusmatch.extend(n_downstream)
                            # print("Match in last Feature")
                            return build_gn_from_features(record_match_foundin, n_upstream_plusmatch, accessionpos)

                elif feature.type != 'CDS' and feature == record.features[-1] and matchreached:
                    n_upstream_plusmatch.extend(n_downstream)  # combine the two lists into one
                    # print("Match in last CDS")
                    return build_gn_from_features(record_match_foundin, n_upstream_plusmatch, accessionpos)


# Main
# print("Creating Clusters for " + args.parsed_out + " as Input file...")

# for each genome in the tab delimited file of matches, find the cluster create genbank file
for each in range(len(genomes)):
    Gene_cluster = find_cluster(genomes[each], accessions[each], num_genomes)

    # Remove special characters from name (add more as needed)
    genbankfile_matchname = Gene_cluster.description.replace(" ", "_").replace(":", "").replace(",", "").replace("|", "_")
    #print(Gene_cluster)
    # print(genbankfile_matchname)
    # print("%i GNs generated. Total=%s" % (each+1, str(len(genomes))))

    if len(genbankfile_matchname) > 150:
        genbankfile_matchname = genbankfile_matchname[:137] + "_trimmed_name"

    # Write GN to args.out
    SeqIO.write(Gene_cluster, path.join(args.out, genbankfile_matchname.replace('/', '_') + '_' + str(accessions[each])
                                        + '_' + str(lines[each].split()[2]) + ".gb"), "genbank")

