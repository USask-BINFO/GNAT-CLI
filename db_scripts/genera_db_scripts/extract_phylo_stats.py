#!/usr/bin/env python3

"""
Script to extract the genera stats from an assembly_summary.txt file
"""

import os
import argparse
import tqdm

def parse_assembly_summary(assembly_summary_file):
    """
    Reads the assembly summary file, extracting the organism name and the gbff name from the ftp link.
    param assembly_summary_file: path to the assembly summary file from NCBI for the gbff directory
    return: list of tuples (organism, gbff_name)
    """
    parsed = []

    if os.path.exists(assembly_summary_file):
        with open(assembly_summary_file, "r") as f:
            lines = f.readlines()[2:]

            with tqdm.tqdm(desc="Reading assembly summary file", unit="lines", total=len(lines)) as pbar:

                for line in lines:
                    line = line.strip().split("\t")
                    organism = line[7]
                    gbff = line[19].split('/')[-1] + "_genomic.gbff.gz"
                    parsed.append((organism, gbff))
                    pbar.update()

    return parsed

def get_db_lines(parsed, gbff_directory):
    """
    For each gbff filename in the parsed assembly summary file, check if the gbff file is in the database's gbff directory
    param parsed: list of tuples (organism, gbff_name) from parse_assembly_summary()
    param gbff_directory: path to the database directory containing the gbff files
    return: N/A. Modifies parsed by removing entries not in the db
    """
    with tqdm.tqdm(desc="Removing lines not in the database", unit="lines", total=len(parsed)) as pbar:
        for i in range(len(parsed) - 1, -1, -1):
            _, gbff = parsed[i]
            if not os.path.exists(os.path.join(gbff_directory, gbff)):
                parsed.pop(i)
            pbar.update()

def extract_genera_stats(parsed):
    """
    Extracts the genera stats from the parsed assembly summary file
    param parsed: list of tuples (organism, gbff_name) from parse_assembly_summary()
    return: dictionary where keys are genera and values are the number of genomes for that genus
    """
    with tqdm.tqdm(desc="Extracting genera stats", unit="lines", total=len(parsed)) as pbar:
        genera_stats = {}

        for organism, gbff in parsed:
            genus = organism.split(" ")[0]
            if genus == "uncultured":
                genus = organism.split(" ")[1]
            if genus in genera_stats:
                genera_stats[genus].append(gbff)
            else:
                genera_stats[genus] = [gbff]
            
            pbar.update()

        sorted_genera_stats = sorted(genera_stats.items(), key=lambda item: len(item[1]), reverse=True)
        del genera_stats

        return sorted_genera_stats
    
def extract_gbff_fasta(genera_gbff_dict, fasta_path):
    """
    For each entry in the dict, creates a new fasta file. Extracts the fasta entries from the main fasta file (fasta_path). Loop through the fasta file once to keep N low. Rewriting the dict to a new dict with gbff name 
    as key and value as the genus name to identify which list it belongs in.
    """

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract genera stats from an assembly_summary.txt file")
    parser.add_argument("assembly_summary_file", help="Path to the assembly_summary.txt file")
    parser.add_argument("gbff_directory", help="Path to the directory containing the gbff files")
    parser.add_argument("output_directory", help="Path to the directory to write the sub-databases to")

    args = parser.parse_args()
    assembly_summary_file = args.assembly_summary_file
    gbff_dir = args.gbff_directory

    parsed = parse_assembly_summary(assembly_summary_file)

    get_db_lines(parsed, gbff_dir)

    genera_stats = extract_genera_stats(parsed)
    del parsed

    # Get the top 20 most frequently-occurring genera in the database
    top_20_genera = []
    while len(top_20_genera) < 20:
        top_20_genera.append(genera_stats.pop(0)) if genera_stats[0][0] != "bacterium" else genera_stats.pop(0)
    del genera_stats

    genera_subdir = os.path.join(args.output_directory, "genera")
    os.makedirs(genera_subdir, exist_ok=True)

    # Write each one to a separate .txt file (each line is a gbff filename)
    for line in top_20_genera:
        genus = line[0]
        gbffs = line[1]
        genus_outdir = os.path.join(genera_subdir, f"{genus}_db")
        os.makedirs(genus_outdir, exist_ok=True)
        with open(os.path.join(genus_outdir, f"{genus}_genomes.txt"), "w") as f:
            for gbff in gbffs:
                f.write(gbff + "\n")
    