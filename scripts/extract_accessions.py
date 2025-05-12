#! /usr/bin/env python3

"""
Script to extract accessions from a given GenBank file and write them to a specified output file.
"""

import argparse
from os import path
from Bio import SeqIO

def extract_accessions(input_file, output_file):
    """
    Extract accessions from a GenBank file and write them to a specified output file.
    
    Args:
        input_file (str): Path to the input GenBank file.
        output_file (str): Path to the output file where accessions will be written.
    """
    if not path.exists(input_file):
        raise FileNotFoundError(f"Input file {input_file} does not exist.")
    
    with open(output_file, "a") as out_f:
        for record in SeqIO.parse(input_file, "genbank"):
            for feature in record.features:
                if 'protein_id' in feature.qualifiers:
                    # Extract the protein ID
                    protein_id = feature.qualifiers['protein_id'][0]
                    # Write the protein ID to the output file
                    out_f.write(f"{protein_id}\n")

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Extract accessions from a GenBank file and write them to a specified output file.")
    parser.add_argument("input_file", help="Path to the input GenBank file.")
    parser.add_argument("output_file", help="Path to the output file where accessions will be written.")
    
    args = parser.parse_args()
    
    # Extract accessions
    extract_accessions(args.input_file, args.output_file)