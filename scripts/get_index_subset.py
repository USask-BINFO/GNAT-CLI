#! /usr/bin/env python3

"""
Script to generate the subset label file from a list of accessions and an indexed label file.
"""

from Index import LabelIndex
import argparse
from os import path

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Generate a subset label file from a list of accessions and an indexed label file.")
    parser.add_argument("index_file", help="Path to the index file (.idx).")
    parser.add_argument("label_file", help="Path to the label file corresponding to the index file.")
    parser.add_argument("accession_file", help="Path to the file containing the list of accessions.")
    parser.add_argument("output_file", help="Path to the output file where the subset label file will be written.")

    args = parser.parse_args()

    # Check the existence of the index file, label file, and accession file
    if not path.exists(args.index_file):
        raise FileNotFoundError(f"Index file {args.index_file} does not exist.")
    if not path.exists(args.label_file):
        raise FileNotFoundError(f"Label file {args.label_file} does not exist.")
    if not path.exists(args.accession_file):
        raise FileNotFoundError(f"Accession file {args.accession_file} does not exist.")
    
    # Verify the index file and label file match
    if not args.label_file == args.index_file[:-4]:
        raise ValueError(f"Label and index file names do not match: {args.label_file} and {args.index_file}.")
    
    # Create an instance of LabelIndex
    index = LabelIndex(args.label_file)

    # Get the header line from the label file
    with open(args.label_file, "r") as label_f:
        header = label_f.readline().strip()

    # Open the output file for writing
    with open(args.output_file, "w") as out_f:
        # Write the header line to the output file
        out_f.write(f"{header}\n")
        # Open the accession file for reading
        with open(args.accession_file, "r") as acc_f:
            for line in acc_f:
                acc = line.strip()

                result = index.search(acc)

                if result is not None:
                    # Write the result to the output file
                    out_f.write(f"{result}\n")
