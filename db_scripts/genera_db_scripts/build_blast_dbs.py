#!/usr/bin/env python3
# Script to run make_dbs_links.sh on each genus subdirectory

import argparse
from tqdm import tqdm
import subprocess
import os
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build a blastdb for each genus subdirectory in the provided directory, and symlinks to the appropriate gbff files from the main database.")
    parser.add_argument("genus_directory", help="Path to the directory containing the genera subdirectories.")
    parser.add_argument("source_gbff", help="Path to the main bacterial gbff directory.")

    args = parser.parse_args()

    genus_dir = args.genus_directory
    gbff_dir = args.source_gbff

    # Verify the provided directories
    if not os.path.isdir(genus_dir):
        print(f"Error: The genus directory '{genus_dir}' does not exist. Exiting.")
        sys.exit(1)
    
    if not os.path.isdir(gbff_dir):
        print(f"Error: The source GBFF directory '{gbff_dir}' does not exist. Exiting.")
        sys.exit(1)
    
    # Get path to make_dbs_links.sh
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_dir, "make_dbs_links.sh")

    # Get the genera in the directory
    subdirs = [dir for dir in os.listdir(genus_dir) if os.path.isdir(os.path.join(genus_dir, dir))]

    if not subdirs:
        print(f"No subdirectories found in '{genus_dir}'. Exiting.")
        sys.exit(1)

    with tqdm(desc="Making genera databases", unit="genera", total=len(subdirs)) as pbar:
        for genus in subdirs:
            genus_path = os.path.join(genus_dir, genus)

            # Call the make_dbs_links.sh script for each genus
            try:
                result = subprocess.run([script_path, '-g', genus_path, '-s', gbff_dir], check=True)
                print(f"Processing genus: {genus}")
                print(result.stdout.decode())
            except subprocess.CalledProcessError as e:
                print(f"Error processing genus '{genus}':", e.stderr.decode())

            pbar.update()