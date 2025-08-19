#!/usr/bin/env python3 

from argparse import ArgumentParser
import os
import gzip
from tqdm import tqdm
from glob import glob

def parse_tax_info(gbff_file):
    """
    Parses out the taxonmic information from the SOURCE and ORGANISM fields of the provided GenBank file with a regular expression.
    :param gbff_file: Path to the gbff.gz file from the gbff/ dir in the database.
    return: (str accession, str source, str organism, str CSV taxonomy, str errors)
    """
    accession = ""
    source = ""
    organism = ""
    taxonomy = []
    errors = ""

    with gzip.open(gbff_file, "rt") as f:
        for line in f:
            l = line.strip().replace(',', ';')
            if not accession and l.startswith("ACCESSION"):
                accession = l.split()[1].strip()
            elif not source and l.startswith("SOURCE"):
                source = l.split("SOURCE", 1)[-1].strip()
            elif not organism and l.startswith("ORGANISM"):
                organism = l.split("ORGANISM", 1)[-1].strip()
            elif organism:
                if line.startswith(' ' * 12):
                    tax_line = l.rstrip(';.').split(';')
                    tax_line = [item.strip() for item in tax_line]

                    if tax_line == [l]:
                        # Organism is multi-line
                        organism += l
                    else:
                        taxonomy += tax_line

                else:
                    break
                # tax_line = line.rstrip(";.").split(";")
                # tax_line = [item.strip() for item in tax_line]
                # if tax_line == [line] or "REFERENCE" in tax_line: # REFERENCE usally comes after ORGANISM section, so check for that as a fallback if it ends in a period or something
                #     break
                # taxonomy += tax_line

    if not (accession and source and organism and taxonomy):
        errors = f"GBFF File '{os.path.basename(gbff_file)}' has at least one missing field."

    if taxonomy:
        taxonomy = ';'.join(taxonomy)
    else:
        taxonomy = ""

    return (accession, source, organism, taxonomy, errors)


if __name__ == "__main__":
    parser = ArgumentParser("Script to build a phylogeny CSV from an old database (before adding it to the pipeline and making it automatic while building db to save computations).")
    parser.add_argument("gbff_dir", help="Path to gbff/ directory of the database to extract phylogeny for.")
    parser.add_argument("-o", "--out_file", help="File to write CSV data to. If not provided defaults to './gbff_taxonomy.csv'", default="./gbff_taxonomy.csv")
    parser.add_argument("-b", "--batch_size", help="Number of files to process before appending to outfile. Defaults to 10,000.", default=10000, type=int)

    args = parser.parse_args()

    if not os.path.isdir(args.gbff_dir):
        raise FileNotFoundError(f"Provided path to gbff directory '{args.gbff_dir}' is not a directory.")
    
    if os.path.exists(args.out_file):
        print(f"Warning: Output file '{args.out_file}' exists. It will be overwritten.")

    with open(args.out_file, 'w') as of:
        of.write("Accession,Source,Organism,Taxonomy,Errors\n")
        results = []
        for gbff_file in tqdm(glob(f"{os.path.join(args.gbff_dir, '*.gbff.gz')}"), desc="Parsing taxonomic data from GBFFs"):
            if results and len(results) % args.batch_size == 0:
                of.writelines(results)
                results = []
            
            data = parse_tax_info(gbff_file)

            try:
                results.append(','.join(data) + '\n')
            except TypeError as t:
                print(data)
                raise t
        
        if len(results) > 0:
            of.writelines(results)