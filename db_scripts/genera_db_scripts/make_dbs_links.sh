#!/usr/bin/env bash

#  - In:
#    - path/to/genus/db/
#      - The path to one of the genus database root directories created by extract_phylo_stats.py
#      - E.g., path/to/databases/genera/{genus}_db/
#    - path/to/source/gbff/
#      - Path to main gbff directory (should be the main bacterial one)
#      - E.g., path/to/databases/bacterial_db/bacterial_genomes/gbff/
#  - Out:
#    - From the provided genus database root directory:
    # - {genus}_db/
    # - {genus}_genomes/
    #     - gbff/
    #     - symlinks to all the gbffs in the file
    #     - {genus}_genomes.txt
    # - {genus}db/
    #     - makeblastdb is run on the symlink genomes in gbff/
    #     - split_pal.py is used to generate the all_pal_names.txt file to make the db queriable with AcrVis

# Parse input
while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--genus_db)
        genus_db="$2"
        shift
        shift
        ;;
        -s|--source_gbff)
        source_gbff="$2"
        shift
        shift
        ;;
        -h|--help)
        echo "    Script to create a blast database for a genus."
        echo "Required Arguments:"
        echo "  [-g][--genus_db]          Path to a particular genus database root directory."
        echo "  [-s][--source_gbff]       Path to the main source gbff directory."
        exit 0
        ;;
        *)
        echo "Invalid argument: $1"
        exit 1
        ;;
    esac
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if required arguments are provided
if [[ -z $genus_db ]] || [[ -z $source_gbff ]]; then
    echo "Required arguments not provided, use [-h][--help] for usage help."
    exit 1
fi

# Check if the provided genus database root directory exists
if [[ ! -d $genus_db ]]; then
    echo "Genus database root directory does not exist."
    exit 1
fi

# Check if the provided source gbff directory exists
if [[ ! -d $source_gbff ]]; then
    echo "Source gbff directory does not exist."
    exit 1
fi

# Check if there's a / at the end of either $genus_db or $source_gbff, if there is, remove it
genus_db="${genus_db%/}"
source_gbff="${source_gbff%/}"

# Create the genus_genomes directory
genus_genomes="$genus_db/$(basename $genus_db | cut -d "_" -f 1)_genomes"
if [[ ! -d $genus_genomes ]]; then
    mkdir $genus_genomes
fi

# Create the genusdb directory
genusdb="$genus_db/$(basename $genus_db | cut -d "_" -f 1)db"
if [[ ! -d $genusdb ]]; then
    mkdir $genusdb
fi

# Create the gbff directory in genus_genomes
genus_gbff="$genus_genomes/gbff"
if [[ ! -d $genus_gbff ]]; then
    mkdir $genus_gbff
fi

# Create symlinks to all the gbff files listed in the {genus}_genomes.txt file
genus_genomes_txt="$genus_db/$(basename $genus_db | cut -d "_" -f 1 )_genomes.txt"
while IFS= read -r line; do
    ln -s $source_gbff/$line $genus_gbff/$line
done < $genus_genomes_txt

# Move the {genus}_genomes.txt file to genus_genomes
mv $genus_genomes_txt $genus_genomes

# Move the all_{genus}_proteins.fasta file to {genus}_genomes/proteins/
genus_proteins="$genus_genomes/proteins"
if [[ ! -d $genus_proteins ]]; then
    mkdir $genus_proteins
fi

fasta_file="$genus_db/all_$( basename $genus_db | cut -d "_" -f 1 )_proteins.fasta"
# Move the fasta file to the proteins directory
if [[ -f $fasta_file ]]; then
    mv "$fasta_file" "$genus_proteins/"
    fasta_file="$genus_proteins/$(basename "$fasta_file")"
else
    echo "Fasta file $fasta_file does not exist."
    exit 1
fi

# Make blastdb in $genusdb off of the fasta file
MM_YYYY=$(date +%m_%Y)
curr_dir=$(pwd)
genus="$(basename $genus_db | cut -d "_" -f 1)"

cd "$genusdb" || exit
makeblastdb -in "$fasta_file" -dbtype prot -title "${genus}"db_$MM_YYYY -out "${genus}"db_$MM_YYYY

# Run split_pal.py to generate the all_pal_names.txt file
cd "$curr_dir" || exit
cd "$genusdb" || exit

cat *.pal | "$script_dir/../split_pal.py"

cd "$curr_dir" || exit
