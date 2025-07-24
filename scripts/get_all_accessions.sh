#! /usr/bin/env bash

# Script to get all accessions from the genbank files in a given job directory
# Usage: get_all_accessions.sh <job_directory>

# Error catching script
source "$(dirname "$0")/utils.sh"

# Help message
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <job_directory>"
    echo "Extracts all accessions from the genbank files in the specified job directory."
    echo "Writes <job_directory>/<job_directory>.accessions."
    exit 0
fi

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <job_directory>"
    exit 1
fi

# Make sure the job directory exists
if [ ! -d "$1" ]; then
    echo "Job directory $1 does not exist."
    exit 1
fi

# Make sure the job directory is an absolute path
curr_dir=$(pwd)
cd "$1" || exit
job_dir=$(pwd)
cd "$curr_dir" || exit

# Make sure the output file is empty
> "$job_dir/$(basename "$job_dir").accessions"

# Loop through all the genbank files in the job directory
for subdir in "$job_dir"/$(basename "$job_dir")_*; do
    # Make sure subdir is a directory
    if [ ! -d "$subdir" ]; then
        echo "Subdirectory $subdir does not exist."
        continue
    fi
    
    # Make our subdir .accessions file
    > "$subdir/$(basename "$subdir").accessions"

    # Loop through each cluster in all_non_redundants
    for clusterdir in "$subdir"/all_non_redundants/*; do
        # Make sure cluster is a directory
        if [ ! -d "$clusterdir" ]; then
            echo "Cluster directory $clusterdir does not exist."
            continue
        fi

        # Loop through each genbank file in clusterdir/genbank/
        for gbfile in "$clusterdir"/genbank_files/*.gb; do
            # Make sure gbfile is a file
            if [ ! -f "$gbfile" ]; then
                echo "Genbank file $gbfile does not exist."
                continue
            fi

            # Extract the accession number from the genbank file
            run_catch extract_accessions.py "$gbfile" "$subdir/$(basename "$subdir").accessions"
        done
    done

    sort -u "$subdir/$(basename "$subdir").accessions" -o "$subdir/$(basename "$subdir").accessions"

    cat "$subdir/$(basename "$subdir").accessions" >> "$job_dir/$(basename "$job_dir").accessions"

done

# Sort the final accessions file
sort -u "$job_dir/$(basename "$job_dir").accessions" -o "$job_dir/$(basename "$job_dir").accessions"

