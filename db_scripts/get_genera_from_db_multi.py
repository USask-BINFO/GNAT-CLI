#!/usr/bin/env python3

import itertools
import multiprocessing as mp
import os
from tqdm import tqdm
import gzip
from Bio import SeqIO

def process_file(filename, out_path):
    genera_in_file = {}
    num_records_in_file = 0

    process_name = mp.current_process().name
    worker_id = int(process_name.split('-')[1])

    output_file = os.path.join(out_path, f"good_{worker_id}.txt")
    unknown_file = os.path.join(out_path, f"unknown_{worker_id}.txt")
    multi_file = os.path.join(out_path, f"multi_{worker_id}.txt")

    with gzip.open(filename, 'rt') as gbff_file:
        for record in SeqIO.parse(gbff_file, "genbank"):
            organism_name = record.annotations.get('organism', 'Unknown unknown')
            genus = organism_name.split()[0]

            if genus != "Unknown":
                if genus in genera_in_file:
                    genera_in_file[genus] += 1
                else:
                    genera_in_file[genus] = 1

            num_records_in_file += 1

    # output_file = f"{output_prefix}_{os.path.basename(filename)}.txt"
    if len(genera_in_file) > 1:
        genus_list = '\t'.join(list(genera_in_file.keys()))
        with open(multi_file, 'a') as multi_genus_file:
            multi_genus_file.write(f"{os.path.basename(filename)}\t{genus_list}\n")
    elif len(genera_in_file) == 0:
        with open(unknown_file, 'a') as no_genus_file:
            no_genus_file.write(f"{os.path.basename(filename)}\n")
    else:
        genus = list(genera_in_file.keys())[0]
        with open(output_file, 'a') as out_file:
            out_file.write(f"{os.path.basename(filename)}\t{genus}\n")

    return num_records_in_file

def process_files_in_batches(directory, batch_size, out_path):
    filenames = (os.path.join(dp, f) for dp, _, fn in os.walk(directory) for f in fn if f.endswith('.gbff.gz'))

    num_processes = min(mp.cpu_count(), 8)
    total_records = 0

    # Set up multiprocessing pool
    with mp.Pool(num_processes) as pool:
        pbar = tqdm(desc="Files processed", unit="file")

        def update_progress(file_count):
            pbar.update(file_count)

        results = []
        while True:
            batch = list(itertools.islice(filenames, batch_size))
            if not batch:
                break

            for filename in batch:
                results.append(pool.apply_async(process_file, args=(filename, out_path)))
            
            # results = pool.starmap(process_file, [(filename, output_prefix) for filename in batch])
            update_progress(len(batch))

        pbar.close()

    # # Aggregate results
    # with open(f"{output_prefix}_combined.txt", 'w') as combined_file:
    #     for filename in os.listdir('.'):
    #         if filename.startswith(output_prefix) and filename.endswith('.txt'):
    #             with open(filename, 'r') as file:
    #                 combined_file.write(file.read())
    total_records = sum(result.get() for result in results)

    return total_records


# Path to the directory containing GBFF.gz files
db_type = "staphylococcus_aureus"

gbff_directory = f'/binfo-nas1/data/Acr/databases/{db_type}_db/{db_type}_genomes/gbff'

# Output file for the gbff name-to-genus mapping
output_directory = f'/binfo-nas1/data/Acr/databases/{db_type}_db/{db_type}_genomes'

num_records = process_files_in_batches(gbff_directory, 100, output_directory)
print(num_records)