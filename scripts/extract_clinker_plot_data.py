#! /usr/bin/env python3

"""
Script to extract plot data (particularly label data) and append it to a summary data file.
"""

from bs4 import BeautifulSoup
import re
import json
from argparse import ArgumentParser
from os import path
from collections import defaultdict
from statistics import mean, median

def read_script_data(filepath):
    """
    Reads the file at filepath and returns the data object in the last <script> element
    :param filepath: Path to clinker plot file
    return: a dict of the data object in the last <script>, or None if something goes wrong
    """
    with open(filepath, "r") as f:
        soup = BeautifulSoup(f.read(), "html.parser")

    scripts = soup.find_all("script")
    last_script = scripts[-1].string

    data_match = re.search(r'\bdata\s*=\s*({.*?})\s*;', last_script, re.DOTALL)

    if data_match:
        json_text = data_match.group(1)
        try:
            data = json.loads(json_text)
            return data
        
        except json.JSONDecodeError as e:
            print("JSON parse error:", e)
    
    else:
        print("No match found")
    return None

def jsonify_summary(filename):
    """
    Returns provided summary JSON file as a JSON
    """
    if path.getsize(filename) == 0:
        return {}
    
    with open(filename, 'r') as f:
        try:
            return json.load(f)
        except json.JSONDecodeError:
            return {} # handle blank file

def to_defaultdict_int(d):
    return defaultdict(int, d)

def to_nested_defaultdict(d):
    return defaultdict(lambda: defaultdict(int), {k: defaultdict(int, v) for k, v in d.items()})

def restore_query_structure(data):
    """
    Casts loaded JSON query structure into the correct types (and fills any gaps with defaults)
    """
    data['gns_per_plot'] = dict(data.get('gns_per_plot', {'avg': 0, 'min': 0, 'max': 0, 'med': 0}))
    data['genes_per_gn'] = dict(data.get('genes_per_gn', {'avg': 0, 'min': 0, 'max': 0, 'med': 0}))

    data['labels_per_gn_unique'] = to_defaultdict_int(data.get('labels_per_gn_unique', {}))
    data['labels_per_gn_unique_co_localized'] = to_defaultdict_int(data.get('labels_per_gn_unique_co_localized', {}))
    data['total_labels'] = to_defaultdict_int(data.get('total_labels', {}))
    data['label_products'] = to_nested_defaultdict(data.get('label_products', {}))
    data['co_localized_label_products'] = to_nested_defaultdict(data.get('co_localized_label_products', {}))

    data['psiblast_per_gn_unique'] = to_defaultdict_int(data.get('psiblast_per_gn_unique', {}))
    data['total_psiblast'] = to_defaultdict_int(data.get('total_psiblast', {}))
    data['psiblast_products'] = to_nested_defaultdict(data.get('psiblast_products', {}))
    data['co_localized_psiblast_products'] = to_nested_defaultdict(data.get('co_localized_psiblast_products', {}))

    data['plot_paths'] = list(data.get('plot_paths', []))

    for field in ['num_plots', 'num_gns', 'num_genes', 'num_labeled', 'num_psiblast', 'num_both']:
        data[field] = int(data.get(field, 0))
    for field in ['labeled_per_plot', 'psiblast_per_plot']:
        data[field] = float(data.get(field, 0.))

    data['gns_per_plot_data'] = [int(x) for x in data.get('gns_per_plot_data', [])]
    data['genes_per_gn_data'] = [int(x) for x in data.get('genes_per_gn_data', [])]

    return data

def build_uid_map(data):
    """
    Build a map of UIDs to Gene dicts
    """
    uid_map = {}
    for cluster in data['clusters']:
        for locus in cluster['loci']:
            for gene in locus['genes']:
                uid = gene['uid']
                if uid in uid_map:
                    print(f"Warning: Gene UID ({uid}) occurs multiple times in data object. Ignorning all but first instance.")
                else:
                    uid_map[uid] = gene

    return uid_map

def build_group_map(data, gene_uid_map=None):
    """
    Build a map of Gene UIDs to a list of the Groups they are in.
    :param data: The data object from read_script_data()
    :param gene_uid_map: The UID : Gene dict map from build_uid_map(). Used to iterate through all the accessions without going through all the clusters again.
    return: A dict map of str UID : [(str group1, int tier_idx), (str group2, int tier_idx), ..., (str groupN, int tier_idx)]
    """
    if not gene_uid_map:
        gene_uid_map = build_uid_map(data)

    uid_group_map = {}

    for gene_uid in gene_uid_map:
        uid_group_map[gene_uid] = []

    for i, tier in enumerate(data['classes']):
        for group in tier['groups']:
            for gene_uid in group['genes']:
                if not re.fullmatch(r"Group \d+", group['baseLabel']):
                    uid_group_map[gene_uid].append((group['baseLabel'], i))

    return uid_group_map

def verify_gene(gene):
    """
    Make sure the gene is a dict and contains the right attributes.
    :param gene: Gene dict from the data object.
    return: True if it's a valid gene, False otherwise.
    """
    if type(gene) != dict:
        return False
    
    expected_keys = ["uid", "label", "names", "start", "end", "strand", "sequence", "translation"]

    if not all(key in gene for key in expected_keys):
        return False
    
    return True

def get_psiblast_match(gene):
    """
    If gene is a PSI-BLAST match, return the match label, else return False
    :param gene: Gene dict from the data object
    return: PSI-BLAST match label if it is a match, else False
    """
    if not verify_gene(gene):
        print(f"Warning: Provided gene is not valid. Returning False.")
        return False

    if "match_very_short" in gene["names"]:
        match_very_short = gene["names"]["match_very_short"]
        if match_very_short is not None and match_very_short.strip() != "":
            return match_very_short
    return False

def extract_plot_data(data, uid_group_map):
    """
    Extract summary data for a given clinker plot
    return: a dict containing the summary data
    """
    # Init dict to store plot data
    plot_data = {
        "plot_name": path.basename(args.clinker_plot),                  # Filename
        "num_gns": 0,                                                   # Number of GNs in the plot
        "num_genes": 0,                                                 # Total number of genes in the plot
        "num_labeled": 0,                                               # Total number of labeled genes (i.e., genes with accessions in the label file)
        "labels_per_gn_unique": defaultdict(int),                       # Counter for unique labels per GN (i.e., iterate by 1 for same label multiple times in a GN)
        "labels_per_gn_unique_co_localized": defaultdict(int),          # Counter for number of times labeled genes are found co-localized
        "total_labels": defaultdict(int),                               # Same counter as labels_per_gn_uniq but counting every duplicate label (even if in same GN)
        "label_products": defaultdict(lambda: defaultdict(int)),        # Map label to its products counters
        "co_localized_products": defaultdict(lambda: defaultdict(int)), # Map label to the products counters of all other genes in the GN
        "num_psiblast": 0,                                              # PSI-BLAST labels (even if overlapping with a label gene)
        "psiblast_labels_per_gn_unique": defaultdict(int),
        "total_psiblast_labels": defaultdict(int),
        "psiblast_products": defaultdict(lambda: defaultdict(int)),
        "psiblast_co_localized_products": defaultdict(lambda: defaultdict(int)),
        "psiblast_label_overlap": defaultdict(lambda: defaultdict(int)), # PSI-BLAST label 1 : [label1_overlap, label2_overlap, ...]
        "num_both": 0,
        "genes_per_gn_data": [],
    }

    for cluster in data['clusters']:
        plot_data['num_gns'] += 1

        gn_labels_seen = []
        psiblast_labels_seen = []

        num_genes_per_gn = 0

        products_excluding_label = defaultdict(int)
        products_excluding_psiblast = defaultdict(int)

        for locus in cluster['loci']:
            for gene in locus['genes']:
                plot_data['num_genes'] += 1
                num_genes_per_gn += 1

                # Gene's label data
                mvs = gene['names'].get('match_very_short')
                alf = gene['names'].get('all_label_file')
                groups = uid_group_map[gene['uid']]
                product = gene['names'].get('product', 'unannotated')

                if not product or product.lower().startswith("uncharacterized") or product.lower().startswith("uncharacterised") or product.lower() == "hypothetical protein":
                    product = "unannotated"

                if len(groups) > 0:
                    for i in range(len(groups)):
                        # Format the group to f'{group_base_label}_{group_tier}'
                        group_label = groups[i][0].split(" (Homologous Group ")[0]
                        group_label_level = f'{group_label}_{groups[i][1]}'

                        # Unique labels per GN counter
                        if group_label_level not in gn_labels_seen:
                            gn_labels_seen.append(group_label_level)

                            # Increment the relevant counters
                            plot_data['labels_per_gn_unique'][group_label_level] += 1

                        plot_data['total_labels'][group_label_level] += 1
                        plot_data['num_labeled'] += 1

                        plot_data['label_products'][group_label_level][product] += 1
                    
                        if mvs and mvs.strip():
                            # Has a PSI-BLAST label and a group label
                            plot_data['psiblast_label_overlap'][mvs][group_label_level] += 1
                            plot_data['num_both'] += 1
                else:
                    products_excluding_label[product] += 1

                if mvs and mvs.strip():
                    # Is a PSI-BLAST match
                    plot_data['num_psiblast'] += 1

                    # Unique PSI-BLAST labels per GN counter
                    if mvs not in psiblast_labels_seen:
                        psiblast_labels_seen.append(mvs)

                        plot_data['psiblast_labels_per_gn_unique'][mvs] += 1
                    
                    plot_data['total_psiblast_labels'][mvs] += 1
                    plot_data['psiblast_products'][mvs][product] += 1
                else:
                    products_excluding_psiblast[product] += 1

        # Now that we've populated gn_labels_seen and psiblast_labels_seen in the GN, if they're longer than 1, update the co-localization dicts
        if len(gn_labels_seen) > 1:
            gn_labels_seen.sort()
            co_str = ';'.join(gn_labels_seen)
            plot_data["labels_per_gn_unique_co_localized"][co_str] += 1

        for label in gn_labels_seen:
            for product in products_excluding_label:
                plot_data['co_localized_products'][label][product] += products_excluding_label[product]

        for psiblast_label in psiblast_labels_seen:
            for product in products_excluding_psiblast:
                plot_data['psiblast_co_localized_products'][psiblast_label][product] += products_excluding_psiblast[product]

        # Add the number of genes in the GN to the data list
        plot_data["genes_per_gn_data"].append(num_genes_per_gn)
    
    return plot_data

def to_jsonable(obj):
    # Turn the dict returned by extract_plot_data into a JSON-able dict
    if isinstance(obj, defaultdict):
        # Recursively convert defaultdict to dict
        inner = {k: to_jsonable(v) for k, v in obj.items()}
        # If it's a counter, sort descending
        if all(isinstance(v, int) for v in inner.values()):
            return dict(sorted(inner.items(), key=lambda x: x[1], reverse=True))
        return inner
    elif isinstance(obj, dict):
        inner = {k: to_jsonable(v) for k, v in obj.items()}
        # Same sorting logic
        if all(isinstance(v, int) for v in inner.values()):
            return dict(sorted(inner.items(), key=lambda x: x[1], reverse=True))
        return inner
    elif isinstance(obj, list):
        return [to_jsonable(item) for item in obj]
    else:
        return obj
    
def merge_int_dicts(dict1, dict2):
    """
    Merges defaultdicts dict 2 into dict 1, returns None
    """
    for key, value in dict2.items():
        dict1[key] += value

def merge_nested_dicts(dict1, dict2):
    """
    Merges defaultdicts dict 2 into dict 1, returns None
    """
    for outer_key, inner_dict in dict2.items():
        for inner_key, value in inner_dict.items():
            dict1[outer_key][inner_key] += value

if __name__ == "__main__":        
    parser = ArgumentParser("Extract clinker plot data, particularly label data, and append it to a summary data file.")

    parser.add_argument("clinker_plot", help="Path to clinker plot to extract data from.")

    file_group = parser.add_mutually_exclusive_group()
    file_group.add_argument("-i", "--input_file", help="Location of an existing summary data file to add to. If not provided, will make a new one in dir(clinker_plot)/summary.json or at [-o --out_location]/summary.json.")
    file_group.add_argument("-o", "--out_location", help="Directory to make output JSON. If blank, will write to dir(clinker_plot)/summary.json. [-i --input] will take precedence.")

    parser.add_argument("-j", "--job_id", help="Job ID that created the clinker plot.")
    parser.add_argument("-q", "--query_name", help="Name of the query that created the clinker plot.")

    args = parser.parse_args()

    # Verify clinker path and extension
    if not path.exists(args.clinker_plot):
        raise FileNotFoundError(f"Error: Path to clinker plot is invalid: {args.clinker_plot}")

    if not args.clinker_plot.endswith(".html"):
        raise ValueError(f"Error: provided clinker file is not an HTML file: {args.clinker_plot}")
    
    clinker_file = path.abspath(args.clinker_plot)

    # Set default location for summary output file
    summary_file = path.join(path.dirname(clinker_file), 'summary.json')

    if args.out_location:
        if not path.isdir(args.out_location):
            print(f"Error: [-o --out_location] is not a valid directory, defaulting to clinker plot directory.")
        else:
            summary_file = path.join(path.abspath(args.out_location), 'summary.json')
    
    if args.input_file:
        if not path.exists(args.input_file):
            print(f"Warning: provided input file not found, creating a new summary file in clinker plot dir.")
        else:
            summary_file = path.abspath(args.input_file)

    job_id = args.job_id if args.job_id else ""
    query_name = args.query_name if args.query_name else ""

    query_structure = {
        "num_plots": 0,
        "num_gns": 0,
        "num_genes": 0,
        "num_labeled": 0,
        "num_psiblast": 0,
        "num_both": 0,
        "gns_per_plot": {'avg': 0, 'min': 0, 'max': 0, 'med': 0},
        "genes_per_gn": {'avg': 0, 'min': 0, 'max': 0, 'med': 0},
        "labels_per_gn_unique": defaultdict(int),
        "labels_per_gn_unique_co_localized": defaultdict(int),
        "total_labels": defaultdict(int),
        "label_products": defaultdict(lambda: defaultdict(int)),
        "co_localized_label_products": defaultdict(lambda: defaultdict(int)),
        "labeled_per_plot": 0.,
        "psiblast_per_gn_unique": defaultdict(int),
        "total_psiblast": defaultdict(int),
        "psiblast_products": defaultdict(lambda: defaultdict(int)),
        "co_localized_psiblast_products": defaultdict(lambda: defaultdict(int)),
        "psiblast_per_plot": 0.,
        "plot_paths": [], 
        "gns_per_plot_data": [],
        "genes_per_gn_data": [],
    }

    main_structure = {
        "job_id": job_id,
        "query_names": {
            query_name: query_structure
        }
    }

    if path.exists(summary_file):
        # Summary file already exists, parse it into main_structure
        main_structure = jsonify_summary(summary_file)

        if main_structure == {}:
            main_structure = {
                "job_id": job_id,
                "query_names": {
                    query_name: query_structure
                }
            }

    # Get the data from the clinker plot
    clinker_data = read_script_data(clinker_file)
    uid_group_map = build_group_map(clinker_data)

    plot_data = to_jsonable(extract_plot_data(clinker_data, uid_group_map))

    # Update the data structure
    if query_name in main_structure["query_names"]:
        # Make sure the dict has the same keys as query_structure
        existing_query_data = main_structure["query_names"][query_name]
        main_structure["query_names"][query_name] = restore_query_structure(existing_query_data)

    else:
        # Add new query entry defaults
        main_structure["query_names"][query_name] = query_structure
    
    # Update query data with plot data
    query_data = main_structure["query_names"][query_name]

    # Update numbers
    query_data["num_plots"] += 1
    query_data["num_gns"] += plot_data["num_gns"]
    query_data["num_genes"] += plot_data["num_genes"]
    query_data["num_labeled"] += plot_data["num_labeled"]
    query_data["num_psiblast"] += plot_data["num_psiblast"]
    query_data["num_both"] += plot_data["num_both"]

    # Update stats
    query_data["gns_per_plot_data"].append(plot_data["num_gns"])
    query_data["genes_per_gn_data"] += plot_data["genes_per_gn_data"]

    gns_vals = query_data["gns_per_plot_data"]
    query_data["gns_per_plot"] = {
        "avg": round(mean(gns_vals), 2),
        "min": min(gns_vals),
        "max": max(gns_vals),
        "med": median(gns_vals)
    }

    genes_vals = query_data["genes_per_gn_data"]
    query_data["genes_per_gn"] = {
        "avg": round(mean(genes_vals), 2),
        "min": min(genes_vals),
        "max": max(genes_vals),
        "med": median(genes_vals)
    }

    # Update counter dicts
    merge_int_dicts(query_data["labels_per_gn_unique"], plot_data["labels_per_gn_unique"])
    merge_int_dicts(query_data["labels_per_gn_unique_co_localized"], plot_data["labels_per_gn_unique_co_localized"])
    merge_int_dicts(query_data["total_labels"], plot_data["total_labels"])
    merge_nested_dicts(query_data["label_products"], plot_data["label_products"])
    merge_nested_dicts(query_data["co_localized_label_products"], plot_data["co_localized_products"])
    merge_int_dicts(query_data["psiblast_per_gn_unique"], plot_data["psiblast_labels_per_gn_unique"])
    merge_int_dicts(query_data["total_psiblast"], plot_data["total_psiblast_labels"])
    merge_nested_dicts(query_data["psiblast_products"], plot_data["psiblast_products"])
    merge_nested_dicts(query_data["co_localized_psiblast_products"], plot_data["psiblast_co_localized_products"])

    # Add relative plot path mapped to (num_genes, num_labeled + num_psiblast - num_both, num_gns)
    if job_id and job_id.strip():
        # If the job_id was provided, try extract the relative path (e.g., the directory == job_id)
        split_path = clinker_file.split(path.sep)
        if job_id in split_path:
            cutoff = 0
            for i, part in enumerate(split_path):
                if part == job_id:
                    cutoff = i
                    break
            clinker_file = path.sep.join(split_path[cutoff:])
        
    path_data = [plot_data["num_genes"], plot_data["num_labeled"] + plot_data["num_psiblast"] - plot_data["num_both"], plot_data["num_gns"]]
    query_data["plot_paths"].append({"path": clinker_file, "counts": path_data})

    # Update summary values that depend on number of plots
    query_data['labeled_per_plot'] = round(query_data['num_labeled'] / query_data['num_plots'], 2)
    query_data['psiblast_per_plot'] = round(query_data['num_psiblast'] / query_data['num_plots'], 2)

        
    # Convert to JSON and write
    with open(summary_file, "w") as f:
        json.dump(to_jsonable(main_structure), f, indent=2)
