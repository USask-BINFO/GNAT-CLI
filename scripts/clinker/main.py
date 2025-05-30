#!/usr/bin/env python3

"""
Command line interface

Cameron Gilchrist
"""

import argparse
import logging
import csv

from pathlib import Path
from collections import defaultdict
from typing import TextIO, Dict, List
from os import path as osPath

from decimal import Decimal

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from LabelFile import LabelFile

#from clinker import __version__, align
#from clinker.plot import plot_clusters, plot_data
#from clinker.classes import find_files, parse_files

import align
from plot import plot_clusters, plot_data
from classes import find_files, parse_files, label_obj_container

__version__ = "0.0.23"

logging.basicConfig(
    level=logging.ERROR,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S"
)
LOG = logging.getLogger(__name__)

def parse_gnat_label_file(label_file_location=None):
    """
    GNAT/AcrVis function to parse the gene labelling file.
    Old format: tsv mapping [0] accession to [1] label with [2] e-value of the match
    New format: csv mapping [0] accession to [1] label with [2] e-value of the match, and additional n*2 items (n=number of class labels), where, e.g., [3] class label name (e.g., Abi2) and [4] is the hex colour to make genes with that label. These additional labels will be used to create n new gene Groups with new colouring and legends. Note this format also has a header line (beginning with #) containing the names of each header (including the name of each class label group (e.g., [3] defense_systems))
    
    :param label_file_location: Location of the label file to parse.

    return: LabelObject (from label_cls.py) containing the parsed label file data. # Rough updating to use LabelFile object instead
    """
    if label_file_location == "" or label_file_location is None:
        LOG.warning("No label file provided. GNAT labelling will not be enabled.")
        return None
    
    try:
        lf_object = LabelFile(label_file_location)
        # lf_object.from_file(label_file_location)
        return lf_object
    

    except FileNotFoundError:
        LOG.warning(f"Label file not found at location: {label_file_location}. GNAT labelling will not be enabled.")
        return None
    except Exception as e:
        LOG.error(f"Unexpected error parsing label file: {e.__class__.__name__}: {e}")
        raise e

def parse_range(string):
    """Extracts the scaffold name, start and end of a scaffold range string.

    Expects the format scaffold:start-stop (e.g. scaffold_1:100-3000).

    Args:
        string (str): Scaffold range string
    Returns:
        scaffold (str): Scaffold name
        start (int): Range start
        end (int): Range end
    Raises:
        ValueError: If range is in invalid format (can't split by '-' or ':')
        TypeError: Range values are not integers
    """
    try:
        scaffold, coordinates = string.split(":")
        start, end = coordinates.split("-")
    except ValueError as e:
        raise ValueError("Expected format scaffold:start-stop (e.g. scaf_1:100-3000)") from e
    if not start.isdigit() or not end.isdigit():
        raise TypeError("Expected range values to be type int")
    return scaffold, int(start), int(end)


def parse_ranges(strings):
    ranges = {}
    for string in strings:
        try:
            scaffold, start, end = parse_range(string)
        except (ValueError, TypeError):
            LOG.exception("Failed to read range, please check it's in the right format")
            raise
        ranges[scaffold] = (start, end)
    return ranges


def parse_gene_functions(fp: TextIO) -> Dict[str, List[str]]:
    """Parses gene functions from a table.

    Gene        Function
    GENE_001    Cytochrome P450
    GENE_002    Methyltransferase
    ...
    """
    functions = defaultdict(list)
    for gene, function in csv.reader(fp):
        functions[function].append(gene)
    return functions


def clinker(
    files,
    session=None,
    identity=0.3,
    delimiter=None,
    decimals=2,
    plot=None,
    output=None,
    force=False,
    no_align=False,
    hide_link_headers=False,
    hide_alignment_headers=False,
    use_file_order=False,
    json_indent=None,
    jobs=None,
    ranges=None,
    matrix_out=None,
    gene_functions=None,
    gnat_label=False
):
    """Entry point for running the script."""
    LOG.info("Starting clinker")

    load_session = session and Path(session).exists()

    # Allow no files, so that user can generate a blank clinker web app
    # and load in previously saved figure data
    if not files:
        LOG.info("No files provided!")
        if plot:
            LOG.info("Opening empty clinker web app...")
            plot_data(
                dict(clusters=[], links=[], groups=[]),
                output=None if plot is True else plot
            )
        return

    # Parse range strings, if any specified
    if ranges:
        ranges = parse_ranges(ranges)

    # Parse any gene functions for grouping, if specified
    if gene_functions:
        gene_functions = parse_gene_functions(gene_functions)

    if load_session:
        LOG.info("Loading session from: %s", session)
        with open(session) as fp:
            globaligner = align.Globaligner.from_json(fp)
        if files:
            paths = find_files(files)
            if not paths:
                LOG.error("No files found")
                raise SystemExit
            LOG.info("Parsing files:")
            clusters = parse_files(paths, ranges=ranges)

            LOG.info("Adding clusters to loaded session and aligning")
            globaligner.add_clusters(*clusters)
            globaligner.align_stored_clusters(cutoff=identity, jobs=jobs)
            groups = globaligner.build_gene_groups_return(functions=gene_functions, use_links=True, hidden=False)
            globaligner.classes.append({"name": "Default", "groups": groups})
            globaligner.groups = globaligner.classes[0]["groups"]
            # globaligner.build_gene_groups(functions=gene_functions)
            load_session = False
    else:
        # Parse files, generate objects
        paths = find_files(files)
        if not paths:
            LOG.error("No files found")
            raise SystemExit
        LOG.info("Parsing files:")
        clusters = parse_files(paths, ranges=ranges, gnat_label=gnat_label)

        # Align all clusters
        if no_align:
            globaligner = align.Globaligner()
            globaligner.add_clusters(*clusters)
            groups = globaligner.build_gene_groups_return(functions=gene_functions, use_links=True, hidden=False)
            globaligner.classes.append({"name": "Default", "groups": groups})
            globaligner.groups = globaligner.classes[0]["groups"]
            # globaligner.build_gene_groups(functions=gene_functions)
        elif len(clusters) == 1:
            globaligner = align.align_clusters(clusters[0], jobs=1)
        else:
            LOG.info("Starting cluster alignments")
            globaligner = align.align_clusters(*clusters, cutoff=identity, jobs=jobs)
            groups = globaligner.build_gene_groups_return(functions=gene_functions, use_links=True, hidden=False)
            globaligner.classes.append({"name": "Default", "groups": groups})
            globaligner.groups = globaligner.classes[0]["groups"]
            # globaligner.build_gene_groups(functions=gene_functions)

    if gnat_label and globaligner:
        label_obj = label_obj_container[0]

        # Add default groups to globalaligner.classes
        if "Default" not in [c["name"] for c in globaligner.classes]:
            globaligner.classes.append({"name": "Default", "groups": globaligner.groups})

        # If the label file has more than 3 columns, we need to build the class (tier) labels
        if label_obj is not None and len(label_obj) > 0 and len(label_obj.header.split(",")) > 3:
            n_label_dict = [{} for _ in range(len(label_obj.header.split(',')) - 3)] # A label dict for each class column {label: [accessions]}
            # n_colour_map = [{} for _ in range((len(label_obj.header.split(',')) - 3) // 2)] # A colour map for each class column {label: [colour(s)]}

            # Go through all the accessions in the plot, and build the label dicts and colour maps (i.e., parse out the label file data for the genes present in the plot)
            for uid, gene in globaligner._genes.items():
                if 'protein_id' in gene.names:
                    # For each of these accessions, get the label from the label file
                    record = label_obj.get_record(gene.names['protein_id'])

                    if record is not None:
                        # Split the record into the labels and the classes
                        split = record.rest.split(",")

                        # Get the label for each class
                        for i in range(2, len(split)):
                            label = split[i]
                            idx = i - 2

                            if label in ["", "None", None]:
                                continue
                            else:
                                # For the label dict, label mapped to accession
                                if label in n_label_dict[idx]:
                                    n_label_dict[idx][label].append(gene.names['protein_id'])
                                else:
                                    n_label_dict[idx][label] = [gene.names['protein_id']]

                                # if colour.strip() in ["", "None", None]:
                                #     colour = ""

                                # if label in n_colour_map[idx]:
                                #     n_colour_map[idx][label].append(colour)
                                # else:
                                #     n_colour_map[idx][label] = [colour]
                
                else:
                    # Assign a placeholder protein_id to the gene
                    # This is a hack to get the re-colouring to work with gene groups with no accessions
                    gene.names['protein_id'] = f"{uid}_placeholder"

            # Done extracting label file data for the genes in the plot

            # Build the default gene functions (based on homology) using the existing groups
            default_groups = globaligner.groups

            n_fn = []
            n_colour_map = []
            n_new_labels = []
            n_default_ordering = []
            n_psiblast_groups = []
            for _ in range(len(label_obj.header.split(',')) - 3):
                default_fn = {}
                default_colours = {}
                default_labels = {}
                default_order = []
                psiblast_groups = []
                for i in range(len(default_groups)):
                    group = default_groups[i]
                    # Get the gene accessions for the group (From the UIDs)
                    default_fn[group.label] = []
                    is_match = False
                    for gene_uid in group.genes:
                        # Get the gene accession from the UID
                        gene_names = globaligner._genes[gene_uid].names
                        if 'protein_id' in gene_names:
                            # If the gene has a protein ID, use that
                            gene_acc = globaligner._genes[gene_uid].names['protein_id']
                        else:
                            gene_acc = gene_uid # Make sure this works TODO

                        if 'psi_blast_match' in gene_names and gene_names['psi_blast_match'] not in [None, 'None'] and gene_names['psi_blast_match'].strip() != "":
                            is_match = True
                        default_fn[group.label].append(gene_acc)
                    psiblast_groups.append(is_match)
                    # default_fn[group.label] = group.genes
                    # print(f"Group {group.label} has the following attributes (and values): {group.__dict__}")
                    default_colours[group.label] = [group.colour]
                    default_labels[group.label] = group.label
                    default_order.append(group.label)
            
                n_fn.append(default_fn)
                n_colour_map.append(default_colours) # Looks like all of these colours are [None], since the groups haven't been built yet?
                # print(f"-----------------------------------------------------------\nDefault Colours: {default_colours}\n-----------------------------------------------------------")
                n_new_labels.append(default_labels)
                n_default_ordering.append(default_order)
                n_psiblast_groups.append(psiblast_groups)

            # Now we have all the data we need to build the new groups (homology extended and non-homologous subclass)
            # 1 - If a gene from the label file is in a homology group, re-name the homology group to the label (e.g., Group 3 --> "Abi2 (Homology Group 3)")
            # 2 - If another gene with a different label from the label file is in the same homology group, prepend it to the label (e.g., Abi2 (Homology Group 3) --> Abi3 | Abi2 (Homology Group 3))
            # 3 - If a gene from the label file isn't in a homology group, create a new group for it in the fn and colour it black
            #   - If the label is already in the fn, add the gene to it (e.g., if 3 already happened)
            # 4 - If genes from the label file share the same label, but are in different homology groups, keep them separate (e.g., Abi2 (Homology Group 3) and Abi2 (Homology Group 4))
            # 5 - Likewise, if two genes from the label file share the same label, but are in different homology groups, and one of the homology groups has another gene with a different label, keep them separate and only prepend the label for the latter one (e.g., Abi2 (Homology Group 3) and Abi3 (Homology Group 4) --> Abi2 | Abi3 (Homology Group 4))

            keys = list(n_fn[0].keys())
            for i in range(len(n_label_dict)):
                label_dict = n_label_dict[i]
                colour_map = n_colour_map[i]
                names = n_new_labels[i]
                default_fn = n_fn[i]
                psiblast_groups = n_psiblast_groups[i]

                acc_to_labels = {v: k for k, vals in label_dict.items() for v in vals}
                acc_to_group = {v: k for k, vals in default_fn.items() for v in vals}

                # Get a list of all the accessions in the label dict
                label_gene_accs = []
                for acc_list in label_dict.values():
                    label_gene_accs += acc_list

                homology_genes_accs = []
                for acc_list in default_fn.values():
                    homology_genes_accs += acc_list

                # Get the intersection of the two lists
                acc_intersection = set(label_gene_accs).intersection(set(homology_genes_accs))
                acc_label_rest = set(label_gene_accs) - set(homology_genes_accs)
                # acc_group_rest = set(homology_genes_accs) - set(label_gene_accs)

                for gene_acc in acc_intersection:
                    # Get the label for the gene
                    label = acc_to_labels[gene_acc]
                    # Get the group for the gene
                    group = acc_to_group[gene_acc]

                    if names[group] == group:
                        # Group name hasn't been updated, overwrite it
                        names[group] = f"{label} (Homologous {group})"
                        n_default_ordering[i][n_default_ordering[i].index(group)] = names[group]
                    else:
                        # Group name has already been updated, check if it has the label already
                        if f"{label} (Homologous " not in names[group] and f"{label} | " not in names[group]:
                            # Group name has been updated with a different label, add this one to it
                            idx = n_default_ordering[i].index(names[group])
                            names[group] = f"{label} | " + names[group]
                            n_default_ordering[i][idx] = names[group]

                
                remove = [] # List of entries to remove from fn (renamed)
                # Now we need to update the fn and colour map with the new labels
                for group in keys:
                    if names[group] == group:
                        # Group name hasn't been updated, colour it grey
                        colour_map[group] = ["#808080"]
                    else:
                        default_fn[names[group]] = default_fn[group]
                        remove.append(group)
                        # Update the colour map with the new label
                        colour_map[names[group]] = colour_map[group]
                
                # Remove the old group names from the fn and colour map
                for group in remove:
                    del default_fn[group]
                    del colour_map[group]

                # Add the rest of the label genes to the fn
                for gene_acc in acc_label_rest:
                    # Get the label for the gene
                    label = acc_to_labels[gene_acc]
                    # Check if the label is already in the fn
                    if label in default_fn:
                        # Label has already been added to the fn, add the gene to it
                        default_fn[label].append(gene_acc)
                    else:
                        # Label hasn't been added to the fn, add it
                        default_fn[label] = [gene_acc]
                        # Add the colour to the colour map (black)
                        colour_map[label] = ["#000000"]

                        # Add the label to the order list
                        n_default_ordering[i].append(label)

                if len(n_default_ordering[i]) > len(psiblast_groups):
                    # Check if the new gene_accs contain PSI-BLAST matches
                    while len(psiblast_groups) < len(n_default_ordering[i]):
                        new_order_label = n_default_ordering[i][len(psiblast_groups)]

                        for gene_acc in default_fn[new_order_label]:
                            # Check if the new group has PSI-BLAST matches, add relevant tag to PSI-BLAST match list
                            record = label_obj.get_record(gene_acc)

                            if record is None or record.rest.split(',')[0].strip() == "":
                                psiblast_groups.append(False)
                            else:
                                psiblast_groups.append(True)
                        
                    n_psiblast_groups[i] = psiblast_groups

                # Homology genes should already all be in there
                
            # Go through the fns/colour_maps and build the groups (in default order)
            for i in range(len(n_fn)):
                fn = n_fn[i]
                colour_map = n_colour_map[i]
                default_order = n_default_ordering[i]
                labeled_groups = []
                psiblast_groups = n_psiblast_groups[i]
                names = n_new_labels[i]

                for j in range(len(default_order)):
                    label = default_order[j]
                    labeled_group = False

                    if label in names:
                        if names[label] != label:
                            # Group name updated, set labeledGroup to true
                            labeled_group = True

                        label = names[label]
                    else:
                        # New group
                        labeled_group = True

                    labeled_groups.append(labeled_group)

                    default_order[j] = label

                    group_name = label

                    if colour_map[group_name] and len(colour_map[group_name]) > 1:
                        LOG.warning(f"Label {group_name} has more than one colour assigned: {colour_map[group_name]}. Using first one.")

                    if not colour_map[group_name]:
                        # New group, pass (it's black)
                        colour_map[group_name] = ["#000000"]

                    colour_map[group_name] = colour_map[group_name][0]

                # Get the class name from the header
                class_name = label_obj.header.split(",")[i + 3]

                # print(labeled_groups)

                new_group = globaligner.build_gene_groups_return(functions=fn, colour_map=colour_map, use_links=False, level=i+1, order=default_order, labeled_groups=labeled_groups, psiblast_groups=psiblast_groups)
                globaligner.classes.append({"name": class_name, "groups": new_group})

    if globaligner.alignments:
        LOG.info("Generating results summary...")
        summary = globaligner.format(
            delimiter=delimiter,
            decimals=decimals,
            link_headers=not hide_link_headers,
            alignment_headers=not hide_alignment_headers,
        )
        if output:
            if (output and Path(output).exists() and not force):
                print(summary)
                LOG.warn("File %s already exists but --force was not specified", output)
            else:
                LOG.info("Writing alignments to: %s", output)
                with open(output, "w") as fp:
                    fp.write(summary)
        else:
            # Not printing summaries to prevent cluttering the terminal if LOG level isn't set to INFO
            if LOG.level == logging.INFO:
                print(summary)
        if matrix_out:
            LOG.info("Writing synteny matrix to: %s", matrix_out)
            matrix = globaligner.format_matrix(normalise=True, as_distance=True)
            with open(matrix_out, "w") as fp:
                fp.write(matrix)
    else:
        LOG.info("No alignments were generated")

    if session and not load_session:
        LOG.info("Saving session to: %s", session)
        with open(session, "w") as fp:
            globaligner.to_json(fp, indent=json_indent)

    # Generate the SVG
    if plot:
        LOG.info("Building clustermap.js visualisation")
        if isinstance(plot, str):
            LOG.info("Writing to: %s", plot)
        plot_clusters(
            globaligner,
            output=None if plot is True else plot,
            use_file_order=use_file_order,
        )

    LOG.info("Done!")
    return globaligner


def get_parser():
    """Creates an ArgumentParser object."""
    parser = argparse.ArgumentParser(
        "clinker",
        description="clinker: Automatic creation of publication-ready"
        " gene cluster comparison figures.\n\n"
        "clinker generates gene cluster comparison figures from GenBank files."
        " It performs pairwise local or global alignments between every sequence"
        " in every unique pair of clusters and generates interactive, to-scale comparison figures"
        " using the clustermap.js library.",
        epilog="Example usage\n-------------\n"
        "Align clusters, plot results and print scores to screen:\n"
        "  $ clinker files/*.gbk\n\n"
        "Only save gene-gene links when identity is over 50%:\n"
        "  $ clinker files/*.gbk -i 0.5\n\n"
        "Save an alignment session for later:\n"
        "  $ clinker files/*.gbk -s session.json\n\n"
        "Save alignments to file, in comma-delimited format, with 4 decimal places:\n"
        "  $ clinker files/*.gbk -o alignments.csv -dl \",\" -dc 4\n\n"
        "Generate visualisation:\n"
        "  $ clinker files/*.gbk -p\n\n"
        "Save visualisation as a static HTML document:\n"
        "  $ clinker files/*.gbk -p plot.html\n\n"
        "Cameron Gilchrist, 2020",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--version", action='version', version=f"clinker v{__version__}")

    inputs = parser.add_argument_group("Input options")
    inputs.add_argument("files", help="Gene cluster GenBank files", nargs="*")
    inputs.add_argument(
        "-r",
        "--ranges",
        help="Scaffold extraction ranges. If a range is specified, only features within"
        " the range will be extracted from the scaffold. Ranges should be formatted"
        " like: scaffold:start-end (e.g. scaffold_1:15000-40000)",
        nargs="+",
    )
    inputs.add_argument(
        "-gf",
        "--gene_functions",
        help="2-column CSV file containing gene functions, used to build gene groups"
        " from same function instead of sequence similarity (e.g. GENE_001,PKS-NRPS).",
        type=argparse.FileType("r")
    )
    inputs.add_argument("-lf", "--labelling_file", default="", help="Location of file to be used for match labelling generated by mod_label_file.py. New argument added for acrvis.")

    alignment = parser.add_argument_group("Alignment options")
    alignment.add_argument(
        "-na",
        "--no_align",
        help="Do not align clusters",
        action="store_true",
    )
    alignment.add_argument(
        "-i",
        "--identity",
        help="Minimum alignment sequence identity [default: 0.3]",
        type=float,
        default=0.3
    )
    alignment.add_argument(
        "-j",
        "--jobs",
        help="Number of alignments to run in parallel (0 to use the number of CPUs) [default: 0]",
        type=int,
        default=0,
    )

    output = parser.add_argument_group("Output options")
    output.add_argument("-s", "--session", help="Path to clinker session")
    output.add_argument("-ji", "--json_indent", type=int, help="Number of spaces to indent JSON [default: none]")
    output.add_argument("-f", "--force", help="Overwrite previous output file", action="store_true")
    output.add_argument("-o", "--output", help="Save alignments to file")
    output.add_argument(
        "-p",
        "--plot",
        nargs="?",
        const=True,
        default=False,
        help="Plot cluster alignments using clustermap.js. If a path is given,"
        " clinker will generate a portable HTML file at that path. Otherwise,"
        " the plot will be served dynamically using Python's HTTP server."
    )
    output.add_argument("-dl", "--delimiter", help="Character to delimit output by [default: human readable]")
    output.add_argument("-dc", "--decimals", help="Number of decimal places in output [default: 2]", default=2)
    output.add_argument(
        "-hl",
        "--hide_link_headers",
        help="Hide alignment column headers",
        action="store_true",
    )
    output.add_argument(
        "-ha",
        "--hide_aln_headers",
        help="Hide alignment cluster name headers",
        action="store_true",
    )
    output.add_argument("-mo", "--matrix_out", help="Save cluster similarity matrix to file")

    viz = parser.add_argument_group("Visualisation options")
    viz.add_argument(
        "-ufo",
        "--use_file_order",
        action="store_true",
        help="Display clusters in order of input files"
    )

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    
    # labelfilelocation[0] = args.labelling_file
    lf = args.labelling_file

    # label_dict = parse_gnat_label_file(lf)
    label_obj = parse_gnat_label_file(lf)
    if label_obj is None:
        gnat_labelling = False
    else:
        gnat_labelling = True
        label_obj_container[0] = label_obj

    clinker(
        args.files,
        session=args.session,
        json_indent=args.json_indent,
        identity=args.identity,
        delimiter=args.delimiter,
        decimals=args.decimals,
        plot=args.plot,
        output=args.output,
        force=args.force,
        no_align=args.no_align,
        hide_link_headers=args.hide_link_headers,
        hide_alignment_headers=args.hide_aln_headers,
        use_file_order=args.use_file_order,
        jobs=args.jobs if args.jobs > 0 else None,
        ranges=args.ranges,
        matrix_out=args.matrix_out,
        gene_functions=args.gene_functions,
        gnat_label=gnat_labelling
    )

if __name__ == "__main__":
    main()


