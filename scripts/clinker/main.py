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

from label_cls import LabelObject

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

    return: LabelObject (from label_cls.py) containing the parsed label file data.
    """
    if label_file_location == "" or label_file_location is None:
        LOG.warning("No label file provided. GNAT labelling will not be enabled.")
        return None
    
    try:
        lf_object = LabelObject()
        lf_object.from_file(label_file_location)
        return lf_object
        # with open(label_file_location, "r") as f:
        #     accessions = {}
        #     is_new = False

        #     for i, line in enumerate(f):

        #         # Skip header line if file is in new format
        #         if i == 0 and line.startswith("#"):
        #             is_new = True
        #             classes = [] # Init a list to store the n class types

        #             for j in range(3, len(line.split(",")), 2):
        #                 classes.append({"label_name": line.split(",")[j], "label_map": {}})
        #             continue

        #         elif i == 0 and not line.startswith("#"):
        #             # Make sure it's in the old format if the first line isn't a header
        #             if not len(line.split()) == 3:
        #                 LOG.warning("Label file does not start with a header line (as expected for a new format file), but also does not contain the expected 3 tab-delimited columns expected of the old format. GNAT labelling will not be enabled.")
        #                 return None


        #         if is_new:
        #             # Parse lines in new format
        #             line = line.strip().split(",")
        #             replace_dupe = False # Flag to replace duplicate accessions' class labels with the lower e-value's ones

        #             accession, label, e_value = line[0], line[1], line[2]

        #             if accessions.get(accession, None) is not None:
        #                 LOG.warning(f"Duplicate accession found in label file: {accession}. The one with the lower E-value will be used.")
        #                 if Decimal(e_value) < Decimal(accessions[accession][1]):
        #                     replace_dupe = True
        #                     old_classes = accessions[accession][2]
        #                     accessions[accession] = [label, e_value, line[3:]]
        #                 else:
        #                     continue
        #             else:
        #                 accessions[accession] = [label, e_value, line[3:]]

        #             if replace_dupe:
        #                 # Remove accession from old class label
        #                 for j in range(0, len(old_classes), 2):
        #                     class_idx = j // 2
        #                     class_label = old_classes[j]

        #                     while accession in classes[class_idx]["label_map"][class_label]["accessions"]:
        #                         classes[class_idx]["label_map"][class_label]["accessions"].remove(accession)

        #             for j in range(3, len(line), 2):
        #                 class_idx = (j - 3) // 2

        #                 class_label = line[j]
        #                 class_colour = line[j+1]

        #                 if class_label not in classes[class_idx]["label_map"]:
        #                     classes[class_idx]["label_map"][class_label] = {"accessions": [accession], "colour": class_colour}

        #                 else:
        #                     # Class label already exists, add accession to it (and make sure the colour is the same)
        #                     if class_colour != classes[class_idx]["label_map"][class_label]["colour"]:
        #                         LOG.warning(f"Class label '{class_label}' has multiple colours in label file (new colour: {class_colour}). Using the first one in the label file: {classes[class_idx]['label_map'][class_label]['colour']}.")
                            
        #                     # Double check that the accession isn't already in the class label (if it is, something is wrong)
        #                     if accession in classes[class_idx]["label_map"][class_label]['accessions']:
        #                         LOG.error(f"Accession {accession} is already in class label {class_label}. This should not be possible.\nTruncated class label dict so far: {str(classes)[:50]}...\nTruncated offending entry: {str(classes[class_idx]['label_map'][class_label]['accessions'])[:100]}...")
        #                         return None
                            
        #                     classes[class_idx]["label_map"][class_label]['accessions'].append(accession)

        #         else:
        #             # Parse lines in old format
        #             accession, label, e_value = line.split()
                    
        #             if accessions.get(accession, None) is not None:
        #                 LOG.warning(f"Duplicate accession found in label file: {accession}. The one with the lower E-value will be used.")
        #                 if Decimal(e_value) < Decimal(accessions[accession][1]):
        #                     accessions[accession] = [label, e_value]
        #                 else:
        #                     continue
        #             else:
        #                 accessions[accession] = [label, e_value]

        #     if is_new:
        #         return {"is_new": is_new, "labels": accessions, "class_labels": classes}
        #     else:
        #         return {"is_new": is_new, "labels": accessions, "class_labels": None}
    

    except FileNotFoundError:
        LOG.warning(f"Label file not found at location: {label_file_location}. GNAT labelling will not be enabled.")
        return None
    except Exception as e:
        LOG.error(f"Unexpected error parsing label file: {e.__class__.__name__}: {e}")
        return None

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
            globaligner.build_gene_groups(functions=gene_functions)
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
            globaligner.build_gene_groups(functions=gene_functions)
        elif len(clusters) == 1:
            globaligner = align.align_clusters(clusters[0], jobs=1)
        else:
            LOG.info("Starting cluster alignments")
            globaligner = align.align_clusters(*clusters, cutoff=identity, jobs=jobs)
            globaligner.build_gene_groups(functions=gene_functions)

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

    # Label Dict structure:
    # { "is_new": bool, 
    #   "labels": {accession: [label, e-value, [split_line[3:] (if is_new)]]}, 
    #   "class_labels": [{"label_name": str, 
    #                     "label_map": {str label: {"accessions": [accession, ...], 
    #                                               "colour": str}}}, 
    #                    ...] or None if not is_new }

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


