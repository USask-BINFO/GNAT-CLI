#! /usr/bin/env python3

import json
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as pex
import plotly.io as pio
import base64
from io import StringIO
import pandas as pd
from collections import Counter, defaultdict

from math import log10
from numpy import nan
import re
from argparse import ArgumentParser
from os import path

# --------------------------------------------
# Helper functions
# --------------------------------------------
def fig_to_svg(fig):
    """
    Convert a matplotlib.pyplot figure to an SVG string.
    :param fig: A figure generated with matplotlib.pyplot (or another library that uses it). E.g., plt.figure()
    return: String SVG representation of the figure.
    """
    f = StringIO()
    fig.savefig(f, format='svg')
    return f.getvalue()

def sort_plots(plot_list):
    """
    Return a sorted plot_list sorted on the 'pct' field.
    :param plot_list: List of dict {'filename', 'path', 'pct'}
    return: plot_list sorted on 'pct'
    """
    return sorted(plot_list, key=lambda plot: plot['score'], reverse=True)

def format_label(label):
    """
    Removes trailing _# from the tier labels.
    :param label: The label to be formatted.
    return: label (excluding _#), tier # (or None if no tier #)
    """
    if "_" not in label:
        return label, None
    
    parts = label.rsplit("_", 1)
    if len(parts) != 2 or not parts[1].isdigit():
        return label, None
    base, tier = parts
    return base, tier

def log_data(data):
    """
    Change counter dict raw counts to log10. If 0, replaces with nan. Works with singly or doubly nested counter dicts.
    :param data: Counter dict {str product: count} or {str label: {str product: count}}
    """
    single_data = {}
    double_data = defaultdict(dict)
    for k, v in data.items():
        if isinstance(v, dict):
            # Nested
            for sub_k, sub_v in v.items():
                if sub_v != 0:
                    double_data[k][sub_k] = log10(sub_v)
                else:
                    double_data[k][sub_k] = nan
        else:
            if v != 0:
                single_data[k] = log10(v)
            else:
                single_data[k] = nan

    if len(single_data) > 1:
        return single_data
    else:
        return double_data

# --------------------------------------------
# Plot functions
# --------------------------------------------
def box_and_whisker_plots(query_data):
    """
    Generate box and whisker plots for the number of genes per GN and the number of GNs per plot. 
    :param query_data: The contents of a given 'query_name' subdict in the summary data JSON (e.g., data['query_names'][query_name_0]).
    return: SVG string box and whisker plots
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    sns.boxplot(y=query_data["gns_per_plot_data"], ax=axes[0], color="#538d6d")
    axes[0].set_title("Genomic Neighbourhoods per Plot")
    axes[0].set_ylabel("Count")
    axes[0].yaxis.get_major_locator().set_params(integer=True)

    sns.boxplot(y=query_data["genes_per_gn_data"], ax=axes[1], color="#b2e2c5")
    axes[1].set_title("Genes per Genomic Neighbourhood")
    axes[1].set_ylabel("Count")
    axes[1].yaxis.get_major_locator().set_params(integer=True)

    return fig_to_svg(fig)

def percent_labelled_pie(query_data):
    """
    Generate pie chart for the number of genes tier-labelled, PSI-BLAST-labelled, and both.
    :param query_data: The contents of a given 'query_name' subdict in the summary data JSON (e.g., data['query_names'][query_name_0]).
    return: SVG string pie chart.
    """
    tier = query_data['num_labeled']
    psiblast = query_data['num_psiblast']
    both = query_data['num_both']
    total_genes = query_data['num_genes']

    only_tier = tier - both
    only_psiblast = psiblast - both
    neither = total_genes - (only_tier + only_psiblast + both)

    labels = ["User-Defined", "Both", "PSI-BLAST", "None"]
    sizes = [only_tier, both, only_psiblast, neither]
    colours = ['#498464', '#105938', '#7db293', '#dddddd']

    # Filter out zero slices
    filtered = [(l, s, c) for l, s, c in zip(labels, sizes, colours) if s > 0]
    if not filtered:
        # Everything is zero, shouldn't happen
        print(f"Warning: When generating pie chart, all slices were filtered out as zero-size. Pie chart will be blank.")
        filtered = [('No Labels Detected', 1, '#dddddd')]

    labels, sizes, colours = zip(*filtered)

    fig, ax = plt.subplots()
    ax.pie(sizes, labels=labels, colors=colours, autopct='%.0f%%', startangle=140)

    ax.set_title("Percentage (%) of Genes Annotated")

    return fig_to_svg(fig)

def bar_chart(label_counts, title):
    """
    Return SVG string for horizontal bar chart with provided label counts.
    :param label_counts: Dict mapping label name to counts
    :param title: Title of the bar chart
    return: SVG string for the bar chart
    """
    fig, ax = plt.subplots(figsize=(8, 4))
    
    labels, values = zip(*label_counts.items())

    trim_tag_labels = []
    for label in labels:
        level = label.split("_")[-1]
        label = "_".join(label.split("_")[:-1])
        trim_tag_labels.append(label)

    ax.barh(trim_tag_labels[::-1], values[::-1], color='#105938')
    ax.set_title(title)
    ax.set_xlabel("Count")
    ax.set_ylabel("Annotation")

    return fig_to_svg(fig)

def co_localized_heatmap(data, title, max_n=None, is_int=True):
    """
    
    """
    # Check if data is one or two layers
    two_layers = False
    for key in data:
        if isinstance(data[key], dict):
            two_layers = True
        break

    if two_layers:
        # Build asymmetric product heatmap
        all_products = set()
        top_products_per_label = defaultdict(dict)

        # Get top N products for each label
        for label, product_dict in data.items():
            sorted_items = sorted(product_dict.items(), key=lambda x: x[1], reverse=True)

            if max_n:
                sorted_items = sorted_items[:max_n]
            
            for product, count in sorted_items:
                top_products_per_label[label][product] = count
                all_products.add(product)
        
        # Build DataFrame
        all_labels = sorted(top_products_per_label.keys())
        all_products = sorted(all_products)

        if not is_int:
            # Log data, floats
            matrix = pd.DataFrame(nan, index=all_labels, columns=all_products)
        else:
            matrix = pd.DataFrame(0, index=all_labels, columns=all_products, dtype=int)

        for label, product_counts in top_products_per_label.items():
            for product, count in product_counts.items():
                matrix.loc[label, product] = count

        width = max(10, 0.5 * len(all_products))
        height = max(8, 0.4 * len(all_labels))
        xlabel = "Product"
        ylabel = "Label"
        bar_label = "Number of Co-Occurrences"
        if is_int:
            bar_label += " (Raw Counts)"
        else:
            bar_label += " (Log10 Scale)"

    else:
        # Build symmetric label heatmap
        label_set = set()
        for combo in data:
            label_set.update(combo.split(';'))
        sorted_labels = sorted(label_set)

        # Init zeros matrix
        if not is_int:
            matrix = pd.DataFrame(nan, index=sorted_labels, columns=sorted_labels)
        else:
            matrix = pd.DataFrame(0, index=sorted_labels, columns=sorted_labels, dtype=int)

        # Dict for label renaming
        name_map = {}
        for label in sorted_labels:
            base, _ = format_label(label)
            name_map[label] = base
        
        # Fill in counts
        for combo, count in data.items():
            labels = combo.split(';')

            for i in labels:
                for j in labels:
                    if pd.isna(matrix.loc[i,j]):
                        matrix.loc[i,j] = 0.0
                    matrix.loc[i,j] += count
        
        # Wipe diagonal
        for label in sorted_labels:
            if not is_int:
                matrix.loc[label, label] = nan
            else:
                matrix.loc[label, label] = 0

        # Rename axes
        matrix.index = [name_map[l] for l in matrix.index]
        matrix.columns = [name_map[l] for l in matrix.columns]

        # Reorder columns and rows
        # matrix = matrix.reindex(index=sorted(matrix.index, reverse=True), columns=sorted(matrix.columns))

        width = max(10, len(matrix.columns) * .5)
        height = max(8, len(matrix.index) * 0.4)
        xlabel = "Label B"
        ylabel = "Label A"
        bar_label = "Number of Co-Occurrences"
        if is_int:
            bar_label += " (Raw Counts)"
        else:
            bar_label += " (Log10 Scale)"

    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(matrix, cmap='Greens', annot=False, ax=ax, cbar_kws={'label': bar_label})
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return fig

def sunburst_plot(label_products, is_psiblast=False):
    """
    Generates a sunburst plot mapping labels to their product annotation counts.
    :param label_products: Dict mapping str label_name to dict str product_annotation: int count
    :param is_psiblast: Bool indicating if labels are from PSI-BLAST, if True doesn't try remove trailing tier # from label
    return: HTML for the sunburst plot
    """
    rows = []

    for label, product_counter in label_products.items():
        for product, count in product_counter.items():
            if not is_psiblast:
                rows.append({"Label": format_label(label)[0], "Product": product, "Count": count})
            else:
                rows.append({"Label": label, "Product": product, "Count": count})

    df = pd.DataFrame(rows)

    fig = pex.sunburst(df, path=["Label", "Product"], values="Count", width=700, height=700)

    # Uncomment to remove labels for sunburst plot (for figure)
    # fig.update_traces(textinfo='none')

    return pio.to_html(fig, include_plotlyjs=False, full_html=False)


# --------------------------------------------
# HTML functions
# --------------------------------------------
def make_query_section(query_name, query_data, num_flanking_genes=15):
    """
    Build the HTML section for a single query. Return the HTML as a string.
    :param query_name: str query name
    :param query_data: The data from the summary JSON for the given query name. (e.g., data['query_names'][query_name_0]).
    :param num_flanking_genes: The max number of flanking genes in each GN. Used to score plots on number of labels. Default=15.
    return: HTML string for the query
    """
    section_string = ""

    # Overview
    section_string += "<h2>Snapshot: Key Metrics</h2>"
    section_string += f"""
    <h3>Totals</h3>
    <pre>
    {query_data['num_plots']} Plots
    {query_data['num_gns']} Genomic Neighbourhoods
    {query_data['num_genes']} Genes
    {query_data['num_psiblast']} Genes Annotated Using PSI-BLAST Results
    {query_data['num_labeled']} Genes Annotated with User-Defined Labels
    {query_data['num_both']} Genes Annotated Using Both User-Defined Labels and PSI-BLAST Results
    </pre>"""

    # Box and whisker plots
    section_string += box_and_whisker_plots(query_data)

    # Pie chart
    section_string += percent_labelled_pie(query_data)

    # Labelling averages
    total_genes = query_data['num_genes']
    total_gns = query_data['num_gns']

    labelled = int(query_data['num_labeled'])
    psiblast = int(query_data['num_psiblast'])
    both = int(query_data['num_both'])

    max_width = 1 + (num_flanking_genes * 2)

    total_labelled = labelled + psiblast - both

    average_labels_per_gene = total_labelled / total_genes if total_genes > 0 else 0
    average_labels_per_gn = total_labelled / total_gns if total_gns > 0 else 0

    section += f"<p><b>Percentage (%) of Genes Annotated: {average_labels_per_gene * 100:.0f}%</b></p>"
    section += f"<p><b>Percentage (%) of Genes Annotated per Genomic Neighbourhood: {average_labels_per_gn * 100:.0f}%</b></p>"

    job_score = total_labelled / (total_gns * max_width)

    # Above average plots
    above_score_plots = []

    for plot_data in query_data['plot_paths']:
        path = plot_data['path']
        total, labelled, gns = plot_data['counts']
        total, labelled, gns = int(total), int(labelled), int(gns)

        if gns == 0:
            # Avoid divide by 0 (shouldn't be possible)
            print(f"WARNING: Plot {path} appears to have 0 genes.")
            continue

        score = labelled / (gns * max_width)

        if score > job_score:
            filename = path.rsplit('/', 1)
            if len(filename) != 2:
                continue
            else:
                above_score_plots.append({"filename": filename[1], "path": path, "score": score})
    
    if len(above_score_plots) > 0:
        above_score_plots = sort_plots(above_score_plots)
        section += f"<p><b>Plots with Above-Average Numbers of Labelled Genes per Maximal Genomic Neighbourhood (Highest to Lowest):</b></p>"
        section += "<ul>"
        for file in above_score_plots:
            section += f"<li><a href=/read_file?path={file['path']}>{file['filename']}</a>"
        section += "</ul>"
    else:
        section += "<p><b>No plots found with above-average numbers of labelled genes.</b></p>"

    # PSI-BLAST Section
    if query_data['num_psiblast'] > 0:
        section += "<h2>PSI-BLAST Annotations Summary</h2>"

        # Bar Charts
        if query_data.get('psiblast_per_gn_unique'):
            section += bar_chart(query_data['psiblast_per_gn_unique'], "Genomic Neighbourhoods Containing PSI-BLAST Annotated Genes")
        
        if query_data.get('total_psiblast'):
            section += bar_chart(query_data['total_psiblast'], "Genes Annotated by PSI-BLAST")

        # Sunburst
        if query_data.get('psiblast_products'):
            title = "Summary of NCBI Product Annotations Linked to PSI-BLAST Labels"
            section += f"<h3>{title}</h3>"
            section += sunburst_plot(query_data['psiblast_products'], psiblast=True)

        # Product heatmap
        if query_data.get('co_localized_psiblast_products'):
            title = "NCBI Product Annotations Linked to PSI-BLAST Labels"
            max_products = 5
            section += f"<h3>{title}</h3>"
            section += f"<p>Retaining only the top {max_products} products from each label.</p>"
            section += co_localized_heatmap(query_data['co_localized_psiblast_products'], title, max_products=max_products, is_int=True)

            title = "Log10 " + title
            section += f"<h3>{title}</h3>"
            section += f"<p>Retaining only the top {max_products} products from each label.</p>"
            logged_data = log_data(query_data['co_localized_psiblast_products'])
            section += co_localized_heatmap(logged_data, title, max_products=max_products, is_int=False)

    # Tier Label Sections



def write_full_report(query_sections, outfile):
    """
    Writes the final HTML report by adding header info and joining provided sections. If len(query_sections) == 1, removes overall header and formats for a single query name.
    :param query_sections: The HTML for each query name.
    :param outfile: The location and name of the file to write to.
    return: None, writes to outfile
    """
    query_names = []

    for section in query_sections:
        query_name = section.split('\n')[0]
        query_name = query_name.split('>')[1]
        query_name = query_name.split('<')[0]

        query_name = query_name.split()[-1]
        query_names.append(query_name)
        

    outstr = """
    <html>
    <head>
        <title>GNAT Summary Report</title>
        <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
    </head>
    <body style="font-family: 'Open Sans', sans-serif; padding: 40px;">
        <h1>"""

    if len(query_names) == 1:
        outstr += f"GNAT Output Summary for {query_names[0]}</h1>"

    else:
        outstr += f"GNAT Output Summary for all Job Queries</h1>"

    for i, query_section in enumerate(query_sections):
        if len(query_names) != 1:
            outstr += f"<h1>Summary for {query_names[i]}</h1>"
        
        outstr += query_section
    
    outstr += """
    </body>
    </html>
    """

    file_exists_flag = False
    newfile = outfile
    # Make sure the outfile doesn't already exist, if it is, add a _(#).html to the filename
    while path.exists(newfile):
        file_exists_flag = True
        no_extension = '.'.join(newfile.split('.')[:-1])

        # Check if there's already a (#)
        num_match = re.search('\(\d+\)$', no_extension)

        if num_match:
            num_match = num_match.group()

        if num_match:
            num = int(num_match.rstrip(')').lstrip('('))
            num += 1

            no_extension = f"{no_extension.rstrip(num_match)}({num}).html"
        
        else:
            no_extension = f"{no_extension}_(0).html"
        
        newfile = no_extension

    if file_exists_flag:
        print(f"Warning: File {outfile} already exists. Writing to {newfile}.")

    with open(newfile, 'w') as f:
        f.write(outstr)

# --------------------------------------------
# Main
# --------------------------------------------
if __name__ == "__main__":    
    parser = ArgumentParser("Generate a visualization file from the GNAT summary JSON.")
    parser.add_argument("json_path", help="Path to summary JSON file. Will write to the same directory if -o | --outdir isn't provided.")
    parser.add_argument("-q", "--query", help="Query name in the JSON summary to generate visualization for. If not provided, generates one visualization file for all query names.")
    parser.add_argument("-o", "--outfile", help="Location and name of the output visualization file. If not provided, outputs 'gnat_summary.html' to the directory containing the summary JSON file.")

    args = parser.parse_args()

    # Check JSON path
    if not path.exists(args.json_path):
        raise FileNotFoundError(f"Error: Summary JSON not found at {args.json_path}. Exiting...")

    json_dir = path.dirname(args.json_path)

    # Check if --query provided
    if args.query and args.query.strip():
        provided_query = args.query.strip() # Check if it's in the summary JSON later
    else:
        provided_query = None

    # Load JSON
    with open(args.json_path) as f:
        summary_json = json.load(f)

    # Check if provided query is in the JSON
    if provided_query not in summary_json['query_names'].keys():
        print(f"Warning: Provided query name {provided_query} not found in the summary JSON. Generating visualization using all query names.")

        provided_query = None

    if provided_query:
        # 