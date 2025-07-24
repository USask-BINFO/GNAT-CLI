#! /usr/bin/env python3

import json
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as pex
import plotly.io as pio
from io import StringIO
import pandas as pd
from collections import Counter, defaultdict

from math import log10
from numpy import nan
import re
from argparse import ArgumentParser
from os import path

# Imports for linear regression
from sklearn.linear_model import LinearRegression
from scipy.stats import norm
from math import sqrt

# --------------------------------------------
# Helper functions
# --------------------------------------------
def fig_to_svg(fig, width_pct='45'):
    """
    Convert a matplotlib.pyplot figure to an SVG string.
    :param fig: A figure generated with matplotlib.pyplot (or another library that uses it). E.g., plt.figure()
    :param width_pct: Desired figure width (in % of page)
    return: String SVG representation of the figure.
    """
    f = StringIO()
    fig.tight_layout()

    # Make font sizes consistent

    fig.savefig(f, format='svg', bbox_inches='tight')
    svg = f.getvalue()

    svg = re.sub(r'width="[^"]+"', '', svg, count=1)
    svg = re.sub(r'height="[^"]+"', '', svg, count=1)

    svg = re.sub(r'^<svg', '<svg style="max-width: 100%; height: auto;"', svg, count=1)

    return f"<div class='hover-div' style='width: {width_pct}%; cursor: pointer;' onclick='downloadSVG(this)'>{svg}</div>"

def side_by_side(*figs):
    """
    Puts the provided SVG/HTML str figs into a <div> arranged side-by-side.
    :param figs: fig strs to put side by side
    return: A str div containing the provided figs arranged side by side
    """
    div_str = "<div style='display: flex; flex-direction: row; justify-content: space-evenly; flex-wrap: wrap;'>"
    for fig in figs:
        div_str += fig
    div_str += "</div>"
    return div_str

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
    
def split_by_tier(data, levels=1, semicolon=False):
    """
    Split a dict or a singly nested dict containing tier labels into subdicts keyed on tier number.
    :data: A dict (or nested dict) containing tiered labels
    :levels: How many levels the dict has (nested = 2). Default = 1 (not nested).
    :semicolon: Boolean flag for 'labels_per_gn_unique_co_localized' since the co-localization map for those is in the key.
    return: A dict containing the mappings of the data dict keyed on tier numbers.
    """
    # Needs to work with: "labels_per_gn_unique" (1)
    #                     "labels_per_gn_unique_co_localized" (1)
    #                     "total_labels" (1)
    #                     "label_products" (2)
    #                     "co_localized_label_products" (2)

    if levels not in [1, 2]:
        # Not allowed
        print(f"Error: 'levels' in 'split_by_tier()' must be one of [1, 2]. Received '{levels}'. Returning None.")
        return None
    
    if data is None:
        print(f"Error: provided data in 'split_by_tier()' is 'None'. Returning None.")
        return None
    
    if not semicolon:
        mapped = {}
        for label in data:
            base, tier = format_label(label)
            if tier not in mapped:
                mapped[tier] = {}
            
            mapped[tier][base] = data[label]
        
        return mapped

    else:
        mapped = {}
        for co_loc_label in data:
            labels = co_loc_label.split(';')

            tiers = {}

            for label in labels:
                base, tier = format_label(label)

                if tier not in tiers:
                    tiers[tier] = []
                
                if base not in tiers[tier]:
                    tiers[tier].append(base)
                
            
            for tier in tiers:
                names = sorted(tiers[tier])

                if len(names) > 1:
                    names_str = ';'.join(names)
                
                    if tier not in mapped:
                        mapped[tier] = {}
                    
                    if names_str not in mapped[tier]:
                        mapped[tier][names_str] = 0

                    mapped[tier][names_str] += data[co_loc_label]
        
        return mapped

def linear_regression(query_data):
    """
    Performs a linear regression on the query plot data (x1: # GNs, x2: total # genes, y: # labelled genes) to predict labelled genes 
    from the number of GNs and the total number of genes in a plot. Then calculates the residuals for the plot data and p-values for 
    those residuals. Returns a pandas dataframe with all the information. Used to find plots with more labelled genes than expected.
    :param query_data: The contents of a given 'query_name' subdict in the summary data JSON (e.g., data['query_names'][query_name_0]).
    return: Tuple containing: (A pandas dataframe with the following headers: 'path', 'total', 'gns', 'labelled', 'residuals', 'z_scores', and 'p_values', A dict containing the intercept, coefficients, and R squared for the model fit)
    """
    total, gns, labelled, path = [], [], [], []

    if query_data.get('plot_paths'):
        for plot_dict in query_data['plot_paths']:
            counts = plot_dict['counts']
            t, l, g = counts[0], counts[1], counts[2]

            total.append(t)
            gns.append(g)
            labelled.append(l)
            path.append(plot_dict['path'])

    df = pd.DataFrame({'path': path, 'total': total, 'gns': gns, 'labelled': labelled})

    model = LinearRegression()

    X, y = df[['total', 'gns']], df.labelled

    model.fit(X, y)

    model_data = {'intercept': model.intercept_, 'coefficients': model.coef_, 'r_square': model.score(X,y)}

    pred = model.predict(X)
    residuals = y - pred

    df['residuals'] = residuals

    # Get the Z score of the residuals (to get the p-values)
    # Z score: X - mean / stdev
    # The mean of the residuals should be approx. == 0 (since there should be an equal number above and below the linear regression line)
    # So Z score: X / stdev

    # stdev is sqrt(var)
    # var is (1/N) * sum((X-mean)^2) ==> (1/N) * sum(X^2)

    squared_res = residuals ** 2
    variance = (1/len(residuals)) * sum(squared_res)
    stdev = sqrt(variance)

    Z = residuals / stdev

    df['z_scores'] = Z
    df['p_values'] = 1 - norm.cdf(Z) # One tailed since we only care about plots with more labelled genes than expected

    return df, model_data

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

    return fig_to_svg(fig, width_pct='22.5')

def bar_chart(label_counts, title):
    """
    Return SVG string for horizontal bar chart with provided label counts.
    :param label_counts: Dict mapping label name to counts
    :param title: Title of the bar chart
    return: SVG string for the bar chart
    """
    fig, ax = plt.subplots(figsize=(8, 4))
    
    labels, values = zip(*label_counts.items())

    # trim_tag_labels = []
    # for label in labels:
    #     level = label.split("_")[-1]
    #     label = "_".join(label.split("_")[:-1])
    #     trim_tag_labels.append(label)

    ax.barh(labels[::-1], values[::-1], color='#105938')
    ax.set_title(title)
    ax.set_xlabel("Count")
    ax.set_ylabel("Annotation")

    return fig_to_svg(fig)

def co_localized_heatmap(data, title, max_n=None, is_int=True, axis_tag=""):
    """
    Generate a heatmap showing label-label or label-product co-localization numbers.
    :param data: A dict or nested dict containing the co-localization data.
    :param title: str title of the plot
    :param max_n: For label-product co-localization, max number of products to keep for each label (the most frequent n).
    :param is_int: True for raw counts, False for Log10 counts
    :param axis_tag: str to prepend to axis labels (e.g., "Tier 1" or "PSI-BLAST")
    return: SVG heatmap str
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
        xlabel = "NCBI Gene Annotations"
        ylabel = f"{axis_tag} Gene Annotations"
        bar_label = "Number of Co-Occurrences"
        if is_int:
            bar_label += " (Raw Counts)"
        else:
            bar_label += " (Log\u2081\u2080 Scale)"

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

        ylabel = xlabel = f"{axis_tag} Gene Annotations"
        bar_label = "Number of Co-Occurrences"
        if is_int:
            bar_label += " (Raw Counts)"
        else:
            bar_label += " (Log\u2081\u2080 Scale)"

    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(matrix, cmap='Greens', annot=False, ax=ax, cbar_kws={'label': bar_label}, square=True)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.yticks(rotation=0)
    return fig_to_svg(fig)

def sunburst_plot(label_products, title, is_psiblast=False, equal_wedges=False):
    """
    Generates a sunburst plot mapping labels to their product annotation counts.
    :param label_products: Dict mapping str label_name to dict str product_annotation: int count
    :param title: str title of the plot
    :param is_psiblast: Bool indicating if labels are from PSI-BLAST, if True doesn't try remove trailing tier # from label
    :param equal_wedges: Bool indicating if wedges should be equally sized (defaults to proportional wedges based on counts)
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

    if equal_wedges:
        fig = pex.sunburst(df, path=["Label", "Product"], width=700, height=700, title=title)
    
    else:
        fig = pex.sunburst(df, path=["Label", "Product"], values="Count", width=700, height=700, title=title)

    # Uncomment to remove labels for sunburst plot (for figure)
    # fig.update_traces(textinfo='none')

    # Make sure titles are consistent
    fig.update_layout(
        title_font=dict(color='black',
                        size=14,
                        family="sans-serif")
    )

    return pio.to_html(fig, include_plotlyjs=False, full_html=False)


# --------------------------------------------
# HTML functions
# --------------------------------------------
def make_query_section(query_data, num_flanking_genes=15):
    """
    Build the HTML section for a single query. Return the HTML as a string.
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

    # Box and whisker plots & pie chart side-by-side
    section_string += side_by_side(box_and_whisker_plots(query_data), percent_labelled_pie(query_data))

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

    section_string += f"<p><b>Percentage (%) of Genes Annotated: {average_labels_per_gene * 100:.0f}%</b></p>"
    section_string += f"<p><b>Average Number of Annotated Genes per Genomic Neighbourhood: {average_labels_per_gn:.2f}</b></p>"

    linear_regression_df, linear_regression_data = linear_regression(query_data)

    # Use a p-value cutoff of 0.15 (there aren't usually huge outliers so 0.05 rarely has a single plot, 0.15 looks like it gives a decent number of plots with interesting labels)
    highly_labelled_plots = linear_regression_df[linear_regression_df['p_values'] <= 0.15].sort_values('p_values', ascending=True)

    if len(highly_labelled_plots) > 0:
        section_string += "<p><b>Plots with More Labelled Genes than Expected (p-value <= 0.15)</b></p>"
        section_string += "<ul>"
        section_string += "<p>p-value: Plot Name</p>"

        for row in highly_labelled_plots.itertuples(index=False):
            section_string += f"<li>{row.p_values:.2f}: <a href=/read_file?path={row.path}>{row.path.split('/')[-1]}</a></li>"

        section_string += "</ul>"
    else:
        section_string += "<p><b>No plots found with more labelled genes than expected.</b></p>"

    # PSI-BLAST Section
    if query_data['num_psiblast'] > 0:
        section_string += "<h2>PSI-BLAST Annotations Summary</h2>"

        # Bar Charts
        bars = []
        if query_data.get('psiblast_per_gn_unique'):
            bars.append(bar_chart(query_data['psiblast_per_gn_unique'], "Total Genomic Neighbourhoods with PSI-BLAST Annotations"))
        
        if query_data.get('total_psiblast'):
            bars.append(bar_chart(query_data['total_psiblast'], "Total Genes with PSI-BLAST Annotations"))

        if bars:
            section_string += side_by_side(*bars)

        # Sunburst
        if query_data.get('psiblast_products'):
            title = "Summary of NCBI Product Annotations Linked to PSI-BLAST Labels"
            # section_string += f"<h3>{title}</h3>"
            section_string += side_by_side(sunburst_plot(query_data['psiblast_products'], title, is_psiblast=True), sunburst_plot(query_data['psiblast_products'], title + " (equal wedges)", is_psiblast=True, equal_wedges=True))

        # Product heatmap
        if query_data.get('co_localized_psiblast_products'):
            hmaps = []
            section_string += "<h3>NCBI Product Annotation Heatmaps</h3>"
            title = "NCBI Product Annotations Linked to PSI-BLAST Labels"
            max_products = 5
            # section_string += f"<h3>{title}</h3>"
            section_string += f"<p>Retaining only the top {max_products} products from each label.</p>"
            hmaps.append(co_localized_heatmap(query_data['co_localized_psiblast_products'], title, max_n=max_products, is_int=True, axis_tag="PSI-BLAST"))

            title = "Log\u2081\u2080 " + title
            # section_string += f"<h3>{title}</h3>"
            logged_data = log_data(query_data['co_localized_psiblast_products'])
            hmaps.append(co_localized_heatmap(logged_data, title, max_n=max_products, is_int=False, axis_tag="PSI-BLAST"))

            section_string += side_by_side(*hmaps)

    # Tier Label Sections
    if query_data.get("num_labeled") and query_data['num_labeled'] > 0:
        section_string += "<h2>User-Defined Annotations Summary by Tier</h2>"

        # Loop over each tier
        tiers = []
        tier_map = {}

        for label in query_data["labels_per_gn_unique"]:
            tier = int(label.split('_')[-1])
            if tier not in tiers:
                tiers.append(tier)
                tier_map[tier] = [label]
            else:
                tier_map[tier].append(label)

        tiers = sorted(tiers)

        tiered_labels_per_gn_unique = split_by_tier(query_data.get('labels_per_gn_unique'))
        tiered_labels_per_gn_unique_co_localized = split_by_tier(query_data.get('labels_per_gn_unique_co_localized'), semicolon=True)
        tiered_total_labels = split_by_tier(query_data.get('total_labels'))
        tiered_label_products = split_by_tier(query_data.get('label_products'), levels=2)
        tiered_co_localized_label_products = split_by_tier(query_data.get('co_localized_label_products'), levels=2)

        for tier in tiers:
            tier = str(tier)
            section_string += f"<h2>Tier {tier} Summary</h2>"

            # Bar Charts
            bars = []
            if tiered_labels_per_gn_unique and tiered_labels_per_gn_unique.get(tier):
                bars.append(bar_chart(tiered_labels_per_gn_unique.get(tier), f"Total Genomic Neighbourhoods with User-Defined Annotations: Tier {tier}"))
            if tiered_total_labels and tiered_total_labels.get(tier):
                bars.append(bar_chart(tiered_total_labels.get(tier), f"Total Genes with User-Defined Annotations: Tier {tier}"))
            if bars:
                section_string += side_by_side(*bars)

            # Sunburst
            if tiered_label_products and tiered_label_products.get(tier):
                section_string += side_by_side(sunburst_plot(tiered_label_products.get(tier), f"Summary of NCBI Product Annotations Linked\n to User-Defined Labels: Tier {tier}", is_psiblast=False), sunburst_plot(tiered_label_products.get(tier), f"Summary of NCBI Product Annotations Linked\n to User-Defined Labels: Tier {tier} (equal wedges)", is_psiblast=False, equal_wedges=True))

            # Label Heatmap
            if tiered_labels_per_gn_unique_co_localized and tiered_labels_per_gn_unique_co_localized.get(tier):
                hmaps = []
                title = f"Co-Occurrence Patterns Among \nTier {tier} User-Annotated Genes"
                hmaps.append(co_localized_heatmap(tiered_labels_per_gn_unique_co_localized.get(tier), title, is_int=True, axis_tag=f"Tier {tier}"))

                logged_data = log_data(tiered_labels_per_gn_unique_co_localized.get(tier))
                
                hmaps.append(co_localized_heatmap(logged_data, "Log\u2081\u2080 " + title, is_int=False, axis_tag=f"Tier {tier}"))

                section_string += side_by_side(*hmaps)

            # Product Heatmap
            if tiered_co_localized_label_products and tiered_co_localized_label_products.get(tier):
                hmaps = []
                section_string += "<h3>NCBI Product Annotation Heatmaps</h3>"
                title = f"NCBI Product Annotations Linked to Tier {tier} User-Defined Labels"
                max_products = 5

                # section_string += f"<h3>{title}</h3>"
                section_string += f"<p>Retaining only the top {max_products} products from each label.</p>"
                hmaps.append(co_localized_heatmap(tiered_co_localized_label_products.get(tier), title, max_n=max_products, is_int=True, axis_tag=f"Tier {tier}"))

                title = "Log\u2081\u2080 " + title
                # section_string += f"<h3>{title}</h3>"
                # section_string += f"<p>Retaining only the top {max_products} products from each label.</p>"
                logged_data = log_data(tiered_co_localized_label_products.get(tier))
                hmaps.append(co_localized_heatmap(logged_data, title, max_n=max_products, is_int=False, axis_tag=f"Tier {tier}"))

                section_string += side_by_side(*hmaps)

    return section_string


def write_full_report(query_sections, query_names, outfile):
    """
    Writes the final HTML report by adding header info and joining provided sections. If len(query_sections) == 1, removes overall header and formats for a single query name.
    :param query_sections: List of the HTML for each query name.
    :param query_names: List of str names of the queries (in same order as query_sections).
    :param outfile: The location and name of the file to write to.
    return: None, writes to outfile
    """

    outstr = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>GNAT Summary Report</title>
        <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
        <style>
            .hover-div {
                border: 2px solid transparent;
            }

            .hover-div:hover {
                border: 2px solid #105938;
            }
        </style>
    </head>
    <body style="font-family: 'Open Sans', sans-serif; padding: 40px;">
        <h1>"""

    if len(query_names) == 1:
        outstr += f"GNAT Output Summary for {query_names[0]}</h1>"

    else:
        outstr += f"GNAT Output Summary for all Job Queries</h1>"

    outstr += "<p>Any figure (other than sunburst plots) may be clicked on to download the SVG. To download the sunburst plot, hover over it and click on the camera button to download a PNG.</p>"

    outstr += "<p>For help interpreting the figures in this report, click <a target='_blank' href='https://gnat-dev.usask.ca/report_help'>here</a>.</p>"

    for i, query_section in enumerate(query_sections):
        if len(query_names) != 1:
            outstr += f"<h1>Summary for {query_names[i]}</h1>"
        
        outstr += query_section
    
    # Add script to allow SVG downloads
    outstr += """
    <script>
        function downloadSVG(container) {
            const svg = container.querySelector('svg');
            const svgData = new XMLSerializer().serializeToString(svg);
            const blob = new Blob([svgData]);
            const element = document.createElement("a");
            element.download = "figure.svg";
            element.href = window.URL.createObjectURL(blob);
            element.click();
            element.remove();
        }
    </script>
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
    parser.add_argument("-o", "--outfile", help="Location and name of the output visualization file. If not provided, outputs 'gnat_summary_report.html' to the directory containing the summary JSON file.")
    parser.add_argument("-n", "--num_flanking_genes", help="Maximum number of flanking genes in each GN. Used to calculate label averages to display the most densly labelled plots. Defaults to 15.", default=15)

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

    flanking = int(args.num_flanking_genes)

    # Parse out location
    if args.outfile:
        outfile = args.outfile
    else:
        # Use JSON dir
        outfile = path.join(json_dir, "gnat_summary_report.html")


    # Load JSON
    with open(args.json_path) as f:
        summary_json = json.load(f)

    # Check if provided query is in the JSON
    if provided_query and provided_query not in summary_json['query_names'].keys():
        print(f"Warning: Provided query name {provided_query} not found in the summary JSON. Generating visualization using all query names.")
        provided_query = None

    if provided_query:
        # If provided query is in the JSON, build the summary for just that query
        query_names = [provided_query]
        query_section = make_query_section(summary_json['query_names'][provided_query], num_flanking_genes=flanking)
        write_full_report([query_section], query_names, outfile)

    else:
        # Build for all queries
        query_sections = []
        query_names = []

        for query_name in summary_json['query_names']:
            query_names.append(query_name)
            query_sections.append(make_query_section(summary_json['query_names'][query_name], num_flanking_genes=flanking))

        write_full_report(query_sections, query_names, outfile)