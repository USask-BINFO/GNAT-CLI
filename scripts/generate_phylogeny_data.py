#! /usr/bin/env python
"""
Script to pull taxonomy data from a list of GNs, and generate phylogenetic trees and pie charts from that data.
"""

import argparse
from Bio import SeqIO
import gzip
from os import path, listdir, makedirs, environ
import re
import sys

from ete3 import TreeStyle, Tree, TextFace, add_face_to_node, NodeStyle

import matplotlib.pyplot as plot

environ["QT_QPA_PLATFORM"] = "offscreen"

# Parse inputs
parser = argparse.ArgumentParser("Pull taxonomy data from a list of GNs and generate phylogenetic tree and pies.")
parser.add_argument("GNs", nargs="+", help="Genomic neighbourhoods (or .gbff.gz files) to generate taxonomic data for. "
                                           "Can also just point to the directory containing the files of interest.")
parser.add_argument("-n", "--name", help="Name of query or gbff source. E.g, AcrIF7, bacterialdb, Aca1. Required.", required=True)
parser.add_argument("-o", "--out", help="Desired location of output. If not specified, outputs to current directory.")
parser.add_argument("-m", "--mode", default="GN", help="Specify file type: GBFF if providing .gbff.gz files. Defaults "
                    "to GN.", choices=["GN", "GBFF"])
parser.add_argument("-dt", "--data_type", help="Expected data source (bacterial, fungal, viral, etc.) for generating " 
                    "phylogenetic data. If 'bacterial', will generate phylogenetic trees and pie charts, else will "
                    "only output a csv containing the data.")

args = parser.parse_args()

if args.out is not None:
    outdir = args.out
else:
    outdir = "."
if outdir[-1] != "/":
    outdir += "/"

if args.data_type is not None:
    phylodir_name = "%s_%s_phylogeny_data" % (args.name, args.data_type)
else:
    phylodir_name = "%s_phylogeny_data" % args.name

if args.data_type == "bacterial":
    if not path.isdir(path.join(outdir, phylodir_name)):
        makedirs(path.join(outdir, phylodir_name))
    outlocation = path.join(outdir, phylodir_name)
else:
    outlocation = outdir

# Check if a directory was provided, and get genome names
genomes = []
if len(args.GNs) == 1:
    if path.isdir(args.GNs[0]):
        for file in listdir(args.GNs[0]):
            if args.mode != "GBFF":
                if '.gb' == file[-3:]:
                    genomes.append(path.join(args.GNs[0], file))
            else:
                if '.gbff.gz' == file[-8:]:
                    genomes.append(path.join(args.GNs[0], file))
    else:
        genomes = args.GNs
else:
    for g in args.GNs:
        genomes.append(g)

# ----------------------------------------------------------------------------------------------------------------------
# Parse taxonomy data

count = 0
tax_list = []
# Parse genomes to generate nested list of taxonomy data
for gn in genomes:
    if args.mode == "GBFF":
        try:
            with gzip.open(gn, "rt") as f:
                count += 1
                recordlist = list(SeqIO.parse(f, "genbank"))
                tax_list.append(recordlist[0].annotations["taxonomy"])
                tax_list[-1].append(recordlist[0].annotations["organism"])
        except Exception as e:
            print("Error: %s when opening file %s" % (str(e), str(gn)))
    else:
        for seq_record in SeqIO.parse(gn, "genbank"):
            count = count + 1
            tax_list.append(seq_record.annotations["taxonomy"])
            tax_list[-1].append(seq_record.annotations["organism"])

# ----------------------------------------------------------------------------------------------------------------------
# Generate taxonomic data csv
taxstr = ""
if args.data_type is not None:
    outfilename = "%s_%s_taxonomy.csv" % (args.name, args.data_type)
else:
    outfilename = "%s_taxonomy.csv" % args.name
for line in tax_list:
    for i in range(len(line)):
        line[i].replace(",", "_")
    taxstr += ",".join(line) + "\n"
with open(path.join(outlocation, outfilename), "w") as f:
    f.write(taxstr)

if args.data_type != "bacterial":   # Exit here to prevent generation of phylo figures if datatype isn't bacterial
    exit()

# ----------------------------------------------------------------------------------------------------------------------
# Only 'bacterial' data_type onward
# Parse lines to avoid crashes in the visualization step (and to make everything look nice)
def remove_erroneous_chars(line):
    """
    Remove characters that would cause problems or look bad from a tax line
    :param line: taxonomy line from tax_list
    return: line with only alphanumeric characters and periods
    """
    for i in range(len(line)):
        result = re.sub('[^a-zA-Z. ]', '', line[i])
        result = re.sub('[_]', ' ', result)
        line[i] = result
    return line


def parse_bacterial(line):
    """
    Parse bacterial taxonomy line
    :param line: taxonomy line from tax_list
    return: Line that is able to be nicely visualized below
    """
    line = remove_erroneous_chars(line)
    for i in range(len(line)):
        if "unclassified" in line[i]:
            line = line[:i+1]
            break
        if i == len(line) - 1:
            l = line[i].split(" ")
            line[i] = " ".join(l[:2])
    return line


bacterial = []
wonky = []
for tax in tax_list:
    # tax = remove_erroneous_chars(tax)
    if re.match("bacteria", tax[0], re.IGNORECASE):
        bacterial.append(parse_bacterial(tax))
        if tax[0] != "Bacteria":
            wonky.append(tax)


# ----------------------------------------------------------------------------------------------------------------------
# TODO: Generate phylogenetic tree
class BasicTree:
    """A simple tree class to be used to build a phylogenetic tree."""
    def __init__(self, val=None, depth=0, parent=None):
        """

        """
        self.val = val
        self.depth = depth
        self.children = []
        self.parent = parent

    def is_root(self):
        if self.depth == 0 and self.parent is None:
            return True
        return False

    def is_empty(self):
        if self.val is None:
            return True
        return False

    def has_children(self):
        if len(self.children) > 0:
            return True
        return False

    def add_child(self, val):
        child = BasicTree(val, self.depth+1, self)
        self.children.append(child)

    def get_children_values(self):
        vals = []
        for child in self.children:
            vals.append(child.val)
        return vals

    def get_parent_val(self):
        if self.is_root():
            return None
        return self.parent.val

    def search(self, val):
        if self.val == val:
            return self
        else:
            if self.has_children():
                for child in self.children:
                    result = child.search(val)
                    if result is not None:
                        return result
            else:
                result = None
        return result

    def to_string_newick(self):
        class Stack:
            def __init__(self, val=None):
                if val is not None:
                    self.__s = [val]
                else:
                    self.__s = []

            def push(self, val):
                self.__s.append(val)

            def pop(self):
                return self.__s.pop()

        newick_string = ""
        stack = Stack()
        if self.has_children():
            stack.push(str(self.val))
            newick_string += "("
            children_strings = []
            for child in self.children:
                children_strings.append(child.to_string_newick())
            for child_string in children_strings:
                newick_string += child_string
                newick_string += ","
            newick_string = newick_string[:-1]
            newick_string += ")" + str(self.val)
        else:
            newick_string += str(self.val)
        if self.is_root():
            return newick_string + ";"
        else:
            return newick_string

# Create the tree
italics = [[]]
trees = []
for bac_line in bacterial:
    for item_index in range(len(bac_line)):
        item = bac_line[item_index]

        if item_index == 0:  # First item in line list (root)
            if len(trees) < 1:
                trees.append(BasicTree(item))
            else:
                found = None
                for t in range(len(trees)):
                    res = trees[t].search(item)
                    if res is not None:
                        found = t
                        break
                if found is None:
                    trees.append(BasicTree(item))
                    italics.append([])

        elif item_index == len(bac_line)-1:  # Last item in line list
            parent = bac_line[item_index-1]
            for t in range(len(trees)):
                res = trees[t].search(parent)
                if res is not None:
                    if bac_line[item_index] not in res.get_children_values():
                        res.add_child(bac_line[item_index])
                        if "p. " not in bac_line[item_index]:
                            if " bacterium" not in bac_line[item_index]:
                                if "unclassified" not in bac_line[item_index]:
                                    italics[t].append(bac_line[item_index])
                                    italics[t].append(bac_line[item_index-1])

        else:  # Neither first nor last item in line list
            parent = bac_line[item_index-1]
            for t in range(len(trees)):
                res = trees[t].search(parent)
                if res is not None:
                    if bac_line[item_index] not in res.get_children_values():
                        res.add_child(bac_line[item_index])

for tree_index in range(len(trees)):
    t = Tree(trees[tree_index].to_string_newick(), format=1)

    ts = TreeStyle()
    ts.show_leaf_name = False

    def my_layout(node):
        if node.name in italics[tree_index]:
            F = TextFace(node.name, tight_text=True, fstyle="italic")
        else:
            F = TextFace(node.name, tight_text=True)
        add_face_to_node(F, node, column=0, position="branch-right")

    ts.layout_fn = my_layout
    ts.branch_vertical_margin = 33
    ts.show_scale = False
    ts.title.add_face(TextFace("Phylogenetic Tree of Matches from "+args.name), column=0)
    ts.scale = 100

    ns = NodeStyle()
    ns["size"] = 0

    for n in t.iter_descendants():
        n.img_style = ns

    t.img_style = ns

    t.render(args.name+"_phylo_tree_"+str(tree_index)+".svg", tree_style=ts)

# ----------------------------------------------------------------------------------------------------------------------
# TODO: Generate Pie Charts
class TaxPie:
    def __init__(self):
        """

        """
        return
