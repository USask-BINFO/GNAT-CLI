#! /usr/bin/env python3

"""
Script to sort the SVG map (now PNG map) in human sort (e.g., lexographically until the cluster number). Also puts 'all_non_redundants' after non-'all_non_redundant' GNs.
This is to make the dashboard display in a better order (and match the order of the directory browser)
"""

import re
from argparse import ArgumentParser
from os import path

def sort_key(map_line):
    non_redundant_after = 2 if "all_non_redundants" in map_line else 1

    match = re.search(r"(.*?_cluster_)(\d+)", map_line)

    if match:
        lex, num = match.groups()
        num = int(num)
    else:
        lex, num = map_line, float("inf")
    
    return (non_redundant_after, lex, num, map_line)

if __name__ == "__main__":
    parser = ArgumentParser("Sort an SVG/PNG map file to make the dashboard display nicer.")
    parser.add_argument("map_file", help="Path to the map file to sort.")

    args = parser.parse_args()

    if not path.exists(args.map_file):
        raise FileNotFoundError(f"Map file not found at {args.map_file}")

    with open(args.map_file, "r") as f:
        lines = [line.strip() for line in f]
    
    sorted_lines = sorted(lines, key=sort_key)

    with open(f"{args.map_file}", "w") as f:
        for line in sorted_lines:
            f.write(f"{line}\n")
