#!/usr/bin/env python3

# This file was used to download the bacterial genomes from the NCBI website
# Only outputs files that are recognized as files
# First arg = link to database
# Link to database should be:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/.../assembly_summary.txt
# where ... is the path to the specific database you want, some examples below:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Pseudomonas_aeruginosa/assembly_summary.txt for Pseudomonas aeruginosa db
# ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt for fungal db

import urllib.request
from os import path
from sys import argv

link = argv[1]

if "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank" != link[:42]:
    raise ValueError("Error: Provided link does not appear to be from the ncbi ftp site " +
    "(ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank). Provided link: '%s'" % link)
if "/assembly_summary.txt" != link[-21:]:
    raise ValueError("Error: Provided link doesn't end in /assembly_summary.txt. This file " +
    "is required to build the database. Provided link: '%s'" % link)

response = urllib.request.urlopen(link)
data = response.read()
text = data.decode("utf-8")

lines = text.split("\n")

lines = lines[2:-1]

for i in range(len(lines)):
    line = lines[i].rstrip()
    organism = line.split("\t")[7]
    #dl_directory = line.split("\t")[19]
    #"""
    try:
        dl_directory = line.split("\t")[19]
    except IndexError:
        string = ""
        for item in line:
            string += str(item)
            string += "\t"
        string += "BIGERROR"
        print(string)
    #"""
    organism = organism.replace(" ", "_")
    dl_name = dl_directory.split("/")[-1]

    filename =  dl_name + "_genomic.gbff.gz"

    url = dl_directory + "/" + filename

    #potentially_existing_filename = os.environ["HOME"] + "/NCBI/bacterial_genomes/gbff/" + filename
    potentially_existing_filename = "/birl2/data/Acr/data/fungal_genomes/gbff/" + filename

    if path.isfile(potentially_existing_filename):
        pass
        #print(url)
        #print(url + "\t" + str(i) + "\t" + lines[i])
    else:
        print(url)
