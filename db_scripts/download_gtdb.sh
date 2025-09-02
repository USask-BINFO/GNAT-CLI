#! /usr/bin/env bash

# Script to download the bacterial and archaeal genomes from the latest release of the Genome Taxonomy Database (GTDB)
# (in the format needed by GNAT)

# Set defaults
curr_dir=$(pwd)
out_location=$(pwd)
script_dir=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            echo "  Script to build a GNAT-compatible database of the latest GTDB genomes (bacterial + archaeal)."
            echo "Potentially Required Arguments:"
            echo "  [-s][--db_scripts]      Path to directory containing db_scripts. Must be specified if "
            echo "                              db_scripts/ not added to PATH."
            echo "Optional Arguments:"
            echo "  [-o][--out_location]    Directory to make the database directory in. Defaults to the current "
            echo "                              directory. Creates the database in a gtdb/ subdirectory at the "
            echo "                              specified directory."
            exit 0
        ;;
        -s|--db_scripts)
            db_scripts="$2"
            shift
            shift
        ;;
        -o|--out_location)
            out_location="$2"
            shift
            shift
        ;;
    esac
done

if [ -z ${db_scripts+x} ]; then
    echo "[-s DB_SCRIPTS] not specified. Assuming db_scripts/ has been added to PATH."
else
    # Make script_dir an absolute path
    script_dir=$(realpath "$db_scripts")
    script_dir="${script_dir%/}/"
fi

# Get the abs path to out location
out_location=$(realpath "$out_location")
out_location="${out_location%/}/"

# Create the gtdb directory in out_location
if [ -d "${out_location}/gtdb" ]; then
    echo "Error: gtdb/ already exists in provided output location, rename or remove it to install gtdb/ in the specified location."
    exit 1
else
    echo "Creating gtdb directory at ${out_location}/gtdb"
    mkdir "${out_location}/gtdb"
fi

cd "${out_location}/gtdb" || exit

# Download the necessary files here
GTDB_latest="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"

bac_tax_file="bac120_taxonomy.tsv"
arc_tax_file="ar53_taxonomy.tsv"

bac_tax_link="${GTDB_latest}${bac_tax_file}"
arc_tax_link="${GTDB_latest}${arc_tax_file}"

NCBI_ftp="ftp://ftp.ncbi.nih.gov/genomes/"

rs_asssum="assembly_summary_refseq.txt"
gb_asssum="assembly_summary_genbank.txt"

refseq_assembly_summary_link="${NCBI_ftp}refseq/${rs_asssum}"
genbank_assembly_summary_link="${NCBI_ftp}genbank/${gb_asssum}"

wget "$bac_tax_link"
wget "$arc_tax_link"

wget "$refseq_assembly_summary_link"
wget "$genbank_assembly_summary_link"

# Ensure the assembly summary files are sorted on the first column (and remove header lines)
rs_sum_sorted="sorted_${rs_asssum}"
gb_sum_sorted="sorted_${gb_asssum}"

tail -n +3 "$rs_asssum" | sort -t$'\t' -k1,1 > "$rs_sum_sorted"
tail -n +3 "$gb_asssum" | sort -t$'\t' -k1,1 > "$gb_sum_sorted"

# Make subdirectories for a bacterial GTDB, archaeal GTDB, and combined GTDB
mkdir gtdb_bacterial
mkdir gtdb_archaeal
mkdir gtdb_combined

touch gtdb_bac_urls.txt
touch gtdb_arc_urls.txt

# Get URLs for bacterial genomes
cut -f1 "$bac_tax_file" | while read -r line; do
    if [ "${line:0:3}" = "RS_" ]; then
        # RefSeq
        gc_acc="${line:3}"
        url=$(look "${gc_acc}" "$rs_sum_sorted" | awk -v id="$gc_acc" -F'\t' '$1==id {print $20}')
        end=$(echo "$url" | rev | cut -d '/' -f1 | rev)
        url="${url}/${end}_genomic.gbff.gz"
        echo "$url" >> gtdb_bac_urls.txt
    elif [ "${line:0:3}" = "GB_" ]; then
        # GenBank
        gc_acc="${line:3}"
        url=$(look "${gc_acc}" "$gb_sum_sorted" | awk -v id="$gc_acc" -F'\t' '$1==id {print $20}')
        end=$(echo "$url" | rev | cut -d '/' -f1 | rev)
        url="${url}/${end}_genomic.gbff.gz"
        echo "$url" >> gtdb_bac_urls.txt
    else
        # Shouldn't happen
        echo "Warning: Line '${line}' does not start with RS_ or GB_. Skipping..."
    fi
done

# Get URLs for archaeal genomes
cut -f1 "$arc_tax_file" | while read -r line; do
    if [ "${line:0:3}" = "RS_" ]; then
        # RefSeq
        gc_acc="${line:3}"
        url=$(look "${gc_acc}" "$rs_sum_sorted" | awk -v id="$gc_acc" -F'\t' '$1==id {print $20}')
        end=$(echo "$url" | rev | cut -d '/' -f1 | rev)
        url="${url}/${end}_genomic.gbff.gz"
        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/898/805/GCF_009898805.1_ASM989880v1
        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/898/805/GCF_009898805.1_ASM989880v1/GCF_009898805.1_ASM989880v1_genomic.gbff.gz

        echo "$url" >> gtdb_arc_urls.txt
    elif [ "${line:0:3}" = "GB_" ]; then
        # GenBank
        gc_acc="${line:3}"
        url=$(look "${gc_acc}" "$gb_sum_sorted" | awk -v id="$gc_acc" -F'\t' '$1==id {print $20}')
        end=$(echo "$url" | rev | cut -d '/' -f1 | rev)
        url="${url}/${end}_genomic.gbff.gz"
        echo "$url" >> gtdb_arc_urls.txt
    else
        # Shouldn't happen
        echo "Warning: Line '${line}' does not start with RS_ or GB_. Skipping..."
    fi
done

# Move URL files to respective directories
mv gtdb_bac_urls.txt gtdb_bacterial
mv gtdb_arc_urls.txt gtdb_archaeal

for subdir in gtdb_bacterial gtdb_archaeal; do
    cd "$subdir"
    # Make a directory to store the downloaded genomes
    mkdir -p "${subdir}_genomes"
    cd "${subdir}_genomes"

    # Download teh genomes
    for url in $(cat ../gtdb*urls.txt); do
        bn=$(basename "$url");
        if [ ! -e "$bn" ]; then
            echo "Downloading URL $url";
            wget -t 1 --timeout=5 "$url";
        else
            echo "Already downloaded $url"
        fi
    done

    # Check for missing genomes
    for genome in $(cat ../gtdb*urls.txt | cut -d "/" -f 11 | sort -u); do
        if [ ! -e "$genome" ]; then
            echo "$genome";
        fi
    done

    # Make a .fasta containing all the new proteins
    echo "Making ${subdir}_proteins.FASTA"
    "${script_dir}"make_gbff_protein_fasta.py . > ../${subdir}_proteins.FASTA

    cd ..

    # Make the subdirectories
    mkdir -p "${subdir}_genomes/proteins"
    mkdir -p "${subdir}_genomes/gbff"

    cat "${subdir}_proteins.FASTA" >> "${subdir}_genomes/proteins/all_"$subdir"_proteins.FASTA"

    # Move new gbffs to main gbff folder
    for gbff in `find "${subdir}_genomes" -name "*.gbff.gz"`; do mv "$gbff" "${subdir}_genomes/gbff"; done

    # Build the BLAST database
    cd "${subdir}_genomes/proteins"
    makeblastdb -in "all_${subdir}_proteins.FASTA" -dbtype prot -title "${subdir}db" -out "${subdir}db"

    mkdir -p ../../${subdir}db

    mv "${subdir}db".* "../../${subdir}db"

    # Generate all_pal_files.txt for PSI-BLAST searches
    cd "../../${subdir}db"

    cat "${subdir}db.pal" | "${script_dir}"split_pal.py

    # Move back to starting dir
    cd "$out_location/gtdb"
done
cd "$curr_dir"
