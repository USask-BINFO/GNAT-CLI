#!/usr/bin/env bash

# Updates or builds a database to be used for the acrvis pipeline
# Creates (or updates) a directory containing a sub-directory with a blast database as well as all of the .gbff files
# that made the blast database.
# May take a very long time to run (month), run inside a screen.
# For argument and usage information, run ./get_ncbi_genomes.sh -h

if [[ $# -lt 1 ]]; then
  echo "No arguments provided, use [-h][--help] for usage help."
  exit 1
fi

# Set default arguments
curr_dir=$(pwd)
out_location=$(pwd)
no_backup="false"
script_dir=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--db_scripts)
      db_scripts="$2"
      shift
      shift
      ;;
    -n|--name)
      name="$2"
      shift
      shift
      ;;
    -l|--link)
      link="$2"
      shift
      shift
      ;;
    -o|--out_location)
      out_location="$2"
      shift
      shift
      ;;
    -e|--existing_db)
      existing_db="$2"
      shift
      shift
      ;;
    -b|--no_backup)
      no_backup="true"
      shift
      ;;
    -h|--help)
      echo "    Script to update/build a blast database."
      echo "Required Arguments:"
      echo "  [-l][--link]              NCBI ftp link to download database from."
      echo "Potentially Required Arguments:"
      echo "  [-s][--db_scripts]        Path to directory containing db_scripts. Must be specified if "
      echo "                              db_scripts/ not added to PATH."
      echo "Optional Arguments:"
      echo "  [-b][--no_backup]         If this flag is specified, and a valid [-e][--existing_db] was "
      echo "                              specified (i.e., a database is being updated), then a backup of"
      echo "                              the old database will not be saved. Default saves a backup."
      echo "  [-e][--existing_db]       Location of existing database directory to update. If not provided,"
      echo "                              will build a new database in [-o][--out_location]."
      echo "  [-h][--help]              Prints this message."
      echo "  [-n][--name]              Name of database to be built (e.g., fungal, bacterial, viral). If "
      echo "                              not provided, defaults to the item before 'assembly_summary.txt'"
      echo "                              (second last item) in the link. Will be ignored if a valid "
      echo "                              [-e][--existing_db] is specified."
      echo "  [-o][--out_location]      Directory to make database directory. Defaults to current directory."
      echo "                              Will be ignored if a valid [-e][--existing_db] is specified."
      echo ""
      echo "Note: Recommended to run this script in a screen, may take a very long time to run, "
      echo "depending on database size. Be sure to have enough space for the database."
      echo ""
      exit 0
      ;;
    *)
      echo "Error: "$1" not recognized as a valid argument. Use [-h][--help] for usage help."
      exit 1
      ;;
  esac
done

echo "Starting..."

# Make sure required argument was provided
if [ -z ${link+x} ]; then
  echo "Must specify a link to download the database from."
  exit 1
fi

if [ -z ${db_scripts+x} ]; then
  echo "[-s DB_SCRIPTS] not specified, assuming db_scripts/ has been added to PATH."
else
  # Make script_dir an absolute path
  script_dir=$(cd "$db_scripts" || exit; pwd; cd "$curr_dir" || exit;)
  script_dir="$script_dir/"
fi

# Get the name to use for the db from the link if not specified
if [ -z ${name+x} ]; then
  name=$(echo "$link" | rev | cut -d "/" -f2 | rev);
fi

# Absolute-ify the path to out_location
out_location=$(cd "$out_location" || exit; pwd; cd "$curr_dir" || exit);

# If creating a new db, create the necessary directories
if [ -z ${existing_db+x} ]; then
  echo "Creating database directory in $out_location"
  echo "Directory name: ${name}_db"
  outdir="$out_location"/"${name}"_db
  mkdir "$outdir"
  existing="null"
else
  existing=$(cd "$existing_db" || exit; pwd; cd "$curr_dir" || exit);
  outdir="$existing"
  echo "Updating database at $existing"
fi

# Date tag for database name
MM_YYYY=$(date +%m_%Y)

# Get list of new ncbi genomes (or all ncbi genomes if building from scratch)
"${script_dir}"get_list_of_new_ncbi_genomes.py "$link" "$existing" | grep -v 'na/na' | sort | uniq > "${outdir}"/new_"${name}"_genomes_urls.txt
# exit 0 # DEBUG
# Make a directory to store the new downloaded genomes
mkdir -p "${outdir}"/new_"${name}"_genomes_$$
cd "${outdir}"/new_"${name}"_genomes_$$ || exit

# Download new genomes
for url in $(cat "${outdir}"/new_"${name}"_genomes_urls.txt); do bn=$(basename "$url"); if [ ! -e "$bn" ]; then echo "Downloading URL $url"; wget -t 1 --timeout=5 $url; else echo "Already downloaded $url"; fi; done

# Check for missing genomes
for genome in $(cat "${outdir}"/new_"${name}"_genomes_urls.txt | cut -d "/" -f 11 | sort | uniq); do if [ ! -e "$genome" ]; then echo "$genome"; fi; done

# Make a .fasta containing all of the new proteins
echo Making new_"${name}"_proteins.FASTA
"${script_dir}"make_gbff_protein_fasta.py . > "${outdir}"/new_"${name}"_proteins.FASTA

echo Updating all_"${name}"_proteins.FASTA

# Make directories in case they don't already exist
mkdir -p "${outdir}"/"${name}"_genomes/proteins
mkdir -p "${outdir}"/"${name}"_genomes/gbff

# Update (or create) master protein file
cat "${outdir}"/new_"${name}"_proteins.FASTA >> "${outdir}"/"${name}"_genomes/proteins/all_"${name}"_proteins.FASTA

# Move newly downloaded gbff files to main gbff folder
for thing in `find . -name "*.gbff.gz"`; do mv "$thing" "${outdir}"/"${name}"_genomes/gbff; done

# Remove unneeded files
rm -f "${outdir}"/new_"${name}"_proteins.FASTA
rmdir new_"${name}"_genomes_*

# Build the BLAST database
cd "${outdir}"/"${name}"_genomes/proteins || exit
makeblastdb -in "${outdir}"/"${name}"_genomes/proteins/all_"${name}"_proteins.FASTA -dbtype prot -title "${name}"db_$MM_YYYY  -out "${name}"db_$MM_YYYY

mkdir -p ${outdir}/${name}db
if [[ "$no_backup" != "true" ]]; then
  # Make a backup directory (if it doesn't exist)
  mkdir -p "${outdir}"/"${name}"db_backup
  # Empty the backup directory (just in case it already exists)
  rm -rf "${outdir}"/"${name}"db_backup/*
  # Move the contents of the main db directory to the backup
  mv "${outdir}"/"${name}"db/* "${outdir}"/"${name}"db_backup
fi

# Move new db into main db folder
mv "${name}"db_"$MM_YYYY".* "${outdir}"/"${name}"db

# Generate all_pal_files.txt for psi-BLAST searches
cd "${outdir}"/"${name}"db || exit

cat "${name}"db_$MM_YYYY.pal | ${script_dir}split_pal.py

# Move back to the directory the database was started in
cd "${curr_dir}" || exit
