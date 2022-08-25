#!/usr/bin/env bash

# Adds new fungal genomes to the existing fungal genome fasta file and /gbff folder. Makes a new blastp database using the updated fasta file.

if [[ $# -lt 1 ]]; then
  echo "No arguments provided, use [-h][--help] for usage help."
  exit 1
fi

curr_dir=$(pwd)
out_location=$(pwd)
no_backup="false"

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
      echo "  [-s][--db_scripts]        Path to directory containing db_scripts."
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

if [ -z ${db_scripts+x} ]; then
  echo "Must specify a path to the directory containing the database scripts."
  exit 1
else
  script_dir=$(cd "$db_scripts" || exit; pwd; cd "$curr_dir" || exit;)
  script_dir="$script_dir/"
fi

if [ -z ${link+x} ]; then
  echo "Must specify a link to download the database from."
  exit 1
fi

if [ -z ${name+x} ]; then
  name=$(echo "$link" | rev | cut -d "/" -f2 | rev);
fi

out_location=$(cd "$out_location" || exit; pwd; cd "$curr_dir" || exit);

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

MM_YYYY=$(date +%m_%Y)

${script_dir}get_list_of_new_ncbi_genomes.py "$link" "$existing" | grep -v 'na/na' | sort | uniq > ${outdir}/new_"${name}"_genomes_urls.txt

mkdir ${outdir}/new_${name}_genomes_$$
cd ${outdir}/new_${name}_genomes_$$

for url in `cat ${outdir}/new_"${name}"_genomes_urls.txt`; do bn=`basename $url`; if [ ! -e $bn ]; then echo "Downloading URL $url"; wget -t 1 --timeout=5 $url; else echo "Already downloaded $url"; fi; done

# Check for missing genomes
for genome in $(cat ${outdir}/new_"${name}"_genomes_urls.txt | cut -d "/" -f 11 | sort | uniq); do if [ ! -e $genome ]; then echo $genome; fi; done

# Check for messed-up GBFF files
# If there are some, they need to be re-downloaded.
# This needs to be done before creating the locus map, otherwise it will give an error.
#../check_gbff.gz_files.py . > ../messed_up_BLAH_gbff_files.txt

echo Making new_"${name}"_proteins.FASTA
${script_dir}make_gbff_protein_fasta.py . > ${outdir}/new_${name}_proteins.FASTA
echo Updating all_${name}_proteins.FASTA
mkdir -p ${outdir}/${name}_genomes/proteins
mkdir -p ${outdir}/${name}_genomes/gbff
cat ${outdir}/new_${name}_proteins.FASTA >> ${outdir}/${name}_genomes/proteins/all_${name}_proteins.FASTA
rm -f ${outdir}/new_${name}_proteins.FASTA
for thing in `find . -name "*.gbff.gz"`; do mv $thing ${outdir}/${name}_genomes/gbff; done
rmdir new_${name}_genomes_*

cd ${outdir}/${name}_genomes/proteins
makeblastdb -in ${outdir}/${name}_genomes/proteins/all_${name}_proteins.FASTA -dbtype prot -title ${name}db_$MM_YYYY  -out ${name}db_$MM_YYYY

mkdir -p ${outdir}/${name}db
if [[ "$no_backup" != "true" ]]; then
  # Backup things
  mkdir -p ${outdir}/${name}db_backup
  rm -rf ${outdir}/${name}db_backup/*
  mv ${outdir}/${name}db/* ${outdir}/${name}db_backup
fi

# Move new db into main db folder
mv ${name}db_$MM_YYYY.* ${outdir}/${name}db

# Generate all_pal_files.txt for psi-BLAST searches
cd ${outdir}/${name}db

cat ${name}db_$MM_YYYY.pal | ${script_dir}split_pal.py

cd ${curr_dir}
# If, in the future, the viral database gets large enough, may need to add split_pal.py into the pipeline. For now, database is small enough that no .pal files are generated.
# Don't do this because it takes a long time and don't think it does anything... Was copy pasted from update_ncbi_bacterial_genomes.sh. Leaving the code in case someone finds it useful in the future.
# Make new locus to filename map
#cd /birl2/users/mdt442/data/bacterial_genomes/gbff
#create_locus_map.py /birl2/users/mdt442/bacterial_genomes/gbff > /birl2/users/mdt442/acr_aca_misc_data/locus_to_filename_map.txt
