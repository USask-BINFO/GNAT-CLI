#! /usr/bin/env bash

# Run psi-BLAST query(ies) against a database. Outputs into current directory as query_name.fasta.all.out
# Note: Will use fasta file basename (e.g., without .fasta) as the name of the query, so make sure that is either the
# name of the protein or the accession of the protein. If an accession is used, version will be cut off so protein name
# is preferred. (e.g. AFQ21949.1 -> AFQ21949)
#
# Two required arguments:
#   $1 must be the path to the database to be queried (MUST CONTAIN all_pal_names.txt)
#   $2 must be the path to the fasta file containing the protein sequence(s) to be used as query(ies)
# Two optional arguments:
#   $3 may be the evalue to be used (default 0.01)
#   $4 may be the number of psi-BLAST iterations (default 4)

EVAL=0.01
NUM_ITER=4
OUT_TAG=""

i=1
while [[ $# -gt 0 ]]; do
  case $i in
    1)
      DBDIR=$1
      shift
      ;;
    2)
      FASTAFILE=$1
      shift
      ;;
    3)
      EVAL=$1
      shift
      ;;
    4)
      NUM_ITER=$1
      shift
      ;;
    5)
      OUT_TAG=$1
      shift
      ;;
  esac
  i=$((i+1))
done

echo "[ $OUT_TAG ] PSI-BLAST Query: $(basename $FASTAFILE)"
touch outfile.out
while read -r db_name;
do
  echo "[ $OUT_TAG ] Querying: $db_name"
  psiblast -query "$FASTAFILE" -db "$DBDIR/$db_name" -num_iterations "$NUM_ITER" -num_threads 4 -outfmt "6 qacc sacc pident qlen slen evalue" -evalue "$EVAL" >> outfile.out;
done < "$DBDIR"/all_pal_names.txt

filename=$(basename "$FASTAFILE")
no_ext="${filename%.*}"

cat outfile.out | grep -v '^$' | grep -v 'CONVERGED' > "$no_ext".out
rm -f outfile.out
