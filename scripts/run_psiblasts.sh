#! /usr/bin/env bash

# Run psi-BLAST query(ies) against a database. Outputs into current directory as query_name.fasta.all.out
# Note: Will use fasta file basename (e.g., without .fasta) as the name of the query, so make sure that is either the
# name of the protein or the accession of the protein. If an accession is used, version will be cut off so protein name
# is preferred. (e.g. AFQ21949.1 -> AFQ21949)
#
# Two required arguments:
#   $1 must be the path to the database to be queried (MUST CONTAIN all_pal_names.txt)
#   $2 must be the path to the directory containing the protein fasta(s) to be used as query(ies)
# Two optional arguments:
#   $3 may be the evalue to be used (default 0.01)
#   $4 may be the number of psi-BLAST iterations (default 4)

EVAL=0.01
NUM_ITER=4

i=1
while [[ $# -gt 0 ]]; do
  case $i in
    1)
      DBDIR=$1
      shift
      ;;
    2)
      FASTADIR=$1
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
  esac
  i=$((i+1))
done

for file in $FASTADIR/*.fasta
do
        echo Psiblast Query: "$file"
        QUERY_NAME=$(basename "$file" | cut -d '.' -f 1)
        echo Query name: "$QUERY_NAME"
        while read -r db_name;
        do
                echo Querying: "$db_name"
                psiblast -query "$file" -db "$DBDIR/$db_name" -num_iterations "$NUM_ITER" -num_threads 4 -outfmt "6 qacc sacc pident qcovhsp length mismatch gapopen qstart qend sstart send qlen slen evalue"  -evalue "$EVAL" -out "$(basename "$file").out"
                cat "$(basename "$file").out" | grep -v '^$' | grep -v 'CONVERGED' >> "$(basename "$file")".all.out
        done < "$DBDIR"/all_pal_names.txt
        rm -f "$(basename "$file")".out
done
