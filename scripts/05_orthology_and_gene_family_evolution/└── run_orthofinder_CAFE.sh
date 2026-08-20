#!/bin/bash

set -euo pipefail

PROTEOME_DIR="proteomes"
THREADS=32

orthofinder \
    -f "$PROTEOME_DIR" \
    -M msa \
    -t "$THREADS" \
    -a "$THREADS"

GENECOUNT="Orthogroups.GeneCount.tsv"
CAFE_INPUT="gene_families_CAFE.tsv"

awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    printf "Desc\tFamily ID"
    for(i=2;i<NF;i++)
        printf "\t%s", $i
    printf "\n"
    next
}
{
    printf "null\t%s", $1
    for(i=2;i<NF;i++)
        printf "\t%s", $i
    printf "\n"
}' "$GENECOUNT" > "$CAFE_INPUT"

FILTERED_CAFE_INPUT="gene_families_CAFE_filtered.tsv"

awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    print
    next
}
{
    max=0
    sum=0

    for(i=3;i<=NF;i++) {

        if($i > max)
            max=$i

        sum += $i
    }

    if(max <= 120 && sum <= 500)
        print
}' "$CAFE_INPUT" > "$FILTERED_CAFE_INPUT"

TREE="species_tree_ultrametric.nwk"

cafe5 \
    -i "$FILTERED_CAFE_INPUT" \
    -t "$TREE" \
    -c "$THREADS" \
    -P 0.05 \
    -o CAFE5_results

echo "Orthology and CAFE5 analyses completed."
