#!/bin/bash

set -euo pipefail

SAMPLE=$1
PROTEINS=$2

THREADS=${THREADS:-32}

INTERPROSCAN=${INTERPROSCAN:-interproscan.sh}

"$INTERPROSCAN" \
    -i "$PROTEINS" \
    -f TSV \
    -appl Pfam \
    --iprlookup \
    --goterms \
    --pathways \
    -cpu "$THREADS" \
    -o "${SAMPLE}_interproscan.tsv"

awk -F'\t' '
BEGIN {
    OFS="\t"
}

NF >= 14 && $14 != "-" && $14 != "" {

    protein=$1

    n=split($14,go,"|")

    for(i=1;i<=n;i++) {

        gsub(/^ +| +$/, "", go[i])

        if(go[i] ~ /^GO:[0-9]+/)
            print protein,go[i]
    }
}
' "${SAMPLE}_interproscan.tsv" \
| sort -u \
> "${SAMPLE}_gene_to_GO.tsv"
