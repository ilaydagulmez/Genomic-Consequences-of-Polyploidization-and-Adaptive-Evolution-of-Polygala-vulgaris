#!/bin/bash

wgd dmd genome_CDS.fa \
    -o dmd

wgd ksd \
    genome_CDS.tsv \
    genome_CDS.fa \
    -o ksd

wgd syn \
    -f transcript \
    -a ID \
    genome_CDS.tsv \
    genome_annotation.gff3 \
    -ks genome_CDS.tsv.ks.tsv \
    --pathiadhore ./i-adhore \
    -o syn

wgd peak \
    --heuristic \
    genome_CDS.tsv.ks.tsv \
    -ap ./syn/iadhore-out/anchorpoints.txt \
    -sm ./syn/iadhore-out/segments.txt \
    -le ./syn/iadhore-out/list_elements.txt \
    -mp ./syn/iadhore-out/multiplicon_pairs.txt \
    -n 1 4 \
    -kc 3 \
    -o wgd_peak
