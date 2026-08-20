#!/bin/bash

SAMPLE=$1
CDS=$2
GFF=$3

mkdir -p "${SAMPLE}_WGD"

wgd dmd \
    "${CDS}" \
    -o "${SAMPLE}_WGD/dmd"

wgd ksd \
    "${SAMPLE}_WGD/dmd/${SAMPLE}_CDS.tsv" \
    "${CDS}" \
    -o "${SAMPLE}_WGD/ksd"
    
wgd syn \
    -f transcript \
    -a ID \
    "${SAMPLE}_WGD/dmd/${SAMPLE}_CDS.tsv" \
    "${GFF}" \
    -ks "${SAMPLE}_WGD/ksd/${SAMPLE}_CDS.tsv.ks.tsv" \
    --pathiadhore ./i-adhore \
    -o "${SAMPLE}_WGD/syn"
    
wgd peak \
    --heuristic \
    "${SAMPLE}_WGD/ksd/${SAMPLE}_CDS.tsv.ks.tsv" \
    -ap "${SAMPLE}_WGD/syn/iadhore-out/anchorpoints.txt" \
    -sm "${SAMPLE}_WGD/syn/iadhore-out/segments.txt" \
    -le "${SAMPLE}_WGD/syn/iadhore-out/list_elements.txt" \
    -mp "${SAMPLE}_WGD/syn/iadhore-out/multiplicon_pairs.txt" \
    -n 1 4 \
    -kc 3 \
    -o "${SAMPLE}_WGD/wgd_peak"

echo "WGD and Ks analysis completed for ${SAMPLE}"
