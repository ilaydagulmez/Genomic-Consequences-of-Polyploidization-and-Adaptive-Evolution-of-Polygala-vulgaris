#!/bin/bash

Helixer.py \
    --fasta-path genome.fasta \
    --lineage land_plant \
    --gff-output-path genome_helixer.gff3

gffread genome_helixer.gff3 \
    -g genome.fasta \
    -x genome_CDS.fa \
    -y genome_protein.fa \
    -w genome_transcripts.fa
