#!/bin/bash

megahit \
    --no-mercy \
    -1 karect_FASTP_1.fastq \
    -2 karect_FASTP_2.fastq
    
SOAPdenovo-fusion \
    -D \
    -s config.txt \
    -p 64 \
    -K 21 \
    -g k21_1 \
    -c ./megahit_out/final.contigs.fa

SOAPdenovo-127mer map \
    -s config.txt \
    -p 64 \
    -g k21_1

SOAPdenovo-127mer scaff \
    -p 64 \
    -g k21_1 \
    -F

redundans.py \
    -v \
    -i karect_FASTP_1.fastq karect_FASTP_2.fastq \
    -f k21_1.scafSeq \
    --nogapclosing \
    -o redundans

# Output: redundans/scaffolds.reduced.fa

# Three successive polishing 

cat karect_FASTP_1.fastq karect_FASTP_2.fastq > karect_FASTP_1_2.fastq

python racon_preprocess.py \
    karect_FASTP_1_2.fastq \
    > final.racon.fastq

minimap2 \
    -ax sr \
    redundans/scaffolds.reduced.fa \
    final.racon.fastq \
    > racon_aln.sam

racon \
    -m 8 \
    -x -6 \
    -g -8 \
    -w 500 \
    -t 64 \
    final.racon.fastq \
    racon_aln.sam \
    redundans/scaffolds.reduced.fa \
    > polished1.fasta

minimap2 \
    -ax sr \
    polished1.fasta \
    final.racon.fastq \
    > 2_racon_aln.sam
    
racon \
    -m 8 \
    -x -6 \
    -g -8 \
    -w 500 \
    -t 64 \
    final.racon.fastq \
    2_racon_aln.sam \
    polished1.fasta \
    > polished2.fasta

minimap2 \
    -ax sr \
    polished2.fasta \
    final.racon.fastq \
    > 3_racon_aln.sam

racon \
    -m 8 \
    -x -6 \
    -g -8 \
    -w 500 \
    -t 64 \
    final.racon.fastq \
    3_racon_aln.sam \
    polished2.fasta \
    > polished3.fasta

ragtag.py scaffold \
    reference.fasta \
    polished3.fasta

# Final assembly: ragtag_output/ragtag.scaffold.fasta
