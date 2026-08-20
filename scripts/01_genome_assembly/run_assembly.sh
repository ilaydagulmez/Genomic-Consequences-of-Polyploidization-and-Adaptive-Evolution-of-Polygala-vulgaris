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
