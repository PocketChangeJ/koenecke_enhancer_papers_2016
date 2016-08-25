#!/bin/bash

FASTQ=$1
BAM=/data/chipseq/bam/new/`basename $FASTQ fastq.gz`bam

bowtie -p 7 -S -m 1 -k 1 -v 2 --best --strata --chunkmbs=512 /data/chipseq/genome/dm3 <(zcat $FASTQ | fastx_trimmer -Q 33 -l 50) | samtools view -F 4 -Sbo $BAM -

