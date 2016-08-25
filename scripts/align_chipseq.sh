#!/bin/bash

FASTQ=$1
BAM=/data/chipseq/bam/`basename $FASTQ fastq.gz`bam

bowtie -p 4 -S -m 1 -k 1 -v 2 --best --strata --chunkmbs=512 /data/chipseq/genome/dm3 <(zcat $FASTQ) | samtools view -F 4 -Sbo $BAM -

