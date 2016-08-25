#!/bin/bash

FASTQ_R1=fastq/$1_read1.fastq.gz
FASTQ_R2=fastq/$1_read2.fastq.gz
BAM=bam/$1.bam

bowtie -p 7 -S -m 1 -k 1 -v 2 -X 1000 \
	--best --strata --chunkmbs=512 \
	/data/chipseq/genome/dm3 \
	-1 <(zcat $FASTQ_R1 | fastx_trimmer -Q 33 -l 25 $FASTQ_R1) \
	-2 <(zcat $FASTQ_R2 | fastx_trimmer -Q 33 -l 25 $FASTQ_R2) | samtools view -F 4 -Sbo $BAM -


