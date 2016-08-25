#!/bin/bash

# 1: fastq

FASTQ=$1
OUTPUT_DIR=tophat/`basename $1 .fastq.gz`
SEGMENT_LENGTH=$((`zcat $FASTQ | head -4 | tail -1 | wc -c`/2))

PATH=/data/software/tophat:$PATH

tophat -p 8 -i 20 -I 5000 --no-coverage-search -z pigz \
       -o $OUTPUT_DIR --segment-length $SEGMENT_LENGTH \
       --transcriptome-index /data/flybase/transcriptome/fb_genes \
       /data/flybase/fb557 \
       $FASTQ

