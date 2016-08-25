#!/bin/bash

# $1: tophat output directory

CQ=/data/software/cufflinks/cuffquant

$CQ -p 3 -o cuffquant/`basename $1` -b /data/flybase/fb557.fa -u /data/flybase/transcriptome/fb_genes.gff $1/accepted_hits.bam

