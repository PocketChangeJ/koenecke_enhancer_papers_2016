#!/bin/bash

CUFFLINKS_PATH=/data/software/cufflinks
CORES=5


$CUFFLINKS_PATH/cuffdiff -p $CORES --use-sample-sheet -o cuffdiff_gd_toll -M /data/flybase/mito_genes.gff /data/flybase/transcriptome/fb_genes.gff sample_sheet_gd_toll.txt

$CUFFLINKS_PATH/cuffdiff -p $CORES --use-sample-sheet -o cuffdiff_gd_rm -M /data/flybase/mito_genes.gff /data/flybase/transcriptome/fb_genes.gff sample_sheet_gd_rm.txt

$CUFFLINKS_PATH/cuffdiff -p $CORES --use-sample-sheet -o cuffdiff_rm_toll -M /data/flybase/mito_genes.gff /data/flybase/transcriptome/fb_genes.gff sample_sheet_rm_toll.txt

