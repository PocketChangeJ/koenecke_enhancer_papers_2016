#!/bin/bash

cat hot_regions.fasta | /data/software/meme/bin/fasta-get-markov > dm3_background.markov

parallel -uj 6 /data/software/meme/bin/fimo --bgfile dm3_background.markov --text motifs.meme {} \| gzip \> \`basename {.}\`_fimo.txt.gz ::: genome/*.fa
