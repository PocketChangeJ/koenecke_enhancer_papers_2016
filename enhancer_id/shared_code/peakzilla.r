library(GenomicRanges)

import.peakzilla <- function(f) {
  f.df <- read.delim(f, stringsAsFactors=FALSE, header=TRUE)
  with(f.df, GRanges(ranges     = IRanges(start=Start, end=End), 
                     seqnames   = X.Chromosome,
                     score      = Score,
                     summit     = Summit,
                     enrichment = FoldEnrichment,
                     fdr        = FDR))
}
