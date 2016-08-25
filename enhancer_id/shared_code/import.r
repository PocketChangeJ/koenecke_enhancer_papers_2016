library(GenomicRanges)
library(magrittr)

import.narrowPeak <- function(filename) {
  np <- read.delim(filename, header=FALSE, stringsAsFactors=FALSE) %>%
        with(GRanges(ranges=IRanges(start=V2+1, end=V3), seqnames=V1, score=V5, signal=V7, pvalue=V8, qvalue=V9, summit=V10))

  np
}
