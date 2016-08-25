library(GenomicRanges)

import.dnase <- function(filename) {
  stopifnot(file.exists(filename))
  dnase.df <- read.delim(filename, stringsAsFactors=FALSE, header=FALSE)
  makeGRangesFromDataFrame(dnase.df, seqnames.field="V2", start.field="V3", end.field="V4")
}

import.macs_peak <- function(filename) {
  stopifnot(file.exists(filename))
  peaks.df <- read.delim(filename, stringsAsFactors=FALSE, header=FALSE)
  with(peaks.df, GRanges(ranges=IRanges(start=V2, end=V3), seqnames=V1, score=V5, name=V4))
}

import.modencode_gff <- function(filename) {
  gr <- import(filename)
  mito.chr <- grep("mito", seqlevels(gr))
  if(length(mito.chr)) seqlevels(gr, force=TRUE) <- seqlevels(gr)[-mito.chr]
  seqlevels(gr) <- paste0("chr", seqlevels(gr))
  filter_chrs(gr)
}

import.peakzilla <- function(f) {
  f.df <- read.delim(f, stringsAsFactors=FALSE, header=TRUE)
  with(f.df, GRanges(ranges     = IRanges(start=Start, end=End), 
                     seqnames   = X.Chromosome,
                     score      = Score,
                     summit     = Summit,
                     enrichment = FoldEnrichment,
                     fdr        = FDR))
}
