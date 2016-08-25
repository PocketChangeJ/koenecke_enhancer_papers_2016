library(GenomicRanges)

flybase_txs <- function() {
  txs <- readRDS("flybase/fb.transcripts.r5.57.rds")
  valid.chrs <- c("chr2L", "chr2LHet", "chr2R", "chr2RHet", "chr3L", "chr3LHet", "chr3R", "chr3RHet", "chr4", "chrX", "chrXHet", "chrYHet")
  txs <- subset(txs, chr %in% valid.chrs)
  txs
}

flybase_txs_granges <- function() {
  genes <- flybase_txs()
  genes.gr <- with(genes, GRanges(ranges     = IRanges(start=start, end=end), 
                                  seqnames   = chr,
                                  strand     = ifelse(strand == 1, "+", "-"), 
                                  fb_tx_id   = fb_tx_id,
                                  fb_gene_id = fb_gene_id, 
                                  fb_symbol  = fb_symbol,
                                  gene_type  = type))
  genes.gr
}

# Assigns the nearest gene by TSS (overlapping gene takes precedence over nearest TSS)
assign_nearest_gene <- function(gr, txs.gr=flybase_txs_granges()) {

  tss.gr <- resize(txs.gr, 1)

  dtn      <- as.data.frame(distanceToNearest(gr, tss.gr, ignore.strand=TRUE))
  dtn.gene <- subset(as.data.frame(distanceToNearest(gr, txs.gr, ignore.strand=TRUE)), distance == 0)

  mcols(gr)$nearest_gene_id <- ""
  mcols(gr)$nearest_gene    <- ""
  mcols(gr)$distance_to_tss <- NA
  mcols(gr)$inside_gene <- ""

  mcols(gr)$nearest_gene_id[dtn$queryHits] <- mcols(tss.gr)$fb_gene_id[dtn$subjectHits]
  mcols(gr)$nearest_gene[dtn$queryHits]    <- mcols(tss.gr)$fb_symbol[dtn$subjectHits]
  mcols(gr)$inside_gene[dtn$queryHits] <- FALSE
  mcols(gr)$distance_to_tss[dtn$queryHits] <- dtn$distance

  mcols(gr)$nearest_gene_id[dtn.gene$queryHits] <- mcols(txs.gr)$fb_gene_id[dtn.gene$subjectHits]
  mcols(gr)$nearest_gene[dtn.gene$queryHits]    <- mcols(txs.gr)$fb_symbol[dtn.gene$subjectHits]
  mcols(gr)$inside_gene[dtn.gene$queryHits] <- TRUE
  mcols(gr)$distance_to_tss[dtn.gene$queryHits] <- dtn$distance[match(dtn.gene$queryHits, dtn$queryHits)]

  gr
}

