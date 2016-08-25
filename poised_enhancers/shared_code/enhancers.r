library(dplyr)
library(magrittr)
library(readr)
library(GenomicRanges)

enhancers.df <- read_csv("data/internal/known_enhancers/TableS1_known_dv_enhancers_reassigned.csv")
enhancers.df$activity %<>% gsub("mesoderm", "m", .) %>%
                           gsub("dorsal ectoderm", "de", .)
names(enhancers.df)[which(names(enhancers.df) == "name")] <- "enhancer_name"

enhancers.gr <- makeGRangesFromDataFrame(enhancers.df, seqnames.field="chr", keep.extra=TRUE)

enrichment_values <- function(tissue, factor, replicate="best", enrichments.df) {
  tissue_ <- tissue
  factor_ <- factor
  replicate_ <- replicate

  if(replicate_ == "best") {
    replicate.df <- enrichments.df %>%
                    filter(tissue == tissue_ & factor == factor_) %>%
                    group_by(replicate) %>%
                    summarize(median_enrichment = median(enrichment)) %>%
                    ungroup %>%
                    summarize(replicate = replicate[which.max(median_enrichment)])
    replicate_ <- replicate.df$replicate[1]
  }
  
  values.df <- subset(enrichments.df, tissue == tissue_ & factor == factor_ & replicate == replicate_)
  values.df[, c("enhancer_name", "tissue", "factor", "replicate", "enrichment")]
}

