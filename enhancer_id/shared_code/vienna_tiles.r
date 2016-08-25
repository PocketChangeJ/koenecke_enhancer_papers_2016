library(GenomicRanges)
library(dplyr)
source("shared_code/stat_tests.r")

vt.df <- readRDS(data_path("stark_enhancers_2014/tile_annotations.df.rds"))
vt.gr <- makeGRangesFromDataFrame(vt.df, seqnames.field="chr", keep.extra=TRUE)
all_tiles.df <- readRDS(data_path("stark_enhancers_2014/all_tiles.df.rds"))
all_tiles.gr <- makeGRangesFromDataFrame(all_tiles.df, seqnames.field="chr", keep.extra=TRUE)

stage_enrichments_for_tiles <- function(tile_ids, group_name, annotations.df, tile_universe) {
  all_stages <- sort(unique(annotations.df$stage))
  
  results.df <- all_stages %>%
                lapply(function(s) {
                  tiles_in_stage <- subset(annotations.df, stage == s)$VTID %>% unique()
                  test.df <- fisher_test_2x2(tile_ids, tiles_in_stage, tile_universe)
                  test.df$stage <- s
                  test.df$group_name <- group_name
                  test.df
                }) %>%
                bind_rows()
}

term_enrichments_for_tiles <- function(tile_ids, group_name, annotations.df, tile_universe) {
  all_terms <- sort(unique(annotations.df$annotation))
  
  results.df <- all_terms %>%
                mclapply(function(term) {
                  tiles_with_term <- subset(annotations.df, annotation == term)$VTID %>% unique()
                  test.df <- fisher_test_2x2(tile_ids, tiles_with_term, tile_universe)
                  test.df$term <- term
                  test.df$group_name <- group_name
                  test.df
                }, mc.cores=6, mc.preschedule=FALSE) %>%
                bind_rows()
}

overlapping_tile_ids <- function(gr, tiles.gr) {
  tiles.gr[countOverlaps(tiles.gr, gr, ignore.strand=TRUE) > 0]$VTID %>% unique()
}
