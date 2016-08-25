library(compiler)

build_overlap_table <- function(motif1_name, motif2_name, gr) {

  motif1.gr <- load_motif(motif1_name)
  motif2.gr <- load_motif(motif2_name)

  ol1 <- findOverlaps(motif1.gr, gr, type="within", ignore.strand=TRUE)
  ol2 <- findOverlaps(motif2.gr, gr, type="within", ignore.strand=TRUE)

  hits1.gr <- motif1.gr[unique(queryHits(ol1))]
  hits2.gr <- motif2.gr[unique(queryHits(ol2))]
  
  strand(hits1.gr) <- "*"
  strand(hits2.gr) <- "*"
  
  hits1.gr <- reduce(hits1.gr)
  hits2.gr <- reduce(hits2.gr)
  
  total_m1 <- length(hits1.gr)
  total_m2 <- length(hits2.gr)
  
  ol <- findOverlaps(query=hits1.gr, subject=hits2.gr)
  
  overlap_m1 <- length(unique(queryHits(ol)))
  overlap_m2 <- length(unique(subjectHits(ol)))
  
  data_frame(motif1 = motif1_name, 
             motif2 = motif2_name,
             m1_total = total_m1,
             m2_total = total_m2,
             overlap_m1  = overlap_m1,
             overlap_m2  = overlap_m2,
             m1_percent = overlap_m1 / total_m1 * 100,
             m2_percent = overlap_m2 / total_m2 * 100)
}
build_overlap_table %<>% cmpfun

overlap_table_for_motifs <- function(motifs, peaks.gr) {
  pairs.df <- t(combn(motifs, 2)) %>%
              as.data.frame(stringsAsFactors=FALSE)
  names(pairs.df) <- c("motif1", "motif2")
  
  message(pn(nrow(pairs.df)), " overlaps to test")
  
  1:nrow(pairs.df) %>%
    checked_mclapply(function(i) {
      if(i %% 1000 == 0) message(i)
      build_overlap_table(pairs.df$motif1[i], pairs.df$motif2[i], peaks.gr)
    }, mc.cores=cores()) %>%
    bind_rows
}

group_by_overlap <- function(sig.df, overlaps.df, percent=50) {
  
  # start with most significant motif
  sig.df <- arrange(sig.df, adj_pv)
  
  results.df <- data_frame(keep_motif="none", child_motif="none")
  
  for(i in unique(sig.df$motif)) {
    if(i %in% results.df$child_motif) next
    overlapping_motifs <- c(subset(overlaps.df, motif1 == i & overlap > percent)$motif2,
                            subset(overlaps.df, motif2 == i & overlap > percent)$motif1) %>%
                          unique
    overlapping_motifs <- overlapping_motifs[!overlapping_motifs %in% results.df$child_motif]
    if(length(overlapping_motifs) > 0) {
      results.df %<>% rbind(data_frame(keep_motif=i, child_motif=overlapping_motifs))
    } else {
      results.df %<>% rbind(data_frame(keep_motif=i, child_motif=i))
    }
  }
  results.df <- subset(results.df, keep_motif != "none")
  
  results.sum <- results.df %>%
                 group_by(keep_motif) %>%
                 summarize(other_names = paste(unique(child_motif), collapse=", "))

  list(df=as.data.frame(results.sum),
       flat=results.df,
       filtered = subset(sig.df, motif %in% results.sum$keep_motif))
}

