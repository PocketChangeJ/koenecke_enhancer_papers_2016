
phastcon_scores_for_motif_instances <- function(motif_name, gr, phastcon.bw) {
  
  motifs.gr <- load_motif(motif_name)

  ol <- findOverlaps(motifs.gr, gr, type="within", ignore.strand=TRUE)

  hits.gr <- motifs.gr[unique(queryHits(ol))]
  strand(hits.gr) <- "*"
  hits.gr <- reduce(hits.gr)
  
  data_frame(motif=motif_name, phastcon_score=regionMeans(hits.gr, phastcon.bw))    
}

phastcon_scores_for_motif_instances_vs_genome <- function(motif_name, gr, phastcon.bw) {
  
  motifs.gr <- load_motif(motif_name)

  if(length(motifs.gr) > 5000) {
    sample.gr <- sample(motifs.gr, 5000)
  } else {
    sample.gr <- motifs.gr
  }

  ol <- findOverlaps(motifs.gr, gr, type="within", ignore.strand=TRUE)

  hits.gr <- motifs.gr[unique(queryHits(ol))]
  strand(hits.gr) <- "*"
  hits.gr <- reduce(hits.gr)
  
  hits.df <- data_frame(group="region", motif=motif_name, phastcon_score=regionMeans(hits.gr, phastcon.bw))
  genome.df <- data_frame(group="genome", motif=motif_name, phastcon_score=regionMeans(sample.gr, phastcon.bw))
  bind_rows(hits.df, genome.df)    
}

