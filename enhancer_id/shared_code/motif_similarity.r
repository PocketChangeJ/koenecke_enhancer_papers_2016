library(PWMEnrich)
library(MotifDb)

add_mirrored_similarities <- function(df) {
  df2 <- df
  names(df2)[1:2] <- rev(names(df2)[1:2])
  rbind(df, df2)
}

cleanup_names <- function(x) {
  gsub("Dmelanogaster-", "", x)
}

convert_to_integer_matrix <- function(m) {
  m <- m * 1000
  mode(m) <- "integer"
  m
}

motif_sim_helper <- function(i, combos.df, motifs) {
  if(i %% 10000 == 0) message(i)
  m1 <- convert_to_integer_matrix(motifs[[combos.df$motif1[i]]])
  m2 <- convert_to_integer_matrix(motifs[[combos.df$motif2[i]]])
  motifSimilarity(m1, m2)
}

similarity_df <- function(motif_names, cores=3) {
  dmel <- query(MotifDb, "Dmelanogaster")
  names(dmel) <- cleanup_names(names(dmel))
  
  combos.df <- as.data.frame(matrix(combn(motif_names, 2, function(x) { matrix(x, nrow=1, ncol=2) } ), ncol=2, byrow=T))
  names(combos.df) <- c("motif1", "motif2")
  combos.df$motif1 <- as.character(combos.df$motif1)
  combos.df$motif2 <- as.character(combos.df$motif2)
  
  message("Total: ", nrow(combos.df))
  combos.df$similarity <- unlist(mclapply(1:nrow(combos.df), motif_sim_helper, combos.df, dmel, mc.cores=cores), use.names=FALSE)
  add_mirrored_similarities(combos.df)
}

