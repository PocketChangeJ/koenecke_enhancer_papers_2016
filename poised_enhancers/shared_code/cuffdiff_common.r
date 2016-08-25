library(dplyr)

source("shared_code/flybase.r")

read_cuffdiff <- function(filename) {
  cuffdiff <- read.delim(filename, stringsAsFactors=FALSE, header=TRUE)
  
  cuffdiff <- cuffdiff[, c("gene_id", "sample_1", "sample_2", "status", "value_1", "value_2", "log2.fold_change.", "q_value", "significant")]
  names(cuffdiff) <- c("fb_gene_id", "sample_1", "sample_2", "status", "RPKM_1", "RPKM_2", "fold_change", "q_value", "significant")

  cuffdiff <- transform(cuffdiff, direction = ifelse(fold_change < 0, "down", "up"))
  cuffdiff <- transform(cuffdiff, significant = ifelse(significant == "yes", TRUE, FALSE))
 
  genes.df <- unique(flybase_txs()[, c("fb_gene_id", "fb_symbol")])
  cuffdiff <- merge(genes.df, cuffdiff, all.y=TRUE)
  cuffdiff
}

read_fpkms <- function(filename, match_str) {
  fpkms <- read.delim(filename, stringsAsFactors=FALSE, header=TRUE)
  gene_id_col <- grep("gene_id", names(fpkms))
  sample_cols <- grep(match_str, names(fpkms))
  fpkms <- fpkms[, c(gene_id_col, sample_cols)]
  names(fpkms)[1] <- "fb_gene_id"
  genes.df <- unique(flybase_txs()[, c("fb_gene_id", "fb_symbol")])
  fpkms <- merge(genes.df, fpkms, all.y=TRUE)
  fpkms
}

summarize_results <- function(cuffdiff.df) {
  summary.df <- cuffdiff.df %>% 
                  group_by(sample_1, sample_2) %>%
                  summarize(up_regulated_count   = length(which(significant & direction == "up")), 
                            down_regulated_count = length(which(significant & direction == "down")))
  summary.df
}



