library(dplyr)
library(tidyr)

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

rnaseq_diff_exp_results <- function() {
  rbind(read_cuffdiff(data_path(file.path("rnaseq", "cuffdiff_gd_toll", "gene_exp.diff"))),
        read_cuffdiff(data_path(file.path("rnaseq", "cuffdiff_gd_rm", "gene_exp.diff"))),
        read_cuffdiff(data_path(file.path("rnaseq", "cuffdiff_rm_toll", "gene_exp.diff")))) %>%
    filter(fb_gene_id != "-")
}

rnaseq_fpkms <- function() {
  fpkms1.df <- read.delim(data_path(file.path("rnaseq", "cuffdiff_gd_toll", "genes.fpkm_tracking")), stringsAsFactors=FALSE) %>%
               select(gene_id, gd_FPKM, toll_FPKM, gd_FPKM) %>%
               gather(tissue, FPKM, toll_FPKM:gd_FPKM) %>%
               mutate(tissue     = gsub("_FPKM", "", tissue),
                      fb_gene_id = gene_id) %>%
               select(-gene_id) %>%
               filter(fb_gene_id != "-")

  fpkms2.df <- read.delim(data_path(file.path("rnaseq", "cuffdiff_gd_rm", "genes.fpkm_tracking")), stringsAsFactors=FALSE) %>%
               select(gene_id, rm_FPKM) %>%
               gather(tissue, FPKM, rm_FPKM) %>%
               mutate(tissue     = gsub("_FPKM", "", tissue),
                      fb_gene_id = gene_id) %>%
               select(-gene_id) %>%
               filter(fb_gene_id != "-")
  
  rbind(fpkms1.df, fpkms2.df)
}

genes_expressed_in_tissue <- function(tissue_names, fpkm_threshold=5) {
  fpkms.df <- rnaseq_fpkms()
  subset(fpkms.df, FPKM > fpkm_threshold & tissue %in% tissue_names)$fb_gene_id %>% unique
}
