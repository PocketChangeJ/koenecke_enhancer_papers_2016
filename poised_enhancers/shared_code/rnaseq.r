source("shared_code/cuffdiff_common.r")

rnaseq_results <- function() {
  rnaseq.df <- read_cuffdiff(data_path("rnaseq/cuffdiff_gd_toll/gene_exp.diff")) %>%
                             filter(sample_1 == "gd" & sample_2 == "toll")

  # Add CG8788 as a duplicate of CG44286, a gene with identical transcripts
  
  row.df <- subset(rnaseq.df, fb_symbol == "CG44286")
  row.df$fb_symbol <- "CG8788"
  row.df$fb_gene_id <- "FBgn0028955"
  bind_rows(rnaseq.df, row.df)
}
