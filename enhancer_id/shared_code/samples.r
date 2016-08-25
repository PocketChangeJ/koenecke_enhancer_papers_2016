library(readr)
library(yaml)

sample_bigwig <- function(sample_name, bigwig="ip", type="raw") {
  
  this_sample <- subset(samples.df, sample == sample_name)
  stopifnot(nrow(this_sample) == 1)
  
  if(bigwig == "ip")
    bigwig_file <- this_sample$bigwig
  else if(bigwig == "wce") 
    bigwig_file <- this_sample$wce_bigwig
  else if(bigwig == "enrichment")
    bigwig_file <- file.path("enrichment", this_sample$enrichment_bigwig)
  else stop("Unknown bigwig type: ", bigwig)
  
  sample_assay <- this_sample$assay
  
  if(type == "rpm") return(data_path(file.path(sample_assay, "bigwigs", "rpm", gsub("\\.bw$", "_rpm.bw", bigwig_file))))
  if(type == "raw") return(data_path(file.path(sample_assay, "bigwigs", bigwig_file)))
  stop("Unrecognized type: ", type)
}

samples.df <- read_csv("data/internal/samples.csv")

