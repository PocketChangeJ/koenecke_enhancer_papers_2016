suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-b", "--bigwig"), 
              type="character",
              default=NA,
              help="Input BigWig file"),
  make_option(c("-r", "--reads"),
              type="integer",
              default=1,
              help="Read count normalization (in millions, default 1)"),
  make_option(c("-s", "--size"),
              type="integer",
              default=150,
              help="Fragment size normalization (default 150 bp)")
)

# OPTION PROCESSING

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$bigwig)) {
  message("No BigWig file specified.")
  q(status=1)
}

suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F))

message("Loading: ", opt$bigwig)
cov <- import(opt$bigwig, as="RleList")
cov <- cov / sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x))))) * opt$reads * 1e6 * opt$size

new_bigwig <- gsub("\\.bw$", "_rpm.bw", basename(opt$bigwig))

message("Writing: ", new_bigwig)
export(cov, new_bigwig)

