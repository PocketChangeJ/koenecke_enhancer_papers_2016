suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-t", "--treatment"), 
              type="character",
              default=NA,
              help="Treatment BigWig file"),
  make_option(c("-c", "--control"), 
              type="character",
              default=NA,
              help="Control BigWig file"),
  make_option(c("-o", "--output"),
              type="character",
              default=NA,
              help="Name of output BigWig file"),
  make_option(c("-w", "--window"),
              type="integer",
              default=101,
              help="Sliding window size"),
  make_option(c("-l", "--transform"),
              type="character",
              default="log2",
              help="Transform to apply to enrichment values (log2 or linear)"),
  make_option(c("--min"),
              type="integer",
              default=-2,
              help="Floor for enrichment values (after transform)"),
  make_option(c("--max"),
              type="integer",
              default=NA,
              help="Ceiling for enrichment values (after transform)"),
  make_option(c("-f", "--filter"),
              type="character",
              default=NA,
              help="Expression to filter out chromosomes")
  )

# OPTION PROCESSING

opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$treatment)) {
  message("No treatment file specified.")
  q(status=1)
}

if(is.na(opt$control)) {
  message("No control file specified.")
  q(status=1)
}

if(is.na(opt$output)) {
  message("No output name specified.")
  q(status=1)
}

if(! opt$transform %in% c("log2", "linear")) {
  message("Transform type must be either 'log2' or 'linear'")
  q(status=1)
}

suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F))

# COMMON FUNCTIONS

total_signal_rlelist <- function(cvg) {
  stopifnot(is(cvg, "RleList"))
  sum(as.numeric(sum(cvg)))
}

filter_chrs <- function(cov, match_expression="H|M|U|_") {
  exclude.chrs <- grep(match_expression, names(cov))
  if(length(exclude.chrs) > 0) cov <- cov[-exclude.chrs]
  cov
}

single_value_rlelist <- function(value, lengths) {
  RleList(lapply(lengths, function(i) Rle(rep(value, times=i))))
}

message("Loading: ", opt$treatment)
ip.cov <- import(opt$treatment, as="RleList")
message("Loading: ", opt$control)
bg.cov <- import(opt$control, as="RleList")

common.chrs <- intersect(names(ip.cov), names(bg.cov))
ip.cov <- ip.cov[common.chrs]
bg.cov <- bg.cov[common.chrs]

if(!is.na(opt$filter)) {
  ip.cov <- filter_chrs(ip.cov, opt$filter)
  bg.cov <- filter_chrs(bg.cov, opt$filter)
}

# START

message("Processing treatment...")
ip.cov <- runsum(ip.cov, opt$window, endrule="constant")
ip.ts  <- total_signal_rlelist(ip.cov)

message("Processing control...")
bg.cov <- runsum(bg.cov, opt$window, endrule="constant")
bg.ts  <- total_signal_rlelist(bg.cov)

message("Calculating enrichments...")

transform_function <- get(ifelse(opt$transform == "log2", "log2", "identity"))
e.cov <- transform_function((ip.cov / ip.ts) / (bg.cov / bg.ts))

if(!is.na(opt$max)) {
  message("Ceiling: ", opt$max)
  e.cov <- pmin(e.cov, single_value_rlelist(opt$max, elementLengths(e.cov)))
}

if(!is.na(opt$min)) {
  message("Floor: ", opt$min)
  e.cov <- pmax(e.cov, single_value_rlelist(opt$min, elementLengths(e.cov)))
}

message("Removing non-finite values...")
e.gr <- as(e.cov, "GRanges")
e.gr <- e.gr[is.finite(e.gr$score)]
e.gr$score <- round(e.gr$score, 2)
e.cov <- coverage(e.gr, weight="score")

message("Saving: ", opt$output)
export(e.cov, opt$output, format="BigWig")
