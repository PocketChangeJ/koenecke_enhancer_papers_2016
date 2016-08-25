suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-f", "--file"), 
              type="character",
              default=NA,
              help="Path of BAM file to process"),
  make_option(c("-e", "--extension"),
              type="character",
              default="native",
              help="Extension length ('native' for no extension, 'auto' to estimate fragment size)")
  )

readBAM <- function(filepath) {
  bam.gr <- granges(readGAlignments(filepath))
  bam.gr
}

write_bigwig <- function(cov, filename) {
  message(id, "Writing bigWig...")
  export(cov, filename)
}

filter_chrs <- function(gr) {
  exclude.chrs <- grep("H|M|U", seqlevels(gr))
  if(length(exclude.chrs) > 0) seqlevels(gr, force=TRUE) <- seqlevels(gr)[-exclude.chrs]
  gr
}

pn <- function(value) {
  prettyNum(value, big.mark=",")
}


# ------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$file)) {
  message("No BAM file specified. Use --help to list available options.")
  q(status=1)
}

suppressPackageStartupMessages(library(Rsamtools, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(chipseq, warn.conflicts=F))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F))
suppressPackageStartupMessages(library(GenomicAlignments))

# used as a prefix for all output messages
id <- "[unknown] "

bam_file      <- opt$file
ext_length    <- opt$extension

var_name <- gsub("\\.bam$", "", bam_file)

id <- paste("[", var_name, "] ", sep="")

bigwig_file   <- paste0(var_name, ".bw")
granges_file <- paste0(var_name, ".granges.rds")

message(id, "Converting BAM to GRanges object:")
message(id, "Input BAM: ", bam_file)
message(id, "Output: ", granges_file)

if(!file.exists(bam_file)) {
	stop("Could not open BAM file: ", bam_file)
}

bam.gr <- readBAM(bam_file)

message(id, "Estimating fragment length...")
estimates <- estimate.mean.fraglen(filter_chrs(bam.gr), method="coverage")
counts <- table(seqnames(bam.gr))
est.frag.size <- as.integer(weighted.mean(estimates, counts[names(estimates)]))

if(ext_length == "native") {
	ext_length <- NULL
} else {
	if(ext_length == "auto")
	  ext_length <- est.frag.size
	else
	  ext_length <- as.integer(ext_length)
}

message(id, "Extension length: ", ifelse(is.null(ext_length), "native", ext_length))

if(!is.null(ext_length)) bam.gr <- trim(resize(bam.gr, ext_length))

message(id, "Saving GRanges object...")
bam.gr <- bam.gr[order(bam.gr)]
saveRDS(bam.gr, file=granges_file)

nothing <- write_bigwig(coverage(bam.gr), bigwig_file)



