library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(lattice)

read_matrix <- function(gr, cov, reverse_reads=FALSE) {
  if(class(cov) == "character") {
    cov <- import.bw(cov, which=gr, as="RleList")
  }
  
  transform_function <- if(reverse_reads) { rev } else { identity }
  o <- order(gr)
  gr <- gr[o]
  rl <- as(gr, "RangesList")
  view <- RleViewsList(rleList=cov[names(rl)], rangesList=rl)
  reads.list <- viewApply(view, function(x) { transform_function(as.numeric(x)) })
  reads.m <- matrix(unlist(sapply(reads.list, as.numeric)), nrow=length(gr), byrow=TRUE)
  reads.m[o, ] <- reads.m
  reads.m
}

standard_metapeak_matrix <- function(regions.gr, sample.cov, upstream=100, downstream=100) {
  
  regions.gr <- resize(regions.gr, width=downstream)
  regions.gr <- trim(resize(regions.gr, width=upstream + width(regions.gr), fix="end"))

  if(class(sample.cov) == "character") {
    stopifnot(file.exists(sample.cov))
    cov_seqlengths <- seqlengths(BigWigFile(sample.cov))
  } else {
    cov_seqlengths <- elementLengths(sample.cov)
  }

  if(anyNA(seqlengths(regions.gr))) {
    seqlengths(regions.gr) <- cov_seqlengths[seqlevels(regions.gr)]
  }

  remove_i <- which(width(regions.gr) != upstream + downstream)
  if(length(remove_i) > 0) {
    warning("Removing ", length(remove_i), " range(s() from `regions.gr` due to chromosome boundaries")
    regions.gr <- regions.gr[-remove_i]
  }
  
  stopifnot(length(regions.gr) > 0)
  stopifnot(length(unique(width(regions.gr))) == 1)
  
  reads <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
  
  i_p <- which(strand(regions.gr) == "+" | strand(regions.gr) == "*")
  i_n <- which(strand(regions.gr) == "-")
  
  if(class(sample.cov) == "character") {
    sample.cov <- import.bw(sample.cov, which=regions.gr, as="RleList")
  }
  
  if(length(i_p) > 0) reads[i_p, ] <- read_matrix(regions.gr[i_p], sample.cov)
  if(length(i_n) > 0) reads[i_n, ] <- read_matrix(regions.gr[i_n], sample.cov, reverse_reads=TRUE)

  reads
}

standard_metapeak <- function(gr, sample.cov, upstream=100, downstream=100, sample_name=NA, smooth=NA, normalize=NA) {
  
  reads <- standard_metapeak_matrix(gr, sample.cov, upstream, downstream)
  
  reads.df <- data.frame(stringsAsFactors=FALSE,
                         tss_distance=(-1 * upstream):(downstream - 1),
                         reads=colMeans(reads, na.rm=TRUE), 
                         sample_name=sample_name)
  
  if(!is.na(normalize)) reads.df$reads <- reads.df$reads / total_signal(sample.cov) * normalize
  
  if(!is.na(smooth)) reads.df$reads <- as.numeric(runmean(Rle(reads.df$reads), k=smooth, endrule="constant"))
  
  reads.df  
}

multisample_standard_metapeak <- function(grl, samples, upstream=100, downstream=100, smooth=NA, cores=1, normalize=NA) {
  if(!is.list(grl)) grl <- list("Unnamed"=grl)

  names(samples) %>%
    mclapply(function(sample_name) {
      names(grl) %>%
      lapply(function(peak_group_name) {
          gr <- grl[[peak_group_name]]
          reads.df <- standard_metapeak(gr, samples[[sample_name]], upstream=upstream, downstream=downstream, smooth=smooth, sample_name=sample_name, normalize=normalize)
          reads.df$peak_group <- peak_group_name
          reads.df
      }) %>%
      bind_rows
    }, mc.cores=cores, mc.preschedule=FALSE) %>%
    bind_rows
}

enrichment_metapeak <- function(gr, sample.cov, bg.cov, upstream=100, downstream=100, sample_name=NA, smooth=NA) {
  
  reads.ip <- standard_metapeak_matrix(gr, sample.cov, upstream, downstream)
  reads.bg <- standard_metapeak_matrix(gr, bg.cov, upstream, downstream)
  
  ts.ip <- total_signal(sample.cov)
  ts.bg <- total_signal(bg.cov)
  
  reads.df <- data.frame(stringsAsFactors=FALSE,
                         tss_distance=(-1 * upstream):(downstream - 1),
                         reads=colMeans(reads.ip, na.rm=TRUE) / ts.ip, 
                         background=colMeans(reads.bg, na.rm=TRUE) / ts.bg,
                         sample_name=sample_name)

  reads.df$enrichment <- reads.df$reads / reads.df$background

  if(!is.na(smooth)) {
    reads.df$reads      <- as.numeric(runmean(Rle(reads.df$reads), k=smooth, endrule="constant"))
    reads.df$background <- as.numeric(runmean(Rle(reads.df$background), k=smooth, endrule="constant"))
    reads.df$enrichment <- reads.df$reads / reads.df$background
  }
  
  reads.df  
}

multisample_enrichment_metapeak <- function(grl, samples, upstream=100, downstream=100, smooth=NA, cores=1) {
  if(!is.list(grl)) grl <- list("Unnamed"=grl)

  names(samples) %>%
  mclapply(function(sample_name) {
    lapply(names(grl), function(peak_group_name) {
      gr <- grl[[peak_group_name]]
      reads.df <- enrichment_metapeak(gr, samples[[sample_name]]$ip, samples[[sample_name]]$wce, 
                                      upstream=upstream, downstream=downstream, smooth=smooth, sample_name=sample_name)
      reads.df$peak_group <- peak_group_name
      reads.df
    }) %>% 
    bind_rows
  }, mc.cores=cores, mc.preschedule=FALSE) %>%
  bind_rows
}

base_frequencies <- function(gr, upstream, downstream, genome=Dmelanogaster) {

  gr <- resize(gr, width=downstream)
  gr <- resize(gr, width=upstream + width(gr), fix="end")

  m <- consensusMatrix(getSeq(genome, gr))[1:4, ]
  m <- m / colSums(m)
  df.bases <- as.data.frame(t(m))
  df.bases$tss_distance <- (-1 * upstream):(downstream - 1)
  df.bases <- melt(df.bases, id.var="tss_distance")
  names(df.bases)[2:3] <- c("base", "frequency")
  df.bases
}


