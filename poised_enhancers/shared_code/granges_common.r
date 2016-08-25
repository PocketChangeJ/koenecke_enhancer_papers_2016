library(GenomicRanges)
library(rtracklayer)

check_coverage_argument <- function(cvg, regions=NULL) {
  if(class(cvg) == "character") {
    stopifnot(file.exists(cvg))
    if(is.null(regions)) {
      cvg <- import(cvg, as="RleList")
    } else {
      cvg <- import(cvg, which=regions, as="RleList")
    }
  } 
  cvg
}

regionApply <- function(regions, cvg, func, ...) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(unlist(
           viewApply(
             Views(cvg, as(regions, "RangesList")),
             function(x) { func(as.numeric(x), ...) },
             simplify=FALSE
           ), use.names=FALSE), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionSums <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewSums(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMeans <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMeans(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewWhichMaxs(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMaxs(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMins <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMins(
             Views(cvg, as(regions, "RangesList"))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionEnrichments <- function(regions, cvg.ip, cvg.bg) {
  ip.ts <- total_signal(cvg.ip)
  bg.ts <- total_signal(cvg.bg)

  cvg.ip <- check_coverage_argument(cvg.ip, regions)
  cvg.bg <- check_coverage_argument(cvg.bg, regions)
  
  return((regionSums(regions, cvg.ip) / ip.ts) / (regionSums(regions, cvg.bg) / bg.ts))
}

total_signal <- function(cov) {
  if(class(cov) == "character") {
    stopifnot(file.exists(cov))
    cache.file <- paste0(cov, ".ts.rds")
    if(file.exists(cache.file)) {
      return(readRDS(cache.file))
    } else {
      cov <- check_coverage_argument(cov)
      ts <- sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x)))))
      saveRDS(ts, file=cache.file)
      return(ts)
    }
  } else {
    return(sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x))))))
  }
}

load_object <- function(filename) {
  if(length(grep("rdata$", tolower(filename))) > 0) {
    o <- get(load(filename))    
  } else {
    o <- faster_readRDS(filename)    
  }  
  updateObject(o)
}

filter_chrs <- function(gr) {
  exclude.chrs <- grep("H|M|U", seqlevels(gr))
  if(length(exclude.chrs) > 0) seqlevels(gr, force=TRUE) <- seqlevels(gr)[-exclude.chrs]
  gr
}
