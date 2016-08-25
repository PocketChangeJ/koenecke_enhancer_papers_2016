library(dplyr)
library(readr)
library(GenomicRanges)

import.dhs <- function(filename) {
  dhs.df <- suppressMessages(read_tsv(filename, col_names=FALSE))
  makeGRangesFromDataFrame(dhs.df, seqnames.field="X2", start.field="X3", end.field="X4")
}

load_dhs <- function() {
  files <- list.files(data_path("dhs"), "bdtnpDnaseAcc", full.names=TRUE)
  stopifnot(length(files) > 0)
  names(files) <- gsub("^bdtnpDnaseAcc(.*)\\.txt\\.gz$", "\\1", basename(files))
  lapply(files, import.dhs)
}

dhs_bigwigs <- function() {
  list(S5  = "dhs/bdtnpDnaseS5R9481.bw",
       S9  = "dhs/bdtnpDnaseS9R9127.bw",
       S10 = "dhs/bdtnpDnaseS10R8816.bw",
       S11 = "dhs/bdtnpDnaseS11R9485.bw",
       S14 = "dhs/bdtnpDnaseS14R9477.bw") %>%
    lapply(data_path)
}
