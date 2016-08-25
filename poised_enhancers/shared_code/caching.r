library(digest)

faster_saveRDS <- function(object, file, cores=4) {
  pigz_out <- pipe(paste0("pigz -p", cores, " > ", file), "wb")
  on.exit(close(pigz_out))
  saveRDS(object, pigz_out)
}

faster_readRDS <- function(file) {
  stopifnot(file.exists(file))
  pigz_in <- pipe(paste0("pigz -dc ", file), "rb")
  on.exit(close(pigz_in))
  readRDS(pigz_in)
}

snap_saveRDS <- function(object, file) {
  snappy_out <- pipe(paste0("snzip > ", file), "wb")
  on.exit(close(snappy_out))
  saveRDS(object, snappy_out)
}

snap_readRDS <- function(file) {
  stopifnot(file.exists(file))
  snappy_in <- pipe(paste0("cat ", file, " | snzip -dc"), "rb")
  on.exit(close(snappy_in))
  readRDS(snappy_in)
}

cache <- function(filename, func, absolute.path=FALSE, cache=TRUE) {
  if(!cache) return(func())
    
  filename <- paste0(filename, ".rds")
  if(!absolute.path) filename <- figure_path(filename)
  hashfile <- paste0(filename, ".hash")
  current_hash <- digest(paste0(as.character(body(func)), collapse="\n"), NULL, ascii=TRUE)
  new_hash_created <- FALSE
  
  if(!file.exists(hashfile)) {
    new_hash_created <- TRUE
    saveRDS(current_hash, file=hashfile)
  }
  
  if(file.exists(filename)) {
    saved_hash <- readRDS(hashfile)
    if(saved_hash != current_hash) {
      warning("For cache file '", filename, "': function has been modified since results were cached.")
      cache_results <- func()
      faster_saveRDS(cache_results, file=filename)
      saveRDS(current_hash, file=hashfile)
      return(cache_results)
    }
    if(new_hash_created) {
      warning("For cache file '", filename, "': function hash was not present.")
    }
    return(faster_readRDS(filename))
  } else {
    cache_results <- func()
    faster_saveRDS(cache_results, file=filename)
    saveRDS(current_hash, file=hashfile)
    return(cache_results)
  }
}

clear_caches <- function() {
  file.remove(list.files(figure_path(), full.names=TRUE))
}
