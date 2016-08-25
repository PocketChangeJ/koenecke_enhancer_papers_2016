
checked_mclapply <- function(...) {
  results <- mclapply(...)
  errors <- which(sapply(results, class) == "try-error")
  if(length(errors) > 0) {
    stop("mclapply() returned errors for these elements: ", paste0(names(results), collapse=", "), "\n", paste0(unlist(results[errors]), collapse="\n"))
  }
  results
}
