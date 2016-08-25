library(magrittr)
library(pander)
library(yaml)

panderOptions("table.split.table", Inf)

options(knitr.project_name = "Repressed Enhancers")

figure_path <- function(filename="") {
  file.path(getOption("knitr.figure_dir"), filename)
}

data_path <- function(append_path=NA) {
  base_path <- yaml.load_file("data/data.yml")$base_path
  if(is.na(append_path))
    return(base_path)
  else
    return(file.path(base_path, append_path))
}

# Format number with commas
pn <- function(i, ...) {
  prettyNum(i, big.mark=",", ...)
}

yesno <- function(booleans) {
  ifelse(booleans == TRUE, "Yes", "No")
}

# Wrap output and code
options(width=80)

# Force knitr to stop evaluation when an error is encountered
knitr::opts_chunk$set(error=FALSE)

# Don't show code blocks by default
knitr::opts_chunk$set(echo=FALSE)

# Don't reformat R code
knitr::opts_chunk$set(tidy=FALSE)

# Set up figure defaults 
knitr::opts_chunk$set(fig.cap="", fig.width=7, fig.height=5, fig.path=figure_path())

# Create output directory if it doesn't exist
if(!file.exists(getOption("knitr.figure_dir"))) dir.create(getOption("knitr.figure_dir"))

source("shared_code/caching.r")
source("shared_code/parallel.r")
