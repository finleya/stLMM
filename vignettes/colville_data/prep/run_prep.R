args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_file <- sub("^--file=", "", file_arg[[1]])
  prep_dir <- dirname(normalizePath(script_file))
} else {
  prep_dir <- getwd()
}

workflow_dir <- dirname(prep_dir)
setwd(workflow_dir)

source(file.path("prep", "FIA_public_plot_data_prep.R"))
source(file.path("prep", "nf_prep.R"))
