#!/usr/bin/Rscript
# BioRemmer v1.0 - Report generator
# Usage: Rscript report.R [working_dir]
#   working_dir: root directory of the BioRemmer run (default: current directory)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) >= 1) args[1] else getwd()

# Resolve the script's own directory to find Biorem_report.Rmd
script_dir <- tryCatch({
  argv <- commandArgs(trailingOnly = FALSE)
  match <- grep("--file=", argv, value = TRUE)
  if (length(match) > 0) {
    normalizePath(dirname(sub("--file=", "", match[1])))
  } else {
    normalizePath(dirname(sys.frames()[[1]]$ofile))
  }
}, error = function(e) normalizePath(getwd()))

rmd_file <- file.path(script_dir, "Biorem_report.Rmd")
if (!file.exists(rmd_file)) {
  stop("Cannot find Biorem_report.Rmd in: ", script_dir)
}

library(rmarkdown)

render(
  input       = rmd_file,
  output_file = file.path(work_dir, "Results", "Biorem_report.html"),
  params      = list(work_dir = work_dir)
)
