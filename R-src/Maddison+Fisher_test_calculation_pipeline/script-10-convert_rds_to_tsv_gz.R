#!/usr/bin/env Rscript

# Convert an RDS result table to TSV or TSV.GZ format.
#
# This script is the final export step of the pipeline. Earlier steps store
# intermediate result tables as RDS files, which are convenient and efficient
# within R but less portable for other tools. This script converts an RDS table
# to a tab-separated text file.
#
# If the output filename ends with ".tar.gz", readr::write_tsv() writes a
# gzip-compressed TSV file. This is useful for large result tables.
#
# Command-line arguments:
#
#   [1] input_rds
#       Input RDS file containing an R data frame or tibble.
#
#   [2] output_tsv
#       Output TSV filename. Use a ".tsv.gz" suffix to write a compressed file.

suppressPackageStartupMessages({
  library(readr)
})

print_usage <- function() {
  cat(
    "Usage:\n",
    "  script-10-convert_rds_to_tsv_gz.R <input_rds> <output_tsv_or_tsv_gz>\n\n",
    "Arguments:\n",
    "  input_rds\n",
    "      Input RDS file containing an R data frame or tibble.\n\n",
    "  output_tsv_or_tsv_gz\n",
    "      Output TSV filename. If the filename ends with .gz, the output is\n",
    "      written as gzip-compressed TSV.\n\n",
    "Examples:\n",
    "  script-10-convert_rds_to_tsv_gz.R result.rds result.tsv\n",
    "  script-10-convert_rds_to_tsv_gz.R result.rds result.tsv.gz\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  print_usage()
  stop("Expected exactly 2 command-line arguments.", call. = FALSE)
}

input_rds <- args[1]
output_tsv <- args[2]

if (!file.exists(input_rds)) {
  stop("Input file does not exist: ", input_rds, call. = FALSE)
}

cat("Input RDS file: ", input_rds, "\n", sep = "")
cat("Output TSV file: ", output_tsv, "\n", sep = "")

result_tbl <- readr::read_rds(input_rds)

if (!is.data.frame(result_tbl)) {
  stop("The input RDS object is not a data frame or tibble.", call. = FALSE)
}

cat("Number of rows: ", nrow(result_tbl), "\n", sep = "")
cat("Number of columns: ", ncol(result_tbl), "\n", sep = "")

readr::write_tsv(result_tbl, output_tsv)

cat("Done.\n")
