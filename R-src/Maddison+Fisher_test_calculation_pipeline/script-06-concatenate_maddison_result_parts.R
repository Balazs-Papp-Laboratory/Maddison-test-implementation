#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-

# Concatenate partial Maddison-test result tables produced by the Java step.
#
# The ger.maddison.DoMaddisonTestForAllPairs Java program can run with multiple
# worker threads. Each worker writes its own partial CSV result file into an
# output directory. This script reads those partial CSV files, concatenates them
# into a single table, and writes the combined result as a compressed RDS file.
#
# By default, the script removes the "table_underscore_separated" column before
# concatenation. That column stores the full two-dimensional simulation
# frequency table for each gene pair. It can be useful for debugging or for
# drawing case-study figures, but it contains large strings and can make the
# output file very large. Removing it keeps the downstream result table smaller
# and easier to load into memory.
#
# Command-line arguments:
#
#   [1] input_dir
#       Directory containing the partial CSV files produced by
#       DoMaddisonTestForAllPairs.java.
#
#   [2] output_rds
#       Output filename for the concatenated compressed RDS table.
#
#   [3] drop_table_underscore_separated
#       Optional logical flag. If TRUE, the "table_underscore_separated" column
#       is removed. If FALSE, it is kept. Default: TRUE.
#
# Accepted values for the third argument:
#   TRUE, true, T, t, 1, yes, y, drop
#   FALSE, false, F, f, 0, no, n, keep

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

print_usage <- function() {
  cat(
    "Usage:\n",
    "  script-06-concatenate_maddison_result_parts.R ",
    "<input_dir> <output_rds> [drop_table_underscore_separated]\n\n",
    "Arguments:\n",
    "  input_dir\n",
    "      Directory containing partial CSV files produced by the Java step.\n\n",
    "  output_rds\n",
    "      Output filename for the concatenated compressed RDS table.\n\n",
    "  drop_table_underscore_separated\n",
    "      Optional. Whether to remove the table_underscore_separated column.\n",
    "      Default: TRUE.\n\n",
    "Examples:\n",
    "  script-06-concatenate_maddison_result_parts.R java-results combined.rds\n",
    "  script-06-concatenate_maddison_result_parts.R java-results combined.rds FALSE\n",
    sep = ""
  )
}

parse_logical_flag <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return(TRUE)
  }
  
  value <- tolower(trimws(x))
  
  if (value %in% c("true", "t", "1", "yes", "y", "drop")) {
    return(TRUE)
  }
  
  if (value %in% c("false", "f", "0", "no", "n", "keep")) {
    return(FALSE)
  }
  
  stop(
    "Invalid value for drop_table_underscore_separated: ", x, "\n",
    "Use TRUE/FALSE, yes/no, 1/0, drop/keep.",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) %in% c(2, 3))) {
  print_usage()
  stop("Expected 2 or 3 command-line arguments.", call. = FALSE)
}

input_dir <- args[1]
output_rds <- args[2]
drop_table_underscore_separated <- if (length(args) == 3) {
  parse_logical_flag(args[3])
} else {
  TRUE
}

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir, call. = FALSE)
}

csv_files <- list.files(
  path = input_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(csv_files) == 0) {
  stop("No CSV files found in input directory: ", input_dir, call. = FALSE)
}

cat("Input directory: ", input_dir, "\n", sep = "")
cat("Output RDS file: ", output_rds, "\n", sep = "")
cat("Drop table_underscore_separated: ", drop_table_underscore_separated, "\n", sep = "")
cat("Number of CSV files to read: ", length(csv_files), "\n", sep = "")

tables <- vector("list", length(csv_files))

for (i in seq_along(csv_files)) {
  input_file <- csv_files[[i]]
  
  cat("Reading file [", i, "/", length(csv_files), "]: ", input_file, "\n", sep = "")
  
  tmp <- readr::read_csv(input_file, show_col_types = FALSE)
  
  if (drop_table_underscore_separated) {
    tmp <- tmp %>%
      select(-any_of("table_underscore_separated"))
  }
  
  tables[[i]] <- tmp
  
  cat("  Table size: ")
  print(object.size(tmp), units = "MB")
}

combined_table <- dplyr::bind_rows(tables)
rm(tables)
gc()

cat("Size of concatenated table: ")
print(object.size(combined_table), units = "MB")

readr::write_rds(combined_table, output_rds, compress = "gz")

cat("Done. Written output file: ", output_rds, "\n", sep = "")
