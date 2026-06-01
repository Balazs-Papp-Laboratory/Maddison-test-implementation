#!/usr/bin/env Rscript

# Rename Maddison result columns to the terminology used in the biological analysis.
#
# This script reads an RDS table produced by the previous Maddison pipeline
# steps, removes intermediate columns that are no longer needed for downstream
# analysis, and renames the remaining columns.
#
# The input table uses a mixture of names originating from the Java implementation,
# the Maddison-test terminology, and earlier project-specific conventions. This
# script converts those names to the terminology used in the biological analysis.
#
# The renaming is rule-based and uses regular expressions. A log table is written
# to document each column transformation and each removed column.
#
# Command-line arguments:
#
#   [1] input_rds
#       Input RDS file containing the Maddison result table.
#
#   [2] output_rds
#       Output RDS file with renamed columns.
#
#   [3] rename_log_csv
#       CSV file documenting the column removals and renamings.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

print_usage <- function() {
  cat(
    "Usage:\n",
    "  script-09-rename_result_columns_for_biological_analysis.R <input_rds> <output_rds> <rename_log_csv>\n\n",
    "Arguments:\n",
    "  input_rds\n",
    "      Input RDS file containing the Maddison result table.\n\n",
    "  output_rds\n",
    "      Output RDS file with renamed columns.\n\n",
    "  rename_log_csv\n",
    "      CSV file documenting removed and renamed columns.\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  print_usage()
  stop("Expected exactly 3 command-line arguments.", call. = FALSE)
}

input_rds <- args[1]
output_rds <- args[2]
rename_log_csv <- args[3]

if (!file.exists(input_rds)) {
  stop("Input file does not exist: ", input_rds, call. = FALSE)
}

cat("Input RDS file: ", input_rds, "\n", sep = "")
cat("Output RDS file: ", output_rds, "\n", sep = "")
cat("Column rename log file: ", rename_log_csv, "\n", sep = "")

result_tbl <- readr::read_rds(input_rds)

cat("Number of rows: ", nrow(result_tbl), "\n", sep = "")
cat("Number of columns before cleanup: ", ncol(result_tbl), "\n", sep = "")

# Intermediate columns that are not needed for downstream biological analysis.
#
# These include raw transition-count columns and odds-ratio component columns.
# The final odds ratios and Maddison/Fisher p-values are kept.
columns_to_remove <- c(
  "G1_absent",
  "G1_gain",
  "G1_loss",
  "G1_present",
  "actual_black_absent",
  "actual_black_gain",
  "actual_black_loss",
  "actual_black_present",
  "G2_gain",
  "G2_absent",
  "odds_A1",
  "odds_B1",
  "odds_A2",
  "odds_B2"
)

columns_to_remove_present <- intersect(names(result_tbl), columns_to_remove)

if (length(columns_to_remove_present) > 0) {
  result_tbl <- result_tbl %>%
    select(-all_of(columns_to_remove_present))
}

# Build the column-renaming table.
#
# The order of the rules matters. First, the original foreground/background
# gene labels are normalized, then contingency-matrix and Maddison-statistic
# names are converted to the biological-analysis terminology.
rename_tbl <- tibble(
  orig = names(result_tbl),
  new_name = orig
) %>%
  mutate(new_name = gsub("cluster1", "Orthogroup_G1", new_name)) %>%
  mutate(new_name = gsub("cluster2", "Orthogroup_G2", new_name)) %>%
  mutate(new_name = gsub("G1", "GTested", new_name)) %>%
  mutate(new_name = gsub("G2", "GBackground", new_name)) %>%
  mutate(
    new_name = gsub(
      "CM_([^_]+)_([^_]+)_([^_]+)_([^_]+)",
      "CM_\\3.\\4_\\1.\\2",
      new_name
    )
  ) %>%
  mutate(
    new_name = gsub(
      "_of_black_",
      "_of_GBackground.present_GTested.",
      new_name
    )
  ) %>%
  mutate(
    new_name = gsub(
      "p_value_([^_]+)_([^_]+)_than_expected",
      "Maddison_p_value_\\1_GTested.\\2_than_expected",
      new_name
    )
  ) %>%
  mutate(
    new_name = gsub(
      "expected_number_of_GBackground\\.([^_]+)_GTested\\.([^_]+)",
      "expected_number_of_GTested.\\2_if_GBackground.\\1",
      new_name
    )
  ) %>%
  mutate(
    new_name = gsub(
      "estimated_median_of_GBackground\\.([^_]+)_GTested\\.([^_]+)",
      "estimated_median_of_GTested.\\2_if_GBackground.\\1",
      new_name
    )
  ) %>%
  mutate(
    new_name = gsub(
      "sample_size_of_simulation",
      "sample_size_of_Maddison_simulation",
      new_name
    )
  ) %>%
  mutate(
    new_name = gsub(
      "black_area_ratio",
      "GBackground.present_ratio",
      new_name
    )
  )

if (anyDuplicated(rename_tbl$new_name)) {
  duplicated_names <- rename_tbl$new_name[duplicated(rename_tbl$new_name)]
  stop(
    "Column renaming would create duplicate column name(s): ",
    paste(unique(duplicated_names), collapse = ", "),
    call. = FALSE
  )
}

stopifnot(identical(names(result_tbl), rename_tbl$orig))

# Create a log table documenting renamed, unchanged, and removed columns.
rename_log_tbl <- rename_tbl %>%
  mutate(
    action = if_else(orig == new_name, "kept", "renamed")
  )

remove_log_tbl <- tibble(
  orig = columns_to_remove_present,
  new_name = NA_character_,
  action = "removed"
)

log_tbl <- bind_rows(rename_log_tbl, remove_log_tbl)

readr::write_csv(log_tbl, rename_log_csv)

names(result_tbl) <- rename_tbl$new_name

cat("Number of columns after cleanup: ", ncol(result_tbl), "\n", sep = "")
cat("Writing renamed result table: ", output_rds, "\n", sep = "")

readr::write_rds(result_tbl, output_rds, compress = "gz")

cat("Done.\n")
