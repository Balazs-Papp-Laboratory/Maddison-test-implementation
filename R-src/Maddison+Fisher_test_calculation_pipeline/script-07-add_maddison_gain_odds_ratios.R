#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-


# Add odds-ratio effect-size measures to Maddison test results.
#
# This script reads the concatenated Maddison result table produced by the
# previous pipeline step and adds two odds-ratio columns. These odds ratios were
# used as descriptive effect-size measures in addition to the simulation-based
# Maddison p-values.
#
# The Maddison p-values quantify statistical significance under a simulation
# null model, but their sensitivity depends on the number of available random
# tree labelings and on the total number of gain/loss events. The odds-ratio
# columns added here provide complementary descriptive measures of effect size.
#
# The two odds-ratio columns describe related but different effect-size
# quantities.
#
# odds_ratio_1 compares the observed black-region gain odds to the gain odds
# expected under the Maddison simulation null model:
#
#   observed black gain odds / expected black gain odds under simulation
#
# It therefore measures gain enrichment relative to the simulation-based
# Maddison null distribution.
#
# odds_ratio_2 is a direct descriptive black-versus-white comparison:
#
#   black-region gain/absent odds / non-black-region gain/absent odds
#
# It compares the gain-to-absent ratio inside the black region with the same
# ratio outside the black region.
#
# Both odds ratios can be useful as effect-size measures, but they have
# different interpretations. This script currently computes odds ratios for
# gain events only, not for loss events.
#
# Command-line arguments:
#
#   [1] input_rds
#       Input RDS file containing the concatenated Maddison result table.
#
#   [2] output_rds
#       Output RDS file with the added odds-ratio columns.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

print_usage <- function() {
  cat(
    "Usage:\n",
    "  script-07-add_maddison_gain_odds_ratios.R <input_rds> <output_rds>\n\n",
    "Arguments:\n",
    "  input_rds\n",
    "      Input RDS file containing the concatenated Maddison result table.\n\n",
    "  output_rds\n",
    "      Output RDS file with the added odds-ratio columns.\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  print_usage()
  stop("Expected exactly 2 command-line arguments.", call. = FALSE)
}

input_rds <- args[1]
output_rds <- args[2]

if (!file.exists(input_rds)) {
  stop("Input file does not exist: ", input_rds, call. = FALSE)
}

cat("Input RDS file: ", input_rds, "\n", sep = "")
cat("Output RDS file: ", output_rds, "\n", sep = "")

maddison_tbl <- readr::read_rds(input_rds)

required_columns <- c(
  "actual_black_gain",
  "actual_black_absent",
  "expected_number_of_black_gains",
  "G1_gain",
  "G1_absent"
)

missing_columns <- setdiff(required_columns, names(maddison_tbl))
if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s) in input table: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

maddison_tbl <- maddison_tbl %>%
  mutate(
    # odds_ratio_1:
    # Observed-to-expected odds ratio for foreground-gain events in the black
    # region relative to the Maddison simulation null model.
    #
    # odds_A1 = observed odds of foreground gains falling inside versus outside
    #           the black region
    # odds_B1 = expected odds of foreground gains falling inside versus outside
    #           the black region under the Maddison simulation null model
    odds_A1 = actual_black_gain / (G1_gain - actual_black_gain),
    odds_B1 = expected_number_of_black_gains /
      (G1_gain - expected_number_of_black_gains),
    odds_ratio_1=case_when(
      abs(expected_number_of_black_gains-actual_black_gain)<1e-10 ~ 1.0,
      abs(G1_gain-expected_number_of_black_gains)<1e-10 ~ 0.0,
      TRUE ~ odds_A1/odds_B1),
    
    # odds_ratio_2:
    # Black-versus-white odds ratio comparing foreground gain events to
    # foreground absent edges.
    #
    # odds_A2 = gain/absent odds inside the black region
    # odds_B2 = gain/absent odds outside the black region
    odds_A2 = actual_black_gain / actual_black_absent,
    odds_B2 = (G1_gain - actual_black_gain) /
      (G1_absent - actual_black_absent),
    odds_ratio_2 = odds_A2 / odds_B2
  )

odds_component_columns <- c("odds_A1", "odds_B1", "odds_A2", "odds_B2")

has_negative_odds_component <- maddison_tbl %>%
  select(all_of(odds_component_columns)) %>%
  summarise(across(everything(), ~ any(.x < 0, na.rm = TRUE))) %>%
  unlist(use.names = TRUE)

if (any(has_negative_odds_component)) {
  stop(
    "Negative value detected in odds component column(s): ",
    paste(names(has_negative_odds_component)[has_negative_odds_component], collapse = ", "),
    call. = FALSE
  )
}

readr::write_rds(maddison_tbl, output_rds, compress = "gz")

cat("Done. Written output file: ", output_rds, "\n", sep = "")
