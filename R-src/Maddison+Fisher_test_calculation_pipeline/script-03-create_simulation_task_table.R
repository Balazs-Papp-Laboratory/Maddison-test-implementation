#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-

library(tidyverse)

print_usage <- function() {
  cat(
    "Usage:\n",
    "  script-03-create_simulation_task_table.R <input_event_count_tsv> <output_task_tsv>\n\n",
    "Arguments:\n",
    "  input_event_count_tsv\n",
    "      TSV file produced by CountActionsOnTree. It must contain at least\n",
    "      the columns 'gene1', 'gain', and 'loss'.\n\n",
    "  output_task_tsv\n",
    "      Output TSV file describing the simulation tasks for\n",
    "      PrepareMaddisonSimulation.java. The output columns are:\n",
    "      gain_plus_loss, n, gains.\n\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  print_usage()
  stop("Expected exactly 2 command-line arguments.", call. = FALSE)
}

input_event_count_tsv <- args[1]
output_task_tsv <- args[2]

cat("Input event-count file: ", input_event_count_tsv, "\n", sep = "")
cat("Output task file: ", output_task_tsv, "\n", sep = "")

tbl1 <- readr::read_tsv(input_event_count_tsv, show_col_types = FALSE)

required_columns <- c("gain", "loss")
missing_columns <- setdiff(required_columns, names(tbl1))
if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s) in input table: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

tasks_tbl <-
  tbl1 %>%
  mutate(gain_plus_loss = gain + loss) %>%
  group_by(gain_plus_loss) %>%
  summarise(
    n = n(),
    gains = paste(sort(unique(gain)), collapse = ","),
    .groups = "drop"
  ) %>%
  select(gain_plus_loss, n, gains) %>%
  arrange(gain_plus_loss)

readr::write_tsv(tasks_tbl, output_task_tsv)

cat("Done.\n")
