#!/usr/bin/env Rscript

# Add Fisher exact test results to Maddison test output.
#
# This script reads an RDS table produced by the previous Maddison pipeline
# steps and adds Fisher exact test statistics based on the 2x2 contingency
# matrix of foreground-gene gain/absent transitions versus background-gene
# present/absent regions.
#
# The Fisher test addresses a question related to the Maddison test: whether
# foreground-gene gain events are enriched when the background gene is present.
# However, unlike the Maddison test, Fisher's exact test does not account for
# phylogenetic dependence or for the constraints imposed by the tree topology.
# It treats the contingency-table counts as if they were independent.
#
# Therefore, the Fisher test is used here as a complementary diagnostic and
# comparison statistic, not as a replacement for the simulation-based Maddison
# test. The Maddison test is the phylogenetically more appropriate test in this
# pipeline, while Fisher's exact test is exact for the simplified contingency
# table.
#
# To improve performance, the script first deduplicates identical 2x2
# contingency matrices. Fisher tests are computed only once for each unique
# matrix, and the results are joined back to the full table.
#
# If the requested number of worker processes is 1, or if the script is running
# on Windows, the Fisher tests are computed sequentially. Parallel execution uses
# fork-based workers and is therefore intended for Unix-like systems such as
# Linux. On Windows, use:
#
#   NN = 1
#
# Command-line arguments:
#
#   [1] input_rds
#       Input RDS file containing Maddison result rows.
#
#   [2] output_rds
#       Output RDS file with added Fisher test columns.
#
#   [3] NN
#       Number of worker processes. Use 1 for sequential execution.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(parallel)
})

print_usage <- function() {
  cat(
    "Usage:\n",
    "  script-08-add_fisher_exact_tests.R <input_rds> <output_rds> <NN>\n\n",
    "Arguments:\n",
    "  input_rds\n",
    "      Input RDS file containing Maddison result rows.\n\n",
    "  output_rds\n",
    "      Output RDS file with added Fisher exact test columns.\n\n",
    "  NN\n",
    "      Number of worker processes. Use 1 for sequential execution.\n",
    "      On Windows, NN is forced to 1 because fork-based parallelism is not available.\n",
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
NN <- suppressWarnings(as.integer(args[3]))

if (is.na(NN) || NN < 1) {
  stop("NN must be a positive integer.", call. = FALSE)
}

if (!file.exists(input_rds)) {
  stop("Input file does not exist: ", input_rds, call. = FALSE)
}

is_windows <- identical(.Platform$OS.type, "windows")

if (is_windows && NN > 1) {
  warning(
    "Fork-based parallel execution is not available on Windows. ",
    "Falling back to sequential execution with NN = 1.",
    call. = FALSE
  )
  NN <- 1
}

cat("Input RDS file: ", input_rds, "\n", sep = "")
cat("Output RDS file: ", output_rds, "\n", sep = "")
cat("Requested worker processes: ", NN, "\n", sep = "")

maddison_tbl <- readr::read_rds(input_rds)

# Remove previously computed Fisher columns, if present.
fisher_column_idx <- grepl("^Fisher_", names(maddison_tbl))

if (any(fisher_column_idx)) {
  cat("Existing Fisher columns found. Removing them before recalculation.\n")
  maddison_tbl <- maddison_tbl[, !fisher_column_idx]
}

# The Fisher test uses the 2x2 contingency matrix:
#
#                  G2_present    G2_absent
#   G1_gain             a             c
#   G1_absent           b             d
#
# where G1 is the foreground/tested gene and G2 is the background gene.
case_tbl <- tidyr::expand_grid(
  G1 = c("gain", "absent"),
  G2 = c("present", "absent")
) %>%
  mutate(name = paste0("CM_G1_", G1, "_G2_", G2))

required_columns <- case_tbl$name
missing_columns <- setdiff(required_columns, names(maddison_tbl))

if (length(missing_columns) > 0) {
  stop(
    "Missing required contingency-matrix column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

# Compute each Fisher test only once for each unique contingency matrix.
fisher_task_tbl <- maddison_tbl %>%
  select(all_of(required_columns)) %>%
  distinct()

cat("Number of result rows: ", nrow(maddison_tbl), "\n", sep = "")
cat("Number of unique Fisher-test matrices: ", nrow(fisher_task_tbl), "\n", sep = "")

fisher_task_tbl <- fisher_task_tbl %>%
  mutate(
    Fisher_odds_ratio = as.numeric(NA),
    Fisher_p_value_both_sided = as.numeric(NA),
    Fisher_p_value_less_G1_gain_than_expected = as.numeric(NA),
    Fisher_p_value_more_G1_gain_than_expected = as.numeric(NA)
  )

run_fisher_tests <- function(task_tbl) {
  if (nrow(task_tbl) == 0) {
    return(task_tbl)
  }
  
  for (i in seq_len(nrow(task_tbl))) {
    if (i == 1 || i %% 500 == 0) {
      cat("[", i, "/", nrow(task_tbl), "]\n", sep = "")
    }
    
    r <- task_tbl[i, ]
    
    mx <- matrix(
      data = c(
        r$CM_G1_gain_G2_present,
        r$CM_G1_absent_G2_present,
        r$CM_G1_gain_G2_absent,
        r$CM_G1_absent_G2_absent
      ),
      nrow = 2,
      ncol = 2,
      dimnames = list(
        c("G1_gain", "G1_absent"),
        c("G2_present", "G2_absent")
      )
    )
    
    fisher_result_two_sided <- fisher.test(x = mx, alternative = "two.sided")
    fisher_result_greater <- fisher.test(x = mx, alternative = "greater")
    fisher_result_less <- fisher.test(x = mx, alternative = "less")
    
    task_tbl$Fisher_p_value_both_sided[i] <- fisher_result_two_sided$p.value
    task_tbl$Fisher_p_value_more_G1_gain_than_expected[i] <- fisher_result_greater$p.value
    task_tbl$Fisher_p_value_less_G1_gain_than_expected[i] <- fisher_result_less$p.value
    task_tbl$Fisher_odds_ratio[i] <- unname(fisher_result_two_sided$estimate)
  }
  
  task_tbl
}

if (NN == 1 || nrow(fisher_task_tbl) <= 1) {
  cat("Running Fisher tests sequentially.\n")
  fisher_task_tbl <- run_fisher_tests(fisher_task_tbl)
} else {
  cat("Running Fisher tests with fork-based parallel workers.\n")
  
  # Split the unique matrices into NN approximately equal subtasks.
  list_of_subtasks <- fisher_task_tbl %>%
    split(sort(rep(seq_len(NN), length.out = nrow(fisher_task_tbl))))
  
  log_file_name <- paste0(
    "tmp-parallel-worker-log-",
    gsub("[ :]", "_", Sys.time()),
    ".txt"
  )
  
  computing_cluster <- parallel::makeForkCluster(
    nnodes = NN,
    outfile = log_file_name
  )
  
  on.exit({
    try(parallel::stopCluster(computing_cluster), silent = TRUE)
  }, add = TRUE)
  
  result_list <- parallel::parLapply(
    cl = computing_cluster,
    X = list_of_subtasks,
    fun = run_fisher_tests
  )
  
  parallel::stopCluster(computing_cluster)
  on.exit(NULL, add = FALSE)
  
  fisher_task_tbl <- dplyr::bind_rows(result_list)
}

# Join Fisher results back to the full Maddison result table.
maddison_tbl <- dplyr::left_join(
  maddison_tbl,
  fisher_task_tbl,
  by = required_columns
)

cat("Writing result to: ", output_rds, "\n", sep = "")
readr::write_rds(maddison_tbl, output_rds, compress = "gz")

cat("Done.\n")
