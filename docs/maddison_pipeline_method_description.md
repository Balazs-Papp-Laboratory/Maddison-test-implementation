# Maddison Simulation Pipeline Method Description

This document describes the main computational steps of the Maddison simulation pipeline implemented in this repository. It focuses on the role of each source file in the workflow, the input and output formats used by the steps, and the main statistical ideas behind the implementation.

## Overview

The pipeline implements a simulation-based workflow related to the concentrated-changes test described by Maddison (1990). The biological use case is to test whether gain or loss events of one binary character are concentrated in a region of a phylogenetic tree defined by another binary character.

In this implementation, the binary characters are gene or orthogroup presence/absence patterns reconstructed on a fixed phylogenetic tree. For each pair of genes or orthogroups:

- the **foreground gene** is the tested gene whose gain and loss events are evaluated;
- the **background gene** defines the region of the tree in which those events are tested for enrichment.

The pipeline combines R and Java code. R scripts prepare input tables, concatenate and post-process results, and export the final table. Java programs perform the tree-based operations and the simulation-heavy parts of the Maddison test.

The main workflow steps are:

1. prepare a binary gene presence/absence matrix in DoubleArrayTree node order;
2. count gain and loss events on the tree for each gene;
3. prepare simulation tasks grouped by total gain/loss counts;
4. pregenerate random tree labelings for the Maddison simulation;
5. run the simulation-based Maddison test for all foreground/background gene pairs;
6. concatenate partial Java output files;
7. add odds-ratio effect-size measures;
8. add Fisher exact test results as auxiliary comparison statistics;
9. rename columns for downstream biological analysis;
10. export the final result table as TSV or TSV.GZ.

## Scope of this implementation

The original Maddison concentrated-changes test allows the focal region of the phylogenetic tree to be any user-defined subset of branches. This implementation covers a more specific use case: the focal, or black, region is derived from a second binary character, referred to here as the background gene.

In the main biological analysis, an edge is treated as part of the black region when the background gene is present at both endpoints of the edge. This corresponds to the `blackNodeConfiguration = "11"` setting in `DoMaddisonTestForAllPairs.java`.

The code also supports other simple definitions of the black region based on the parent-to-child transition state of the background gene, for example background-gene gain edges or loss edges. However, the current workflow is still narrower than the most general formulation of the Maddison test, because the black region is not supplied as an arbitrary user-defined branch subset.

The software should therefore be understood as an independent implementation of a simulation-based workflow for a specific Maddison-test use case, rather than a complete general-purpose implementation of all possible concentrated-changes test configurations.

## Input data

The pipeline expects two main biological inputs:

1. a phylogenetic tree in Newick format;
2. a gene or orthogroup copy-number table whose rows correspond to genes or orthogroups and whose columns correspond to tree nodes.

In the biological project for which this pipeline was developed, the copy-number table was produced upstream with the **Count** software. Count is a phylogenetic gene-content analysis program for homolog family sizes, phylogenetic profiles, and other numerical copy-number-like characters along a phylogeny. In this document, phrases such as "Count output" or "Count-derived table" refer to that upstream software step, not merely to the general act of counting copies.

The input table may contain copy-number values:

```text
0 = gene absent
1 = one copy present
2, 3, ... = multiple copies present
```

For the Maddison test implemented here, copy numbers are converted to binary presence/absence values:

```text
0 = gene absent
1 = gene present
```

Any positive copy number is treated as presence. The pipeline therefore does not model copy-number expansion or reduction events; it uses only binary gain/loss dynamics.

The tree and the binary matrix must use matching node names. Later Java steps require that the matrix columns are ordered according to the node-index order used by the Java `DoubleArrayTree` representation.

---

## `script-01-prepare_DAT_ordered_binary_gene_presence_matrix.R`

This script prepares the binary gene presence/absence matrix used by later steps of the pipeline.

It reads the extended gene-family count table produced by the preceding Count pipeline. In that table, each row corresponds to a gene or orthogroup, and the tree-node columns contain copy-number values. The first column is expected to be named `name` and contains the gene or orthogroup identifiers. The table may also contain the summary columns `Gains`, `Losses`, `Expansions`, and `Reductions`; these columns are not used by the Maddison test and are removed.

The script converts all positive copy-number values to `1`, leaving zero values as `0`. The resulting matrix is binary:

```text
0 = gene absent
1 = gene present
```

The script also reorders the node columns to match the node-index order of the Java `DoubleArrayTree` representation built from the input Newick tree. This ordering is required by the later Java and R steps. In the pipeline, the `DAT_ordered` suffix in filenames refers to this DoubleArrayTree node order.

---

## `ger.tree.actions.CountActionsOnTree.java`

`ger.tree.actions.CountActionsOnTree` is the second main pipeline step. It reads the Newick tree and the DAT-ordered binary gene presence/absence matrix created by `script-01-prepare_DAT_ordered_binary_gene_presence_matrix.R`.

The program expects four command-line arguments:

```text
[1] path to the Newick tree file
[2] name of the root node of the subtree to analyse; use "root" for the full tree
[3] path to the DAT-ordered binary gene presence/absence matrix
[4] output TSV filename
```

The current pipeline typically analyses the full tree, so the second argument is usually `root`.

For each gene or orthogroup, the program counts the number of tree edges belonging to each of the four possible parent-to-child binary transition categories:

```text
0 -> 0  absent
0 -> 1  gain
1 -> 0  loss
1 -> 1  present
```

The output is a tab-separated table with one row per gene or orthogroup and the columns:

```text
gene1, absent, gain, loss, present
```

The program assumes that the input matrix contains columns for all tree nodes in the analysed tree or subtree, and that the matrix column names match the tree node names.

---

## `script-03-create_simulation_task_table.R` and `PrepareMaddisonSimulation.java`

The Maddison simulation test requires many random binary tree labelings. A random tree labeling is a simulated gene presence/absence pattern on the fixed phylogenetic tree, assigning `0` or `1` to every tree node.

The simulation does not use arbitrary random labelings. For a given foreground gene, the total numbers of gains and losses observed on the full tree are fixed. The simulation therefore needs random tree labelings with the same total `(gain, loss)` pair as the foreground gene.

This preparation step is split between two programs:

- `script-03-create_simulation_task_table.R`
- `PrepareMaddisonSimulation.java`

The R script determines which `(gain, loss)` pairs will be needed. The Java program generates and stores random tree labelings matching those requirements.

### `script-03-create_simulation_task_table.R`

The input of `script-03-create_simulation_task_table.R` is the table produced by `CountActionsOnTree.java`. For each gene or orthogroup, this table contains the number of gain and loss events on the tree.

The script computes:

```text
gain_plus_loss = gain + loss
```

It then groups genes by the `gain_plus_loss` value. For each group, it records which gain counts occur for that total number of changes. Since the loss count can be reconstructed as:

```text
loss = gain_plus_loss - gain
```

it is sufficient to store the total number of changes and the list of required gain counts.

The output is a task table, typically named `tasks.tsv`, with the columns:

```text
gain_plus_loss
n
gains
```

Column meanings:

- `gain_plus_loss`: total number of gain and loss events in the group;
- `n`: number of genes with this total number of changes; this is mainly informative or diagnostic;
- `gains`: comma-separated list of gain counts required for this `gain_plus_loss` value.

Each row of `tasks.tsv` describes one simulation-generation task for the Java program.

### Why tasks are grouped by `gain + loss`

The expensive part of the simulation preparation is finding random tree labelings that match a given `(gain, loss)` pair.

The implemented random-labeling generator can efficiently fix the total number of changes:

```text
total_changes = gain + loss
```

However, it cannot directly and efficiently force the exact split between gains and losses. The Java program therefore generates many random labelings with a fixed total number of changes and then counts how many gains and losses were actually produced. It keeps only those labelings whose gain/loss split matches one of the required pairs.

Grouping by `gain + loss` improves efficiency. If several required `(gain, loss)` pairs share the same total number of changes, one block of generated random labelings can contribute samples to several target categories. Only labelings whose gain/loss split is not needed are discarded.

### `PrepareMaddisonSimulation.java`

`PrepareMaddisonSimulation.java` reads the task table produced by `script-03-create_simulation_task_table.R`. Each row defines one simulation task with a fixed `gain_plus_loss` value and a list of accepted gain counts.

The program generates random tree labelings for each task. These pregenerated samples are later used by `DoMaddisonTestForAllPairs.java` to compute the simulation-based Maddison p-values.

The program expects ten command-line arguments:

```text
[1] numOfWorkers
    Number of parallel worker threads used to generate simulation samples.

[2] treeNwkFilename
    Path to the input phylogenetic tree in Newick format.

[3] rootNodeName
    Name of the root node of the subtree to analyse. Use "root" for the full tree.

[4] iterationsPerBlock
    Number of random tree labelings generated in one block for a given gain_plus_loss value.

[5] minIterations
    Minimum number of accepted samples to collect for each required gain/loss pair.

[6] maxIterations
    Maximum number of samples to store for each gain/loss pair. Additional hits are not saved.

[7] maxTrial
    Hard upper limit on the number of random labelings attempted for a task.

[8] taskFileName
    Task table produced by script-03-create_simulation_task_table.R.

[9] outDirName
    Output directory for pregenerated simulation samples and the broker file.

[10] logFileName
    Log file used to report task-broker progress during long runs.
```

### Runtime control for sample search

`PrepareMaddisonSimulation.java` uses an accept/reject strategy. It generates random labelings and keeps only those whose gain/loss split matches a required target pair.

For some `(gain, loss)` pairs, the probability of finding a matching random labeling can be very small. In extreme cases, the program may find no suitable samples, or too few samples within a reasonable time. The parameters `iterationsPerBlock`, `minIterations`, `maxIterations`, and `maxTrial` control this process.

The `maxTrial` parameter is a safety limit that ensures that a task finishes after a bounded number of attempts. The trade-off is that rare or difficult gain/loss pairs may end up with too few simulation samples. In that case, the downstream Maddison test may be unreliable or unavailable for the affected foreground genes.

### Generating a random tree labeling

For a fixed tree and a fixed total number of changes, the program randomly selects `total_changes` non-root nodes. Each selected node represents an edge where a state change occurs between the parent and the child.

A selected change position produces one of the two state changes:

```text
0 -> 1  gain
1 -> 0  loss
```

The type of change is not selected directly. It is determined by the state of the parent node, which depends on the previous changes along the path from the root.

The root is not selected as a change position because it has no parent edge. Its state is handled separately: the program considers both possible root states, `0` and `1`, and propagates the binary states through the tree for both cases. This is important because changing the root state swaps the interpretation of downstream changes as gains or losses.

Non-selected child nodes inherit the state of their parent. Selected child nodes switch state relative to their parent. After the full tree labeling has been generated, the program counts the resulting gains and losses. The labeling is saved only if the resulting gain count appears in the `gains` list for the corresponding task row. Otherwise, the labeling is discarded.

### Why direct sampling for a fixed `(gain, loss)` pair is not used

The difficulty is that gain and loss types cannot be chosen independently at the selected nodes. Whether a selected change is a gain or a loss depends on the parent state, and the parent state depends on the previous changes on the tree.

Therefore, it is not straightforward to generate uniformly random tree labelings with a pre-specified exact number of gains and losses. The implemented solution uses the accept/reject strategy described above: it fixes the total number of changes, generates random change positions, propagates states over the tree, counts the resulting gains and losses, and keeps the labeling only if it matches a required target pair.

### Pregenerated samples and restartability

Simulation-sample generation is one of the slowest parts of the pipeline. Therefore, `PrepareMaddisonSimulation.java` pregenerates samples and stores them for later use rather than generating them during every p-value calculation.

The generated samples are written to ZIP files in the output directory. The directory also contains a broker file:

```text
_BROKER_FILE.txt
```

This file records progress for each `gain_plus_loss` task, including trial counts, observed gain-count distributions, and the ZIP files containing accepted samples.

The broker file has two roles:

1. it indexes pregenerated samples for the later Maddison test;
2. it allows long sample-generation runs to be interrupted and resumed.

When the program is restarted, it reconstructs the broker state from the existing broker file and continues from the previously generated samples.

---

## `DoMaddisonTestForAllPairs.java`

`ger.maddison.DoMaddisonTestForAllPairs` performs the main simulation-based Maddison test calculation.

For every `(foreground gene, background gene)` pair, it tests whether gain or loss events of the foreground gene are enriched in the black region defined by the background gene.

The program uses as input:

- the phylogenetic tree in Newick format;
- the DAT-ordered binary gene presence/absence matrix;
- the directory of pregenerated random tree labelings;
- an optional setting controlling the definition of the black region.

The output is a CSV result table. When several worker threads are used, each worker writes a separate partial CSV file. A later R step concatenates these partial files.

### Terminology

The output table uses historical column names:

- `cluster1` = foreground gene, i.e. the tested gene;
- `cluster2` = background gene, i.e. the gene defining the black region.

In the biological analysis, these may correspond to orthogroups rather than individual genes. In the computation, they are treated as binary presence/absence characters.

### Command-line arguments

The program expects six or seven command-line arguments:

```text
[1] numOfWorkers
    Number of parallel worker threads.

[2] treeNwkFilename
    Path to the input phylogenetic tree in Newick format.

[3] rootNodeName
    Name of the root node of the subtree to analyse. Use "root" for the full tree.

[4] genePresenceMatrixFilename
    Path to the DAT-ordered binary gene presence/absence matrix.

[5] dirnameOfPrecalculatedSamples
    Directory containing pregenerated random tree labelings.

[6] outDirName
    Output directory where the result CSV files are written.

[7] blackNodeConfiguration
    Optional parameter controlling which background-gene transitions define the black region.
```

### Input matrix and tree

The program reads the binary gene presence/absence table into an `IntMatrix`. Rows correspond to genes or orthogroups, and columns correspond to tree nodes.

The matrix values are binary:

```text
0 = gene absent
1 = gene present
```

The program expects the matrix columns to be in the same order as the Java `DoubleArrayTree` node order. It checks that tree node names and matrix column names match and that they appear in the same order.

Along each parent-to-child edge, four transitions are possible:

```text
0 -> 0  absent
0 -> 1  gain
1 -> 0  loss
1 -> 1  present
```

These categories are represented internally in the order:

```text
ABSENT, GAIN, LOSS, PRESENT
```

### Black-region definition

The Maddison test asks whether foreground-gene gain/loss events are concentrated in a selected region of the tree. In this implementation, the black region is derived from the background gene.

For each background gene, the program precomputes which tree edges belong to the black region. Events are biologically associated with edges, but the implementation represents an edge by the index of its child node. The black region is therefore stored as a `boolean[]` mask over node indices.

The `blackNodeConfiguration` parameter specifies which background-gene parent-to-child transitions are included in the black region:

```text
00  background gene absent at both endpoints
01  background-gene gain
10  background-gene loss
11  background gene present at both endpoints
```

Categories can be combined, for example:

```text
01+11+10
```

In the main biological use case, the configuration is:

```text
blackNodeConfiguration = "11"
```

This means that an edge is black when the background gene is present at both endpoints.

### Parallelization

The program splits foreground genes across worker threads. Each worker processes its assigned foreground genes against all background genes and writes a separate CSV result file. A later R step concatenates the partial result files.

### Main computation

For each foreground gene, the program:

1. reads the gene's binary presence/absence pattern on the tree;
2. counts its total `absent`, `gain`, `loss`, and `present` transitions;
3. loads pregenerated random labelings with the same total gain and loss counts;
4. iterates over all background genes;
5. retrieves the precomputed black-region mask for the background gene;
6. counts the observed foreground-gene gains and losses in the black region;
7. estimates the null distribution of black-region gains and losses from the simulated labelings;
8. calculates p-values, expected values, and medians;
9. writes one output row for the foreground/background gene pair.

### Simulation null model

For a foreground gene, let:

```text
G = total number of gains
L = total number of losses
```

The program uses pregenerated random tree labelings with the same total `(G, L)` counts. For a given background gene, each simulated labeling is evaluated against the background-gene black region.

For each simulated labeling, the program counts:

```text
x = number of gains in the black region
y = number of losses in the black region
```

This gives a Monte Carlo estimate of the two-dimensional discrete distribution:

```text
P(x, y | G, L, black region)
```

The frequency table is stored in the output column:

```text
table_underscore_separated
```

### One-dimensional distributions and p-values

From the two-dimensional frequency table, the program derives one-dimensional distributions for black-region gains and losses.

It then estimates:

- the expected number of black-region gains;
- the expected number of black-region losses;
- the median number of black-region gains;
- the median number of black-region losses.

The observed foreground-gene gain/loss counts in the black region are compared to these simulated distributions. The program calculates four one-sided p-values:

```text
p_value_less_gains_than_expected
p_value_more_gains_than_expected
p_value_less_loss_than_expected
p_value_more_loss_than_expected
```

For example, `p_value_more_gains_than_expected` measures how unusual it is, under the simulation null model, to observe at least as many foreground-gene gain events in the black region as in the real data.

### Contingency-matrix output columns

The program also computes a 4 x 4 contingency matrix describing the joint transition categories of the foreground and background genes on tree edges.

Foreground-gene categories:

```text
G1_absent
G1_gain
G1_loss
G1_present
```

Background-gene categories:

```text
G2_absent
G2_gain
G2_loss
G2_present
```

The output columns use the prefix `CM_`, for example:

```text
CM_G1_absent_G2_absent
CM_G1_gain_G2_absent
CM_G1_loss_G2_present
CM_G1_present_G2_present
```

These columns are not the main input to the Maddison p-value calculation. They are auxiliary outputs useful for diagnostics, effect-size summaries, and the later Fisher exact test step.

### Important output columns

The output CSV contains one row per `(cluster1, cluster2)` pair. Important columns include:

```text
cluster1
    Foreground gene or orthogroup.

cluster2
    Background gene or orthogroup.

black_area_ratio
    Proportion of the analysed tree region marked as black.

sample_size_of_simulation
    Number of simulation samples loaded for the foreground gene's total (G, L) value.
    If this number is too small, the resulting p-values should be treated as unreliable.

G1_absent, G1_gain, G1_loss, G1_present
    Foreground-gene transition counts on the full tree.

actual_black_absent, actual_black_gain, actual_black_loss, actual_black_present
    Foreground-gene transition counts restricted to the black region.

CM_...
    4 x 4 transition contingency matrix between foreground and background genes.

expected_number_of_black_gains
    Expected number of black-region gains estimated from the simulation distribution.

expected_number_of_black_loss
    Expected number of black-region losses estimated from the simulation distribution.

estimated_median_of_black_gains
    Estimated median number of black-region gains.

estimated_median_of_black_loss
    Estimated median number of black-region losses.

p_value_less_gains_than_expected
p_value_more_gains_than_expected
p_value_less_loss_than_expected
p_value_more_loss_than_expected
    One-sided Maddison simulation p-values.

table_underscore_separated
    String representation of the two-dimensional simulation frequency table.
```

If no pregenerated simulation samples are available for a foreground gene's total `(G, L)` count, `sample_size_of_simulation` is set to `0`, the frequency-table column is empty, and p-values and summary statistics are written as `NaN`.

---

## `script-06-concatenate_maddison_result_parts.R`

This script follows `DoMaddisonTestForAllPairs.java`.

When the Java program runs in parallel, each worker thread writes a separate partial CSV file. This script reads those partial CSV files, concatenates them row-wise, and saves the combined result as a compressed RDS file.

The script expects two or three command-line arguments:

```text
[1] input_dir
    Directory containing the partial CSV files produced by the Java step.

[2] output_rds
    Output RDS filename for the concatenated result table.

[3] drop_table_underscore_separated
    Optional logical flag. If TRUE, the table_underscore_separated column is removed.
    Default: TRUE.
```

The `table_underscore_separated` column contains a string representation of the two-dimensional Maddison simulation frequency table. It can be useful for debugging or for drawing case-study figures, but it can substantially increase file size and memory use in large runs. Therefore, the script removes this column by default.

This script performs no new statistical calculation. It only concatenates partial CSV files and optionally removes a large debug/case-study column.

---

## `script-07-add_maddison_gain_odds_ratios.R`

This script adds two odds-ratio effect-size measures to the Maddison result table. Both are descriptive measures of gain enrichment, but they have different interpretations.

### `odds_ratio_1`

`odds_ratio_1` compares the observed black-region gain odds to the gain odds expected under the Maddison simulation null model:

```text
odds_ratio_1 =
    observed black gain odds /
    expected black gain odds under Maddison simulation
```

It therefore describes gain enrichment in the black region relative to the simulation-based Maddison null model.

### `odds_ratio_2`

`odds_ratio_2` is a direct black-versus-white descriptive comparison:

```text
odds_ratio_2 =
    black-region gain/absent odds /
    non-black-region gain/absent odds
```

It compares the foreground gain-to-absent ratio inside the black region with the same ratio outside the black region.

Both odds ratios can be useful as effect-size summaries, but they do not measure the same quantity:

- `odds_ratio_1`: gain enrichment relative to the Maddison simulation null model;
- `odds_ratio_2`: direct black-versus-white gain/absent comparison.

The script currently computes odds ratios for gain events only, not for loss events.

---

## `script-08-add_fisher_exact_tests.R`

This script augments the Maddison result table with Fisher exact test p-values and a Fisher odds ratio.

It uses four columns from the 4 x 4 transition contingency matrix to build a 2 x 2 table:

```text
CM_G1_gain_G2_present
CM_G1_absent_G2_present
CM_G1_gain_G2_absent
CM_G1_absent_G2_absent
```

The resulting Fisher-test matrix is:

```text
                 G2_present    G2_absent
G1_gain              a             c
G1_absent            b             d
```

Here:

- `G1` is the foreground/tested gene;
- `G2` is the background gene;
- `G1_gain` represents foreground-gene gain events;
- `G1_absent` represents edges where the foreground gene is absent;
- `G2_present` corresponds to the background-gene present region;
- `G2_absent` corresponds to the non-black region.

The Fisher test evaluates whether foreground-gene gain events are enriched in the region where the background gene is present.

### Relationship to the Maddison test

The Fisher test addresses a biological question related to the Maddison test, but it does not account for the phylogenetic tree structure. It treats the entries of the 2 x 2 contingency matrix as if they were independent observations. This is not fully appropriate for phylogenetic data because the possible gain/loss events and their distribution depend on the tree topology and ancestral states.

Therefore, in this pipeline, the Fisher test is used as a complementary diagnostic and comparison statistic. It is not a replacement for the simulation-based Maddison test, which is the phylogenetically more appropriate test in this workflow.

### Fisher output columns

For each unique 2 x 2 contingency matrix, the script runs Fisher tests with:

- `alternative = "two.sided"`
- `alternative = "greater"`
- `alternative = "less"`

It adds the following columns:

```text
Fisher_odds_ratio
Fisher_p_value_both_sided
Fisher_p_value_less_G1_gain_than_expected
Fisher_p_value_more_G1_gain_than_expected
```

The `Fisher_p_value_more_G1_gain_than_expected` column corresponds to the enrichment direction where foreground-gene gains are more frequent in the background-present region.

### Deduplication and parallelization

Fisher tests can be slow when the number of gene pairs is large. The script therefore first deduplicates identical 2 x 2 contingency matrices. Fisher tests are computed only once for each unique matrix, and the results are joined back to the full result table.

The script can run in parallel. The third command-line argument controls the number of worker processes:

```text
[1] input_rds
    Input RDS file.

[2] output_rds
    Output RDS file.

[3] NN
    Number of worker processes.
```

If `NN = 1`, the script runs sequentially. On Unix/Linux systems, `NN > 1` enables fork-based parallel execution. On Windows, fork-based parallelism is not available, so the script falls back to single-threaded execution.

---

## `script-09-rename_result_columns_for_biological_analysis.R`

This script converts the column names of the Maddison result table to the terminology used in the downstream biological analysis.

The input is an RDS table containing the Maddison test results, odds-ratio columns, and Fisher test columns. The script removes some intermediate columns that are not needed later and renames the remaining columns using regular-expression rules.

The script expects three command-line arguments:

```text
[1] input_rds
    Input RDS file.

[2] output_rds
    Output RDS file with renamed columns.

[3] rename_log_csv
    CSV log file documenting column removals and renamings.
```

The script does not perform any new statistical calculation. It only modifies the table structure and column names to make the results easier to interpret in the biological analysis.

The rename log is useful because it records how each original Maddison/Java-style column name was converted to the biological-analysis column name, and which columns were removed.

---

## `script-10-convert_rds_to_tsv_gz.R`

This script is the final export step of the pipeline.

Earlier R steps store result tables in RDS format, which is efficient and convenient within R but less portable for non-R tools. This script converts an RDS table to a tab-separated text file.

The script expects two command-line arguments:

```text
[1] input_rds
    Input RDS file.

[2] output_tsv_or_tsv_gz
    Output TSV filename. If the filename ends with .gz, the output is written as gzip-compressed TSV.
```

The script does not modify the data and performs no statistical calculation. It only converts the file format from R-specific RDS to the more general TSV or TSV.GZ format.
