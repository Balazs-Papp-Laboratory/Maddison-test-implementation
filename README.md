
# Maddison concentrated-changes test Nextflow pipeline

This repository contains a Nextflow implementation of a simulation-based version of the Maddison concentrated-changes test for binary characters on a phylogenetic tree.

The method tests whether gains or losses of a binary character are concentrated in a specified region of a phylogenetic tree more often than expected under a random placement model. In this implementation, the tested events are gene gain/loss events, and the tree region of interest is defined by the reconstructed presence or absence of a second, background gene.

The pipeline was designed for batch processing of many gene-pair tests on large phylogenetic trees.

## Background

The implemented method is based on:

Maddison, W. P. 1990.  
*A method for testing the correlated evolution of two binary characters: are gains or losses concentrated on certain branches of a phylogenetic tree?*  
Evolution 44(3):539–557.

The original exact recursive calculation becomes computationally impractical for large trees and many changes. This pipeline therefore uses a simulation-based approximation. To make large batch analyses feasible, random tree labelings can be precomputed and reused for many gene-pair tests.

## Repository structure

```text
.
├── main.nf
├── nextflow.config
├── compile_java.sh
├── README.md
├── java-src/
│   └── ger/
│       ├── maddison/
│       └── tree/
├── java-bin/
│   └── generated after Java compilation
├── R-src/
│   └── Maddison+Fisher_test_calculation_pipeline/
└── test-data/
    └── small/
        ├── for_count.nwk
        └── count_result_FAMILY_extended.tsv
```

Main files:

- `main.nf`: main Nextflow workflow.
- `nextflow.config`: default pipeline configuration.
- `compile_java.sh`: compiles the Java source files into `java-bin/`.
- `java-src/`: Java source code for tree representation and Maddison simulation.
- `R-src/`: R helper scripts used by the workflow.
- `test-data/small/`: small bundled test dataset for a quick test run.

## Requirements

The pipeline requires:

- Nextflow
- Java Development Kit, including `javac`
- Singularity or Apptainer, if using the default container configuration
- Bash

By default, R code is expected to run inside the container configured in `main.nf`.

The current container image is set in `main.nf`:

```nextflow
container = "mesti90/ecoli_hgt_all_in:1.1"
```

If you want to use a different container or run without containers, edit `main.nf` and/or `nextflow.config` accordingly.

## Installation

Clone the repository:

```bash
git clone <repository-url>
cd <repository-directory>
```

Compile the Java sources:

```bash
chmod +x compile_java.sh
./compile_java.sh
```

This creates the `java-bin/` directory containing the compiled Java classes. The Nextflow workflow uses this directory through:

```nextflow
params.java_bin_dir = "$projectDir/java-bin"
```

## Quick test run

A small test dataset is included in:

```text
test-data/small/
```

After compiling the Java code, run:

```bash
nextflow run main.nf
```

By default, `nextflow.config` uses:

```text
test-data/small/for_count.nwk
test-data/small/count_result_FAMILY_extended.tsv
```

The default output directory is configured by `params.outdir` in `nextflow.config`.

## Running on your own data

Use the `--tree` and `--count_family` parameters:

```bash
nextflow run main.nf \
  --tree /path/to/tree.nwk \
  --count_family /path/to/count_result_FAMILY_extended.tsv \
  --outdir /path/to/results
```

The Java classes must be compiled before running the workflow:

```bash
./compile_java.sh
```

## Input files

### 1. Phylogenetic tree

The tree must be provided in Newick format:

```text
/path/to/tree.nwk
```

The terminal labels in the tree must match the sample columns in the gene presence/absence table.

### 2. Gene presence/absence count table

The second input is a tab-separated table such as:

```text
count_result_FAMILY_extended.tsv
```

The table is expected to contain:

- a `name` column identifying the gene or gene family,
- one column per tree tip / genome / terminal taxon,
- optional summary columns such as:
  - `Gains`
  - `Losses`
  - `Expansions`
  - `Reductions`

The R preprocessing step removes the summary columns and converts all positive values in the sample columns to `1`, so the working matrix is binary:

```text
0 = absent
1 = present
```

The tree tip names and the matrix column names must match exactly.

## Output

The workflow writes intermediate and final results to the output directory specified by:

```nextflow
params.outdir
```

or by the command-line option:

```bash
--outdir /path/to/results
```

Important output files include:

- binary gene presence/absence matrices ordered according to the tree,
- counts of gain/loss events on the tree,
- Maddison simulation task lists,
- Maddison simulation results,
- concatenated result tables,
- result tables extended with odds ratios,
- result tables extended with Fisher exact tests,
- final renamed and compressed TSV output.

The exact output filenames are defined in `main.nf`.

## Main workflow steps

The pipeline performs the following major steps:

1. Read the Newick tree and gene presence/absence table.
2. Reorder the gene matrix according to the Double-Array-Tree representation.
3. Count gain and loss events for each gene on the tree.
4. Create a task list for gene-pair Maddison tests.
5. Prepare simulation input for each relevant gain/loss configuration.
6. Run Maddison simulation tests for all selected gene pairs.
7. Concatenate Java result tables.
8. Add odds ratios.
9. Add Fisher exact test results.
10. Rename result columns.
11. Convert the final result table to compressed TSV format.

## Configuration

The most important parameters are defined in `nextflow.config`:

```nextflow
params {
    r_script_dir = "$projectDir/R-src"
    java_src_dir = "$projectDir/java-src"
    java_bin_dir = "$projectDir/java-bin"

    tree = "$projectDir/test-data/small/for_count.nwk"
    count_family = "$projectDir/test-data/small/count_result_FAMILY_extended.tsv"

    outdir = "$launchDir/test-result"
}
```

For files that are part of the repository, such as scripts and test data, `$projectDir` is used.

For user-selected output locations, `$launchDir` or an explicit `--outdir` path can be used.

## Resource settings

The default resource settings are defined in `nextflow.config`, for example:

```nextflow
process.cpus = 15
process.memory = "300 GB"
```

These values may be too high for a local laptop run. For small tests or local debugging, reduce them, for example:

```nextflow
process.cpus = 2
process.memory = "8 GB"
```

For HPC execution, configure the appropriate Nextflow executor and resource settings in `nextflow.config` or in a separate profile.

## Containers

The workflow is configured to use Singularity/Apptainer by default:

```nextflow
singularity.enabled = true
singularity.autoMounts = true
singularity.cacheDir = "$launchDir/containers"
singularity.runOptions = "--bind /"
```

If you do not want to use Singularity/Apptainer, disable it in `nextflow.config` and ensure that all required R packages and Java dependencies are available in your local environment.

## Java compilation

The Java source code is not compiled automatically by the workflow. Compile it manually before running Nextflow:

```bash
./compile_java.sh
```

The script runs:

```bash
mkdir -p java-bin
javac -d java-bin $(find java-src -name "*.java")
```

If you modify any Java source file, rerun:

```bash
./compile_java.sh
```

before restarting the workflow.

## Re-running the workflow

Nextflow can reuse completed tasks with:

```bash
nextflow run main.nf -resume
```

Use this after interrupted runs or after changing only downstream parts of the workflow.

## Notes on the current implementation

This pipeline currently focuses on the case where the black region of the tree is defined by presence of the background gene.

In Maddison terminology, the pipeline tests whether gain/loss events of the tested gene are enriched in the tree region defined by the state of another binary character.

Some legacy scripts may still be present in the repository but are not called by the active workflow in `main.nf`.

## Scope of this implementation

This repository implements a simulation-based version of the Maddison concentrated-changes test for a specific use case.

In the original Maddison concentrated-changes test, the focal region of the phylogenetic tree can be any user-defined subset of branches. In the current implementation, this region is derived from a second binary character, referred to here as the background gene. By default, branches are treated as part of the focal region when the background gene is present at both endpoints of the branch.

Thus, the current software implements a restricted form of the more general Maddison test. A future version is planned to support arbitrary user-defined branch subsets.


## Troubleshooting

### Java classes not found

If you see an error similar to:

```text
Error: Could not find or load main class ger...
```

compile the Java sources:

```bash
./compile_java.sh
```

and check that `java-bin/` exists.

### Tree labels do not match table columns

The preprocessing step requires that the tree tip labels and the sample columns in the count table match. If they do not, the workflow will stop with an error.

Check:

- spelling,
- underscores versus spaces,
- missing samples,
- extra samples in the count table,
- extra tips in the tree.

### Not enough memory

Reduce the test size or increase the memory setting in `nextflow.config`:

```nextflow
process.memory = "300 GB"
```

For local testing, use a much smaller value only if the test dataset is small.

### Container problems

If Singularity/Apptainer cannot pull or run the configured container, either fix the container configuration or disable Singularity and install the required software manually.

## Citation

If you use this pipeline, cite the original Maddison concentrated-changes test:

Maddison, W. P. 1990.  
*A method for testing the correlated evolution of two binary characters: are gains or losses concentrated on certain branches of a phylogenetic tree?*  
Evolution 44(3):539–557.


## Relationship to Maddison (1990) and MacClade

This software is an independent implementation of a simulation-based workflow inspired by the concentrated-changes test described by Maddison (1990). It is designed for the specific use case analysed in the associated biological manuscript, where the focal tree region is derived from the reconstructed presence of a background gene.

The implementation does not include or derive from MacClade source code. MacClade and the original concentrated-changes test should be cited separately where appropriate.

## License

his software is released under the MIT License. See the [LICENSE](LICENSE) file for details.

You are free to use, modify, and redistribute this software under the terms of the MIT License. If you use the software, or results generated by it, in scientific work, please cite the repository and the associated manuscript.


## Contact

Maintainer: Gergely Fekete : 
BRC /
Institute of Biochemistry /
Synthetic and Systems Biology Unit /
Balázs Papp Laboratory

```text
https://sysbiol.brc.hu//papp-balazs-lab-index.html
```
