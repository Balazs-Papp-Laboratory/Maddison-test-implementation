#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline version
pipeline_version = "v0.1"

// Container versions
container = "mesti90/ecoli_hgt_all_in:1.1"


workflow.onComplete {
    def mystr = "===========\n\nMaddison concentrated-changes pipeline started at: $workflow.start.\n\nInput tree: $params.tree\nInput count table: $params.count_family\n\nCommand invoked:\n$workflow.commandLine.\n\nCompleted at: $workflow.complete.\nTotal duration: $workflow.duration.\n========\n\n"
    def file = new File('finished_workflow.info')
    file << mystr
    println mystr
}

workflow.onError {
    def mystr = "===========\n\nMaddison concentrated-changes pipeline started at: $workflow.start.\n\nInput tree: $params.tree\nInput count table: $params.count_family\n\nCommand invoked:\n$workflow.commandLine.\n\nError message:\n$workflow.errorMessage.\n=========\n\n"
    def file = new File('workflow.error')
    file << mystr
    println mystr
}


// testing input parameters 

if (!params.tree) {
    error "Missing required parameter: --tree or params.tree in nextflow.config"
}

if (!params.count_family) {
    error "Missing required parameter: --count_family or params.count_family in nextflow.config"
}


if (!file(params.tree).exists()) {
    error "Tree file does not exist: ${params.tree}"
}

if (!file(params.count_family).exists()) {
    error "Count-family file does not exist: ${params.count_family}"
}


// Maddison processes

process prepare_DAT_ordered_binary_gene_presence_matrix {
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path tree
	path count_family

	output:
	path "gene_presence_mx-DAT_ordered.tsv"
	path "gene_presence_mx-DAT_ordered.rds"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-01-prepare_DAT_ordered_binary_gene_presence_matrix.R \
 	$tree\
 	$count_family\
 	gene_presence_mx-DAT_ordered.tsv\
 	gene_presence_mx-DAT_ordered.rds\
 	${params.java_bin_dir}\
 	-Xmx2g
	"""
}

process maddison_CountActionsOnTree {
    container "$container"
	storeDir "${params.outdir}"

	input:
	path tree
	path gene_presence_mx_DAT_ordered_tsv

	output:
	path "event_cnt.tsv"

	script:
	"""
	java -cp ${params.java_bin_dir}  -Xmx3g \
	ger.tree.actions.CountActionsOnTree \
   	$tree \
	root \
	$gene_presence_mx_DAT_ordered_tsv \
	"event_cnt.tsv"
	"""
}

process maddison_create_tasklist {
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path event_cnt_tsv
	
	output:
	path "tasks.tsv"
	
	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-03-create_simulation_task_table.R \
 	$event_cnt_tsv\
 	"tasks.tsv"
	"""
}

process maddison_PrepareMaddisonSimulation{
    container "$container"
	storeDir "${params.outdir}"

	input:
	path tree
	path maddison_tasklist
	path gene_presence_mx_DAT_ordered_tsv

	output:
	path "dir_of_pregenerated_tree_samples"
	path "sample_generator_log.txt"
	
	script:
	"""
	java -cp ${params.java_bin_dir} -Xmx${task.cpus}g \
	ger.maddison.PrepareMaddisonSimulation \
	${task.cpus} \
    $tree\
    root \
    500000 2000 250000 50000000 \
    $maddison_tasklist \
    dir_of_pregenerated_tree_samples \
    sample_generator_log.txt
    """
}

// The workflow branches here into alternative definitions of the black region.
// The following processes are duplicated with different parameters.

// BLACK 11
process maddison_DoMaddisonTestForAllPairs_BGpresent_black {
    container "$container"
	storeDir "${params.outdir}"

	input:
	path tree
	path gene_presence_mx_DAT_ordered_tsv
	path dir_of_pregenerated_tree_samples

	output:
	path "Maddison_result_BGpresent_black/"

	script:
	"""
	java -cp ${params.java_bin_dir} -Xmx${task.cpus}g \
	ger.maddison.DoMaddisonTestForAllPairs \
	${task.cpus} \
    $tree\
    root \
    $gene_presence_mx_DAT_ordered_tsv \
	$dir_of_pregenerated_tree_samples \
	"Maddison_result_BGpresent_black/" \
    "11"
    """
}

process maddison_concat_result_tables_BGpresent_black{
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path Maddison_result_BGpresent_black

	output:
	path "concated_reduced_Maddison_result_BGpresent_black.rds"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-16-concat_Maddison_result_tables_of_java.R \
 	$Maddison_result_BGpresent_black \
 	concated_reduced_Maddison_result_BGpresent_black.rds
	"""
}

process maddison_add_odds_ratios_BGpresent_black{
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path maddison_result_file_BGpresent_black
	
	output:
	path "maddison_OR_BGpresent_black.rds"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-19A-extend_java_result_add_odds_ratios.R \
 	$maddison_result_file_BGpresent_black \
 	"maddison_OR_BGpresent_black.rds"
	"""
}

process maddison_add_fisher_BGpresent_black {
    container "$container"
    containerOptions "--no-home"
    storeDir "${params.outdir}"

    input:
    path maddison_or_file_BGpresent_black

    output:
    path "maddison_fisher_BGpresent_black.rds"

    script:
    """
    Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-19B-extend_java_result_add_Fisher_tests.R \
    $maddison_or_file_BGpresent_black \
    "maddison_fisher_BGpresent_black.rds" \
    ${task.cpus}
    """
}

process maddison_rename_columns_BGpresent_black{
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path maddison_result_file_BGpresent_black
	
	
	output:
	path "Maddison+Fisher-result_BGpresent_black.rds"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-19C1-rename_columns_only.R\
 	$maddison_result_file_BGpresent_black \
 	"Maddison+Fisher-result_BGpresent_black.rds"\
 	"column_rename_log_BGpresent_black.csv"
	"""
}


process maddison_convert_result_file_to_tsv_gz_BGpresent_black{
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path maddison_result_file_rds_BGpresent_black
	
	
	output:
	path "Maddison+Fisher-result_BGpresent_black.tsv.gz"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-90-convert_rds_to_tsv_gz.R\
 	$maddison_result_file_rds_BGpresent_black \
 	"Maddison+Fisher-result_BGpresent_black.tsv.gz"
	"""
}







// BLACK 00
process maddison_DoMaddisonTestForAllPairs_BGabsent_black {
    container "$container"
	storeDir "${params.outdir}"

	input:
	path tree
	path gene_presence_mx_DAT_ordered_tsv
	path dir_of_pregenerated_tree_samples

	output:
	path "Maddison_result_BGabsent_black/"

	script:
	"""
	java -cp ${params.java_bin_dir} -Xmx${task.cpus}g \
	ger.maddison.DoMaddisonTestForAllPairs \
	${task.cpus} \
    $tree\
    root \
    $gene_presence_mx_DAT_ordered_tsv \
	$dir_of_pregenerated_tree_samples \
	"Maddison_result_BGabsent_black/" \
    "00"
    """
}

process maddison_concat_result_tables_BGabsent_black{
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path Maddison_result_BGabsent_black

	output:
	path "concated_reduced_Maddison_result_BGabsent_black.rds"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-16-concat_Maddison_result_tables_of_java.R \
 	$Maddison_result_BGabsent_black \
 	concated_reduced_Maddison_result_BGabsent_black.rds
	"""
}

process maddison_add_odds_ratios_BGabsent_black{
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path maddison_result_file_BGabsent_black
	
	output:
	path "maddison_OR_BGabsent_black.rds"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-19A-extend_java_result_add_odds_ratios.R \
 	$maddison_result_file_BGabsent_black \
 	"maddison_OR_BGabsent_black.rds"
	"""
}

process maddison_add_fisher_BGabsent_black {
    container "$container"
    containerOptions "--no-home"
    storeDir "${params.outdir}"

    input:
    path maddison_or_file_BGabsent_black

    output:
    path "maddison_fisher_BGabsent_black.rds"

    script:
    """
    Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-19B-extend_java_result_add_Fisher_tests.R \
    $maddison_or_file_BGabsent_black \
    "maddison_fisher_BGabsent_black.rds" \
    ${task.cpus}
    """
}

process maddison_rename_columns_BGabsent_black{
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path maddison_result_file_BGabsent_black
	
	
	output:
	path "Maddison+Fisher-result_BGabsent_black.rds"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-19C1-rename_columns_only.R\
 	$maddison_result_file_BGabsent_black \
 	"Maddison+Fisher-result_BGabsent_black.rds"\
 	"column_rename_log_BGabsent_black.csv"
	"""
}


process maddison_convert_result_file_to_tsv_gz_BGabsent_black{
	container "$container"
    containerOptions "--no-home"
	storeDir "${params.outdir}"

	input:
	path maddison_result_file_rds_BGabsent_black
	
	
	output:
	path "Maddison+Fisher-result_BGabsent_black.tsv.gz"

	script:
	"""
	Rscript --vanilla ${params.r_script_dir}/Maddison+Fisher_test_calculation_pipeline/script-90-convert_rds_to_tsv_gz.R\
 	$maddison_result_file_rds_BGabsent_black \
 	"Maddison+Fisher-result_BGabsent_black.tsv.gz"
	"""
}





//------------------------------------------------------------------


workflow {
    //Maddison
	prepare_DAT_ordered_binary_gene_presence_matrix(params.tree, params.count_family)
	maddison_CountActionsOnTree(params.tree, prepare_DAT_ordered_binary_gene_presence_matrix.out[0]) | maddison_create_tasklist
	maddison_PrepareMaddisonSimulation(params.tree, maddison_create_tasklist.out, prepare_DAT_ordered_binary_gene_presence_matrix.out[0])
	maddison_DoMaddisonTestForAllPairs_BGpresent_black(params.tree,prepare_DAT_ordered_binary_gene_presence_matrix.out[0],maddison_PrepareMaddisonSimulation.out[0]) | maddison_concat_result_tables_BGpresent_black | maddison_add_odds_ratios_BGpresent_black | maddison_add_fisher_BGpresent_black | maddison_rename_columns_BGpresent_black | maddison_convert_result_file_to_tsv_gz_BGpresent_black
}
