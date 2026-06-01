#!/usr/bin/env bash
set -x

Rscript --vanilla \
  R-src/Maddison+Fisher_test_calculation_pipeline/script-01-prepare_DAT_ordered_binary_gene_presence_matrix.R  \
  test-data/small/tree.nwk \
  test-data/small/count_result_FAMILY_extended.tsv \
  test-data/small/gene_presence_mx-DAT_ordered.tsv \
  test-data/small/gene_presence_mx-DAT_ordered.rds \
  ./java-bin \
  -Xmx2g


java -cp java-bin  -Xmx1g \
	ger.tree.actions.CountActionsOnTree \
   	test-data/small/tree.nwk\
	root \
	test-data/small/gene_presence_mx-DAT_ordered.tsv \
	test-data/small/event_cnt.tsv



Rscript --vanilla R-src/Maddison+Fisher_test_calculation_pipeline/script-03-create_simulation_task_table.R \
 	test-data/small/event_cnt.tsv \
 	test-data/small/tasks.tsv


java -cp java-bin -Xmx1g \
	ger.maddison.PrepareMaddisonSimulation \
	8 \
    test-data/small/tree.nwk\
    root \
    5000 20 2500 50000 \
    test-data/small/tasks.tsv \
    test-data/small/dir_of_pregenerated_tree_samples \
    test-data/small/sample_generator_log.txt


java -cp java-bin -Xmx1g \
	ger.maddison.DoMaddisonTestForAllPairs \
	3 \
    test-data/small/tree.nwk\
    root \
    test-data/small/gene_presence_mx-DAT_ordered.tsv \
	test-data/small/dir_of_pregenerated_tree_samples \
	test-data/small/maddison_result_dir/ \
    "11"


Rscript --vanilla R-src/Maddison+Fisher_test_calculation_pipeline/script-06-concatenate_maddison_result_parts.R \
 	test-data/small/maddison_result_dir/ \
 	test-data/small/concated_madison_result_file.rds \
 	"keep" 


# this step is not needed, only for debug
Rscript --vanilla R-src/Maddison+Fisher_test_calculation_pipeline/script-10-convert_rds_to_tsv_gz.R\
 	test-data/small/concated_madison_result_file.rds \
 	test-data/small/concated_madison_result_file.tsv


Rscript --vanilla R-src/Maddison+Fisher_test_calculation_pipeline/script-07-add_maddison_gain_odds_ratios.R \
 	test-data/small/concated_madison_result_file.rds \
 	test-data/small/maddison_result_ext1.rds


Rscript --vanilla R-src/Maddison+Fisher_test_calculation_pipeline/script-08-add_fisher_exact_tests.R \
    test-data/small/maddison_result_ext1.rds \
    test-data/small/maddison_result_ext2.rds \
    3


Rscript --vanilla R-src/Maddison+Fisher_test_calculation_pipeline/script-09-rename_result_columns_for_biological_analysis.R \
 	test-data/small/maddison_result_ext2.rds \
 	test-data/small/maddison_result_renamed_cols.rds \
 	test-data/small/rename_colnames_log.csv 


Rscript --vanilla R-src/Maddison+Fisher_test_calculation_pipeline/script-10-convert_rds_to_tsv_gz.R\
 	test-data/small/maddison_result_renamed_cols.rds \
 	test-data/small/maddison_result_renamed_cols.tsv

