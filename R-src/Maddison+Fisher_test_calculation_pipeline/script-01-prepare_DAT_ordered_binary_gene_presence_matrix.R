#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-


# Prepare a binary gene presence/absence matrix in DoubleArrayTree node order.
#
# This script reads the extended gene-family count table produced by the Count
# pipeline. The input table contains one row per gene or orthogroup and one
# column per node of the phylogenetic tree. The first column, named "name",
# contains the gene or orthogroup identifiers. Four additional summary columns
# ("Gains", "Losses", "Expansions", and "Reductions") may also be present; these
# columns are not used by the Maddison test and are removed.
#
# The node columns may contain copy-number values: 0 means that the gene is
# absent at the corresponding tree node, 1 means that it is present in one copy,
# and values greater than 1 indicate multiple copies. For the Maddison
# concentrated-changes test, copy number is ignored, so all positive values are
# converted to 1. The resulting matrix is therefore binary:
#
#   0 = gene absent
#   1 = gene present
#
# The script also reorders the matrix columns to match the node order used by
# the Java DoubleArrayTree representation of the input Newick tree. This order is
# required by downstream Java and R steps in the analysis pipeline. The output
# filenames use the "DAT_ordered" suffix to indicate this DoubleArrayTree node
# order.





args = commandArgs(trailingOnly=TRUE)
if(length(args)==6)
  {
     tree_filename<-args[[1]]
     input_count_table_filename<-args[[2]]
     output_filename_tsv<-args[[3]]
     output_filename_rds<-args[[4]]
     java_classpath<-args[[5]] # java_classpath="bin"
     java_VM_parameters<-args[[6]] #java_VM_parameters ="-Xmx2g"
     
     
     
     # tree_filename <- "test-data/small/for_count.nwk"
     # input_count_table_filename <- "test-data/small/count_result_FAMILY_extended.tsv"
     # java_classpath <- "java-bin"
     # output_filename_tsv <- "test-data/small/gene_presence_mx-DAT_ordered.tsv"
     # output_filename_rds <- "test-data/small/gene_presence_mx-DAT_ordered.rds"
     # java_VM_parameters ="-Xmx2g"
     

     stopifnot(file.exists(tree_filename))
     stopifnot(file.exists(input_count_table_filename))
     stopifnot(endsWith(output_filename_tsv, "tsv"))
     stopifnot(endsWith(output_filename_rds, "rds"))
  
     dir.create(dirname(output_filename_tsv),recursive = TRUE,showWarnings = FALSE)
     dir.create(dirname(output_filename_rds),recursive = TRUE,showWarnings = FALSE)

  }else{
    cat("Expected 6 command-line arguments:\n",
        "  (1) Newick tree file\n",
        "  (2) input count table produced by the Count pipeline\n",
        "  (3) output TSV filename\n",
        "  (4) output RDS filename\n",
        "  (5) Java classpath, e.g. java-bin\n",
        "  (6) Java VM parameters, e.g. -Xmx2g\n")
  stop()
}

library(tidyverse)
library(tictoc)



########################x
library(rJava)
rJava::.jinit(parameters =java_VM_parameters, classpath =java_classpath, force.init = TRUE) # initialize the java VM

###############################

###############################
tree_J_obj <- .jnew("ger/tree/DoubleArrayTree",.jnew("java/io/File", tree_filename))

names_on_tree <-tree_J_obj$getNames()

cat("Node cnt in the tree: ", length(names_on_tree))
#####



gene_presence_mx<-read_tsv(input_count_table_filename)

gene_presence_mx<-gene_presence_mx %>% select(-any_of(c("Gains","Losses","Expansions","Reductions")))

colnames(gene_presence_mx)<-gsub(" ","_",colnames(gene_presence_mx))

tmp<-setdiff(names_on_tree ,names(gene_presence_mx))
cat("Labels present on tree and missing from FAMILY file:", paste(sprintf("\"%s\"",tmp), collapse = ", "),"\n")
tmp<-setdiff(names(gene_presence_mx),names_on_tree )
cat("Labels present in the FAMILY file and missing from the tree:", paste(sprintf("\"%s\"",tmp), collapse = ", "),"\n")


stopifnot(length(names_on_tree)+1==ncol(gene_presence_mx) )
stopifnot(setdiff(names_on_tree,names(gene_presence_mx)) %>% length()==0)
extra_columns <- setdiff(names(gene_presence_mx), names_on_tree)
stopifnot(identical(extra_columns, "name"))

gene_presence_mx<-gene_presence_mx %>% select(name, all_of(names_on_tree))
gene_names<-gene_presence_mx$name
gene_presence_mx<-as.matrix(gene_presence_mx %>% select(-name))
row.names(gene_presence_mx)<-gene_names


stopifnot( all(colnames(gene_presence_mx)==names_on_tree) )

gene_presence_mx[gene_presence_mx>0]<-1
stopifnot(all(as.vector(gene_presence_mx) %in% c(0, 1)))
#######################################################

write.table(gene_presence_mx,file = output_filename_tsv, quote=FALSE, sep='\t', col.names = NA)
write_rds(gene_presence_mx,output_filename_rds, compress = "gz")
cat("Done.\n")
