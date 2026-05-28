#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-

rm(list = ls())
library(tidyverse)


#####################################
## CLI paracssori parameterek atvetele

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  cat("WARNING: Nem kaptam paracssori parametert!\n\n")
  stop()
}else if(length(args)==2){
  file1<-args[1]
  file_out<-args[2]
}else{
  stop("ERROR: 2 parancssori parameter kell! (1) input.rds ,outpot.csv.gz")
}

tbl1<-read_rds(file1)
write_tsv(tbl1, file_out)

cat("Done.")
