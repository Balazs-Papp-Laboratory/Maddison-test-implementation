#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-


## A java-ban kiszamolt Maddison teszt eredménytablazatat finomitom
## atnevezem az oszlopokat, es mindket gen alapjan hozzatesze m a metaadatokat

rm(list = ls())
library(tidyverse)


#####################################
## CLI paracssori parameterek atvetele

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  cat("WARNING: Nem kaptam paracssori parametert!\n\n")
  stop()
  
}else if(length(args)==3){
  file1<-args[1]
  file_out<-args[2]
  log_file_out<-args[3]
}else{
  stop("ERROR: 3 parancssori parameter kell! (1) Maddison result tablazat, (2) output file (3) log file neve")
}


cat("Fo tablazat:",file1,"\n")
cat("Output file: ",file_out,"\n")
cat("colum_rename_log_file: ",log_file_out,"\n")

###################################
### betoltes



tbl1<-read_rds(file1)

nrow_tbl1<-nrow(tbl1)

###########x
# tbl1 %>%  filter(odds_ratio_1<0)
# tbl1 %>%  filter(odds_ratio_2<0)
###################################
# torlom az oszlopokat amik nem kellenek

columns_to_remove<-str_split("G1_absent G1_gain G1_loss G1_present actual_black_absent actual_black_gain actual_black_loss actual_black_present G2_gain G2_absent odds_A1 odds_B1 odds_A2 odds_B2",pattern = " ")[[1]]
columns_to_remove<-intersect(colnames(tbl1),columns_to_remove)

tbl1<-tbl1 %>% select(-all_of(columns_to_remove))

###################################
## oszlop átnevezések
#
# Eszter kérése alapjan (2022.05.30 megbeszelees)

rename_tbl<-tibble(orig=names(tbl1), new_name=orig)
rename_tbl<-rename_tbl %>% 
  mutate(new_name=gsub("cluster1", "Orthogroup_G1", new_name)) %>%
  mutate(new_name=gsub("cluster2", "Orthogroup_G2", new_name)) %>% 
  mutate(new_name=gsub("G1", "GTested", new_name)) %>%
  mutate(new_name=gsub("G2", "GBackground", new_name)) %>% 
  #mutate(new_name=gsub("CM_([^_]+)_([^_]+)_([^_]+)_([^_]+)", "CM_\\1.\\2_\\3.\\4", new_name)) %>% 
  mutate(new_name=gsub("CM_([^_]+)_([^_]+)_([^_]+)_([^_]+)", "CM_\\3.\\4_\\1.\\2", new_name)) %>% 
  mutate(new_name=gsub("_of_black_", "_of_GBacground.present_GTested.", new_name)) %>% 
  mutate(new_name=gsub("p_value_([^_]+)_([^_]+)_then_expected", "Maddison_p_value_\\1_GTested.\\2_then_expected", new_name)) %>%
  mutate(new_name=gsub("expedted_number_of_GBacground.([^_]+)_GTested.([^_]+)", "expedted_number_of_GTested.\\2_if_GBackground.\\1", new_name)) %>%
  mutate(new_name=gsub("estimated_median_of_GBacground.([^_]+)_GTested.([^_]+)", "estimated_median_of_GTested.\\2_if_GBackground.\\1", new_name)) %>% 
  mutate(new_name=gsub("sample_size_of_simulation", "sample_size_of_Maddison_simulation", new_name)) %>% 
  mutate(new_name=gsub("black_area_ratio", "GBackground.present_ratio", new_name))
# %>% 
#   mutate(new_name=gsub("gains", "gain", new_name)) 
# %>% 
#   View()

log_tbl<-bind_rows(rename_tbl, tibble(orig=columns_to_remove, new_name=NA))
write_csv(log_tbl, log_file_out)

stopifnot(all(names(tbl1)==rename_tbl$orig))
names(tbl1)<-rename_tbl$new_name


# cat("writeing :", file_out,"\n")
# if( grepl("\\.rds$", file_out) )
# {
#   write_rds(tbl1,file_out, compress = "gz")
# }
# if(grepl("\\.tsv$", file_out) || grepl("\\.tsv\\.gz$", file_out) )
# {
#   write_tsv(tbl1,file_out)
# }
# if(grepl("\\.csv$", file_out) || grepl("\\.csv\\.gz$", file_out) )
# {
#   write_csv(tbl1,file_out)
# }

write_rds(tbl1, file_out, compress = "gz")


cat("Done.")
