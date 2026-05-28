#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-


## A java-ban kiszamolt Maddison teszt eredményeit bővítem 
## ugy hogy kiszamito 2 féle odds ratiot
## Az eredmenytugyanabba a tablazatba es fileba irja vissza,
## csak uj oszlopokat tesz hozza

## UPDATE - 2023-02-22
## A script tovabbra is a teljes fajlt exportalja, de mas neven, mint az input

library(tidyverse)

#####################################
## CLI paracssori parameterek atvetele

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  cat("WARNING: Nem kaptam paracssori parametert, ezert test adatokkal futok!\n")
  #file1<-"data-processed/Maddison-resut-java-reduced-TEST.rds"
  #file1<-"data-processed/subtrees/subtree-root/Maddison-result-java-reduced.rds"
  #file1<-"data-processed/subtrees/subtree-Node9179/Maddison-result-java-reduced.rds"
   file1<-"data-processed2/ST131_fasttree_2022.10.19/ST131_fasttree/subtrees/subtree-root/Maddison-result-java-reduced.rds"
  stop()
}else if(length(args)==2){
  file1<-args[1]
  out_file1<-args[2]
}else{
  stop("ERROR: 2 parancssori parameter kell: input file , output file")
}

###################################
# betoltes

cat("loading file: ",file1,"\n")
tbl_java<-read_rds(file1)


# tbl_java %>%  filter(odds_ratio_1<0)
# tbl_java %>%  filter(odds_ratio_2<0)


######################################
## ODDS RATIO-k kiszamitasa


tbl_java<-tbl_java %>%
 # select(-G2_gain, G2_absent) %>% 
  mutate(
  # G2_gain=CM_G1_absent_G2_gain+CM_G1_gain_G2_gain+CM_G1_loss_G2_gain+CM_G1_present_G2_gain,
  # G2_absent=CM_G1_absent_G2_absent+CM_G1_gain_G2_absent+CM_G1_loss_G2_absent+CM_G1_present_G2_absent,
  
  odds_A1=actual_black_gain/(G1_gain-actual_black_gain), # fekete gain-ek osztva a feher gain-ekkel 
  odds_B1=expedted_number_of_black_gains/(G1_gain-expedted_number_of_black_gains), # ujgyanez várható értékekkel
  odds_ratio_1=case_when(
    abs(expedted_number_of_black_gains-actual_black_gain)<1e-10 ~ 1.0,
    abs(G1_gain-expedted_number_of_black_gains)<1e-10 ~ 0.0,
    TRUE ~ odds_A1/odds_B1),
  # odds_ratio_1 =odds_A1/odds_B1,
  # ha a fekete gain-ek száma nagy, akkor az odds_ratio_1 nagy
  # azaz  ha Gbackground=presence növeli a gain-ek számát azt a  odds_ratio_1>1 jelzi
  
  
  odds_A2=actual_black_gain/actual_black_absent, # fekete gain/absent
  odds_B2=(G1_gain-actual_black_gain)/(G1_absent-actual_black_absent), # feher gain/absent arany
  odds_ratio_2 =odds_A2/odds_B2
  # ha a fekete gain-ek száma nagy, akkor az odds_ratio_2 nagy
  # azaz  ha Gbackground=presence növeli a gain-ek számát az absent-ekhez képest a  odds_ratio_2>1 jelzi
  
    
)
##############################x
# tbl_java %>% filter(odds_ratio_2<0) %>% View()
# tbl_java %>% filter(odds_ratio_1<0)


########################################
## consistency test

tmp<-tbl_java %>%  filter(odds_A1 <0 | odds_B1<0 | odds_A2<0 | odds_B2<0 )
tmp1<-tmp %>%  filter((odds_A1 <0 | odds_B1<0) & !(odds_ratio_1==1 | odds_ratio_1==0))  
tmp2<-tmp %>%  filter(odds_A2 <0 | odds_B2<0 ) 
stopifnot(nrow(tmp1)==0)
stopifnot(nrow(tmp2)==0)

# Ez csak a hibakereseshez szukseges kód
# 
# tmp1  %>%  View()
# tmp2  %>%  View()
# # odds_B2=(G2_gain-actual_black_gain)/(G2_absent-actual_black_absent)
# tmp2 %>% select(odds_B2,G2_gain,actual_black_gain, G2_absent,actual_black_absent)
# 
# # odds_A1=actual_black_gain/(G1_gain-actual_black_gain),
# # odds_B1_=expedted_number_of_black_gains/(G1_gain-expedted_number_of_black_gains), 
# tmp1 %>% select(odds_A1,actual_black_gain,G1_gain,odds_B1_,expedted_number_of_black_gains) 
# 
# 
# x<-tmp1[1,c(-1,-2)]
# t(x)
# #odds_A2=
# x$actual_black_gain/x$actual_black_absent # fekete gain/absent
# 
# CM_names<-grep("CM" ,names(x) , value = TRUE)
# tbl1<-tibble(v=unlist(x[CM_names]), lab1=gsub("^CM_([^_]+)_([^_]+)_([^_]+)_([^_]+)","\\1_\\2",CM_names), lab2=gsub("^CM_([^_]+)_([^_]+)_([^_]+)_([^_]+)","\\3_\\4",CM_names))
# tbl1 %>%  ggplot(mapping = aes(x=lab1, y=lab2, fill=v, label=v))+geom_tile() +geom_text(color="gray")+scale_fill_continuous(trans="log")

#######################################x
### eredmeny mentese

# out_file1<-file1
# out_file2<-gsub("\\.rds",".csv.gz",file1)

# make output filename different from input filename
#out_file1 <- "maddison_OR.rds"

cat("Writeing result to the files :\n  ",out_file1, "\n  ", sep="")

write_rds(tbl_java, out_file1, compress = "gz")
#write_csv(tbl_java,gsub("\\.rds",".csv.gz",file1))

cat("Done.")
