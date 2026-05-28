#!/usr/bin/env Rscript
# -*- coding: UTF-8 -*-


# A ger.tree.actions.CountActionPairsOnTree java program párhuzamosan számítja ki
# a (gain/loss/presence/absence)×(gain/loss/presence/absence) eseménypárokat az osszes gén párra.
# Az eredményét sok kis táblázatbe teszi el egy könytárban.
# Minden Worker thread külön fileba ir, mert így egyszerubb volt megoldani.. 
# Ez az R script ragasztja ossze a darabokat egyetlen nagy tablazatta

library(tidyverse)


args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  cat("WARNING: Nem kaptam paracssori parametert, ezert test adatokkal futok!\n")
  stop()
}else if(length(args)==2){
  dir1<-args[1]
  out_file1<-args[2]
}else{
  stop("ERROR: 2 parancssori parameter kell: input_dir , output_file")
}


args = commandArgs(trailingOnly=TRUE)
dir1<-args[1]
#dir1<-"data-processed/Maddison-resut-java/"
if(substring(dir1, nchar(dir1))!="/"){dir1<-paste0(dir1,"/")}
if(!dir.exists(dir1))
{
  cat("Directory not exist, or is not a directory.: '",dir1,"'\n", sep="")
}

list1<-list();
filelist1<-list.files(dir1)
for(filename in filelist1){
  #DEBUG filename <- filelist1[[1]]
  cat("loading file: ", filename,"\n")
  
  f1<-paste0(dir1,filename)
  
  tmp<-read_csv(f1)
  tmp<- tmp %>%  select(-table_undersore_separated)
  list1[[length(list1)+1]]<-tmp
  
  print(object.size(tmp),units = "MB")
}
tbl1<-do.call(bind_rows,list1)
rm(list1)
cat("Size of big concatenated table:\n")
print(object.size(tbl1),units = "MB")

gc()

# write_csv(tbl1,gsub("/$",".csv.gz",dir1))
# write_csv(tbl1,gsub("/$",".csv",dir1))
# write_rds(tbl1,gsub("/$",".rds",dir1), compress = "gz")
#############

#tbl1 %>% distinct(cluster1) %>%  arrange(cluster1)


tbl2<-tbl1 #%>% select(-table_undersore_separated)
#write_csv(tbl2,gsub("/$","-reduced.csv.gz",dir1))
# write_csv(tbl2,gsub("/$","-reduced.csv",dir1))
write_rds(tbl2,out_file1, compress = "gz")
# 
# paste(names(tbl1), collapse = ", ")
# 
# tbl3<-tbl1 %>% select(cluster1, cluster2, black_area_ratio, sample_size_of_simulation,  G1_gain, G1_loss, actual_black_gain, actual_black_loss, expedted_number_of_black_gains, expedted_number_of_black_loss, estimated_median_of_black_gains, estimated_median_of_black_loss, p_value_less_gains_then_expected, p_value_more_gains_then_expected, p_value_less_loss_then_expected, p_value_more_loss_then_expected)
# write_csv(tbl3,gsub("/$","-minimal.csv.gz",dir1))
# write_csv(tbl3,gsub("/$","-minimal.csv",dir1))
# write_rds(tbl3,gsub("/$","-minimal.rds",dir1), compress = "gz")

