
## A java-ban kiszamolt Maddison teszt eredményeit bővítem 
## ugy hogy Fishert szetteketvegzek. Uj oszlopkent beteszi
## a p-értékeket (mindkét oldalit + a kétoldalit)
## és az Fisher_odds_ratio-t
## az eredmenyt ugyanabba a fileba írja vissza
##
## ket CLI parametere: az rds file, és a hasznalni kivant magok szama

library(tidyverse)
#library(parallel)

#####################################
## CLI paracssori parameterek atvetele

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  cat("WARNING: Nem kaptam paracssori parametert, ezert test adatokkal futok!\n")
  file1<-"data-processed2/ST131_fasttree_2022.10.19/ST131_fasttree/subtrees/subtree-root/Maddison-result-java-reduced.rds"
  NN<-3
  stop()
}else if(length(args)==3){
  file1<-args[1]
  out_file1<-args[2]
  NN<-as.integer(args[3])  # number of parralel worker threads
}else{
  stop("ERROR: 3 parancssori parameter kell: input_file.rds, output_file.rds, number_of_threads")
}

###################################
# betoltes

cat("Num of paralell workers:", NN,"\n")
cat("loading file: ",file1,"\n")
tbl_java<-read_rds(file1)



### ###################################################x
### az ODDS RATIO-kat hozzaadtuk
### johetnek a Fisher tesztek


idx1<-grepl("Fisher_",names(tbl_java))

if(any(idx1)){
  cat ("WARNING: The table contains Fisher test results yet.\n I remove the columns and recalculate them.\n")
  tbl_java<-tbl_java[!idx1]
}
# write_rds(tbl_java, file1, compress = "gz")
# write_csv(tbl_java,gsub("\\.rds",".csv.gz",file1))

case_tbl<-expand_grid(G1=c("gain","absent"), G2=c("present","absent")) %>% 
  mutate(name=paste0("CM_G1_",G1,"_G2_",G2))


# A Fisher teszt elvégzésére csak erre a 4 oszlopra van szukseg.
# van egy csomo sor, ami e négy oszlop szempontjabol duplikatum;
# mindegyikbol csak egyet veszek
tbl_Fisher_test<-tbl_java %>% select(all_of(case_tbl$name)) %>% distinct()
stopifnot(ncol(tbl_Fisher_test)==4) # ugye mind a negy oszlop letezik?

# letrehozok 4 ures oszlopot, ahova majd az eredmenyt rakom
tbl_Fisher_test<-tbl_Fisher_test %>%
  add_column(
    Fisher_odds_ratio=as.numeric(NA),
    Fisher_p_value_both_sided=as.numeric(NA) ,
    Fisher_p_value_less_G1_gain_then_expected=as.numeric(NA),
    Fisher_p_value_more_G1_gain_then_expected=as.numeric(NA))

worker_function<-function(tbl_Fisher_test)
{
  for( i in 1:nrow(tbl_Fisher_test)){
    if(i==1 || i%%500==0){cat("[",i,"/",nrow(tbl_Fisher_test),"]\n", sep = "")}
    
    r <- tbl_Fisher_test[i,]
    mx <- matrix(data=c(r$CM_G1_gain_G2_present,r$CM_G1_absent_G2_present,r$CM_G1_gain_G2_absent,r$CM_G1_absent_G2_absent),
      nrow = 2, ncol = 2,dimnames = list(c("G1_gain","G1_absent"),c("G2_present","G2_absent")))
    
    #  Ez olyan test adat, hogy a G2~present esetén nagon megnő a G1 gain-ek gyakorisága
    # ilyenkor a Fisher teszt által visszaadott odds_ratio>1, 
    # és az alternative = "greater" esetén kapunk szignifikáns p-értéket
    # mx=matrix(data=c(100,10,30,30),   nrow = 2, ncol = 2,dimnames = list(c("G1_gain","G1_absent"),c("G2_present","G2_absent")))
    
    
    
    fisher_result_2s<-fisher.test(x = mx,alternative = "two.sided")
    fisher_result_greather<-  fisher.test(x = mx,alternative = "greater")
    fisher_result_less<-fisher.test(x = mx,alternative = "less")
    tbl_Fisher_test$Fisher_p_value_both_sided[i]<-fisher_result_2s$p.value
    tbl_Fisher_test$Fisher_p_value_more_G1_gain_then_expected[i]<-fisher_result_greather$p.value
    tbl_Fisher_test$Fisher_p_value_less_G1_gain_then_expected[i]<-fisher_result_less$p.value
    tbl_Fisher_test$Fisher_odds_ratio[i]<-fisher_result_2s$estimate
  }
  return(tbl_Fisher_test)
}


log_file_name<-paste0("tmp-pararlell-worker-log-",gsub("[ :]","_",Sys.time()),".txt")
computing_cluster<-parallel::makeForkCluster(nnodes = NN, outfile=log_file_name)

# szetszedm a tablazatot annyi reszre, ahány paralell worker thread lesz
list_of_subtasks<-tbl_Fisher_test %>%  split(sort(rep(1:NN, length.out=nrow(tbl_Fisher_test))))
rm(tbl_Fisher_test)
# gc()
res<-parallel::parLapply(computing_cluster,X = list_of_subtasks,fun = worker_function)

parallel::stopCluster(computing_cluster)

tbl_Fisher_test<-do.call(bind_rows,res)

# a Fisher teszt eredmenyeit visszateszem a fo tablazatba
tbl_java<-full_join(tbl_java,tbl_Fisher_test, by=case_tbl$name)

#######################################x
### eredmeny mentese

# out_file1<-"maddison_fisher.rds"
# out_file2<-gsub("\\.rds",".csv.gz",file1)

cat("Writeing result to the files :\n  ",out_file1, "\n", sep="")

write_rds(tbl_java, out_file1, compress = "gz")
#write_csv(tbl_java,gsub("\\.rds",".csv.gz",file1))

cat("Done.")
