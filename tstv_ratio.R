#ratio of transitions to transversions for the virus of interest
#run as Rscript tstv_ratio.R $virusstem.virus.analyses/$virusstem.snp.sums.tsv $virusstem
library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

##
snp.sums=read.table(args[1],header = FALSE,col.names = c("sum","snp_type"))
virusstem<-args[2]
#for testing
#snp.sums=read.table("nonKallithea_viruses_diversity/Linvill_Road.snp.sums.tsv",header = FALSE,col.names = c("sum","snp_type"))

#counting the number of transitions vs transversions 
ts<-snp.sums %>% 
  group_by(snp_type) %>% 
  filter(snp_type == "AG" | snp_type == "GA" | snp_type == "TC" | snp_type == "CT") 
tv<-snp.sums %>% 
  group_by(snp_type) %>% 
  filter(snp_type == "AC" | snp_type == "CA" | snp_type == "TG" | snp_type == "GT" | snp_type == "GC" | snp_type == "CG" | snp_type == "AT" | snp_type == "TA") 

if (as.numeric(sum(tv$sum)) != 0) {
  tstv<-(as.numeric(sum(ts$sum)))/(as.numeric(sum(tv$sum)))
  tstv<-round(tstv)
} else {
  tstv<-(as.numeric(sum(ts$sum)))
}

if (str_detect(virusstem,"Vesanto")){
  
  seg<-gsub('^.*_S','S',virusstem)
  write.table(tstv, file = paste("Vesanto.virus.remapping.",seg,"/",virusstem,".tstv.tsv",sep = ""),append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE)
  
}else{
  
  write.table(tstv, file = paste(virusstem,".virus.analyses/",virusstem,".tstv.tsv",sep = ""),append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE)
  
}

