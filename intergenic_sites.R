###################################################
## Number of intergenic sites in a virus genome? ##
###################################################

##Written by Megan A. Wallace
##07-09/2020

#Script to calculate the number of intergenic sites in a virus genome
#run as Rscript intergenic.sites.R $virusid $refseqlength $virusstem.virus.analyses/$virusid.edited.cds.gtf

library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

##for shell

#setting the virusid and virusstem from the shell variables
virusid<-args[1]
virusstem<-gsub('[A-Z]{2}[0-9]{6}_','',virusid)
virusstem<-gsub('_virus','',virusstem)

#reading in the length of the reference genome 
refseqlength=as.numeric(as.character(args[2]))

#reading in the cds regions from the gtf file, and then summing them and subtracting from the refseqlength to get the number of intergenic sites
read.table(args[3], header = FALSE, col.names = c("seqname","source","feature","start","end","score","strand","frame","id"), sep = "\t")->cds_positions
#making a vector of all the cds positions and then filtering out the overlapping ones 
total_cds_seq<-vector()
for (i in 1:length(cds_positions$start)){
  start<-as.numeric(cds_positions$start[i])
  end<-as.numeric(cds_positions$end[i])
  cds_seq<-seq(start,end,by = 1)
  total_cds_seq<-c(total_cds_seq,cds_seq)
}

unique(total_cds_seq)->nonoverlapping_cds_positions
length(nonoverlapping_cds_positions)->cds_length

refseqlength-cds_length->intergenic_sites

#exporting the result 
if (str_detect(virusstem,"Vesanto")){
  
  seg<-gsub('^.*_S','S',virusstem)
  
  write.table(intergenic_sites,file = paste("Vesanto.virus.remapping.",seg,"/",virusstem,".intergenic.sites.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
} else {
  
  write.table(intergenic_sites,file = paste(virusstem,".virus.analyses/",virusstem,".intergenic.sites.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
}
