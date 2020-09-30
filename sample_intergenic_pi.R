######################################
## Local intergenic pi calculations ##
######################################

##Written by Megan A. Wallace
##07-09/2020

#Script to calculate per sample pi for the intergenic regions of a specific virus, and create a table of results for the virus with values for each sample 
#Run as Rscript sample_intergenic_pi.R bam.list.txt $virusid $refseqlength cds.gtf

library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

##for shell
#reading in the length of the reference genome 
refseqlength=as.numeric(as.character(args[3]))

#reading in the cds regions from the gtf file, and then summing them and subtracting from the refseqlength to get the number of intergenic sites
read.table(args[4], header = FALSE, col.names = c("seqname","source","feature","start","end","score","strand","frame","id"), sep = "\t")->cds_positions
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

refseqlength-cds_length->no_intergenic_sites

#setting the virusstem from the shell variables
virusid<-args[2]
virusstem<-gsub('[A-Z]{2}[0-9]{6}_','',virusid)
virusstem<-gsub('_virus','',virusstem)

#reading in the list of mpileups to retrieve the samples for the intergenic nucl div tables
sample_list=read.table(args[1],header = FALSE, sep = "\t")

if (str_detect(virusid,"Vesanto")){
  
  seg<-gsub('^.*_S','S',virusid)
  
  length(sample_list$V1)->no_samples
  
  #initialising data frame for by gene and total results
  virus_per_sample_intergenic_pi_results<-data.frame(seqname = as.character(rep.int(virusid,times = no_samples)), sample = character(length = no_samples), intergenic_sites = as.numeric(rep.int(no_intergenic_sites,times = no_samples)), no_snps = numeric(length = no_samples), pi = numeric(length = no_samples), stringsAsFactors = FALSE)
  
  for (i in 1:no_samples) {
    
    pop<-gsub(paste(".",virusstem,".bwa.noTIRs_InDel.bam",sep = ""),'',as.character(sample_list$V1[i]))
    pop<-gsub(paste("Vesanto.virus.remapping.",seg,"/",sep = ""),'',pop)
    
    virus_per_sample_intergenic_pi_results$sample[i]<-pop
    
    if (file.exists(paste("Vesanto.virus.remapping.",seg,"/",pop,".",virusstem,".intergenic.per.site.nucl.div.tsv",sep = ""))){
      
      nucl.div.table<-read.table(file = paste("Vesanto.virus.remapping.",seg,"/",pop,".",virusstem,".intergenic.per.site.nucl.div.tsv",sep = ""), sep = "\t", header = TRUE)
      
      #PI
      total_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0])
      
      total_nucl_div <- sum(nucl.div.table$site_nucl_div_maf)
      
      pi<-total_nucl_div/no_intergenic_sites
      
      virus_per_sample_intergenic_pi_results$no_snps[i]<-total_snps
      virus_per_sample_intergenic_pi_results$pi[i]<-pi
      
    } else {
      
      virus_per_sample_intergenic_pi_results$no_snps[i]<-0
      virus_per_sample_intergenic_pi_results$pi[i]<-0
      
    }
    
  }
  
  write.table(virus_per_sample_intergenic_pi_results,file = paste("Vesanto.virus.remapping.",seg,"/",virusid,".per.sample.intergenic.pi.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
} else {
  
  length(sample_list$V1)->no_samples
  
  #initialising data frame for by gene and total results
  virus_per_sample_intergenic_pi_results<-data.frame(seqname = as.character(rep.int(virusid,times = no_samples)), sample = character(length = no_samples), intergenic_sites = as.numeric(rep.int(no_intergenic_sites,times = no_samples)), no_snps = numeric(length = no_samples), pi = numeric(length = no_samples), stringsAsFactors = FALSE)
  
  for (i in 1:no_samples) {
    
    pop<-gsub(paste(".",virusstem,".bwa.noTIRs_InDel.bam",sep = ""),'',as.character(sample_list$V1[i]))
    pop<-gsub(paste(virusstem,".virus.analyses/",sep = ""),'',pop)
    
    virus_per_sample_intergenic_pi_results$sample[i]<-pop
    
    if (file.exists(paste(virusstem,".virus.analyses/",pop,".",virusstem,".intergenic.per.site.nucl.div.tsv",sep = ""))) {
      nucl.div.table<-read.table(file = paste(virusstem,".virus.analyses/",pop,".",virusstem,".intergenic.per.site.nucl.div.tsv",sep = ""), sep = "\t",header=TRUE)
      #nucl.div.table<-read.table(file = "nonKallithea_viruses_diversity/2014_07.Kallithea.intergenic.per.site.nucl.div.tsv", sep = "\t", header = TRUE)
      
      #PI
      total_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0])
      
      total_nucl_div <- sum(nucl.div.table$site_nucl_div_maf)
      
      pi<-total_nucl_div/no_intergenic_sites
      
      virus_per_sample_intergenic_pi_results$no_snps[i]<-total_snps
      virus_per_sample_intergenic_pi_results$pi[i]<-pi
      
    } else {
      
      virus_per_sample_intergenic_pi_results$no_snps[i]<-0
      virus_per_sample_intergenic_pi_results$pi[i]<-0
      
    }
    
  }
  
  write.table(virus_per_sample_intergenic_pi_results,file = paste(virusstem,".virus.analyses/",virusid,".per.sample.intergenic.pi.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
}
