##################################
## Local piA & piS calculations ##
##################################

#script to calculate whole genome piA and piS for each sample of a specific virus from per site nucleotide diversity tables organised into cds regions, and the average no of syn and nsyn sites per gene
#run as Rscript sample_piA_piS.R $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt $virusstem.virus.analyses/$virusstem.genewise.syn-nsyn.sites.tsv $virusid

library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

##for shell
#reading in the no. of syn and nsyn sites
syn.nsyn.sites=read.table(args[2],header = TRUE, sep = "\t")
total_nsyn_sites <- sum(syn.nsyn.sites$nsyn_sites)
total_syn_sites <- sum(syn.nsyn.sites$syn_sites)

#for testing
#syn.nsyn.sites=read.table(file = "nonKallithea_viruses_diversity/Vesanto_UA_Kan_2016_57_S01.genewise.syn-nsyn.sites.tsv",header = TRUE, sep = "\t")
#nucl.div.table=read.table(file = "nonKallithea_viruses_diversity/Vesanto_UA_Kan_2016_57_S01.25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv",header = TRUE,sep = "\t")
#virusid<-"Vesanto_virus_UA_Kan_2016_57_S01"

#setting the virusstem from the shell variables
virusid<-args[3]
virusstem<-gsub('[A-Z]{2}[0-9]{6}_','',virusid)
virusstem<-gsub('_virus','',virusstem)

#reading in the list of mpileups to retrieve the samples for the genewise nucl div tables
sample_list=read.table(args[1],header = FALSE, sep = "\t")

if (str_detect(virusid,"Vesanto")){
  
  seg<-gsub('^.*_S','S',virusid)
  
  length(sample_list$V1)->no_samples
  
  #initialising data frame for results
  virus_per_sample_piA_piS_results<-data.frame(seqname = as.character(rep.int(virusid,times = no_samples)), sample = character(length = no_samples), nsyn_sites = as.numeric(rep.int(total_nsyn_sites,times = no_samples)), syn_sites = as.numeric(rep.int(total_syn_sites,times = no_samples)),nsyn_snps = numeric(length = no_samples), syn_snps = numeric(length = no_samples), piA = numeric(length = no_samples), piS = numeric(length = no_samples), piA_over_piS = numeric(length = no_samples), stringsAsFactors = FALSE)
  
  for (i in 1:no_samples) {
    
    pop<-gsub(paste(".",virusstem,".bwa.noTIRs_InDel.bam",sep = ""),'',as.character(sample_list$V1[i]))
    pop<-gsub(paste("Vesanto.virus.remapping.",seg,"/",sep = ""),'',pop)
    
    virus_per_sample_piA_piS_results$sample[i]<-pop
    
    if (file.exists(paste("Vesanto.virus.remapping.",seg,"/",pop,".",virusstem,".syn-nsyn.per.site.nucl.div.tsv",sep = ""))) {
      
      nucl.div.table<-read.table(file = paste("Vesanto.virus.remapping.",seg,"/",pop,".",virusstem,".syn-nsyn.per.site.nucl.div.tsv",sep = ""), sep = "\t", header = TRUE)
      
      #WHOLE GENOME PIA AND PIS
      total_nsyn_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0 & nucl.div.table$type=="non-syn"])
      total_syn_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0 & nucl.div.table$type=="syn"])
      
      total_syn_nucl_div <- sum(nucl.div.table$site_nucl_div_maf[nucl.div.table$type=="syn"])
      total_nsyn_nucl_div <- sum(nucl.div.table$site_nucl_div_maf[nucl.div.table$type=="non-syn"])
      
      piA<-total_nsyn_nucl_div/total_nsyn_sites
      piS<-total_syn_nucl_div/total_syn_sites
      
      piA_over_piS<-piA/piS
      
      virus_per_sample_piA_piS_results$nsyn_snps[i]<-total_nsyn_snps
      virus_per_sample_piA_piS_results$syn_snps[i]<-total_syn_snps
      virus_per_sample_piA_piS_results$piA[i]<-piA
      virus_per_sample_piA_piS_results$piS[i]<-piS
      virus_per_sample_piA_piS_results$piA_over_piS[i]<-piA_over_piS
    
      } else {
      
      virus_per_sample_piA_piS_results$nsyn_snps[i]<-0
      virus_per_sample_piA_piS_results$syn_snps[i]<-0
      virus_per_sample_piA_piS_results$piA[i]<-0
      virus_per_sample_piA_piS_results$piS[i]<-0
      virus_per_sample_piA_piS_results$piA_over_piS[i]<-0
      
    }
    
  }
  
  write.table(virus_per_sample_piA_piS_results,file = paste("Vesanto.virus.remapping.",seg,"/",virusid,".per.sample.piA_piS.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
  } else {
  
  length(sample_list$V1)->no_samples
  
  #initialising data frame for results
  virus_per_sample_piA_piS_results<-data.frame(seqname = as.character(rep.int(virusid,times = no_samples)), sample = character(length = no_samples), nsyn_sites = as.numeric(rep.int(total_nsyn_sites,times = no_samples)), syn_sites = as.numeric(rep.int(total_syn_sites,times = no_samples)),nsyn_snps = numeric(length = no_samples), syn_snps = numeric(length = no_samples), piA = numeric(length = no_samples), piS = numeric(length = no_samples), piA_over_piS = numeric(length = no_samples), stringsAsFactors = FALSE)
  
  for (i in 1:no_samples) {
    
    pop<-gsub(paste(".",virusstem,".bwa.noTIRs_InDel.bam",sep = ""),'',as.character(sample_list$V1[i]))
    pop<-gsub(paste(virusstem,".virus.analyses/",sep = ""),'',pop)
    
    virus_per_sample_piA_piS_results$sample[i]<-pop
    
    if (file.exists(paste(virusstem,".virus.analyses/",pop,".",virusstem,".syn-nsyn.per.site.nucl.div.tsv",sep = ""))){
      nucl.div.table<-read.table(file = paste(virusstem,".virus.analyses/",pop,".",virusstem,".syn-nsyn.per.site.nucl.div.tsv",sep = ""), sep = "\t", header = TRUE)
      
      #WHOLE GENOME PIA AND PIS
      total_nsyn_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0 & nucl.div.table$type=="non-syn"])
      total_syn_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0 & nucl.div.table$type=="syn"])
      
      total_syn_nucl_div <- sum(nucl.div.table$site_nucl_div_maf[nucl.div.table$type=="syn"])
      total_nsyn_nucl_div <- sum(nucl.div.table$site_nucl_div_maf[nucl.div.table$type=="non-syn"])
      
      piA<-total_nsyn_nucl_div/total_nsyn_sites
      piS<-total_syn_nucl_div/total_syn_sites
      
      piA_over_piS<-piA/piS
      
      virus_per_sample_piA_piS_results$nsyn_snps[i]<-total_nsyn_snps
      virus_per_sample_piA_piS_results$syn_snps[i]<-total_syn_snps
      virus_per_sample_piA_piS_results$piA[i]<-piA
      virus_per_sample_piA_piS_results$piS[i]<-piS
      virus_per_sample_piA_piS_results$piA_over_piS[i]<-piA_over_piS
      
    } else {
      
      virus_per_sample_piA_piS_results$nsyn_snps[i]<-0
      virus_per_sample_piA_piS_results$syn_snps[i]<-0
      virus_per_sample_piA_piS_results$piA[i]<-0
      virus_per_sample_piA_piS_results$piS[i]<-0
      virus_per_sample_piA_piS_results$piA_over_piS[i]<-0
    }
    
  }
  
  write.table(virus_per_sample_piA_piS_results,file = paste(virusstem,".virus.analyses/",virusid,".per.sample.piA_piS.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
}
