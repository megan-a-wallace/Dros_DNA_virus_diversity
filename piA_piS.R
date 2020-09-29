#####################
## Total piA & piS ##
#####################

##Written by Megan A. Wallace
##07-09/2020

#script to calculate whole genome and per gene piA and piS from per site nucleotide diversity tables organised into cds regions, and the average no of syn and nsyn sites per gene
#run as Rscript piA_piS.R $virusstem.virus.analyses/$virusstem.genewise.syn-nsyn.sites.tsv $virusstem.virus.analyses/$virusstem.25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv $virusid

library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

##for shell
#reading in the no. of syn and nsyn sites
syn.nsyn.sites=read.table(args[1],header = TRUE, sep = "\t")
#reading in the genewise nucl div table 
nucl.div.table=read.table(args[2],header = TRUE, sep = "\t")
#setting the virusstem from the shell variables
virusid<-args[3]
virusstem<-gsub('[A-Z]{2}[0-9]{6}_','',virusid)
virusstem<-gsub('_virus','',virusstem)

#for testing
#syn.nsyn.sites=read.table(file = "nonKallithea_viruses_diversity/Kallithea.genewise.syn-nsyn.sites.tsv",header = TRUE, sep = "\t")
#nucl.div.table=read.table(file = "nonKallithea_viruses_diversity/Kallithea.25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv",header = TRUE,sep = "\t")
#virusid<-"KX130344_Kallithea_virus"

length(syn.nsyn.sites$cds)->no_cds

#WHOLE GENOME PIA AND PIS
total_nsyn_sites <- sum(syn.nsyn.sites$nsyn_sites)
total_syn_sites <- sum(syn.nsyn.sites$syn_sites)

total_nsyn_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0 & nucl.div.table$type=="non-syn"])
total_syn_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0 & nucl.div.table$type=="syn"])

total_syn_nucl_div <- sum(nucl.div.table$site_nucl_div_maf[nucl.div.table$type=="syn"])
total_nsyn_nucl_div <- sum(nucl.div.table$site_nucl_div_maf[nucl.div.table$type=="non-syn"])

piA<-total_nsyn_nucl_div/total_nsyn_sites
piS<-total_syn_nucl_div/total_syn_sites

piA_over_piS<-piA/piS

#initialising data frame for by gene and total results
piA_piS_results<-data.frame(seqname = as.character(rep.int(virusid,times = no_cds+1)), cds = character(length = no_cds + 1), nsyn_sites = numeric(length = no_cds + 1), syn_sites = numeric(length = no_cds + 1),nsyn_snps = numeric(length = no_cds + 1), syn_snps = numeric(length = no_cds + 1), piA = numeric(length = no_cds + 1), piS = numeric(length = no_cds + 1), piA_over_piS = numeric(length = no_cds + 1), stringsAsFactors = FALSE)

#adding the whole genome values to the bottom of the table
piA_piS_results$cds[no_cds+1]<-"whole_genome"
piA_piS_results$nsyn_sites[no_cds+1]<-total_nsyn_sites
piA_piS_results$syn_sites[no_cds+1]<-total_syn_sites
piA_piS_results$nsyn_snps[no_cds+1]<-total_nsyn_snps
piA_piS_results$syn_snps[no_cds+1]<-total_syn_snps
piA_piS_results$piA[no_cds+1]<-piA
piA_piS_results$piS[no_cds+1]<-piS
piA_piS_results$piA_over_piS[no_cds+1]<-piA_over_piS

for (i in 1:no_cds) {
  
  cds_no_sites<-syn.nsyn.sites[i,]
  cds_name<-as.character(cds_no_sites$cds)
  cds_name->piA_piS_results$cds[i]
  
  nsyn_sites<-cds_no_sites$nsyn_sites
  syn_sites<-cds_no_sites$syn_sites
  piA_piS_results$nsyn_sites[i]<-nsyn_sites
  piA_piS_results$syn_sites[i]<-syn_sites
  
  nsyn_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0 & nucl.div.table$type=="non-syn" & nucl.div.table$cds == cds_name ])
  syn_snps <- length(nucl.div.table$position[nucl.div.table$site_nucl_div_maf>0 & nucl.div.table$type=="syn" & nucl.div.table$cds == cds_name])
  piA_piS_results$nsyn_snps[i]<-nsyn_snps
  piA_piS_results$syn_snps[i]<-syn_snps
  
  syn_nucl_div <- sum(nucl.div.table$site_nucl_div_maf[nucl.div.table$type=="syn" & nucl.div.table$cds == cds_name])
  nsyn_nucl_div <- sum(nucl.div.table$site_nucl_div_maf[nucl.div.table$type=="non-syn" & nucl.div.table$cds == cds_name])
  
  piA<-nsyn_nucl_div/nsyn_sites
  piS<-syn_nucl_div/syn_sites
  
  piA_over_piS<-piA/piS
  
  piA->piA_piS_results$piA[i]
  piS->piA_piS_results$piS[i]
  
  piA_over_piS->piA_piS_results$piA_over_piS[i]
  
}

if (str_detect(virusid,"Vesanto")){
  
  seg<-gsub('^.*_S','S',virusid)
  
  write.table(piA_piS_results,file = paste("Vesanto.virus.remapping.",seg,"/",virusid,".genewise.wholepop.piA_piS.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
} else {
  
  write.table(piA_piS_results,file = paste(virusstem,".virus.analyses/",virusid,".genewise.wholepop.piA_piS.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
}
