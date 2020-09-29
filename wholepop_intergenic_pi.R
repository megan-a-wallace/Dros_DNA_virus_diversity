#########################
## Total intergenic pi ##
#########################

##Written by Megan A. Wallace
##07-09/2020

#Script to calculate per site nucleotide diversity in intergenic regions, from snp files generated from whole population merged mpileup files, and then calculate intergenic pi from this table
#Run as Rscript wholepop_intergenic_pi.R snp.file $virusid $refseqlength $cds.gtf

library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

##for shell

#setting the virusid and virusstem from the shell variables
virusid<-args[2]
virusstem<-gsub('[A-Z]{2}[0-9]{6}_','',virusid)
virusstem<-gsub('_virus','',virusstem)

#setting minimum allele freq
maf<-0.01

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

snp_input_table=read.table(args[1], header = FALSE, col.names = c("seqname","position","ref","cov","A","T","C","G","N"), sep = "\t")

no_snps<-length(snp_input_table$position)

snp_div<-data.frame(seqname = snp_input_table$seqname, position = snp_input_table$position, ref = snp_input_table$ref, cov = snp_input_table$cov, A_count = snp_input_table$A, T_count = snp_input_table$T, C_count = snp_input_table$C, G_count = snp_input_table$G, A_freq = numeric(length = no_snps), A_freq_maf = numeric(length = no_snps), T_freq = numeric(length = no_snps), T_freq_maf = numeric(length = no_snps), C_freq = numeric(length = no_snps), C_freq_maf = numeric(length = no_snps), G_freq = numeric(length = no_snps), G_freq_maf = numeric(length = no_snps), site_nucl_div = numeric(length = no_snps), site_nucl_div_maf = numeric(length = no_snps) )

for (k in 1:no_snps){
  
  #calculating allele frequencies and adding them to the data table of snps
  #where one of the minor allele frequencies is less than 1%, I've changed it to 0 in the freq_maf col, so that you have the option to calculate nucleotide diversity only considering alleles with 1% or higher frequency
  
  snp_div$cov[k]->pos_cov
  A_count<-snp_div$A_count[k]
  A_freq<-A_count/pos_cov
  A_count_maf<-ifelse(A_freq < maf, 0 , A_count )
  
  snp_div$A_freq[k]<-A_freq
  snp_div$A_freq_maf[k]<-ifelse(A_freq < maf, 0 , A_freq )
  
  T_count<-snp_div$T_count[k]
  T_freq<-T_count/pos_cov
  T_count_maf<-ifelse(T_freq < maf, 0, T_count)
  
  snp_div$T_freq[k]<-T_freq
  snp_div$T_freq_maf[k]<-ifelse(T_freq < maf, 0, T_freq )
  
  C_count<-snp_div$C_count[k]
  C_freq<-C_count/pos_cov
  C_count_maf<-ifelse(C_freq < maf, 0, C_count )
  
  snp_div$C_freq[k]<-C_freq
  snp_div$C_freq_maf[k]<-ifelse(C_freq < maf, 0, C_freq )
  
  G_count<-snp_div$G_count[k]
  G_freq<-G_count/pos_cov
  G_count_maf<-ifelse(G_freq < maf, 0, G_count )
  
  snp_div$G_freq[k]<-G_freq
  snp_div$G_freq_maf[k]<-ifelse(G_freq < maf, 0, G_freq )
  
  #now using the counts to calculate the nucleotide diversity for that site
  
  total_count<-pos_cov
  total_count_maf<-sum(A_count_maf,T_count_maf,C_count_maf,G_count_maf)
  
  snp_div$site_nucl_div[k]<-sum( (A_count*(total_count-A_count))/(total_count*(total_count-1)),(T_count*(total_count-T_count))/(total_count*(total_count-1)),(C_count*(total_count-C_count))/(total_count*(total_count-1)),(G_count*(total_count-G_count))/(total_count*(total_count-1)))
  
  snp_div$site_nucl_div_maf[k]<-sum( (A_count_maf*(total_count_maf-A_count_maf))/(total_count_maf*(total_count_maf-1)),(T_count_maf*(total_count_maf-T_count_maf))/(total_count_maf*(total_count_maf-1)),(C_count_maf*(total_count_maf-C_count_maf))/(total_count_maf*(total_count_maf-1)),(G_count_maf*(total_count_maf-G_count_maf))/(total_count_maf*(total_count_maf-1)))
  
}

#exporting the nucl div table
if (str_detect(virusstem,"Vesanto")){
  
  seg<-gsub('^.*_S','S',virusstem)
  
  write.table(snp_div,file = paste("Vesanto.virus.remapping.",seg,"/",virusstem,".wholepop.intergenic.per.site.nucl.div.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
} else {
  
  write.table(snp_div,file = paste(virusstem,".virus.analyses/",virusstem,".wholepop.intergenic.per.site.nucl.div.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
}

#Now using the nucleotide diversity table to calculate intergenic pi
total_snps <- length(snp_div$position[snp_div$site_nucl_div_maf>0])

total_nucl_div <- sum(snp_div$site_nucl_div_maf)

pi<-total_nucl_div/no_intergenic_sites

#creating data frame for results
wholepop_intergenic_pi_results<-data.frame(seqname = as.character(rep.int(virusid,times = 1)), intergenic_sites = as.numeric(rep.int(no_intergenic_sites,times = 1)), no_snps = as.numeric(rep.int(total_snps,times = 1)), pi = as.numeric(rep.int(pi,times = 1)), stringsAsFactors = FALSE)

#exporting the dataframe
if (str_detect(virusstem,"Vesanto")){
  
  seg<-gsub('^.*_S','S',virusstem)
  
  write.table(wholepop_intergenic_pi_results,file = paste("Vesanto.virus.remapping.",seg,"/",virusstem,".wholepop.intergenic.pi.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
} else {
  
  write.table(wholepop_intergenic_pi_results,file = paste(virusstem,".virus.analyses/",virusstem,".wholepop.intergenic.pi.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
}
