#########################################
## Local per-site nucleotide diversity ##
#########################################

##Written by Megan A. Wallace
##07-09/2020

#script to calculate per site nucleotide diversity tables from snp.tsv files of a single sample
#the type of sample can be cds or intergenic 
#run as Rscript sample_per_site_nucl_div.R snp.file type_of_file $virusid $popstem

library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

##for shell
#reading in the type of snps (intergenic or cds)
type<-args[2]

#setting the virusid and virusstem from the shell variables
virusid<-args[3]
virusstem<-gsub('[A-Z]{2}[0-9]{6}_','',virusid)
virusstem<-gsub('_virus','',virusstem)
#setting the sample id from the shell variables
pop<-args[4]

#setting minimum allele freq
maf<-0.01

#reading in the snp table
if (str_detect(type,"cds")){
  snp_input_table=read.table(args[1], header = FALSE, col.names = c("seqname","position","ref","cov","A","T","C","G","type","strand","codon_change","aa_change"), sep = "\t")
  
  no_snps<-length(snp_input_table$position)
  
  snp_div<-data.frame(seqname = snp_input_table$seqname, sample = as.character(rep.int(pop,times = no_snps)), position = snp_input_table$position, ref = snp_input_table$ref, cov = snp_input_table$cov, A_count = snp_input_table$A, T_count = snp_input_table$T, C_count = snp_input_table$C, G_count = snp_input_table$G, type = snp_input_table$type, strand = snp_input_table$strand, codon_change = snp_input_table$codon_change, AA_change = snp_input_table$aa_change, A_freq = numeric(length = no_snps), A_freq_maf = numeric(length = no_snps), T_freq = numeric(length = no_snps), T_freq_maf = numeric(length = no_snps), C_freq = numeric(length = no_snps), C_freq_maf = numeric(length = no_snps), G_freq = numeric(length = no_snps), G_freq_maf = numeric(length = no_snps), site_nucl_div = numeric(length = no_snps), site_nucl_div_maf = numeric(length = no_snps) )

}else{
  
  snp_input_table=read.table(args[1], header = FALSE, col.names = c("seqname","position","ref","cov","A","T","C","G","N"), sep = "\t")
  
  no_snps<-length(snp_input_table$position)
  
  snp_div<-data.frame(seqname = snp_input_table$seqname, sample = as.character(rep.int(pop,times = no_snps)), position = snp_input_table$position, ref = snp_input_table$ref, cov = snp_input_table$cov, A_count = snp_input_table$A, T_count = snp_input_table$T, C_count = snp_input_table$C, G_count = snp_input_table$G, A_freq = numeric(length = no_snps), A_freq_maf = numeric(length = no_snps), T_freq = numeric(length = no_snps), T_freq_maf = numeric(length = no_snps), C_freq = numeric(length = no_snps), C_freq_maf = numeric(length = no_snps), G_freq = numeric(length = no_snps), G_freq_maf = numeric(length = no_snps), site_nucl_div = numeric(length = no_snps), site_nucl_div_maf = numeric(length = no_snps) )
  
}

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

if (str_detect(virusstem,"Vesanto") & str_detect(type,"cds")){
  
  seg<-gsub('^.*_S','S',virusstem)
  
  write.table(snp_div,file = paste("Vesanto.virus.remapping.",seg,"/",pop,".",virusstem,".syn-nsyn.per.site.nucl.div.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
  } else if (str_detect(virusstem,"Vesanto") & str_detect(type,"intergenic")) {
  
    seg<-gsub('^.*_S','S',virusstem)
    
    write.table(snp_div,file = paste("Vesanto.virus.remapping.",seg,"/",pop,".",virusstem,".intergenic.per.site.nucl.div.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
    
  } else if (str_detect(type,"cds")) {
    
    write.table(snp_div,file = paste(virusstem,".virus.analyses/",pop,".",virusstem,".syn-nsyn.per.site.nucl.div.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
    
  } else {
    
    write.table(snp_div,file = paste(virusstem,".virus.analyses/",pop,".",virusstem,".intergenic.per.site.nucl.div.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
    
}
