###Per site nucleotide diversity for total population merged bam mpileup files 

#Aim : to calculate per site nucleotide diversity for each virus in the dataset, which can then be used to calculate piA and piS

###packages
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(devtools)

#list of viruses to calculate piA and piS for 
read.table("gt5samples.virus.list.txt",sep="\t",header=FALSE)->virus_list

virus_list<-sapply(virus_list$V1, as.character)
n_viruses<-as.numeric(length(virus_list))

#LOOP TO SPLIT SNP FILES BY GENE ANNOTATIONS AND CALCULATE PER SITE PI

#setting minimum allele frequency
maf<-0.01

for (i in 1:n_viruses){
  
  #virusid<-virus_list[34]
  virusid<-virus_list[i]
  virusstem<-gsub('[A-Z]{2}[0-9]{6}_','',virusid)
  virusstem<-gsub('_virus','',virusstem)
  
  if (str_detect(virusstem,"Vesanto")){
    
    #for vesanto
    seg<-gsub('^.*_S','S',virusstem)
    
    #import the gtf cds annotation file for the virus/segment variant
    read.table(file = paste(virusstem,".virus.analyses/",virusstem,".popoolation1.gtf", sep = ""), header = FALSE, col.names = c("seqname","source","feature","start","end","score","strand","frame","id"), sep = "\t")->cds_positions
    #read.table(file = paste("nonKallithea_viruses_diversity/",virusstem,".popoolation1.gtf", sep = ""), header = FALSE, col.names = c("seqname","source","feature","start","end","score","strand","frame","id"), sep = "\t")->cds_positions
    
    no_cds<-length(cds_positions$seqname)
    
    #import the snps for the virus/segment variant
    read.table(file = paste("Vesanto.virus.remapping.",seg,"/",virusstem,".25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.tsv", sep = ""), header = FALSE, col.names = c("seqname","position","ref","cov","A","T","C","G","type","strand","codon_change","aa_change"), sep = "\t")->whole_virus_snps
    #read.table(file = paste("nonKallithea_viruses_diversity/",virusstem,".25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.tsv", sep = ""), header = FALSE, col.names = c("seqname","position","ref","cov","A","T","C","G","type","strand","codon_change","aa_change"), sep = "\t")->whole_virus_snps
    
    #initiating the per virus/segment results table for per site nucl div
    virus_nuc_div_table<-data.frame(seqname = character(length = 0), cds = character(length = 0), position = numeric(length = 0), ref = character(length = 0), cov = numeric(length = 0), A_count = numeric(length = 0) , T_count = numeric(length = 0), C_count = numeric(length = 0), G_count = numeric(length = 0), type = character(length = 0), strand = character(length = 0) , codon_change = character(length = 0), AA_change = character(length = 0), A_freq = numeric(length = 0), A_freq_maf = numeric(length = 0), T_freq = numeric(length = 0), T_freq_maf = numeric(length = 0), C_freq = numeric(length = 0), C_freq_maf = numeric(length = 0), G_freq = numeric(length = 0), G_freq_maf = numeric(length = 0), site_nucl_div = numeric(length = 0), site_nucl_div_maf = numeric(length = 0))
    
   for (j in 1:no_cds) {
     
     #creating a snp data file for each cds in the virus/segment
     cds_info<-cds_positions[j,]
     
     cds_name<-gsub('ID=|;','',cds_info$id)
     cds_snps<-as.data.frame(whole_virus_snps[whole_virus_snps$position>=cds_info$start & whole_virus_snps$position<=cds_info$end ,])
     
     if (nrow(cds_snps)==0){
       
       print(paste(virusstem,"virus has no snps in the cds",cds_name, sep = " "))
       
     } else {
     
     cds<-rep.int(cds_name,times = length(cds_snps$position))
     cds_snps<-cbind(cds,cds_snps)
     
     #Calculating allele freqs, and making any positions where the allele is at a freq of less than 1% = zero
     #ie. only looking at minor alleles at at least a freq of 1%
     
     cds_no_snps<-length(cds_snps$position)
     
     cds_snp_div<-data.frame(seqname = cds_snps$seqname, cds = cds_snps$cds, position = cds_snps$position, ref = cds_snps$ref, cov = cds_snps$cov, A_count = cds_snps$A, T_count = cds_snps$T, C_count = cds_snps$C, G_count = cds_snps$G, type = cds_snps$type, strand = cds_snps$strand, codon_change = cds_snps$codon_change, AA_change = cds_snps$aa_change, A_freq = numeric(length = cds_no_snps), A_freq_maf = numeric(length = cds_no_snps), T_freq = numeric(length = cds_no_snps), T_freq_maf = numeric(length = cds_no_snps), C_freq = numeric(length = cds_no_snps), C_freq_maf = numeric(length = cds_no_snps), G_freq = numeric(length = cds_no_snps), G_freq_maf = numeric(length = cds_no_snps), site_nucl_div = numeric(length = cds_no_snps), site_nucl_div_maf = numeric(length = cds_no_snps) )
     
     for (k in 1:cds_no_snps){
       
       #calculating allele frequencies and adding them to the data table of snps
       #where one of the minor allele frequencies is less than 1%, I've changed it to 0 in the freq_maf col, so that you have the option to calculate nucleotide diversity only considering alleles with 1% or higher frequency
       
       cds_snp_div$cov[k]->pos_cov
       A_count<-cds_snp_div$A_count[k]
       A_freq<-A_count/pos_cov
       A_count_maf<-ifelse(A_freq < maf, 0 , A_count )
       
       cds_snp_div$A_freq[k]<-A_freq
       cds_snp_div$A_freq_maf[k]<-ifelse(A_freq < maf, 0 , A_freq )
       
       T_count<-cds_snp_div$T_count[k]
       T_freq<-T_count/pos_cov
       T_count_maf<-ifelse(T_freq < maf, 0, T_count)
       
       cds_snp_div$T_freq[k]<-T_freq
       cds_snp_div$T_freq_maf[k]<-ifelse(T_freq < maf, 0, T_freq )
       
       C_count<-cds_snp_div$C_count[k]
       C_freq<-C_count/pos_cov
       C_count_maf<-ifelse(C_freq < maf, 0, C_count )
       
       cds_snp_div$C_freq[k]<-C_freq
       cds_snp_div$C_freq_maf[k]<-ifelse(C_freq < maf, 0, C_freq )
       
       G_count<-cds_snp_div$G_count[k]
       G_freq<-G_count/pos_cov
       G_count_maf<-ifelse(G_freq < maf, 0, G_count )
       
       cds_snp_div$G_freq[k]<-G_freq
       cds_snp_div$G_freq_maf[k]<-ifelse(G_freq < maf, 0, G_freq )
       
       #now using the counts to calculate the nucleotide diversity for that site
       
       total_count<-pos_cov
       total_count_maf<-sum(A_count_maf,T_count_maf,C_count_maf,G_count_maf)
       
       ##calculating nucleotide diversity for the whole sequence with the formula below
       
       ##For a monomorphic site the minor allele freq is 0, so nucleotide diversity is just 0 at that site, and you don't need to add it into the sum
       ##for a triallelic/biallelic site
       #take a list of counts (Ji values) including the major allele
       ##L = total no. of sites 
       ##sum_i = sum over this for each of the count values
       ##Ji = J sub i eg. the count of i allele
       ##sum( sum_i( Ji (n-Ji) ) / (n(n-1)) ) ) / L
       ##So, for a triallelic site w counts of 230, 0, 780, 100, (ie. cov is 1110) the nucleotide diversity for that site would be;
       ##sum_i( Ji (n-Ji) ) / (n(n-1)) )
       ##sum( (230*(1110-230))/(1110*(1110-1)),(0*(1110-0))/(1110*(1110-1)),(780*(1110-780))/(1110*(1110-1)),(100*(1110-100))/(1110*(1110-1)))
       
       cds_snp_div$site_nucl_div[k]<-sum( (A_count*(total_count-A_count))/(total_count*(total_count-1)),(T_count*(total_count-T_count))/(total_count*(total_count-1)),(C_count*(total_count-C_count))/(total_count*(total_count-1)),(G_count*(total_count-G_count))/(total_count*(total_count-1)))
       
       cds_snp_div$site_nucl_div_maf[k]<-sum( (A_count_maf*(total_count_maf-A_count_maf))/(total_count_maf*(total_count_maf-1)),(T_count_maf*(total_count_maf-T_count_maf))/(total_count_maf*(total_count_maf-1)),(C_count_maf*(total_count_maf-C_count_maf))/(total_count_maf*(total_count_maf-1)),(G_count_maf*(total_count_maf-G_count_maf))/(total_count_maf*(total_count_maf-1)))
       
     }
   
     #adding the cds diversity table to the whole virus nucl diversity table 
     virus_nuc_div_table<-rbind(virus_nuc_div_table,cds_snp_div)
     
     }
   
  }
    
    #writing the whole virus nucl diversity table to a file
    write.table(virus_nuc_div_table,file = paste("Vesanto.virus.remapping.",seg,"/",virusstem,".25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv",sep = ""),sep = "\t", dec = ".", row.names = FALSE)
    #write.table(virus_nuc_div_table,file = paste("nonKallithea_viruses_diversity/",virusstem,".25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv",sep = ""),sep = "\t", dec = ".", row.names = FALSE)
  
  } else {
    
    #import the gtf cds annotation file for the virus/segment variant
    read.table(file = paste(virusstem,".virus.analyses/",virusid,".edited.cds.gtf", sep = ""), header = FALSE, col.names = c("seqname","source","feature","start","end","score","strand","frame","id"), sep = "\t")->cds_positions
    #read.table(file = paste("nonKallithea_viruses_diversity/",virusid,".edited.cds.gtf", sep = ""), header = FALSE, col.names = c("seqname","source","feature","start","end","score","strand","frame","id"), sep = "\t")->cds_positions
    
    no_cds<-length(cds_positions$seqname)
    
    #import the snps for the virus/segment variant
    read.table(file = paste(virusstem,".virus.analyses/",virusstem,".25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.tsv", sep = ""), header = FALSE, col.names = c("seqname","position","ref","cov","A","T","C","G","type","strand","codon_change","aa_change"), sep = "\t")->whole_virus_snps
    #read.table(file = paste("nonKallithea_viruses_diversity/",virusstem,".25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.tsv", sep = ""), header = FALSE, col.names = c("seqname","position","ref","cov","A","T","C","G","type","strand","codon_change","aa_change"), sep = "\t")->whole_virus_snps
    
    #initiating the per virus/segment results table for per site nucl div
    virus_nuc_div_table<-data.frame(seqname = character(length = 0), cds = character(length = 0), position = numeric(length = 0), ref = character(length = 0), cov = numeric(length = 0), A_count = numeric(length = 0) , T_count = numeric(length = 0), C_count = numeric(length = 0), G_count = numeric(length = 0), type = character(length = 0), strand = character(length = 0) , codon_change = character(length = 0), AA_change = character(length = 0), A_freq = numeric(length = 0), A_freq_maf = numeric(length = 0), T_freq = numeric(length = 0), T_freq_maf = numeric(length = 0), C_freq = numeric(length = 0), C_freq_maf = numeric(length = 0), G_freq = numeric(length = 0), G_freq_maf = numeric(length = 0), site_nucl_div = numeric(length = 0), site_nucl_div_maf = numeric(length = 0))
    
    for (j in 1:no_cds) {
      
      #creating a snp data file for each cds in the virus/segment
      cds_info<-cds_positions[j,]
      cds_name<-gsub('gene_id |; |transcript_id.*','',cds_info$id)
      cds_snps<-as.data.frame(whole_virus_snps[whole_virus_snps$position>=cds_info$start & whole_virus_snps$position<=cds_info$end ,])
      
      if (nrow(cds_snps)==0){
        
        print(paste(virusstem,"virus has no snps in the cds",cds_name, sep = " "))
        
      } else {
        
        cds<-rep.int(cds_name,times = length(cds_snps$position))
        cds_snps<-cbind(cds,cds_snps)
        
        #Calculating allele freqs, and making any positions where the allele is at a freq of less than 1% = zero
        #ie. only looking at minor alleles at at least a freq of 1%
        
        cds_no_snps<-length(cds_snps$position)
        
        cds_snp_div<-data.frame(seqname = cds_snps$seqname, cds = cds_snps$cds, position = cds_snps$position, ref = cds_snps$ref, cov = cds_snps$cov, A_count = cds_snps$A, T_count = cds_snps$T, C_count = cds_snps$C, G_count = cds_snps$G, type = cds_snps$type, strand = cds_snps$strand, codon_change = cds_snps$codon_change, AA_change = cds_snps$aa_change, A_freq = numeric(length = cds_no_snps), A_freq_maf = numeric(length = cds_no_snps), T_freq = numeric(length = cds_no_snps), T_freq_maf = numeric(length = cds_no_snps), C_freq = numeric(length = cds_no_snps), C_freq_maf = numeric(length = cds_no_snps), G_freq = numeric(length = cds_no_snps), G_freq_maf = numeric(length = cds_no_snps), site_nucl_div = numeric(length = cds_no_snps), site_nucl_div_maf = numeric(length = cds_no_snps) )
        
        for (k in 1:cds_no_snps){
          
          #calculating allele frequencies and adding them to the data table of snps
          #where one of the minor allele frequencies is less than 1%, I've changed it to 0 in the freq_maf col, so that you have the option to calculate nucleotide diversity only considering alleles with 1% or higher frequency
          
          cds_snp_div$cov[k]->pos_cov
          A_count<-cds_snp_div$A_count[k]
          A_freq<-A_count/pos_cov
          A_count_maf<-ifelse(A_freq < maf, 0.000 , A_count )
          
          cds_snp_div$A_freq[k]<-A_freq
          cds_snp_div$A_freq_maf[k]<-ifelse(A_freq < maf, 0.000 , A_freq )
          
          T_count<-cds_snp_div$T_count[k]
          T_freq<-T_count/pos_cov
          T_count_maf<-ifelse(T_freq < maf, 0.000, T_count)
          
          cds_snp_div$T_freq[k]<-T_freq
          cds_snp_div$T_freq_maf[k]<-ifelse(T_freq < maf, 0.000, T_freq )
          
          C_count<-cds_snp_div$C_count[k]
          C_freq<-C_count/pos_cov
          C_count_maf<-ifelse(C_freq < maf, 0.000, C_count )
          
          cds_snp_div$C_freq[k]<-C_freq
          cds_snp_div$C_freq_maf[k]<-ifelse(C_freq < maf, 0.000, C_freq )
          
          G_count<-cds_snp_div$G_count[k]
          G_freq<-G_count/pos_cov
          G_count_maf<-ifelse(G_freq < maf, 0.000, G_count )
          
          cds_snp_div$G_freq[k]<-G_freq
          cds_snp_div$G_freq_maf[k]<-ifelse(G_freq < maf, 0.000, G_freq )
          
          #now using the counts to calculate the nucleotide diversity for that site
          
          total_count<-pos_cov
          total_count_maf<-sum(A_count_maf,T_count_maf,C_count_maf,G_count_maf)
          
          cds_snp_div$site_nucl_div[k]<-sum( (A_count*(total_count-A_count))/(total_count*(total_count-1)),(T_count*(total_count-T_count))/(total_count*(total_count-1)),(C_count*(total_count-C_count))/(total_count*(total_count-1)),(G_count*(total_count-G_count))/(total_count*(total_count-1)))
          
          cds_snp_div$site_nucl_div_maf[k]<-sum( (A_count_maf*(total_count_maf-A_count_maf))/(total_count_maf*(total_count_maf-1)),(T_count_maf*(total_count_maf-T_count_maf))/(total_count_maf*(total_count_maf-1)),(C_count_maf*(total_count_maf-C_count_maf))/(total_count_maf*(total_count_maf-1)),(G_count_maf*(total_count_maf-G_count_maf))/(total_count_maf*(total_count_maf-1)))
          
        }
        
        #adding the cds diversity table to the whole virus nucl diversity table 
        virus_nuc_div_table<-rbind(virus_nuc_div_table,cds_snp_div)
        
      }
      
    }
  
    #writing the whole virus nucl diversity table to a file
    write.table(virus_nuc_div_table,file = paste(virusstem,".virus.analyses/",virusstem,".25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv",sep = ""),sep = "\t", dec = ".", row.names = FALSE)
    #write.table(virus_nuc_div_table,file = paste("nonKallithea_viruses_diversity/",virusstem,".25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv",sep = ""),sep = "\t", dec = ".", row.names = FALSE) 
      
  }
  
}
