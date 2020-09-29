#counting the number of synonymous and non-synonymous sites for each cds for the virus of interest
#run as Rscript syn_nsyn_sites.R $virusstem.virus.analyses/$virusstem.cds.regions.fasta $virusstem.virus.analyses/snl.$virusstem.txt $virusstem
library(dplyr)
library(stringr)
library(seqinr)
args = commandArgs(trailingOnly=TRUE)

##for shell
#reading in the snl table
snl.table=read.table(args[2],header = FALSE,col.names = c("codon","nsyn_length"),skip = 5,sep = ":")
snl.table$syn_length<-(3-snl.table$nsyn_length)
#setting the virusstem from the shell variables
virusstem<-args[3]
#reading in the fasta file of cds sequences
cds.seqs=read.fasta(args[1],forceDNAtolower = FALSE, strip.desc = TRUE)

#for testing
#snl.table=read.table("nonKallithea_viruses_diversity/snl.Kallithea.txt",header = FALSE,col.names = c("codon","nsyn_length"),skip = 5,sep = ":")
#snl.table$syn_length<-(3-snl.table$nsyn_length)
#virusstem<-"Kallithea"
#cds.seqs=read.fasta("nonKallithea_viruses_diversity/Kallithea.cds.regions.fasta", forceDNAtolower = FALSE, strip.desc = TRUE)

length(cds.seqs)->no_cds

#initialising data frame for results
syn_nsyn_sites_per_cds<-data.frame(cds = character(length = no_cds), nsyn_sites = numeric(length = no_cds), syn_sites = numeric(length = no_cds), stringsAsFactors = FALSE)

for (i in 1:length(cds.seqs))
{
  #isolating the DNA sequence of the cds of interest
  cds.seqs[[i]]->cds_seq
  
  #the name of the cds of interest
  gsub("\\(.*\\)","",as.character(attr(cds.seqs[[i]],"name")))->cds_name
  
  #creating a matrix of codons for the cds
  length(cds_seq)->cds_length
  codon<-matrix(data = cds_seq, nrow = cds_length/3, ncol = 3, byrow = TRUE)
  codon<-paste(codon[,1], codon[,2], codon[,3], sep = "")
  
  #counting the number of each codon in the cds, and appending the syn and nsyn length for that codon to the row
  codon_freqs<-as.data.frame(table(codon))
  codon_freqs<-merge(codon_freqs,snl.table,by.x = "codon",by.y = "codon",all.x = TRUE,all.y = FALSE)
  #Calculating the no of avg. syn and non-syn sites for that cds
  nsyn_sites<-sum(codon_freqs$Freq*codon_freqs$nsyn_length)
  syn_sites<-sum(codon_freqs$Freq*codon_freqs$syn_length)
  
  #adding to the data table
  syn_nsyn_sites_per_cds$cds[i]<-cds_name
  syn_nsyn_sites_per_cds$nsyn_sites[i]<-nsyn_sites
  syn_nsyn_sites_per_cds$syn_sites[i]<-syn_sites
}

if (str_detect(virusstem,"Vesanto")){
  
  seg<-gsub('^.*_S','S',virusstem)

  write.table(syn_nsyn_sites_per_cds,file = paste("Vesanto.virus.remapping.",seg,"/",virusstem,".genewise.syn-nsyn.sites.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)

} else {
  
  write.table(syn_nsyn_sites_per_cds,file = paste(virusstem,".virus.analyses/",virusstem,".genewise.syn-nsyn.sites.tsv",sep = ""), sep = "\t", dec = ".", row.names = FALSE)
  
}

