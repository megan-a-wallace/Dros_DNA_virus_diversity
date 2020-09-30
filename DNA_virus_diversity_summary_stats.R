####################################################
## Calculating Drosophila DNA virus summary stats ##
####################################################

##Written by Megan A. Wallace
##07-09/2020

#Script to calculate diversity summary stats for Kallithea, Vesanto and Linvill Road virus samples from European Drosophila
##Creates table of results for summary stats, combined table of genewise piA & piS, and for Vesanto, generates the SNP totals 

require(dplyr); require(stringr); require(seqinr)

##importing total population data for all viruses and haplotypes
total_piA_piS_data<-read.table(file = "nonKallithea_viruses_diversity/all_viruses_whole_genome_piA_piS.tsv",col.names = c("seqname","type","nsyn_sites","syn_sites","nsyn_snps","syn_snps","piA","piS","piA_over_piS"))
total_intergenic_pi_data<-read.table(file = "nonKallithea_viruses_diversity/all_viruses_intergenic_pi.tsv",col.names = c("seqname","type","intergenic_sites","intergenic_snps","intergenic_pi"))

#and total genewise diversity stats
Kallithea_genewise_data<-read.table(file = "nonKallithea_viruses_diversity/KX130344_Kallithea_virus.genewise.wholepop.piA_piS.tsv",header = TRUE)
Linvill_Road_genewise_data<-read.table(file = "nonKallithea_viruses_diversity/KX648536_Linvill_Road_virus.genewise.wholepop.piA_piS.tsv",header = TRUE)
#Made a combined table of the per gene piA and piS for all vesanto segments in diversity calculations.sh script
Vesanto_segs_genewise_data<-read.table(file = "nonKallithea_viruses_diversity/Vesanto.segs.genewise.piA.piS.results.tsv", col.names = c("seqname","cds","nsyn_sites","syn_sites","nsyn_snps","syn_snps","piA","piS","piA_over_piS"))

###################################
## Kallithea virus summary stats ##
###################################

#importing Kallithea virus local (per sample) diversity statistics
Kallithea_intergenic_data<-read.table(file = "nonKallithea_viruses_diversity/KX130344_Kallithea_virus.per.sample.intergenic.pi.tsv", header = TRUE)
Kallithea_piA_piS_data<-read.table(file = "nonKallithea_viruses_diversity/KX130344_Kallithea_virus.per.sample.piA_piS.tsv",header = TRUE)

Kallithea_total_intergenic_pi<-signif(total_intergenic_pi_data$intergenic_pi[total_intergenic_pi_data$seqname=="KX130344_Kallithea_virus"],digits = 4)
Kallithea_total_piS<-signif(total_piA_piS_data$piS[total_piA_piS_data$seqname=="KX130344_Kallithea_virus"],digits = 4)
Kallithea_total_piA_over_piS<-signif(total_piA_piS_data$piA_over_piS[total_piA_piS_data$seqname=="KX130344_Kallithea_virus"],digits = 4)

Kallithea_mean_local_piS<-signif(mean(Kallithea_piA_piS_data$piS),digits = 4)
Kallithea_mean_local_intergenic_pi<-signif(mean(Kallithea_intergenic_data$pi),digits = 4)
Kallithea_mean_local_piA_over_piS<-signif(mean(Kallithea_piA_piS_data$piA)/mean(Kallithea_piA_piS_data$piS),digits = 4)

######################################
## Linvill Road virus summary stats ##
######################################

#Linvill Road virus summary stats
Linvill_Road_intergenic_data<-read.table(file = "nonKallithea_viruses_diversity/KX648536_Linvill_Road_virus.per.sample.intergenic.pi.tsv", header = TRUE)
Linvill_Road_piA_piS_data<-read.table(file = "nonKallithea_viruses_diversity/KX648536_Linvill_Road_virus.per.sample.piA_piS.tsv",header = TRUE)

Linvill_Road_total_intergenic_pi<-signif(total_intergenic_pi_data$intergenic_pi[total_intergenic_pi_data$seqname=="KX648536_Linvill_Road_virus"],digits = 4)
Linvill_Road_total_piS<-signif(total_piA_piS_data$piS[total_piA_piS_data$seqname=="KX648536_Linvill_Road_virus"],digits = 4)
Linvill_Road_total_piA_over_piS<-signif(total_piA_piS_data$piA_over_piS[total_piA_piS_data$seqname=="KX648536_Linvill_Road_virus"],digits = 4)

Linvill_Road_mean_local_piS<-signif(mean(Linvill_Road_piA_piS_data$piS),digits = 4)
Linvill_Road_mean_local_intergenic_pi<-signif(mean(Linvill_Road_intergenic_data$pi),digits = 4)
Linvill_Road_mean_local_piA_over_piS<-signif(mean(Linvill_Road_piA_piS_data$piA)/mean(Linvill_Road_piA_piS_data$piS),digits = 4)

###########################################################
## Vesanto virus summary stats (all segs and haplotypes) ##
###########################################################

filter(total_intergenic_pi_data, grepl('Vesanto', seqname))->Vesanto_intergenic_data
filter(total_piA_piS_data, grepl('Vesanto', seqname))->Vesanto_piA_piS_data
#Mean across segment variants for total piS
Vesanto_segs_total_piS<-signif(mean(Vesanto_piA_piS_data$piS),digits = 4)
#Mean across segment variants for total intergenic pi
Vesanto_segs_total_intergenic_pi<-signif(mean(Vesanto_intergenic_data$intergenic_pi),digits = 4)
#Mean across segment variants for total piA/piS (ie. mean total piA across segs/ mean total piS across segs)
Vesanto_segs_total_piA_over_piS<-signif(mean(Vesanto_piA_piS_data$piA)/mean(Vesanto_piA_piS_data$piS),digits = 4)

#mean local values for all Vesanto segs 
Vesanto_segs_intergenic_data<-read.table(file = "nonKallithea_viruses_diversity/Vesanto.segs.per.sample.intergenic.pi.results.tsv", col.names = c("seqname","sample","intergenic_sites","no_snps","pi") )
Vesanto_segs_piA_piS_data<-read.table(file = "nonKallithea_viruses_diversity/Vesanto.segs.per.sample.piA.piS.results.tsv", col.names = c("seqname","sample","nsyn_sites","syn_sites","nsyn_snps","syn_snps","piA","piS","piA_over_piS"))
#mean local piS across all vesanto segments 
Vesanto_segs_mean_local_piS<-signif(mean(Vesanto_segs_piA_piS_data$piS), digits = 4)
#mean local intergenic pi across all vesanto segments 
Vesanto_segs_mean_local_intergenic_pi<-signif(mean(Vesanto_segs_intergenic_data$pi), digits = 4)
Vesanto_segs_mean_local_piA_over_piS<-signif(mean(Vesanto_segs_piA_piS_data$piA)/mean(Vesanto_segs_piA_piS_data$piS), digits = 4)

##########################################
## PRINTING SUMMARY STATS RESULTS TABLE ##
##########################################

Dros_DNA_virus_summary_stats_table<-data.frame(virus = as.character(c("Kallithea","Linvill_Road","Vesanto")), total_piS = as.numeric(c(Kallithea_total_piS,Linvill_Road_total_piS,Vesanto_segs_total_piS)), total_intergenic_pi = as.numeric(c(Kallithea_total_intergenic_pi,Linvill_Road_total_intergenic_pi,Vesanto_segs_total_intergenic_pi)), total_piA_over_piS = as.numeric(c(Kallithea_total_piA_over_piS,Linvill_Road_total_piA_over_piS,Vesanto_segs_total_piA_over_piS)), mean_local_piS = as.numeric(c(Kallithea_mean_local_piS,Linvill_Road_mean_local_piS,Vesanto_segs_mean_local_piS)), mean_local_intergenic_pi = as.numeric(c(Kallithea_mean_local_intergenic_pi, Linvill_Road_mean_local_intergenic_pi, Vesanto_segs_mean_local_intergenic_pi)), mean_local_piA_over_piS = as.numeric(c(Kallithea_mean_local_piA_over_piS, Linvill_Road_mean_local_piA_over_piS, Vesanto_segs_mean_local_piA_over_piS)))

write.table(Dros_DNA_virus_summary_stats_table, file = "nonKallithea_viruses_diversity/Dros_DNA_virus_summary_stats_table.tsv", sep = "\t", row.names = FALSE)

#########################################
## Total genewise piA & piS data table ## 
#########################################

total_cds<-sum(length(Kallithea_genewise_data$cds)-1,length(Linvill_Road_genewise_data$cds)-1,length(Vesanto_segs_genewise_data$cds))

genewise_piA_piS_combined_virus_table<-data.frame(virus = as.character(c(rep.int("Kallithea virus",95),rep.int("Linvill Road virus",2),rep.int("Vesanto virus",60))), contig_id = character(length = total_cds), cds = character(length = total_cds), avg_nsyn_sites = numeric(length = total_cds), avg_syn_sites = numeric(length = total_cds), nsyn_snps = numeric(length = total_cds), syn_snps = numeric(length = total_cds), piA = numeric(length = total_cds), piS = numeric(length = total_cds), piA_over_piS = numeric(length = total_cds))

genewise_piA_piS_combined_virus_table[1:length(Kallithea_genewise_data$cds)-1,2:10]<-Kallithea_genewise_data[1:95,]
genewise_piA_piS_combined_virus_table[96:97,2:10]<-Linvill_Road_genewise_data[1:2,]
genewise_piA_piS_combined_virus_table[98:total_cds,2:10]<-Vesanto_segs_genewise_data

write.table(genewise_piA_piS_combined_virus_table, file = "nonKallithea_viruses_diversity/genewise_piA_piS_combined_virus_table.tsv", sep = "\t", row.names = FALSE)

################################
## Vesanto segment SNP totals ##
################################

#importing data table 
read.table(file = "nonKallithea_viruses_diversity/Vesanto.segs.combined.SNP.numbers.tsv", header = FALSE, sep = "\t")->vesanto_local_snp_table
colnames(vesanto_local_snp_table)<-c("seqname", "no_samples","total_SNPs_across_samples","unique_SNPs_across_samples","single_sample_snps","two_sample_snps","gr_than_two_sample_snps")
sum(vesanto_local_snp_table$no_samples)->total_vesanto_samples
sum(vesanto_local_snp_table$unique_SNPs_across_samples)->vesanto_total_unique_local_snps
sum(vesanto_local_snp_table$single_sample_snps)->vesanto_total_local_single_sample_snps
sum(vesanto_local_snp_table$two_sample_snps)->vesanto_local_two_sample_snps
sum(vesanto_local_snp_table$gr_than_two_sample_snps)->vesanto_local_gr_two_sample_snps

sum(Vesanto_intergenic_data$intergenic_snps)->total_vesanto_intergenic_snps
sum(Vesanto_piA_piS_data$nsyn_snps)->total_vesanto_nsyn_snps
sum(Vesanto_segs_piA_piS_data$syn_snps)->total_vesanto_syn_snps
total_vesanto_snps<-sum(total_vesanto_intergenic_snps,total_vesanto_nsyn_snps,total_vesanto_syn_snps)
