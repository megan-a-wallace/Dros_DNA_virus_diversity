############################################
## Kallithea virus genomic diversity plot ##
############################################

##Author : Megan Wallace
##September 2020
##Figure showing variation in total piS and intergenic pi across the Kallithea genome, above variation in the % of samples showing InDel support

#setwd("C:/Users/User/Dropbox/PhD - 1st Year/popgen_drosEU")
setwd("C:/Users/s1667991/Dropbox/PhD - 1st Year/popgen_drosEU")

require(tidyverse); require(evobiR); require(seqinr)

virusid="KX130344_Kallithea_virus"
virusstem<-gsub('[A-Z]{2}[0-9]{6}_|_virus','',virusid)

#####################################
## Importing the genome info (gff) ##
#####################################

##Importing genome data (location of CDS) and list of samples
refseqlength=152388
read.table(file = "nonKallithea_viruses_diversity/KX130344_Kallithea_virus.edited.cds.gtf", col.names = c("seqname","source","feature","start","end","score","strand","frame","id"), sep = "\t")->Kallithea_cds_positions
read.table(file = "nonKallithea_viruses_diversity/Kallithea.mpileup.25.bam.list.txt", sep = "\t")->Kallithea_sample_list
gsub('Kallithea.virus.analyses/|.Kallithea.bwa.noTIRs_InDel.bam','',Kallithea_sample_list$V1)->Kallithea_sample_list

##############################
## Importing the InDel data ##
##############################

###Importing InDel data 
##Combined table of indellic locations across all Kallithea positive samples in a binary format, coded as support for a gap location if at least 5 reads support it in that sample.
indel.positions=read.table(file = "nonKallithea_viruses_diversity/Kallithea.25.combined.indel.positions.txt", sep = "\t", col.names = c("virus","sample","position","indel_support"))

##############################
## Preparing the InDel Data ##
##############################

#Labelling each position in the indel incidence data as cds or intergenic 
indel.positions$cds<-character(length = length(indel.positions$position))
indel.positions$cds<-ifelse(sapply(indel.positions$position, function(p) any(Kallithea_cds_positions$start <= p & Kallithea_cds_positions$end >= p)),"cds","ig")

#creating a dataframe showing the number of samples with indel support
indel.positions$sample<-as.factor(indel.positions$sample)
no_samples<-as.numeric(nlevels(indel.positions$sample))

indel_sums <- indel.positions %>%
  group_by(position) %>%
  summarise(sum_indel = sum(indel_support))
indel_sums$cds<-character(length = length(indel_sums$position))
indel_sums$cds<-ifelse(sapply(indel_sums$position, function(p) any(Kallithea_cds_positions$start <= p & Kallithea_cds_positions$end >= p)),"cds","ig")

indel_window<-10
indel_perc <- cbind(indel_sums,indel_sums$sum_indel/no_samples)
sliding_indel_perc<-SlidingWindow("mean",indel_perc[,4],indel_window,1)
sliding_indel_pos<-seq(median(seq(1,indel_window,by=1)),refseqlength-(indel_window-median(seq(1,indel_window,by=1))),by = 1)
sliding_indel_data<-as.data.frame(cbind(sliding_indel_pos,sliding_indel_perc))
colnames(sliding_indel_data)<-c("pos","perc")

#and preparing the data which will become the plotted polygons for intergenic and coding regions
#adding a col into the sliding perc data frame to colour the line by CDS or non-CDS status
cds_col<-"#E69F00"
ig_col<-"darkblue"
sliding_indel_data$col<-character(length = length(sliding_indel_data$pos))
sliding_indel_data$col<-ifelse(sapply(sliding_indel_data$pos, function(p) any(Kallithea_cds_positions$start <= p & Kallithea_cds_positions$end >= p)),cds_col,ig_col)

#setting up the data so theres a separate percentage col for ig regions and cds regions
sliding_indel_data_polygon<-as.data.frame(cbind(sliding_indel_data$pos,sliding_indel_data$col))
colnames(sliding_indel_data_polygon)<-c("pos","col")
sliding_indel_data_polygon$perc_intergenic<-sliding_indel_data$perc
sliding_indel_data_polygon$perc_cds<-sliding_indel_data$perc

sliding_indel_data_polygon$perc_intergenic[sliding_indel_data_polygon$col==cds_col]<-0
sliding_indel_data_polygon$perc_cds[sliding_indel_data_polygon$col==ig_col]<-0

###########################################
## Inspecting the distribution of InDels ##
###########################################

#Adding a binary col to the sums data frame (to indicate that an InDel is found at a position at least once), then using this to investigate whether indel incidence is more likely in intergenic or coding regions of the genome
indel_sums<-data.frame(indel_sums,binary_support=numeric(length = length(indel_sums$position)))
indel_sums$binary_support[indel_sums$sum_indel>=1]<-1
indel_sums$no_samples=as.numeric(rep.int(no_samples,times = length(indel_sums$binary_support)))

no_intergenic_indels<-as.numeric(length(indel_sums$binary_support[indel_sums$binary_support==1 & indel_sums$cds=="ig"]))
no_cds_indels<-as.numeric(length(indel_sums$binary_support[indel_sums$binary_support==1 & indel_sums$cds=="cds"]))

no_cds_non_indels<-as.numeric(length(indel_sums$binary_support[indel_sums$binary_support==0 & indel_sums$cds=="cds"]))
no_intergenic_non_indels<-as.numeric(length(indel_sums$binary_support[indel_sums$binary_support==0 & indel_sums$cds=="ig"]))

#Now using chi squared test to examine whether rows and columns of the contingency table are statistically independent or not - eg. is the distribution of indels and non-indels in the genome independent of the distribution of cds and intergenic sites

#Our data table is the no_cds indels, no_integenic indels, no_cds non indels, no_intergenic non indels
chisq_test_data<-matrix(data = c(no_intergenic_indels,no_cds_indels,no_intergenic_non_indels,no_cds_non_indels), nrow = 2, ncol = 2, byrow = TRUE)
rownames(chisq_test_data)<-c("indels","non-indels")
colnames(chisq_test_data)<-c("intergenic","cds")
Kallithea_indel_chisq<-chisq.test(chisq_test_data)

#to look at the residuals 
Kallithea_indel_chisq$residuals

#Chi-square test for independence on the number of indels found in cds and intergenic regions, compared to the expected numbers based on the number of total cds and intergenic sites found a strong positive association between intergenic regions and finding indels (X-squared = 3236, df = 1, p-value < 2.2e-16) 

######################################
## Importing the piS and intpi data ##
######################################

##Importing total piS and intergenic pi data 
# intergenic and syn-nsyn per position nucleotide diversity files for total/merged population 
read.table(file = "nonKallithea_viruses_diversity/Kallithea.25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv", sep = "\t", header = TRUE)->Kallithea.per.position.syn.nsyn.nucl.div
read.table(file = "nonKallithea_viruses_diversity/Kallithea.wholepop.intergenic.per.site.nucl.div.tsv", sep = "\t", header = TRUE)->Kallithea.per.position.intergenic.nucl.div

######################################
## Preparing the piS and intpi data ##
######################################

#making data frame of nucleotide diversity fr the total/merged population
data.frame(position = as.numeric(seq(1,refseqlength,by = 1)), type = character(length = refseqlength), nucl_div = numeric(length = refseqlength))->Kallithea.nucl.div.data
Kallithea.nucl.div.data$type<-ifelse(sapply(Kallithea.nucl.div.data$position, function(p) any(Kallithea_cds_positions$start <= p & Kallithea_cds_positions$end >= p)),"cds","ig")
#putting the nucleotide diversity from intergenic sites into the combined data frame
Kallithea.nucl.div.data$nucl_div[Kallithea.per.position.intergenic.nucl.div$position[Kallithea.per.position.intergenic.nucl.div$site_nucl_div_maf>0]]<-Kallithea.per.position.intergenic.nucl.div$site_nucl_div_maf[Kallithea.per.position.intergenic.nucl.div$site_nucl_div_maf>0]
#putting the nucleotide diversity from the synonymous coding sites into the combined data frame 
Kallithea.nucl.div.data$nucl_div[Kallithea.per.position.syn.nsyn.nucl.div$position[Kallithea.per.position.syn.nsyn.nucl.div$site_nucl_div_maf>0 & Kallithea.per.position.syn.nsyn.nucl.div$type=="syn"]]<-Kallithea.per.position.syn.nsyn.nucl.div$site_nucl_div_maf[Kallithea.per.position.syn.nsyn.nucl.div$site_nucl_div_maf>0 & Kallithea.per.position.syn.nsyn.nucl.div$type=="syn"]

#Figuring out the number of synonymous sites per codon for the denominator of the sliding window
#importing the non-synonymous length table for Kallithea virus (generated using popoolation)
snl.table=read.table("nonKallithea_viruses_diversity/snl.Kallithea.txt",header = FALSE,col.names = c("codon","nsyn_length"),skip = 5,sep = ":")
snl.table$syn_length<-(3-snl.table$nsyn_length)
snl.table<-snl.table[,c(1,3)]#removing the nsyn length values
cds.seqs=read.fasta("nonKallithea_viruses_diversity/Kallithea.cds.regions.fasta", forceDNAtolower = FALSE, strip.desc = TRUE)
ref.seq=read.fasta("nonKallithea_viruses_diversity/KX130344_Kallithea_virus.fasta",forceDNAtolower = FALSE, strip.desc = TRUE)
refseq_mat<-matrix(data = ref.seq$KX130344_Kallithea_virus, nrow = refseqlength, ncol = 1, byrow = TRUE)

#initialising data frame for results
syn_sites_per_codon<-data.frame(position = seq(1,refseqlength,1), cds = rep.int("intergenic",times = refseqlength), sequence = refseq_mat, codon = character(length = refseqlength), int_syn_sites = rep.int(1,times = refseqlength), stringsAsFactors = FALSE)

#Adding the names of the genes/intergenic to the data frame from the gtf gene positions
for (k in 1:length(Kallithea_cds_positions$start)){
  
  Kallithea_cds_positions$start[k]->start
  Kallithea_cds_positions$end[k]->end
  gsub("gene_id ","",Kallithea_cds_positions$id[k])->name
  gsub("; transcript_id.*","",name)->name
  
  syn_sites_per_codon$cds<-ifelse(sapply(syn_sites_per_codon$position, function(p) any(start <= p & end >= p)),name,syn_sites_per_codon$cds)
  
}

#adding the codons that each position is a part of into the table
refseq_codons<-matrix(data = refseq_mat, nrow = refseqlength/3, ncol = 3, byrow = TRUE)
refseq_codons<-paste(refseq_codons[,1], refseq_codons[,2], refseq_codons[,3], sep = "")
codon_vec<-vector()

for (m in 1:length(refseq_codons)) {
  rep.int(refseq_codons[m],times = 3)->three_codons
  codon_vec<-c(codon_vec,three_codons)
}
syn_sites_per_codon$codon<-codon_vec

#adding the synonymous length of each coding codon into the table
syn_sites_per_codon$int_syn_sites[syn_sites_per_codon$cds!="intergenic"]<-snl.table$syn_length[match(syn_sites_per_codon$codon,snl.table$codon)]
#and now dividing the coding site syn lengths by 3 so its per site rather than per codon
syn_sites_per_codon$int_syn_sites_per_position<-1
syn_sites_per_codon$int_syn_sites_per_position[syn_sites_per_codon$cds!="intergenic"]<-syn_sites_per_codon$int_syn_sites[syn_sites_per_codon$cds!="intergenic"]/3

####And now calculating a slding window of piS/intergenic pi along the Kallithea genome
#edited a function I found on this guthub so that instead of a mean, it would divide by the sum of number of syn and intergenic sites
#http://coleoguy.blogspot.com/2014/04/sliding-window-analysis.html
#takes a vector of nucleotide diversities, and vector of per position intergenic or synonymous site numbers
slideFunct <- function(nucl_div, sites, window, step){
  total <- length(nucl_div)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- sum(nucl_div[spots[i]:(spots[i]+window)])/sum(sites[spots[i]:(spots[i]+window)])
  }
  return(result)
}

##Creating two sliding window data series, one with a windw size of 1000, and another with 5000
window<-1000
step<-200
sliding_nucl_div<-slideFunct(Kallithea.nucl.div.data$nucl_div,syn_sites_per_codon$int_syn_sites_per_position,window = window,step = step)
sliding_nucl_div_pos<-seq(median(seq(1,window,by=1)),refseqlength-(window-median(seq(1,window,by=1))),by = step)
sliding_nucl_div_data<-as.data.frame(cbind(sliding_nucl_div_pos,sliding_nucl_div))
colnames(sliding_nucl_div_data)<-c("pos","pi")

##And a second line on the plot with a 5kb window, to show larger area patterns 
window<-5000
step<-1000
sliding_nucl_div_2<-slideFunct(Kallithea.nucl.div.data$nucl_div,syn_sites_per_codon$int_syn_sites_per_position,window = window,step = step)
sliding_nucl_div_pos_2<-seq(median(seq(1,window,by=1)),refseqlength-(window-median(seq(1,window,by=1))),by = step)

sliding_nucl_div_data_2<-as.data.frame(cbind(sliding_nucl_div_pos_2,sliding_nucl_div_2))
colnames(sliding_nucl_div_data_2)<-c("pos","pi")

###############################
## And now creating the plot ##
###############################

###sliding piS and intergenic pi
################################

#making some fake data to set up the plot
runif(refseqlength,min(sliding_nucl_div_data$pi),max(sliding_nucl_div_data$pi))->nucl_div_fake 
seq(from=1,to=refseqlength,by=1)->pos_fake
fake_data_nucl_div<-data.frame(pos_fake,nucl_div_fake)

par(mar = c(3, 5.5, 0.5, 0.5), # change the margins
    lwd = 1.5,# increase the line thickness
    cex.axis = 1.2, # increase default axis label size
    cex.lab = 1.3)

par(mfrow=c(2,1))#creating two plots, one on top of the other

#plotting the fake data
plot(fake_data_nucl_div$pos_fake,fake_data_nucl_div$nucl_div_fake,axes = FALSE, main="", xlab="",
     ylab="", pch = '')

#Adding grey boxes to the plot to indicate non-coding regions
#rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"),.)
rect(0, 0, 152388, 0.01, border = NA, col = "azure2")
#cds in white 
rect_xmin<-Kallithea_cds_positions$start
rect_xmax<-Kallithea_cds_positions$end
rect(rect_xmin, 0,rect_xmax, 0.01, border = NA, col = "white")

#adding y axis
axis(2, at = seq(0,0.008,by = 0.002), tick = TRUE, pos = par("usr")[1]+4500,labels = FALSE)
text(x = par("usr")[1]+2000,
     y = seq(0,0.008,by = 0.002),
     labels = seq(0,0.008,by = 0.002),
     adj = 1,
     xpd = NA,
     cex = 1.2)
mtext(side = 2, line = 2.75, "pi at intergenic & \nsynonymous sites", cex = 1.3)

#adding x axis
vec<-c(0,15,30,45,60,75,90,105,120,135,150,153)
axis(1, at = vec*1000, tick = TRUE, pos = par("usr")[3]+0.0001, labels = FALSE)
text(x = vec*1000,
     y = par("usr")[3]-0.0004,
     labels = c(vec[1:11],""),
     cex = 1.2,
     xpd = NA)

#plotting the real data
lines(sliding_nucl_div_data$pos,sliding_nucl_div_data$pi,lwd = 1.95,col="darkblue")

#plotting the real data with a larger sliding window
lines(sliding_nucl_div_data_2$pos,sliding_nucl_div_data_2$pi,lwd = 3.7,col="orange")

#Adding legend to the plot
legend_lines_coords_ig<-matrix(data = c(vec[1]*1000+600,0.0062,vec[2]*1000+600,0.0062) ,nrow = 2,ncol = 2,byrow = TRUE)
legend_lines_coords_cds<-matrix(data = c(vec[1]*1000+600,0.0058,vec[2]*1000+600,0.0058) ,nrow = 2,ncol = 2,byrow = TRUE)
lines(legend_lines_coords_ig, col=ig_col, lwd=1.95)
lines(legend_lines_coords_cds, col=cds_col, lwd=3.7)

legend(vec[2]*1000-2000,0.0066,legend=c("1000 bp window","5000 bp window"),plot=T,bty="n",cex = 0.9, y.intersp = 1, adj = c(0,0.5))

####And the second, lower panel : distribution of indels along the Kallithea virus genome 

#making some fake data to set up plot
runif(refseqlength,0,100)->perc_fake 
seq(1,refseqlength,by=1)->pos_fake_indel
fake_data_perc<-data.frame(pos_fake_indel,perc_fake)

#plotting the fake data
plot(fake_data_perc$pos_fake_indel,fake_data_perc$perc_fake,axes = FALSE, main="", xlab="",
     ylab="", pch = '')

#Adding grey boxes to the plot to indicate non-coding regions
rect(0, 0, 152388, 100, border = NA, col = "azure2")
#cds in white 
rect_xmin<-Kallithea_cds_positions$start
rect_xmax<-Kallithea_cds_positions$end
rect(rect_xmin, 0,rect_xmax, 100, border = NA, col = "white")

#adding y axis
axis(2, at = seq(0,100,by = 20), tick = TRUE, pos = par("usr")[1]+4500,labels = FALSE)
text(x = par("usr")[1]+2000,
     y = seq(0,100,by = 20),
     labels = seq(0,100,by = 20),
     ## Rotate the labels by 35 degrees.
     xpd = NA,
     adj = 1,
     cex = 1.2)
mtext(side = 2, line = 2.75, "% of samples \nwith gap support", cex = 1.3)

#adding x axis
axis(1, at = vec*1000, tick = TRUE, pos = par("usr")[3]+2, labels = FALSE)
text(x = vec*1000,
     y = par("usr")[3]-6,
     labels = c(vec[1:11],""),
     cex = 1.2,
     xpd = NA)
mtext(side = 1, line = 1.75, "Position in reference genome (Kb)", cex = 1.3)

#plotting the real data
#with polygons
#plotting two polygons, one for the intergenic regions and one for the CDS
polygon(sliding_indel_data_polygon$pos,sliding_indel_data_polygon$perc_cds*100,col = cds_col, border = cds_col, lwd = 1.8)
polygon(sliding_indel_data_polygon$pos,sliding_indel_data_polygon$perc_intergenic*100,col = ig_col, border = ig_col, lwd = 1.8)

#Adding legend to the plot
legend_lines_coords_ig<-matrix(data = c(vec[1]*1000+600,94,vec[2]*1000+600,94) ,nrow = 2,ncol = 2,byrow = TRUE)
legend_lines_coords_cds<-matrix(data = c(vec[1]*1000+600,86,vec[2]*1000+600,86) ,nrow = 2,ncol = 2,byrow = TRUE)
lines(legend_lines_coords_ig, col=ig_col, lwd=2.5)
lines(legend_lines_coords_cds, col=cds_col, lwd=2.5)

legend(vec[2]*1000-2000,par("usr")[4]-1.75,legend=c("Intergenic regions","CDS"),plot=T,bty="n",cex = 0.9, y.intersp = 1.15, adj = c(0,0.5))

##moved plots a little closer together in inkscape and added pi greek symbol to axis