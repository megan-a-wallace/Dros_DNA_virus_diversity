# Dros_DNA_virus_diversity

A collection of scripts to analyse Drosophila DNA viruses, found in the collections and sequencing conducted by the DrosEU consortium

These script require piccard tools, popoolation, and popoolation2 to be installed (also uses gatk3 and bedops in conda environments)

Scripts should be run in the same folder as competitively mapped bam files against Drosophila viruses, with the pre-processing script run first. 

The pre-processing script takes the two included fasta files (virusgenomes.fasta & virusgenomes_nonTIRs.fasta) as input, and generates the bam files, and other files needed for further analyses (excluding the genome gff and gtf files). It will also generate a final list of analysis ready bam files ($virus.final.sample.list.txt). 

After running the pre-processing script, run the diversity calculations shell script (in the same folder as $virus.final.sample.list.txt) to calculate piA, piS and intergenic pi for both local samples and total populations. 

This calculations script also requires the gtf files for each virus or segment(.popoolation1.gtf for all, plus either .prokka.edited2.cds.gtf or .edited.cds.gtf), which need to be added the .virus.analyses folders created in the pre-processing script. 

This script will call the following R scripts;
 - tstv_ratio.R
 - syn_nsyn_sites.R
 - per_site_nucl_div.R
 - piA_piS.R
 - wholepop_intergenic_pi.R
 - intergenic_sites.R
 - sample_per_site_nucl.div
 - sample_piA_piS.R
 - sample_intergenic_pi.R
  
The results can be summarised into tables using the DNA_virus_diversity_summary_stats.R file, which takes as input the final tables of piA and piS generated in the calculations script.

The combined_Kalltihea_genome_diversity.R script uses some of the calculations script output files to create a figure showing total piS and intergenic pi across the Kallithea virus genome, and the distribution of indel incidence across the virus genome. 

The DNA_virus_prevalence_INLA_models_and_plots.R code takes the Supplementary_File_S4_justGLMM.txt file with relative mapped reads for each virus, and builds models of spatial and temporal variation in prevalence, with relevant plots of the model output and spatial fields
