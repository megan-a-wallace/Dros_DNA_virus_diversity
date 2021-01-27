# Dros_DNA_virus_diversity

A collection of scripts to analyse Drosophila DNA viruses, found in the collections and sequencing conducted by the DrosEU consortium

Scripts should be run in the same folder as competitively mapped bam files against Drosophila viruses and contaminants, with the pre-processing script run first. 

The pre-processing script takes the two included fasta files (virusgenomes.fasta & virusgenomes_nonTIRs.fasta) as input, and generates all further files needed for the analyses (excluding virus genome gtf and gff files)

After running the pre-processing script, run the diversity calculations shell script to calculate piA, piS and intergenic pi for both local samples and total populations.

The combined Kallithea virus diversity script uses some of these output files to create a figure showing total piS and intergenic pi across the Kallithea virus genome, and the distribution of indel incidence across the virus genome. 

The INLA model code takes the input file with relative mapped reads for each virus, and builds models of spatial and temporal variation in prevalence, with relevant plots of the model output and spatial fields
