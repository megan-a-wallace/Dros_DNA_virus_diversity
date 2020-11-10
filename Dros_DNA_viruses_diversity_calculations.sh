#############################################
## Drosophila DNA virus Diversity analyses ##
#############################################

## Written by Megan A. Wallace in 07-09/2020

####Diversity analyses on the DNA viruses found in the DrosEU dataset 
##This script should be run in the same folder $virus.final.sample.list.txt files which have been generated from the preprocessing script, 
#which should have been run before this to generate analysis ready bam files for mpileup in a .virus.analyses file for each virus / segment of interest
#This script will also have created and indexed a fasta file for each of the viruses of interest in this folder (without TIRs), this needs created otherwise
#run as ./Dros_DNA_viruses_diversity_calculations.sh 

for v in $(cat gt5samples.virus.list.txt)
do
	##Setting up the id and file stem for each virus
	virusid=$(echo $v)

	isseg=$(echo $v | grep '^Vesanto')
	if [ -n "$isseg" ]
	then
		virusstem=$(echo $v | sed 's/\(_virus\)//')
	else
		virusstem=$(echo $v | sed -r 's/^[A-Z]{2}[0-9]{6}_//' | sed 's/virus.*/virus/' | sed 's/\(_virus\|_Nudivirus\)//')
	fi
	
	##########################################
	## Generating list of bams for analysis ##
	##########################################

	##For each virus, taking the final.sample.list.txt and excluding the 2015_01 and 02 files if they are in it, as they aren't from Europe
	##then creating two versions of the list, the original (in which all samples have at least a read depth of 10 across 95% of the viral genome) and another
	# with a higher threshold for read depth of 25 to see if this makes any difference qualititatively
	
	refseqlength=$(grep -A1 "$virusid" linear.virus.noTIRs.fasta | awk '/^>/{getline; getline; print}' | awk '{print length($1)}')
	#calculating the bp consitituting 95% of the genome length, and rounding value up to whole number using awk
	perc="0.95"
	bptarget=$(awk -v perc="${perc}" -v refseqlength="$refseqlength" 'BEGIN{bpt=(perc*refseqlength); i=int(bpt); print (bpt-i<0.5)?i:i+1}' )
	
	#First for the Vesanto bams 
	echo making bam lists for the $virusstem mpileup files 
	
	if [ -n "$isseg" ]
	then
		seg=$(echo $virusstem | sed 's/^.*_S/S/g' )
		
		#renaming the original bam lists from the competitive mapping so they aren't mixed up
		LIST10="$virusstem.virus.analyses/$virusstem.compmap.mpileup.10.bam.list.txt"
		if [ -f "$LIST10" ]
		then 
			echo $virusstem compmap 10 fold coverage bam list already exists
		else
			mv $virusstem.virus.analyses/$virusstem.mpileup.10.bam.list.txt $virusstem.virus.analyses/$virusstem.compmap.mpileup.10.bam.list.txt
		fi
		
		LIST25="$virusstem.virus.analyses/$virusstem.compmap.mpileup.25.bam.list.txt"
		if [ -f "$LIST25" ]
		then 
			echo $virusstem compmap 25 fold coverage bam list already exists
		else
			mv $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt $virusstem.virus.analyses/$virusstem.compmap.mpileup.25.bam.list.txt
		fi
		
		#And now making the list of bam files for analysis with only one Vesanto haplotype per segment per sample
		for s in $(cat $virusstem.final.sample.list.txt)
		do 
			popstem=$(echo $s)

			#removing 2015_01 and 02 from sample list and creating bam list w samples with rd 10 across 95% of the virus genome 
			isnonEU=$(echo $s | grep -E '^2015_01|^2015_02' )
			if [ -n "$isnonEU" ]
			then
				echo sample $popstem is not from europe - removed from further analyses
			else
				echo Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs_InDel.bam >> Vesanto.virus.remapping.$seg/$virusstem.mpileup.10.bam.list.txt
			fi
			
			#now creating bam list with just samples where rd is 25 or greater across 95% of the virus genome
			#using bedtools to measure per position read depth and count the positions where read depth is > twenty four
			depthcount=$(genomeCoverageBed -d -ibam Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs_InDel.bam | awk '$3 > 24'| grep -w "$virusid" -c)
			
			echo in the $popstem sample $depthcount bp of the $virusstem virus genome have greater than twenty four fold read depth, $bptarget bp are needed for analyses
			
			if [ -n "$isnonEU" ]
			then
				echo sample $popstem is not from europe - removed from further analyses
			elif [ $depthcount -ge $bptarget ]
			then 
				echo sample $popstem has at least twenty five fold read depth across ninety five percent of the $virusstem virus genome
				echo Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs_InDel.bam >> Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt
			else
				echo sample $popstem does not have twenty five fold $virusstem read depth across ninety five percent of the genome 
			fi	
		done 
		
	else
		
		#and for the non-Vesanto bams...
		
		for s in $(cat $virusstem.final.sample.list.txt)
		do 
			popstem=$(echo $s)

			#removing 2015_01 and 02 from sample list and creating bam list w samples with rd 10 across 95% of the virus genome 
			isnonEU=$(echo $s | grep -E '^2015_01|^2015_02' )
			if [ -n "$isnonEU" ]
			then
				echo sample $popstem is not from europe - removed from further analyses
			else
				echo $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs_InDel.bam >> $virusstem.virus.analyses/$virusstem.mpileup.10.bam.list.txt
			fi
			
			#now creating bam list with just samples where rd is 25 or greater across 95% of the virus genome
			#using bedtools to measure per position read depth and count the positions where read depth is > twenty four
			depthcount=$(genomeCoverageBed -d -ibam $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs_InDel.bam | awk '$3 > 24'| grep -w "$virusid" -c)
			echo in the $popstem sample $depthcount bp of the $virusstem virus genome have greater than twenty four fold read depth, $bptarget bp are needed for analyses
			if [ -n "$isnonEU" ]
			then
				echo sample $popstem is not from europe - removed from further analyses
			elif [ $depthcount -ge $bptarget ]
			then 
				echo sample $popstem has at least twenty five fold read depth across ninety five percent of the $virusstem virus genome
				echo $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs_InDel.bam >> $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt
			else
				echo sample $popstem does not have twenty five fold $virusstem read depth across ninety five percent of the genome 
			fi	
		done 
	fi
	
	
	################################
	##							  ##
	## TOTAL POPULATION DIVERSITY ##
	##							  ##
	################################
	
	
	#####################################################################################################################
	##  Merging the bam files in the 25.mpileup lists to make a cross-sample merged bam of each whole virus population ##
	#####################################################################################################################
	
	#Taking the bam files for merging before they are realigned around indels, and then doing this on the whole merged bam 
	
	#altering the bam lists so the non_realigned bams are selected
	if [ -n "$isseg" ]
	then 
		
		#for Vesanto
		seg=$(echo $virusstem | sed 's/^.*_S/S/g' )
		cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt | sed 's/noTIRs_InDel/noTIRs/' > Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list_norealignment.txt 
		
		#cycling through these single sample bam files and reducing the no. of reads used so that none contribute more than an average coverage of target fold to the analysis 
		
		#first calculating the no. of 151 bp reads needed for avergage targetx coverage across the genome, and rounding value up to whole number using awk
		cov="500"
		readtarget=$(awk -v cov="${cov}" -v refseqlength="$refseqlength" 'BEGIN{rt=((cov*refseqlength)/151); i=int(rt); print (rt-i<0.5)?i:i+1 }')
		#readtarget rounded up so its divisible by 2...to keep read pairs intact
		if [ $(($readtarget%2 )) -eq 0 ] 
		then 
			readtarget_rp=$(echo $readtarget)
		else 
			readtarget_rp=$(($readtarget + 1 ))
		fi
					
		echo $readtarget_rp reads needed per bam file for $cov fold $virusstem virus coverage, reference sequence is $refseqlength bp long 
		
		#cycling through the bams and counting the number of reads, asking if they are greater than the read target and then if they are, downsampling the bam to the read target 
		for b in $(cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list_norealignment.txt)
		do 
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/Vesanto.virus.remapping.$seg\///" | sed "s/.$virusstem.bwa.noTIRs.bam//")
			
			no_reads=$(samtools view -@ 6 -c $bam)
			
			if [ $no_reads -gt $readtarget ]
			then 
				echo $popstem sample of $virusstem virus has greater than $cov fold coverage, downsampling the bam file
				#calculating the proportion of reads to keep 
				prop_keep=$( awk -v readtarget="$readtarget" -v no_reads="$no_reads" 'BEGIN { pk=(readtarget/no_reads); print pk }') 
				
				echo downsampling the $popstem sample of $virusstem virus so that $prop_keep proportion of the reads remain 
				samtools view -@ 6 -b -h -s $prop_keep $bam > Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.$cov.fc.bam
				
				#and adding the name of the downsampled bam file to a new list of bam files for merging 
				echo Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.$cov.fc.bam >> Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list_norealignment.$cov.fc.txt
				
			else 
				echo $popstem sample of $virusstem virus has less than $cov fold coverage, the bam file is not downsampled 
				samtools view -@ 6 -b -h $bam > Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.$cov.fc.bam
				
				#and adding the name of the downsampled bam file to a new list of bam files for merging 
				echo Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.$cov.fc.bam >> Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list_norealignment.$cov.fc.txt
			fi
		done 
		
		#now merging the un-realigned bam files in the list
		echo merging the bams for $virusstem
		samtools merge -@ 8 -b Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list_norealignment.$cov.fc.txt Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bam
		
		#indexing the merged bam
		echo indexing the merged bam file for $virusstem virus
		samtools index -@ 6 Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bam Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bai
		
		#and now running indel realignment on the merged bam file
		echo realigning around indels in the merged bam for $virusstem virus
		
		source activate gatk3_env
		gatk3 -T RealignerTargetCreator -R /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -I Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bam -o Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc.indeltarget.list
		##Now realigning around these indels
		gatk3 -T IndelRealigner -R /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -I Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bam -targetIntervals Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc.indeltarget.list -o Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc_InDel.bam
		##this also creates a new bam.bai file
		source deactivate
		
		#creating an mpileup file from the merged bam 
		echo creating an mpileup file from the merged bam file for $virusstem but with an increased base quality filter of 40
		samtools mpileup -B -f /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc_InDel.bam -q 30 -Q 40 > Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup
		
		##Calculating the max read depth in the merged bam mpileup file to use in the allele freq script
		maxreaddepth_fc=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup | awk -v max=0 '{if($4>max) {max=$4}} END {print max}')
		echo the max read depth of $virusstem virus in the merged and downsampled sample is $maxreaddepth_fc 
		
		#calculating the read depth cut off for a 1% freq minor allele
		avgreaddepth_fc=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup | awk '{total += $4; count++ } END {ard=(total/count); i=int(ard); print (ard-i<0.5)?i:i+1 }')
		oneperc_rd_fc=$(( avgreaddepth_fc / 100 ))
		
		#calculating the median readdepth in the mpileup file, for potential downsampling
		medianreaddepth_fc=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup | cut -f4 | sort -n | awk '{ count[NR] = $1} END {if (NR % 2){print count[(NR + 1) / 2]} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2} }')
		#and variance in read depth to evaluate whether downsampling is needed 
		variancereaddepth_fc=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup | cut -f4 | awk 'BEGIN {sum=0;summ=0} {arr[NR]=$1; sum+=$1} END {average=sum/NR; for (i=1;i<=NR;i++) {summ+=(a[i]-average)^2} print summ }')
		
		#making sure that none of the mafs are less than 0 for Vesanto bams
		if (( $oneperc_rd_fc < 1 ))
		then
			oneperc_rd_fc="2"
		else
			oneperc_rd_fc=$(echo $oneperc_rd_fc )
		fi
		
		echo the average read depth of $virusstem virus in the merged and downsampled mpileup is $avgreaddepth_fc, the median read depth is $medianreaddepth_fc, variance is $variancereaddepth_fc and the average count of a minor allele for one percent frequency is $oneperc_rd_fc
		
		#downsampling the mpileup file to the median readdepth 
		echo creating an mpileup file with a depth limit of the median depth, $medianreaddepth_fc from the merged and downsampled bam file for $virusstem but with an increased base quality filter of 40
		samtools mpileup -B -f /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.noTIRs.$cov.fc_InDel.bam -q 30 -Q 40 -d $medianreaddepth_fc > Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup
		
		#recalculating stats
		##Calculating the max read depth in the merged bam mpileup file to use in the allele freq script
		maxreaddepth_fc=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup | awk -v max=0 '{if($4>max) {max=$4}} END {print max}')
		#calculating the read depth cut off for a 1% freq minor allele
		avgreaddepth_fc=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup | awk '{total += $4; count++ } END {ard=(total/count); i=int(ard); print (ard-i<0.5)?i:i+1 }')
		oneperc_rd_fc=$(( avgreaddepth_fc / 100 ))
		
		#making sure that none of the mafs are less than 0 for Vesanto bams
		if (( $oneperc_rd_fc < 1 ))
		then
			oneperc_rd_fc="2"
		else
			oneperc_rd_fc=$(echo $oneperc_rd_fc )
		fi
		
		#calculating the median readdepth in the mpileup file
		medianreaddepth_fc=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup | cut -f4 | sort -n | awk '{ count[NR] = $1} END {if (NR % 2){print count[(NR + 1) / 2]} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2} }')
		
		echo the max read depth of the $virusstem two level downsampled mpileup file is now $maxreaddepth_fc, average read depth is $avgreaddepth_fc and the average count of a minor allele for one percent frequency is $oneperc_rd_fc
		
		#Now using popoolation 2 to generate a synchronised mpileup file from the single merged mpileup file 
		echo generating a synchronised mpileup file for $virusstem from a single merged and two level downsampled bam
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/mpileup2sync.pl --input Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup --output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.sync --fastq-type sanger --min-qual 40
		
		echo identifying indellic regions in the $virusstem merged and two level downsampled mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup --output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.gtf --min-count 7
		
		#and saving these regions as a bed file
		echo converting the indellic regions identified in the $virusstem virus merged and downsampled bam mpileup files to bed format for IGV
		source activate bedops_env
		gtf2bed < Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.gtf > Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.bed
		source deactivate
		
		echo masking identified indellic regions in the $virusstem virus merged and downsampled bam synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.gtf --input Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.sync --output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync
		
		echo calculating allele frequencies and pairwise SNP differences from the $virusstem virus merged bam synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/snp-frequency-diff.pl --input Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync --output-prefix Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync --min-count 2 --min-coverage 5 --max-coverage $maxreaddepth_fc
		
		#and now with min allele count so that no SNPs are identified with <1% freq
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/snp-frequency-diff.pl --input Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync --output-prefix Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync_oneperc_maf --min-count $oneperc_rd_fc --min-coverage 5 --max-coverage $maxreaddepth_fc
		
		#############################################
		## Total population PiA / PiS calculations ##
		#############################################
		
		##Calculating the TSTV ratio from the SNPs for each virus 
		cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync_oneperc_maf_rc | awk 'NR >1 {print $0}' | cut -f 8,9 | sed 's/\t//' | sort  | uniq -c > Vesanto.virus.remapping.$seg/$virusstem.$cov.fc.snp.sums.tsv
		
		Rscript tstv_ratio.R Vesanto.virus.remapping.$seg/$virusstem.$cov.fc.snp.sums.tsv $virusstem
		tstv=$(cat Vesanto.virus.remapping.$seg/$virusstem.tstv.tsv)
		
		echo the tstv ratio for $virusstem virus is $tstv
		
		echo creating non-synonymous length table for $virusstem virus
		perl popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/create-syn-nonsynmatrix.pl --codon-table popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/codon-table.txt --output Vesanto.virus.remapping.$seg/snl.$virusstem.txt --transversion-penalty $tstv
		
		#Filtering the reference genome by the cds regions in order to count the number of each codon in the CDS 
		echo saving the fasta sequences from each cds region in $virusstem virus 
		source activate bedops_env 
		awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' $virusstem.virus.analyses/$virusid.prokka.edited2.cds.gtf | gtf2bed - > $virusstem.virus.analyses/$virusid.prokka.edited2.cds.gtf.bed
		source deactivate 
		bedtools getfasta -name -s -fi $virusid.noTIRs.fasta -bed $virusstem.virus.analyses/$virusid.prokka.edited2.cds.gtf.bed | sed 's/ /_/g' > Vesanto.virus.remapping.$seg/$virusstem.cds.regions.fasta
		
		#using an Rscript to calculate the avg. no of syn and nsyn sites in each gene in each virus
		Rscript syn_nsyn_sites.R Vesanto.virus.remapping.$seg/$virusstem.cds.regions.fasta Vesanto.virus.remapping.$seg/snl.$virusstem.txt $virusstem 
		
		##Identifying syn and non-syn SNPs in the whole population mpileup files - in order to manually calculate piA and piS 
		#calculating syn-nonsyn pi with a sliding window approach and outputting the SNPs detected
		#gtf annotation file should be in the gff3 style
		
		#using popoolation1 to filter indellic regions out of the mpileup file rather than the sync file 
		echo filtering indellic regions out of the merged whole population mpileup file for $virusstem virus
		perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup --gtf Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.gtf --output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.mpileup
		
		#Doesn't actually much matter what values you put in for pool size etc here as I'm not going to use the pi values, just the syn nonsyn statuses
		echo calculating nsyn syn pi across the whole genome with a sliding window approach for merged bam sample of $virusstem virus 
		perl popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --fastq-type sanger --gtf $virusstem.virus.analyses/$virusstem.popoolation1.gtf --pileup Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.mpileup --codon-table popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table Vesanto.virus.remapping.$seg/snl.$virusstem.txt --output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.sliding.wholegenome.syn-nsyn.pi.tsv --snp-output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.sliding.wholegenome.syn-nsyn.snps.tsv --pool-size 40 --min-qual 40 --min-count 1 --min-coverage 5 --max-coverage $maxreaddepth_fc --window-size 30 --step-size 24
		
		#because I used a sliding window - removing any duplicate snps which were in more than one window, and removing the window headers and empty lines 
		cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.sliding.wholegenome.syn-nsyn.snps.tsv | grep -v '^>' | grep "\S" | uniq > Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.wholegenome.syn-nsyn.snps.tsv
		
	else 
		
		#for Kallithea and Linvil road
		cat $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt | sed 's/noTIRs_InDel/noTIRs/' > $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list_norealignment.txt 
		
		#cycling through these single sample bam files and reducing the no. of reads used so that none contribute more than an average coverage of target fold to the analysis 
		
		#first calculating the no. of 151 bp reads needed for avergage targetx coverage across the genome, and rounding value up to whole number using awk
		cov="500"
		readtarget=$(awk -v cov="${cov}" -v refseqlength="$refseqlength" 'BEGIN{rt=((cov*refseqlength)/151); i=int(rt); print (rt-i<0.5)?i:i+1 }')
		
		if [ $(($readtarget%2 )) -eq 0 ] 
		then 
			readtarget_rp=$(echo $readtarget)
		else 
			readtarget_rp=$(($readtarget + 1 ))
		fi
					
		echo $readtarget_rp reads needed per bam file for $cov fold $virusstem virus coverage, reference sequence is $refseqlength bp long 
		
		#cycling through the bams and counting the number of reads, asking if they are greater than the read target and then if they are, downsampling the bam to the read target 
		for b in $(cat $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list_norealignment.txt)
		do 
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/$virusstem.virus.analyses\///" | sed "s/.$virusstem.bwa.noTIRs.bam//")
			
			no_reads=$(samtools view -@ 6 -c $bam)
			
			if [ $no_reads -gt $readtarget ]
			then 
				echo $popstem sample of $virusstem virus has greater than $cov fold coverage, downsampling the bam file
				#calculating the proportion of reads to keep 
				prop_keep=$( awk -v readtarget="$readtarget" -v no_reads="$no_reads" 'BEGIN { pk=(readtarget/no_reads); print pk }') 
				
				echo downsampling the $popstem sample of $virusstem virus so that $prop_keep proportion of the reads remain 
				samtools view -@ 6 -b -h -s $prop_keep $bam > $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.$cov.fc.bam
				
				#and adding the name of the downsampled bam file to a new list of bam files for merging 
				echo $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.$cov.fc.bam >> $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list_norealignment.$cov.fc.txt
				
			else 
				echo $popstem sample of $virusstem virus has less than $cov fold coverage, the bam file is not downsampled 
				samtools view -@ 6 -b -h $bam > $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.$cov.fc.bam
				
				#and adding the name of the downsampled bam file to a new list of bam files for merging 
				echo $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.$cov.fc.bam >> $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list_norealignment.$cov.fc.txt
			fi
		done 
		
		#now merging the un-realigned bam files in the list
		echo merging the bams for $virusstem
		samtools merge -@ 8 -b $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list_norealignment.$cov.fc.txt $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bam
		
		#indexing the merged bam 
		echo indexing the merged bam file for $virusstem virus
		samtools index -@ 6 $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bam $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bai
		
		#and now running indel realignment on the merged bam file
		echo realigning around indels in the merged bam for $virusstem virus 
		
		source activate gatk3_env
		gatk3 -T RealignerTargetCreator -R /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -I $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bam -o $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc.indeltarget.list
		##Now realigning around these indels
		gatk3 -T IndelRealigner -R /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -I $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc.bam -targetIntervals $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc.indeltarget.list -o $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc_InDel.bam
		##this also creates a new bam.bai file
		source deactivate
		
		#creating an mpileup file from the merged bam 
		echo creating an mpileup file from the merged bam file for $virusstem but with an increased base quality filter of 40
		samtools mpileup -B -f /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc_InDel.bam -q 30 -Q 40 > $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup
		
		##Calculating the max read depth in the merged bam mpileup file to use in the allele freq script
		maxreaddepth_fc=$(cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup | awk -v max=0 '{if($4>max) {max=$4}} END {print max}')
		echo the max read depth of $virusstem virus in the merged and downsampled sample is $maxreaddepth_fc 
		
		#calculating the read depth cut off for a 1% freq minor allele
		avgreaddepth_fc=$(cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup | awk '{total += $4; count++ } END {ard=(total/count); i=int(ard); print (ard-i<0.5)?i:i+1 }')
		oneperc_rd_fc=$(( avgreaddepth_fc / 100 ))
		
		#calculating the median readdepth in the mpileup file, for potential downsampling
		medianreaddepth_fc=$(cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup | cut -f4 | sort -n | awk '{ count[NR] = $1} END {if (NR % 2){print count[(NR + 1) / 2]} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2} }')
		and variance in read depth to evaluate whether downsampling is needed 
		variancereaddepth_fc=$(cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup | cut -f4 | awk 'BEGIN {sum=0;summ=0} {arr[NR]=$1; sum+=$1} END {average=sum/NR; for (i=1;i<=NR;i++) {summ+=(a[i]-average)^2} print summ }')
		
		echo the average read depth of $virusstem virus in the merged and downsampled mpileup is $avgreaddepth_fc, the median read depth is $medianreaddepth_fc, variance is $variancereaddepth_fc and the average count of a minor allele for one percent frequency is $oneperc_rd_fc
		
		#downsampling the mpileup file to the median readdepth 
		echo creating an mpileup file with a depth limit of the median depth, $medianreaddepth_fc from the merged and downsampled bam file for $virusstem but with an increased base quality filter of 40
		samtools mpileup -B -f /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta $virusstem.virus.analyses/$virusstem.25.merged.bwa.noTIRs.$cov.fc_InDel.bam -q 30 -Q 40 -d $medianreaddepth_fc > $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup
		
		#recalculating stats
		##Calculating the max read depth in the merged bam mpileup file to use in the allele freq script
		maxreaddepth_fc=$(cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup | awk -v max=0 '{if($4>max) {max=$4}} END {print max}')
		#calculating the read depth cut off for a 1% freq minor allele
		avgreaddepth_fc=$(cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup | awk '{total += $4; count++ } END {ard=(total/count); i=int(ard); print (ard-i<0.5)?i:i+1 }')
		oneperc_rd_fc=$(( avgreaddepth_fc / 100 ))
		#calculating the median readdepth in the mpileup file
		medianreaddepth_fc=$(cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup | cut -f4 | sort -n | awk '{ count[NR] = $1} END {if (NR % 2){print count[(NR + 1) / 2]} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2} }')
		
		#echo the max read depth of the $virusstem downsampled mpileup file is now $maxreaddepth, average read depth is $avgreaddepth and the average count of a minor allele for one percent frequency is $oneperc_rd
		echo the max read depth of the $virusstem two level downsampled mpileup file is now $maxreaddepth_fc, average read depth is $avgreaddepth_fc and the average count of a minor allele for one percent frequency is $oneperc_rd_fc
		
		#Now using popoolation 2 to generate a synchronised mpileup file from the single merged mpileup file 
		echo generating a synchronised mpileup file for $virusstem from a single merged and two level downsampled bam
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/mpileup2sync.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.sync --fastq-type sanger --min-qual 40
		
		echo identifying indellic regions in the $virusstem merged and two level downsampled mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.gtf --min-count 7
		
		#and saving these regions as a bed file
		echo converting the indellic regions identified in the $virusstem virus merged and downsampled bam mpileup files to bed format for IGV
		source activate bedops_env
		gtf2bed < $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.gtf > $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.bed
		source deactivate
		
		echo masking identified indellic regions in the $virusstem virus merged and downsampled bam synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.gtf --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.sync --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync
		
		echo organising the merged bam synchronised mpileup files for $virusstem virus by genes
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/create-genewise-sync.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync --gtf $virusstem.virus.analyses/$virusid.edited.cds.gtf --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.genewise.sync
		
		echo creating a merged bam sync file with only non-CDS regions in it for $virusstem virus 
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf $virusstem.virus.analyses/$virusid.edited.cds.gtf --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.intergenic.sync
		
		echo calculating allele frequencies and pairwise SNP differences from the $virusstem virus merged bam synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/snp-frequency-diff.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync --output-prefix $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync --min-count 2 --min-coverage 5 --max-coverage $maxreaddepth_fc
		
		echo calculating allele frequencies and pairwise SNP differences from the $virusstem virus merged bam genewise synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/snp-frequency-diff.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.genewise.sync --output-prefix $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.genewise.sync --min-count 2 --min-coverage 5 --max-coverage $maxreaddepth_fc
		
		echo calculating allele frequencies and pairwise SNP differences from the $virusstem virus merged bam intergenic synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/snp-frequency-diff.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.intergenic.sync --output-prefix $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.intergenic.sync --min-count 2 --min-coverage 5 --max-coverage $maxreaddepth_fc
		
		#and now with min allele count so that no SNPs are identified with <1% freq
		echo calculating allele frequencies and pairwise SNP differences from the $virusstem virus merged bam synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/snp-frequency-diff.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync --output-prefix $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync_oneperc_maf --min-count $oneperc_rd_fc --min-coverage 5 --max-coverage $maxreaddepth_fc
		
		echo calculating allele frequencies and pairwise SNP differences from the $virusstem virus merged bam genewise synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/snp-frequency-diff.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.genewise.sync --output-prefix $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.genewise.sync_oneperc_maf --min-count $oneperc_rd_fc --min-coverage 5 --max-coverage $maxreaddepth_fc
		
		echo calculating allele frequencies and pairwise SNP differences from the $virusstem virus merged bam intergenic synchronised mpileup file
		perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/snp-frequency-diff.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.intergenic.sync --output-prefix $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.intergenic.sync_oneperc_maf --min-count $oneperc_rd_fc --min-coverage 5 --max-coverage $maxreaddepth_fc
		
		#############################################
		## Whole population PiA / PiS calculations ##
		#############################################
		
		##Calculating the TSTV ratio from the SNPs for each virus 
		cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indelfiltered.sync_oneperc_maf_rc | awk 'NR >1 {print $0}' | cut -f 8,9 | sed 's/\t//' | sort  | uniq -c > $virusstem.virus.analyses/$virusstem.$cov.fc.snp.sums.tsv
		
		Rscript tstv_ratio.R $virusstem.virus.analyses/$virusstem.$cov.fc.snp.sums.tsv $virusstem
		tstv=$(cat $virusstem.virus.analyses/$virusstem.tstv.tsv)
		
		echo the tstv ratio for $virusstem virus is $tstv
		
		echo creating non-synonymous length table for $virusstem virus
		perl popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/create-syn-nonsynmatrix.pl --codon-table popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/codon-table.txt --output $virusstem.virus.analyses/snl.$virusstem.txt --transversion-penalty $tstv
		
		#Filtering the reference genome by the cds regions in order to count the number of each codon in the CDS 
		echo saving the fasta sequences from each cds region in $virusstem virus 
		source activate bedops_env 
		awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' $virusstem.virus.analyses/$virusid.edited.cds.gtf | gtf2bed - > $virusstem.virus.analyses/$virusid.edited.cds.gtf.bed
		source deactivate 
		bedtools getfasta -name -s -fi $virusid.noTIRs.fasta -bed $virusstem.virus.analyses/$virusid.edited.cds.gtf.bed > $virusstem.virus.analyses/$virusstem.cds.regions.fasta
		
		#using an Rscript to calculate the avg. no of syn and nsyn sites in each gene in each virus
		Rscript syn_nsyn_sites.R $virusstem.virus.analyses/$virusstem.cds.regions.fasta $virusstem.virus.analyses/snl.$virusstem.txt $virusstem 
		
		##Identifying syn and non-syn SNPs in the whole population mpileup files - in order to manually calculate piA and piS 
		#calculating syn-nonsyn pi with a sliding window approach and outputting the SNPs detected
		#gtf annotation file should be in the gff3 style
		
		#using popoolation1 to filter indellic regions out of the mpileup file rather than the sync file 
		echo filtering indellic regions out of the merged whole population mpileup file for $virusstem virus
		perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup --gtf $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.mpileup.indel.regions.gtf --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.mpileup
		
		#Doesn't actually much matter what values you put in for pool size etc here as I'm not going to use the pi values, just the syn nonsyn statuses - I've put a slightly overlapping sliding window to try to speed up the process of outputting the snps
		echo calculating nsyn syn pi across the whole genome with a sliding window approach for merged bam sample of $virusstem virus 
		perl popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --fastq-type sanger --gtf $virusstem.virus.analyses/$virusstem.popoolation1.gtf --pileup $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.mpileup --codon-table popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table $virusstem.virus.analyses/snl.$virusstem.txt --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.sliding.wholegenome.syn-nsyn.pi.tsv --snp-output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.sliding.wholegenome.syn-nsyn.snps.tsv --pool-size 40 --min-qual 40 --min-count 1 --min-coverage 5 --max-coverage $maxreaddepth_fc --window-size 15 --step-size 12 --dissable-corrections
		
		#because I used a sliding window - removing any duplicate snps which were in more than one window, and removing the window headers and empty lines 
		cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.sliding.wholegenome.syn-nsyn.snps.tsv | grep -v '^>' | grep "\S" | uniq > $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.wholegenome.syn-nsyn.snps.tsv
	fi
done

#Running the per_site_nucl_div.R rscript which calculates site wise nucleotide diversity for all viruses/segments with CDS regions labelled
#With the output of this script, you can calculate piA and piS for the whole genome, and per gene, using the no. of syn and nsyn sites in each gene 
Rscript per_site_nucl_div.R

for v in $(cat gt5samples.virus.list.txt)
do
	cov="500"
	
	##Setting up the id and file stem for each virus
	virusid=$(echo $v)
	
	isseg=$(echo $v | grep '^Vesanto')
	if [ -n "$isseg" ]
	then
		virusstem=$(echo $v | sed 's/\(_virus\)//')
	else
		virusstem=$(echo $v | sed -r 's/^[A-Z]{2}[0-9]{6}_//' | sed 's/virus.*/virus/' | sed 's/\(_virus\|_Nudivirus\)//')
	fi
	
	refseqlength=$(grep -A1 "$virusid" linear.virus.noTIRs.fasta | awk '/^>/{getline; getline; print}' | awk '{print length($1)}')
	
	#First for the Vesanto bams
	if [ -n "$isseg" ]
	then
		seg=$(echo $virusstem | sed 's/^.*_S/S/g' )
	
		#calculating total population piA and piS
		Rscript piA_piS.R Vesanto.virus.remapping.$seg/$virusstem.genewise.syn-nsyn.sites.tsv Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv $virusid
		
		#Combining the genewise Vesanto piA piS files 
		cat Vesanto.virus.remapping.$seg/$virusid.genewise.wholepop.piA_piS.tsv | awk 'NR >1 {print $0}' | grep -v 'whole_genome' >> Vesanto.segs.genewise.piA.piS.results.tsv 
		
		type=cds
		
		no_cds_nsyn_snps=$(cat Vesanto.virus.remapping.$seg/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 5 )
		no_cds_syn_snps=$(cat Vesanto.virus.remapping.$seg/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 6 )
		no_nsyn_sites=$(cat Vesanto.virus.remapping.$seg/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 3 )
		no_syn_sites=$(cat Vesanto.virus.remapping.$seg/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 4 )
		
		if [ $no_cds_nsyn_snps == 0 ] && [ $no_cds_syn_snps == 0 ]
		then 
			echo the whole population merged sample of $virusstem virus contains no snps at more than 0.01 allele frequency in cds
			echo -e "$virusstem \t $type" >> wholepop_no_snp_samples.txt

			echo -e "$virusid \t $type \t $no_nsyn_sites \t $no_syn_sites \t $no_cds_nsyn_snps \t $no_cds_syn_snps \t 0 \t 0 \t 0 " >> all_viruses_whole_genome_piA_piS.tsv
			
		else 
			echo the whole population merged sample of $virusstem virus contains some snps at more than 0.01 allele frequency in cds
			
			#and echo the whole genome data to a combined table of whole genome piA and piS
			piA=$(cat Vesanto.virus.remapping.$seg/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 7 )
			piS=$(cat Vesanto.virus.remapping.$seg/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 8 )
			piA_over_piS=$(cat Vesanto.virus.remapping.$seg/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 9 )
			
			echo -e "$virusid \t $type \t $no_nsyn_sites \t $no_syn_sites \t $no_cds_nsyn_snps \t $no_cds_syn_snps \t $piA \t $piS \t $piA_over_piS " >> all_viruses_whole_genome_piA_piS.tsv
			
		fi
		
		#Making an mpileup file with just the intergenic regions by masking the cds regions
		echo making an mpileup file with just the intergenic regions for $virusstem
		#adding an empty transcript id field to the gtf file to make sure it works
		awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' $virusstem.virus.analyses/$virusid.prokka.edited2.cds.gtf > $virusstem.virus.analyses/$virusid.prokka.edited3.cds.gtf
		perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.mpileup --gtf $virusstem.virus.analyses/$virusid.prokka.edited3.cds.gtf --output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.mpileup
		
		#and generating the intergenic snps using a sliding window for the intergenic regions 
		echo generating a list of the snps in intergenic regions of $virusstem 
		perl popoolation_1.2.2/popoolation_1.2.2/Variance-sliding.pl --measure pi --input Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.mpileup --output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.sliding.pi.tsv --snp-output Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.sliding.pi.snps.tsv --fastq-type sanger --pool-size 40 --min-count 1 --min-coverage 5 --max-coverage $maxreaddepth_fc --window-size 15 --step-size 12 --min-qual 40 --dissable-corrections
			
		#removing duplicate rows from snp file as sliding window was used
		cat Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.sliding.pi.snps.tsv | grep -v '^>' | grep "\S" | uniq > Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.intergenic.snps.tsv
		
		#Running this script to generate nucleotide diversity tables for each of the wholepop intergenic snp files, and then calculate intergenic pi from them
		type=intergenic
			
		if [ -s Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.intergenic.snps.tsv ]
		then
			echo the whole population merged sample $virusstem virus contains some snps in intergenic regions			
			Rscript wholepop_intergenic_pi.R Vesanto.virus.remapping.$seg/$virusstem.25.merged.bwa.$cov.fc.intergenic.snps.tsv $virusid $refseqlength $virusstem.virus.analyses/$virusid.prokka.edited2.cds.gtf
			
			#if no snps are found at an allele frequency greater than 0.01, then also add the sample to the no_snp list 
			no_intergenic_snps=$(cat Vesanto.virus.remapping.$seg/$virusstem.wholepop.intergenic.pi.tsv | grep $virusid | cut -f 3 )
			
			if [ $no_intergenic_snps == 0 ]
			then 
				echo the whole population merged sample of $virusstem virus contains no snps at more than 0.01 allele frequency in intergenic regions
				echo -e "$virusstem \t $type" >> wholepop_no_snp_samples.txt
				
				#and echo the data to a combined table of intergenic pi
				no_intergenic_sites=$(cat Vesanto.virus.remapping.$seg/$virusstem.wholepop.intergenic.pi.tsv | grep $virusid | cut -f 2 )
				echo -e "$virusid \t $type \t $no_intergenic_sites \t $no_intergenic_snps \t 0" >> all_viruses_intergenic_pi.tsv
				
			else 
				echo the whole population merged sample of $virusstem virus contains snps at at least 0.01 allele frequency in intergenic regions 
				
				#and echo the data to a combined table of intergenic pi
				no_intergenic_sites=$(cat Vesanto.virus.remapping.$seg/$virusstem.wholepop.intergenic.pi.tsv | grep $virusid | cut -f 2 )
				intergenic_pi=$(cat Vesanto.virus.remapping.$seg/$virusstem.wholepop.intergenic.pi.tsv | grep $virusid | cut -f 4 )
				echo -e "$virusid \t $type \t $no_intergenic_sites \t $no_intergenic_snps \t $intergenic_pi" >> all_viruses_intergenic_pi.tsv
			fi
		else
			echo the whole population merged sample of $virusstem virus contains no snps in intergenic regions	
			echo -e "$virusstem \t $type" >> wholepop_no_snp_samples.txt
			
			Rscript intergenic_sites.R $virusid $refseqlength $virusstem.virus.analyses/$virusid.prokka.edited2.cds.gtf
			no_intergenic_sites=$(cat Vesanto.virus.remapping.$seg/$virusstem.intergenic.sites.tsv | grep -v 'x' )
			
			echo -e "$virusid \t $type \t $no_intergenic_sites \t 0 \t 0" >> all_viruses_intergenic_pi.tsv
		fi
		
	else 
		#For Kallithea and Linvill Road viruses
		
		#calculating total population piA and piS
		Rscript piA_piS.R $virusstem.virus.analyses/$virusstem.genewise.syn-nsyn.sites.tsv $virusstem.virus.analyses/$virusstem.25.merged.bwa.500.fc.wholegenome.syn-nsyn.snps.nulc.div.tsv $virusid
		
		#Adding the virus whole genome sample to the list of whole genome samples without snps if there are no syn and nsyn snps in the whole genome pi table at greater than 0.01 allele freq 
		
		type=cds
		
		no_cds_nsyn_snps=$(cat $virusstem.virus.analyses/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 5 )
		no_cds_syn_snps=$(cat $virusstem.virus.analyses/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 6 )
		no_nsyn_sites=$(cat $virusstem.virus.analyses/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 3 )
		no_syn_sites=$(cat $virusstem.virus.analyses/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 4 )
		
		if [ $no_cds_nsyn_snps == 0 ] && [ $no_cds_syn_snps == 0 ]
		then 
			echo the whole population merged sample of $virusstem virus contains no snps at more than 0.01 allele frequency in cds
			echo -e "$virusstem \t $type" >> wholepop_no_snp_samples.txt
			
			echo -e "$virusid \t $type \t $no_nsyn_sites \t $no_syn_sites \t $no_cds_nsyn_snps \t $no_cds_syn_snps \t 0 \t 0 \t 0 " >> all_viruses_whole_genome_piA_piS.tsv
			
		else 
			echo the whole population merged sample of $virusstem virus contains some snps at more than 0.01 allele frequency in cds
			
			#and echo the whole genome data to a combined table of whole genome piA and piS
			piA=$(cat $virusstem.virus.analyses/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 7 )
			piS=$(cat $virusstem.virus.analyses/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 8 )
			piA_over_piS=$(cat $virusstem.virus.analyses/$virusid.genewise.wholepop.piA_piS.tsv | grep 'whole_genome' | cut -f 9 )
			
			echo -e "$virusid \t $type \t $no_nsyn_sites \t $no_syn_sites \t $no_cds_nsyn_snps \t $no_cds_syn_snps \t $piA \t $piS \t $piA_over_piS " >> all_viruses_whole_genome_piA_piS.tsv
		fi
		
		#Making an mpileup file with just the intergenic regions by masking the cds regions
		echo making an mpileup file with just the intergenic regions for $virusstem
		perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.mpileup --gtf $virusstem.virus.analyses/$virusid.edited.cds.gtf --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.mpileup
		
		#and generating the intergenic snps using a sliding window for the intergenic regions 
		echo generating a list of the snps in intergenic regions of $virusstem 
		perl popoolation_1.2.2/popoolation_1.2.2/Variance-sliding.pl --measure pi --input $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.mpileup --output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.sliding.pi.tsv --snp-output $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.sliding.pi.snps.tsv --fastq-type sanger --pool-size 40 --min-count 1 --min-coverage 5 --max-coverage $maxreaddepth_fc --window-size 15 --step-size 12 --min-qual 40 --dissable-corrections
			
		#removing duplicate rows from snp file as sliding window was used
		cat $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.indelfiltered.intergenic.sliding.pi.snps.tsv | grep -v '^>' | grep "\S" | uniq > $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.intergenic.snps.tsv
		
		#Running this script to generate nucleotide diversity tables for each of the wholepop intergenic snp files, and then calculate intergenic pi from them
		type=intergenic
			
		if [ -s $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.intergenic.snps.tsv ]
		then
			echo the whole population merged sample $virusstem virus contains some snps in intergenic regions			
			Rscript wholepop_intergenic_pi.R $virusstem.virus.analyses/$virusstem.25.merged.bwa.$cov.fc.intergenic.snps.tsv $virusid $refseqlength $virusstem.virus.analyses/$virusid.edited.cds.gtf
			
			#if no snps are found at an allele frequency greater than 0.01, then also add the sample to the no_snp list 
			no_intergenic_snps=$(cat $virusstem.virus.analyses/$virusstem.wholepop.intergenic.pi.tsv | grep $virusid | cut -f 3 )
			
			if [ $no_intergenic_snps == 0 ]
			then 
				echo the whole population merged sample of $virusstem virus contains no snps at more than 0.01 allele frequency in intergenic regions 
				echo -e "$virusstem \t $type" >> wholepop_no_snp_samples.txt
				
				no_intergenic_sites=$(cat $virusstem.virus.analyses/$virusstem.wholepop.intergenic.pi.tsv | grep $virusid | cut -f 2 )
	
				echo -e "$virusid \t $type \t $no_intergenic_sites \t $no_intergenic_snps \t 0" >> all_viruses_intergenic_pi.tsv
				
			else 
				echo the whole population merged sample of $virusstem virus contains snps at at least 0.01 allele frequency in intergenic regions
				
				#and echo the data to a combined table of intergenic pi
				no_intergenic_sites=$(cat $virusstem.virus.analyses/$virusstem.wholepop.intergenic.pi.tsv | grep $virusid | cut -f 2 )
				intergenic_pi=$(cat $virusstem.virus.analyses/$virusstem.wholepop.intergenic.pi.tsv | grep $virusid | cut -f 4 )
				echo -e "$virusid \t $type \t $no_intergenic_sites \t $no_intergenic_snps \t $intergenic_pi" >> all_viruses_intergenic_pi.tsv
			fi
		else
			echo the whole population merged sample of $virusstem virus contains no snps in intergenic regions	
			echo -e "$virusstem \t $type" >> wholepop_no_snp_samples.txt
			
			Rscript intergenic_sites.R $virusid $refseqlength $virusstem.virus.analyses/$virusid.edited.cds.gtf
			no_intergenic_sites=$(cat $virusstem.virus.analyses/$virusstem.intergenic.sites.tsv | grep -v 'x' )
			
			echo -e "$virusid \t $type \t $no_intergenic_sites \t 0 \t 0" >> all_viruses_intergenic_pi.tsv
		fi	
	fi
	
	##############################################
	## Distribution of indels around the genome ##
	##############################################
	
	#Generating single sample mpileup files for each of the bams in the 25 fold coverage bam lists 
	
	if [ -n "$isseg" ]
	then
		seg=$(echo $virusstem | sed 's/^.*_S/S/g' )
		echo generating mpileups and identifying indellic regions in the $virusstem mpileup files
		for b in $(cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.10.bam.list.txt)
		do 
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/Vesanto.virus.remapping.$seg\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			#creating a pileup formatted file for each bam 
			echo creating a pileup file for the $popstem bam file mapped to $virusstem virus
			samtools mpileup -B -f /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -q 30 -Q 40 -d 500 $bam > Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup
			
			#identifying indels in the mpileup file, only counting if at least 5 reads support the indel and with window set to 1 so no window around it
			echo indentifying indels in the $popstem bam mapped to $virusstem virus
			perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup --output Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup.just.indels.gtf --min-count 5 --indel-window 1
			
			#making a bed file to look at in igv
			source activate bedops_env
			gtf2bed < Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup.just.indels.gtf > Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup.just.indels.bed
			source deactivate
			
			#Then converting the gtf file to a file with the virus, population, position and whether an indel (supported by at least 5 reads) is there as binary
			echo creating txt file of indel positions for $popstem sample of $virusstem virus
			#this creates a file w a 0 at each position in the reference sequence 
			paste <(awk -v virusstem="$virusstem" -v refseqlength="$refseqlength" 'BEGIN{ORS="\n"; for(c=0;c<refseqlength;c++) print virusstem}') <(awk -v popstem="$popstem" -v refseqlength="$refseqlength" 'BEGIN{ORS="\n"; for(c=0;c<refseqlength;c++) print popstem}') <(seq 1 $refseqlength) <(awk -v refseqlength="$refseqlength" 'BEGIN{ORS="\n"; for(c=0;c<refseqlength;c++) print '0'}') > Vesanto.virus.remapping.$seg/$popstem.$virusstem.indel.positions.txt
			#And this creates a binary file w 1 or 0 at each position to indicate whether an indel is supported there
			while read -r virus source feature start end x y z gene_id transcript_id
			do
				cat Vesanto.virus.remapping.$seg/$popstem.$virusstem.indel.positions.txt | awk -v start="$start" -v end="$end" -F'\t' 'BEGIN{OFS="\t";} {if ($3 >= start && $3 <= end) print $1,$2,$3,"1"; else print $0}' > Vesanto.virus.remapping.$seg/$popstem.$virusstem.indel.positions.tmp && mv Vesanto.virus.remapping.$seg/$popstem.$virusstem.indel.positions.tmp Vesanto.virus.remapping.$seg/$popstem.$virusstem.indel.positions.txt
			done < Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup.just.indels.gtf
		
			echoing the name of each indel position file into a list - this time for the 10 cov bams
			echo Vesanto.virus.remapping.$seg/$popstem.$virusstem.indel.positions.txt >> Vesanto.virus.remapping.$seg/$virusstem.10.indel.positions.list.txt
		done
		
		#making another list for the 25 cov bams 
		for b in $(cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt)
		do
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/Vesanto.virus.remapping.$seg\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			echo Vesanto.virus.remapping.$seg/$popstem.$virusstem.indel.positions.txt >> Vesanto.virus.remapping.$seg/$virusstem.25.indel.positions.list.txt
		done
		
		#Generating a file for each virus/variant with the indel regions for all populations
		#make array of the indel position text files for the virus
		
		indel_positions_25=($(cat Vesanto.virus.remapping.$seg/$virusstem.25.indel.positions.list.txt))
		
		#combining the indel position files for each virus into one dataframe
		echo combining the indel position text files for $virusstem virus
		cat $(echo ${indel_positions_25[@]}) > Vesanto.virus.remapping.$seg/$virusstem.25.combined.indel.positions.txt
		
		#examine how many positions of the reference genome have an indel supported by at least one sample
		no_indels=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.combined.indel.positions.txt | awk '{ if($4==1) print $0 }' | wc -l )
		no_indel_positions=$(cat Vesanto.virus.remapping.$seg/$virusstem.25.combined.indel.positions.txt | awk '{ if($4==1) print $0 }' | sort -u -k3,3 | wc -l )
		no_samples=$(cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt | wc -l )
		
		echo $virusstem virus has $no_indels indels supported across all samples, on $no_indel_positions positions of the virus genome, which is $refseqlength bp long 
		
		echo -e "$virusstem \t $refseqlength \t $no_samples \t $no_indels \t $no_indel_positions" >> all_viruses_indel_data.tsv 
		
	else
	
		echo generating mpileups and identifying indellic regions in the $virusstem mpileup files
		for b in $(cat $virusstem.virus.analyses/$virusstem.mpileup.10.bam.list.txt)
		do 
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/$virusstem.virus.analyses\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			#creating a pileup formatted file for each bam 
			echo creating a pileup file for the $popstem bam file mapped to $virusstem virus
			samtools mpileup -B -f /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -q 30 -Q 40 -d 500 $bam > $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup
			
			#identifying indels in the mpileup file, only counting if at least 5 reads support the indel and with window set to 1 so no window around it
			echo indentifying indels in the $popstem bam mapped to $virusstem virus
			perl /mnt/drive3-6tb/Kallithea_Diversity/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup --output $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup.just.indels.gtf --min-count 5 --indel-window 1
			
			#making a bed file to look at in igv
			source activate bedops_env
			gtf2bed < $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup.just.indels.gtf > $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup.just.indels.bed
			source deactivate
			
			#Then converting the gtf file to a file with the virus, population, position and whether an indel (supported by at least 5 reads) is there as binary
			echo creating txt file of indel positions for $popstem sample of $virusstem virus
			#this creates a file w a 0 at each position in the reference sequence
			paste <(awk -v virusstem="$virusstem" -v refseqlength="$refseqlength" 'BEGIN{ORS="\n"; for(c=0;c<refseqlength;c++) print virusstem}') <(awk -v popstem="$popstem" -v refseqlength="$refseqlength" 'BEGIN{ORS="\n"; for(c=0;c<refseqlength;c++) print popstem}') <(seq 1 $refseqlength) <(awk -v refseqlength="$refseqlength" 'BEGIN{ORS="\n"; for(c=0;c<refseqlength;c++) print '0'}') > $virusstem.virus.analyses/$popstem.$virusstem.indel.positions.txt
			
			#And this creates a binary file w 1 or 0 at each position to indicate whether an indel is supported there
			while read -r virus source feature start end x y z gene_id transcript_id
			do
				cat $virusstem.virus.analyses/$popstem.$virusstem.indel.positions.txt | awk -v start="$start" -v end="$end" -F'\t' 'BEGIN{OFS="\t";} {if ($3 >= start && $3 <= end) print $1,$2,$3,"1"; else print $0}' > $virusstem.virus.analyses/$popstem.$virusstem.indel.positions.tmp && mv $virusstem.virus.analyses/$popstem.$virusstem.indel.positions.tmp $virusstem.virus.analyses/$popstem.$virusstem.indel.positions.txt
			done < $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup.just.indels.gtf
		
			#echoing the name of each indel position file into a list - this time for the 10 cov bams
			echo $virusstem.virus.analyses/$popstem.$virusstem.indel.positions.txt >> $virusstem.virus.analyses/$virusstem.10.indel.positions.list.txt
		
		done
		
		#making another list for the 25 cov bams 
		for b in $(cat $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt)
		do
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/$virusstem.virus.analyses\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			echo $virusstem.virus.analyses/$popstem.$virusstem.indel.positions.txt >> $virusstem.virus.analyses/$virusstem.25.indel.positions.list.txt
		done
		
		#Generating a file for each virus/variant with the indel regions for all populations
		make array of the indel position text files for the virus
		
		indel_positions_25=($(cat $virusstem.virus.analyses/$virusstem.25.indel.positions.list.txt))
		
		#combining the indel position files for each virus into one dataframe
		echo combining the indel position text files for $virusstem virus
		cat $(echo ${indel_positions_25[@]}) > $virusstem.virus.analyses/$virusstem.25.combined.indel.positions.txt 
		
		#examine how many positions of the reference genome have an indel supported by at least one sample
		no_indels=$(cat $virusstem.virus.analyses/$virusstem.25.combined.indel.positions.txt | awk '{ if($4==1) print $0 }' | wc -l )
		no_indel_positions=$(cat $virusstem.virus.analyses/$virusstem.25.combined.indel.positions.txt | awk '{ if($4==1) print $0 }' | sort -u -k3,3 | wc -l )
		no_samples=$(cat $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt | wc -l )
		
		echo $virusstem virus has $no_indels indels supported across all samples, on $no_indel_positions positions of the virus genome, which is $refseqlength bp long 
		
		echo -e "$virusstem \t $refseqlength \t $no_samples \t $no_indels \t $no_indel_positions" >> all_viruses_indel_data.tsv 
		
	fi
		
	################################################################
	## Local (per sample) piA, piS and intergenic pi calculations ##
	################################################################
	
	if [ -n "$isseg" ]
	then
		
		#for Vesanto segments 
		seg=$(echo $virusstem | sed 's/^.*_S/S/g' )
		
		for b in $(cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt)
		do
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/Vesanto.virus.remapping.$seg\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			
			echo identifying indels in the pileup file for the $popstem mpileup file mapped to $virusstem virus
			perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup --output Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup.indel.regions.gtf --min-count 7
			
			echo masking indels in the pileup file for the $popstem mpileup file mapped to $virusstem virus
			#mask these indellic regions
			perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup --gtf Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.mpileup.indel.regions.gtf --output Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.mpileup
			
			#Making an mpileup file with just the intergenic regions by masking the cds regions
			echo making an mpileup file with just the intergenic regions for $popstem sample of $virusstem
			#adding an empty transcript id field to the gtf file to make sure it works
			awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id "";"; }' $virusstem.virus.analyses/$virusid.prokka.edited2.cds.gtf > $virusstem.virus.analyses/$virusid.prokka.edited3.cds.gtf
			perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.mpileup --gtf $virusstem.virus.analyses/$virusid.prokka.edited3.cds.gtf --output Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.intergenic.mpileup
			
			#calculating the nsyn and syn sites in each virus/segment using the sliding window syn-nsyn pi script
			echo calculating nsyn syn pi with a sliding window approach for $popstem sample of $virusstem virus 
			perl popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --fastq-type sanger --gtf $virusstem.virus.analyses/$virusstem.popoolation1.gtf --pileup Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.mpileup --codon-table popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table Vesanto.virus.remapping.$seg/snl.$virusstem.txt  --output Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.sliding.syn-nsyn.pi.tsv --snp-output Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.sliding.syn-nsyn.snps.tsv --pool-size 40 --min-qual 40 --min-count 1 --min-coverage 5 --max-coverage 550 --window-size 30 --step-size 24 --dissable-corrections
			
			#removing duplicate rows from snp file as sliding window was used
			cat Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.sliding.syn-nsyn.snps.tsv | grep -v '^>' | grep "\S" | uniq > Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.syn-nsyn.snps.tsv
			
			#running this Rscript to calculate per site nucleotide diversity for each snp in the sample
			type=cds
			
			if [ -s Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.syn-nsyn.snps.tsv ]
			then
				echo the $popstem sample of $virusstem virus contains some snps in coding regions			
				Rscript sample_per_site_nucl_div.R Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.syn-nsyn.snps.tsv $type $virusid $popstem
			else
				echo the $popstem sample of $virusstem virus contains no snps in coding regions	
				echo -e "$virusstem \t $popstem \t $type" >> no_snp_samples.txt
			fi
			
			#output from this will be the snps and thier per site nucleotide diversity
			
			#Calculating pi across the intergenic regions - to get the intergenic snps
			echo calculating pi across these intergenic regions for the $popstem sample of $virusstem virus 
			perl popoolation_1.2.2/popoolation_1.2.2/Variance-sliding.pl --measure pi --input Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.intergenic.mpileup --output Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.intergenic.sliding.pi.tsv --snp-output Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.intergenic.sliding.pi.snps.tsv --fastq-type sanger --pool-size 40 --min-count 1 --min-coverage 5 --max-coverage 550 --window-size 30 --step-size 24 --min-qual 40 --dissable-corrections
			
			#removing duplicate rows from snp file as sliding window was used
			cat Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.intergenic.sliding.pi.snps.tsv | grep -v '^>' | grep "\S" | uniq > Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.intergenic.snps.tsv
			
			#running this Rscript to calculate per site nucleotide diversity for each snp in the sample
			type=intergenic
			
			if [ -s Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.intergenic.snps.tsv ]
			then
				echo the $popstem sample of $virusstem virus contains some snps in coding regions			
				Rscript sample_per_site_nucl_div.R Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.indelfiltered.intergenic.snps.tsv $type $virusid $popstem
			else
				echo the $popstem sample of $virusstem virus contains no snps in intergenic regions	
			fi
			#output from this will be the snps and thier per site nucleotide diversity
			
		done
		
		#Now running this R script to cycle through the per sample nucleotide diversity files for each virus and calculate piA and piS
		#run as Rscript sample_piA_piS.R mpileup.list syn-nsyn.sites.tsv $virusid
		Rscript sample_piA_piS.R Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt Vesanto.virus.remapping.$seg/$virusstem.genewise.syn-nsyn.sites.tsv $virusid
		
		#And this one to calculate intergenic pi 
		Rscript sample_intergenic_pi.R Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt $virusid $refseqlength $virusstem.virus.analyses/$virusid.prokka.edited2.cds.gtf
		
		#Combining the Vesanto piA piS and intergenic results files for all segments
		cat Vesanto.virus.remapping.$seg/$virusid.per.sample.piA_piS.tsv | awk 'NR >1 {print $0}' >> Vesanto.segs.per.sample.piA.piS.results.tsv 
		cat Vesanto.virus.remapping.$seg/$virusid.per.sample.intergenic.pi.tsv | awk 'NR >1 {print $0}' >> Vesanto.segs.per.sample.intergenic.pi.results.tsv 
		
		#And now adding samples with no snps into a file for potential haplotype analysis
		for b in $(cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt)
		do
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/Vesanto.virus.remapping.$seg\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			
			no_nsyn_snps=$(cat Vesanto.virus.remapping.$seg/$virusid.per.sample.piA_piS.tsv | grep $popstem | cut -f 5)
			no_syn_snps=$(cat Vesanto.virus.remapping.$seg/$virusid.per.sample.piA_piS.tsv | grep $popstem | cut -f 6)
			no_intergenic_snps=$(cat Vesanto.virus.remapping.$seg/$virusid.per.sample.intergenic.pi.tsv | grep $popstem | cut -f 4 )
			
			if [[ $no_nsyn_snps == 0 ]] && [[ $no_syn_snps == 0 ]] && [[ $no_intergenic_snps == 0 ]]
			then
				echo -e "$virusid \t $popstem" >> no_snp_samples.tsv
			else 
				echo the sample $popstem of $virusstem virus contains some snps at greater than 0.01 allele frequency
			fi
		done 
		
		##Counting unique SNPs across the per sample bam files 
		for b in $(cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt)
		do 
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/Vesanto.virus.remapping.$seg\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			
			#adding intergenic into the type column of the SNPs so that it can be combined with the syn-nonsyn SNPs
			if [ -f Vesanto.virus.remapping.$seg/$popstem.$virusstem.intergenic.per.site.nucl.div.tsv ]
			then 
				while read -r seqname sample position ref cov A_count T_count C_count G_count A_freq A_freq_maf T_freq T_freq_maf C_freq C_freq_maf G_freq G_freq_maf site_nucl_div site_nucl_div_maf
				do
					echo -e "$seqname \t $sample \t $position \t $ref \t $cov \t $A_count \t $T_count \t $C_count \t $G_count \t intergenic \t NA \t NA \t NA \t $A_freq \t $A_freq_maf \t $T_freq \t $T_freq_maf \t $C_freq \t $C_freq_maf \t $G_freq \t $G_freq_maf \t $site_nucl_div \t $site_nucl_div_maf" >> Vesanto.virus.remapping.$seg/$popstem.$virusstem.intergenic.per.site.nucl.div_for_combination.tsv
				done < Vesanto.virus.remapping.$seg/$popstem.$virusstem.intergenic.per.site.nucl.div.tsv
			else
				echo the sample $popstem of $virusstem virus contains no nucl div file for intergenic regions
			fi
			
			#reading the intergenic nucl div file for each sample and filtering out lines where the nucl_div_maf = 0
			if [ -f Vesanto.virus.remapping.$seg/$popstem.$virusstem.intergenic.per.site.nucl.div_for_combination.tsv ]
			then 
				cat Vesanto.virus.remapping.$seg/$popstem.$virusstem.intergenic.per.site.nucl.div_for_combination.tsv | awk 'NR >1 {print $0}' | awk '{if ($23 > 0) print $0}' >> Vesanto.virus.remapping.$seg/$virusstem.combined.per.sample.maf1.nucl.div.tsv
			else
				echo the sample $popstem of $virusstem virus contains no nucl div file for intergenic regions
			fi
			
			#reading the syn-nonsyn nucl div file for each sample and filtering out lines where the nucl_div_maf = 0
			if [ -f Vesanto.virus.remapping.$seg/$popstem.$virusstem.syn-nsyn.per.site.nucl.div.tsv ]
			then 
				cat Vesanto.virus.remapping.$seg/$popstem.$virusstem.syn-nsyn.per.site.nucl.div.tsv | awk 'NR >1 {print $0}' | awk '{if ($23 > 0) print $0}' >> Vesanto.virus.remapping.$seg/$virusstem.combined.per.sample.maf1.nucl.div.tsv
			else
				echo the sample $popstem of $virusstem virus contains no nucl div file for coding regions
			fi
		done
		
		#counting the number of SNPs across all samples
		total_SNPs_across_samples=$(cat Vesanto.virus.remapping.$seg/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | wc -l )
		 
		#counting the number of unique position SNPs across all samples
		unique_SNPs_across_samples=$(cat Vesanto.virus.remapping.$seg/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq | wc -l )
		
		#counting the number of SNPs which only appear in 1 sample
		single_sample_SNPs=$(cat Vesanto.virus.remapping.$seg/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq -c | sort -n -k 1,1 | awk '{if($1<2)print $0}' | wc -l )
		
		#the number that only appear in 2
		two_sample_SNPs=$(cat Vesanto.virus.remapping.$seg/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq -c | sort -n -k 1,1 | awk '{if($1==2)print $0}' | wc -l )
		
		#the number that appear in more than 2
		more_than_two_sample_SNPs=$(cat Vesanto.virus.remapping.$seg/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq -c | sort -n -k 1,1 | awk '{if($1>2)print $0}' | wc -l )
		
		#the number of samples analysed for that segment variant
		no_samples=$( cat Vesanto.virus.remapping.$seg/$virusstem.mpileup.25.bam.list.txt | wc -l )
		
		#saving a tsv file in order to plot the no. of samples which each SNP is found in 
		cat Vesanto.virus.remapping.$seg/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq -c | sort -n -k 1,1 >> $Vesanto.virus.remapping.$seg/$virusstem.per.sample.snp.frequencies.tsv
		
		#printing stats for each virus
		echo there are $total_SNPs_across_samples SNPs across all the per sample bam files for $virusstem virus , $unique_SNPs_across_samples of which are from unique positions. $single_sample_SNPs are from only one sample, $two_sample_SNPs are from two samples, and $more_than_two_sample_SNPs appear in more than two samples
		
		#making a table of the vesanto segs and segment variants, and thier SNPs
		echo -e "$virusstem \t $no_samples \t $total_SNPs_across_samples \t $unique_SNPs_across_samples \t $single_sample_SNPs \t $two_sample_SNPs \t $more_than_two_sample_SNPs" >> Vesanto.segs.combined.SNP.numbers.tsv
		
	else
		
		#For Kallithea and Linvill Road
		for b in $(cat $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt)
		do 
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/$virusstem.virus.analyses\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			
			echo identifying indels in the pileup file for the $popstem mpileup file mapped to $virusstem virus
			#identify indels in each of the mpileup files
			perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup --output $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup.indel.regions.gtf --min-count 7
			
			echo masking indels in the pileup file for the $popstem mpileup file mapped to $virusstem virus
			#mask these indellic regions
			perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup --gtf $virusstem.virus.analyses/$popstem.$virusstem.bwa.mpileup.indel.regions.gtf --output $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.mpileup
			
			#Making an mpileup file with just the intergenic regions by masking the cds regions
			echo making an mpileup file with just the intergenic regions for $popstem sample of $virusstem
			perl popoolation_1.2.2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.mpileup --gtf $virusstem.virus.analyses/$virusid.edited.cds.gtf --output $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.intergenic.mpileup
			
			#calculating the nsyn and syn sites in each virus/segment using the sliding window syn-nsyn pi script
			echo calculating nsyn syn pi with a sliding window approach for $popstem sample of $virusstem virus 
			perl popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --fastq-type sanger --gtf $virusstem.virus.analyses/$virusstem.popoolation1.gtf --pileup $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.mpileup --codon-table popoolation_1.2.2/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table $virusstem.virus.analyses/snl.$virusstem.txt --output $virusstem.virus.analyses/$popstem.$virusstem.bwa.sliding.syn-nsyn.pi.tsv --snp-output $virusstem.virus.analyses/$popstem.$virusstem.bwa.sliding.syn-nsyn.snps.tsv --pool-size 40 --min-qual 40 --min-count 1 --min-coverage 5 --max-coverage 550 --window-size 30 --step-size 24
			
			#removing duplicate rows from snp file as sliding window was used
			cat $virusstem.virus.analyses/$popstem.$virusstem.bwa.sliding.syn-nsyn.snps.tsv | grep -v '^>' | grep "\S" | uniq > $virusstem.virus.analyses/$popstem.$virusstem.bwa.syn-nsyn.snps.tsv
			
			#running this Rscript to calculate per site nucleotide diversity for each snp in the sample
			#run as Rscript sample_per_site_nucl_div.R snp.file type_of_file $virusid $popstem
			
			type=cds
			
			if [ -s $virusstem.virus.analyses/$popstem.$virusstem.bwa.syn-nsyn.snps.tsv ]
			then
				echo the $popstem sample of $virusstem virus contains some snps in coding regions			
				Rscript sample_per_site_nucl_div.R $virusstem.virus.analyses/$popstem.$virusstem.bwa.syn-nsyn.snps.tsv $type $virusid $popstem
			else
				echo the $popstem sample of $virusstem virus contains no snps in coding regions	
				echo -e "$virusstem \t $popstem \t $type" >> no_snp_samples.txt
			fi
			
			#output from this will be the snps and thier per site nucleotide diversity
			
			#Calculating pi across the intergenic regions - to get the intergenic snps
			echo calculating pi across these intergenic regions for the $popstem sample of $virusstem virus 
			perl popoolation_1.2.2/popoolation_1.2.2/Variance-sliding.pl --measure pi --input $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.intergenic.mpileup --output $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.intergenic.sliding.pi.tsv --snp-output $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.intergenic.sliding.pi.snps.tsv --fastq-type sanger --pool-size 40 --min-count 1 --min-coverage 5 --max-coverage 550 --window-size 30 --step-size 24 --min-qual 40
			
			#removing duplicate rows from snp file as sliding window was used
			cat $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.intergenic.sliding.pi.snps.tsv | grep -v '^>' | grep "\S" | uniq > $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.intergenic.snps.tsv
			
			#running this Rscript to calculate per site nucleotide diversity for each snp in the sample
			#run as Rscript sample_per_site_nucl_div.R snp.file type_of_file $virusid $popstem
			
			type=intergenic
			
			if [ -s $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.intergenic.snps.tsv ]
			then
				echo the $popstem sample of $virusstem virus contains some snps in intergenic regions			
				Rscript sample_per_site_nucl_div.R $virusstem.virus.analyses/$popstem.$virusstem.bwa.indelfiltered.intergenic.snps.tsv $type $virusid $popstem
			else
				echo the $popstem sample of $virusstem virus contains no snps in intergenic regions	
			fi
			
			#output from this will be the snps and thier per site nucleotide diversity

		done
		
		#Now running this R script to cycle through the per sample nucleotide diversity files for each virus and calculate piA and piS
		#run as Rscript sample_piA_piS.R mpileup.list syn-nsyn.sites.tsv $virusid
		Rscript sample_piA_piS.R $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt $virusstem.virus.analyses/$virusstem.genewise.syn-nsyn.sites.tsv $virusid
		
		#And this one to calculate intergenic pi 
		Rscript sample_intergenic_pi.R  $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt $virusid $refseqlength $virusstem.virus.analyses/$virusid.edited.cds.gtf
		
		#And now adding samples with no snps into a file for potential haplotype analysis
		for b in $(cat $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt)
		do
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/$virusstem.virus.analyses\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			
			no_nsyn_snps=$(cat $virusstem.virus.analyses/$virusid.per.sample.piA_piS.tsv | grep $popstem | cut -f 5)
			no_syn_snps=$(cat $virusstem.virus.analyses/$virusid.per.sample.piA_piS.tsv | grep $popstem | cut -f 6)
			no_intergenic_snps=$(cat $virusstem.virus.analyses/$virusid.per.sample.intergenic.pi.tsv | grep $popstem | cut -f 4 )
			
			if [[ $no_nsyn_snps == 0 ]] && [[ $no_syn_snps == 0 ]] && [[ $no_intergenic_snps == 0 ]]
			then
				echo -e "$virusid \t $popstem" >> no_snp_samples.tsv
			else 
				echo the sample $popstem of $virusstem virus contains some snps at greater than 0.01 allele frequency
			fi
		done 
		
		##Counting unique SNPs across the per sample bam files 
		for b in $(cat $virusstem.virus.analyses/$virusstem.mpileup.25.bam.list.txt)
		do 
			bam=$(echo $b)
			popstem=$(echo $bam | sed "s/$virusstem.virus.analyses\///" | sed "s/.$virusstem.bwa.noTIRs_InDel.bam//")
			
			#adding intergenic into the type column of the SNPs so that it can be combined with the syn-nonsyn SNPs
			if [ -f $virusstem.virus.analyses/$popstem.$virusstem.intergenic.per.site.nucl.div.tsv ]
			then 
				while read -r seqname sample position ref cov A_count T_count C_count G_count A_freq A_freq_maf T_freq T_freq_maf C_freq C_freq_maf G_freq G_freq_maf site_nucl_div site_nucl_div_maf
				do
					echo -e "$seqname \t $sample \t $position \t $ref \t $cov \t $A_count \t $T_count \t $C_count \t $G_count \t intergenic \t NA \t NA \t NA \t $A_freq \t $A_freq_maf \t $T_freq \t $T_freq_maf \t $C_freq \t $C_freq_maf \t $G_freq \t $G_freq_maf \t $site_nucl_div \t $site_nucl_div_maf" >> $virusstem.virus.analyses/$popstem.$virusstem.intergenic.per.site.nucl.div_for_combination.tsv
				done < $virusstem.virus.analyses/$popstem.$virusstem.intergenic.per.site.nucl.div.tsv
			else
				echo the sample $popstem of $virusstem virus contains no nucl div file for intergenic regions
			fi
			
			#reading the intergenic nucl div file for each sample and filtering out lines where the nucl_div_maf = 0
			if [ -f $virusstem.virus.analyses/$popstem.$virusstem.intergenic.per.site.nucl.div_for_combination.tsv ]
			then 
				cat $virusstem.virus.analyses/$popstem.$virusstem.intergenic.per.site.nucl.div_for_combination.tsv | awk 'NR >1 {print $0}' | awk '{if ($23 > 0) print $0}' >> $virusstem.virus.analyses/$virusstem.combined.per.sample.maf1.nucl.div.tsv
			else
				echo the sample $popstem of $virusstem virus contains no nucl div file for intergenic regions
			fi
			
			#reading the syn-nonsyn nucl div file for each sample and filtering out lines where the nucl_div_maf = 0
			if [ -f $virusstem.virus.analyses/$popstem.$virusstem.syn-nsyn.per.site.nucl.div.tsv ]
			then 
				cat $virusstem.virus.analyses/$popstem.$virusstem.syn-nsyn.per.site.nucl.div.tsv | awk 'NR >1 {print $0}' | awk '{if ($23 > 0) print $0}' >> $virusstem.virus.analyses/$virusstem.combined.per.sample.maf1.nucl.div.tsv
			else
				echo the sample $popstem of $virusstem virus contains no nucl div file for coding regions
			fi
		done
		
		#counting the number of SNPs across all samples
		total_SNPs_across_samples=$(cat $virusstem.virus.analyses/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | wc -l )
		 
		#counting the number of unique position SNPs across all samples
		unique_SNPs_across_samples=$(cat $virusstem.virus.analyses/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq | wc -l )
		
		#counting the number of SNPs which only appear in 1 sample
		single_sample_SNPs=$(cat $virusstem.virus.analyses/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq -c | sort -n -k 1,1 | awk '{if($1<2)print $0}' | wc -l )
		
		#the number that only appear in 2
		two_sample_SNPs=$(cat $virusstem.virus.analyses/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq -c | sort -n -k 1,1 | awk '{if($1==2)print $0}' | wc -l )
		
		#the number that appear in more than 2
		more_than_two_sample_SNPs=$(cat $virusstem.virus.analyses/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq -c | sort -n -k 1,1 | awk '{if($1>2)print $0}' | wc -l )
		
		#saving a tsv file in order to plot the no. of samples which each SNP is found in 
		cat $virusstem.virus.analyses/$virusstem.combined.per.sample.maf1.nucl.div.tsv | sort -n -k 3,3 | cut -f 3 | sort -n | uniq -c | sort -n -k 1,1 >> $virusstem.virus.analyses/$virusstem.per.sample.snp.frequencies.tsv
		
		#printing stats for each virus
		echo there are $total_SNPs_across_samples SNPs across all the per sample bam files for $virusstem virus , $unique_SNPs_across_samples of which are from unique positions. $single_sample_SNPs are from only one sample, $two_sample_SNPs are from two samples, and $more_than_two_sample_SNPs appear in more than two samples
	fi
done 
