###########################################################
## Pre-processing bam files for virus diversity analyses ##
###########################################################

###Written by Megan A. Wallace in 07-09/2020 

##This script should be run in the same folder as all of the virusmap.bam files, which contain reads which have been competitively mapped to multiple viruses from each population
##Also needs to be a fasta file containing the reference genomes for the viruses of interest in the folder
##and a second fasta file with the non-TIR sequences, must have same ids as the viruses in the first fasta
##input = .fasta file containing the viruses and segments of viruses you want to assess data coverage for and prep fastq files for analysis
#run as ./Dros_DNA_viruses_preprocessing.sh DrosEU_DNA_viruses_franalysis.fasta DrosEU_DNA_viruses_franalysis_nonTIR.fasta

##Taking the .fasta of virus genomes
virusfasta="$1"
#and the fasta of virus genomes without TIRs
virusfastanoTIRs="$2"
#Creating a linearised version of the fastas for grepping etc. 
sed -e 's/\(^>.*$\)/#\1#/' $virusfasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > linear.virus.fasta
sed -e 's/\(^>.*$\)/#\1#/' $virusfastanoTIRs | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > linear.virus.noTIRs.fasta
#extracting all the identifier names from the fastas - creating a list of virus genomes or segments to potentially analyse
perl -ne 'if(/^>(\S+)/){print "$1\n"}' $virusfasta > virus.list.txt
echo created lists of viruses to analyse
###Creating a list of population samples to analyse from the .virusmap.bam files in the current folder
ls *.virusmap.bam | sed 's/.virusmap.bam$//g' > all.sample.list.txt
echo created sample list of all populations

##Assesing the number of reads mapping to each of the non-segmented viruses, and excluding viruses from the analysis with <6 samples with sufficient read nos.
for v in $(cat virus.list.txt)
do
	#assigning the variable that will be used to filter the bam files by virus
	virusid=$(echo $v)
	
	#will be used for file naming - creating longer stems with more detials for the segments of Vesanto
	isseg=$(echo $v | grep '^Vesanto')
	if [ -n "$isseg" ]
	then
		virusstem=$(echo $v | sed 's/\(_virus\)//' )
	else
		virusstem=$(echo $v | sed -r 's/^[A-Z]{2}[0-9]{6}_//' | sed 's/virus.*/virus/' | sed 's/\(_virus\|_Nudivirus\)//' )
	fi
	
	#fetching the virus sequence of interest with header in single line form, then retreiving the sequence, and measuring its length 
	refseqlength=$(grep -A1 "$virusid" linear.virus.fasta | awk '/^>/{getline; getline; print}' | awk '{print length($1)}')
	#calculating the no. of 151 bp reads needed for avergage 10x coverage across 95% of the genome, and rounding value up to whole number using awk
	perc="9.5"
	readtarget=$(awk -v perc="${perc}" -v refseqlength="$refseqlength" 'BEGIN{rt=((perc*refseqlength)/151); i=int(rt); print (rt-i<0.5)?i:i+1 }')
	echo $readtarget reads needed per bam file for sufficient $virusstem virus coverage, reference sequence is $refseqlength bp long 
	
	echo Now creating list of samples with sufficient $virusstem virus coverage
	touch $virusstem.sample.list.txt
	
	for s in $(cat all.sample.list.txt)
	do
		popstem=$(echo $s)
		#counting the number of reads mapping to the focal virus in the bam file, no mapping quality filter, but only including proper read pairs
		readcount=$(samtools view -b -@ 6 -c -f 2 $popstem.virusmap.bam $virusid)
		if [ $readcount -ge $readtarget ] 
		then
			echo $popstem >> $virusstem.sample.list.txt
		fi
	done	
	
	#storing the number of samples with enough virus coverage to analyse
	virussampleno=$( wc -l < $virusstem.sample.list.txt )
	echo $virussampleno samples have sufficient $virusstem virus coverage for analysis
	
	#creating a new list of virus ids with only those non-segmented viruses which have >5 samples with sufficient coverage, and segmented viruses with >0 samples with sufficient coverage
	if [ $virussampleno -gt 5 ]
		then
			echo $virusid >> gt5samples.virus.list.txt	
		else
			if [ -n "$isseg" ] && [ $virussampleno -gt 0 ]
				then
					echo $virusid >> gt5samples.virus.list.txt
				else
					echo too few samples have sufficient $virusstem virus coverage for analysis 
					echo too few samples have sufficient $virusstem virus coverage for analysis >> virus.preprocessing.log.txt
			fi
	fi
done

echo Now preparing the viruses with sufficient samples for analysis

for a in $(cat gt5samples.virus.list.txt)
do	
	virusid=$(echo $a)
	
	isseg=$(echo $a | grep '^Vesanto')
	if [ -n "$isseg" ]
	then
		virusstem=$(echo $a | sed 's/\(_virus\)//')
	else
		virusstem=$(echo $a | sed -r 's/^[A-Z]{2}[0-9]{6}_//' | sed 's/virus.*/virus/' | sed 's/\(_virus\|_Nudivirus\)//')
	fi
	
	#checking if there's a directory for focal virus analyses, and if not, creating one
	DIR="$virusstem.virus.analyses"
	if [ -d "$DIR" ] 
	then
		echo Directory "$DIR" already exists
	else
		mkdir "$DIR"
		echo created directory for $virusstem virus analyses
	fi
	
	##Now filtering the competitively mapped reads for those mapping to the focal virus only
	echo filtering the bam files with sufficient reads for $virusstem virus
	for s in $(cat $virusstem.sample.list.txt)
	do 
		popstem=$(echo $s)
		samtools view -b -@ 6 -f 2 $popstem.virusmap.bam $virusid > $virusstem.virus.analyses/$popstem.$virusstem.bam 
		samtools index $virusstem.virus.analyses/$popstem.$virusstem.bam
	done 
	
	##Now converting the filtered bam files to fastq, zipping them, splitting them into lanes, and then quality trimming them , before remapping to the focal virus genome 
	##Samtools option - so that I can use mutliple threads as bazam was having some heap size issues...creates 2 file, one w the F reads and one w the R reads, discarding singletons
	for s in $(cat $virusstem.sample.list.txt)
	do
		popstem=$(echo $s)
		echo converting $popstem.$virusstem.bam to fastq
		samtools sort -@8 -n $virusstem.virus.analyses/$popstem.$virusstem.bam | samtools bam2fq -@8 -c6 -1 $virusstem.virus.analyses/$popstem.$virusstem.1.fastq.gz -2 $virusstem.virus.analyses/$popstem.$virusstem.2.fastq.gz -0 /dev/null -s /dev/null -
		
		#for each sample, split the fastq files into lanes + flowcells, and then trim based on quality using trimadapt
		#Splitting the fastq files from the paired end reads into lanes using awk
		awk -v popstem="$popstem" -v virusstem="$virusstem" 'BEGIN {FS = ":"} {flowcell=$3 ; lane=$4 ; print > virusstem".virus.analyses/"popstem"."virusstem"."flowcell"."lane".1.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > virusstem".virus.analyses/"popstem"."virusstem"."flowcell"."lane".1.fastq"}}' <(gzip -dc $virusstem.virus.analyses/$popstem.$virusstem.1.fastq.gz)
		awk -v popstem="$popstem" -v virusstem="$virusstem" 'BEGIN {FS = ":"} {flowcell=$3 ; lane=$4 ; print > virusstem".virus.analyses/"popstem"."virusstem"."flowcell"."lane".2.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > virusstem".virus.analyses/"popstem"."virusstem"."flowcell"."lane".2.fastq"}}' <(gzip -dc $virusstem.virus.analyses/$popstem.$virusstem.2.fastq.gz)
		
		#trimming based on quality with cutadapt 
		#first making a list of the unique flowcell and lane combos for this sample in the virus folder 
		ls $virusstem.virus.analyses/$popstem.$virusstem.*.*.*.fastq | sed -r "s/$virusstem.virus.analyses\/$popstem.$virusstem.//" | sed 's/.1.fastq//' | sed 's/.2.fastq//' | uniq > $virusstem.virus.analyses/$popstem.$virusstem.lanelist.txt
		
		##if the samples are from 2014, a Nextseq trim is used to account for two colour chemistry, otherwise a normal -q 18 is used 
		is2014=$(echo $s | grep '^2014_')
		if [ -n "$is2014" ]
		then
			for l in $(cat $virusstem.virus.analyses/$popstem.$virusstem.lanelist.txt)
			do
				lane=$(echo $l)
				#zipping the fastq files from each of the lanes 
				gzip $virusstem.virus.analyses/$popstem.$virusstem.$lane.1.fastq 
				gzip $virusstem.virus.analyses/$popstem.$virusstem.$lane.2.fastq 
				#trimming with cutadapt
				cutadapt -j 4 --nextseq-trim=18 --minimum-length 75 --pair-filter=any -o $virusstem.virus.analyses/trimmed18_$popstem.$virusstem.$lane.1.fastq.gz -p $virusstem.virus.analyses/trimmed18_$popstem.$virusstem.$lane.2.fastq.gz $virusstem.virus.analyses/$popstem.$virusstem.$lane.1.fastq.gz $virusstem.virus.analyses/$popstem.$virusstem.$lane.2.fastq.gz
			done
		else
			for l in $(cat $virusstem.virus.analyses/$popstem.$virusstem.lanelist.txt)
			do
				lane=$(echo $l)
				#zipping the fastq files from each of the lanes
				gzip $virusstem.virus.analyses/$popstem.$virusstem.$lane.1.fastq 
				gzip $virusstem.virus.analyses/$popstem.$virusstem.$lane.2.fastq 
				#trimming with cutadapt
				cutadapt -j 4 -q 18 --minimum-length 75 --pair-filter=any -o $virusstem.virus.analyses/trimmed18_$popstem.$virusstem.$lane.1.fastq.gz -p $virusstem.virus.analyses/trimmed18_$popstem.$virusstem.$lane.2.fastq.gz $virusstem.virus.analyses/$popstem.$virusstem.$lane.1.fastq.gz $virusstem.virus.analyses/$popstem.$virusstem.$lane.2.fastq.gz
			done 
		fi
	done 
	
	#Now remapping the samples to the virus genomes using bwa -mem, and using a fasta of the Vesanto segments with the TIRs removed to re-map to 
	echo Creating fasta index and dictionary for $virusid taking the sequences from the non-TIR fasta 
	#extracting the fasta file in wrapped form for the virus of interest, will be useful for mapping later 
	grep -A1 "$virusid" linear.virus.noTIRs.fasta | fasta_formatter -w 60 > $virusid.noTIRs.fasta
	#also indexing the single virus fasta file 
	samtools faidx $virusid.noTIRs.fasta
	#and making a .dict file of the .fasta file using picard, for later indel realignment 
	java -jar /localdisk/home/s1667991/picard-2.22.8/picard.jar CreateSequenceDictionary R=$virusid.noTIRs.fasta O=$virusid.noTIRs.dict VERBOSITY=ERROR QUIET=true &> /dev/null
	
	#making a bwa index of the fasta
	#checking if theres a directory for focal virus bwa index, and if not, creating one
	DIR="$virusid.noTIRs.bwa.index"
	if [ -d "$DIR" ] 
	then
		echo Directory "$DIR" already exists
	else
		mkdir "$DIR"
		echo created directory for $virusstem virus without TIRs bwa index 
	fi
	
	bwa index -p $virusid.noTIRs.bwa.index/$virusid.noTIRs.bwa.index.fa $virusid.noTIRs.fasta
	
	#now mapping these reads to the virus reference genome without TIRs using bwa mem
	#added read group at this point too
	for s in $(cat $virusstem.sample.list.txt)
	do
		popstem=$(echo $s)
		#for library I've just used the $popstem as each sample had its own lib
		#ID and PU now includes the sample for freebayes sample identification - so that multiple samples dont have the same ID if they were on the same lane
		#filtering for only MAPQ >30, no unmapped reads, or reads where the mate is unmapped, and no non-primary alignments
		echo mapping $popstem reads to the $virusstem virus reference genome with no TIRs and adding read groups before sorting and quality filtering the bam using samtools 
		for l in $(cat $virusstem.virus.analyses/$popstem.$virusstem.lanelist.txt)
		do
			lane=$(echo $l)
			bwa mem -M -t 12 -v 1 -R $(echo "@RG\tID:$lane.$popstem\tSM:$popstem\tLB:$popstem\tPL:ILLUMINA\tPU:$lane.$popstem") /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.bwa.index/$virusid.noTIRs.bwa.index.fa $virusstem.virus.analyses/trimmed18_$popstem.$virusstem.$lane.1.fastq.gz $virusstem.virus.analyses/trimmed18_$popstem.$virusstem.$lane.2.fastq.gz | samtools view -b -@ 6 -q 30 -u -F 4 -F 8 -F 0x100 - | samtools sort -n -@ 6 - | samtools fixmate -cm -@ 6 - - | samtools sort -@ 6 -m 6G - > $virusstem.virus.analyses/$popstem.$virusstem.$lane.bwa.noTIRs.bam 
			#removing the untrimmed fastq.gz files to save space 
			rm $virusstem.virus.analyses/$popstem.$virusstem.$lane.1.fastq.gz 
			rm $virusstem.virus.analyses/$popstem.$virusstem.$lane.2.fastq.gz 
		done
		
	#####################
	#removing duplicates#
	#####################
	
		echo removing duplicates and merging the $virusstem virus bam files from $popstem 
		#making an array of the bam files to be merged + dedupped, for each sample
		ls $virusstem.virus.analyses/$popstem.$virusstem.*.*.bwa.noTIRs.bam > $virusstem.virus.analyses/$popstem.$virusstem.noTIRs.bam.list.txt
		bam_array=($(sed -e 's/^/INPUT=/' $virusstem.virus.analyses/$popstem.$virusstem.noTIRs.bam.list.txt))
		#For the 2014 samples, using the default optical duplicate pixel distance of 100, but for 2015 + 2016 using 2500 as they were run on the HiSeq
		is2014=$(echo $s | grep '^2014_')
		if [ -n "$is2014" ]
		then
			java -jar /localdisk/home/s1667991/picard-2.22.8/picard.jar MarkDuplicates $(echo ${bam_array[@]}) REMOVE_DUPLICATES=true OUTPUT=$virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.dedup.bam METRICS_FILE=$virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.dedup.metrics.txt VALIDATION_STRINGENCY=SILENT
		else
			java -jar /localdisk/home/s1667991/picard-2.22.8/picard.jar MarkDuplicates $(echo ${bam_array[@]}) OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 REMOVE_DUPLICATES=true OUTPUT=$virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.dedup.bam METRICS_FILE=$virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.dedup.metrics.txt VALIDATION_STRINGENCY=SILENT
		fi
	
		##preparing the bam file for final assesment of coverage, and indel realignment by removing any remaining unpaired reads, coordinate sorting and indexing the bam file 
		echo doing a final removal of any unpaired reads and sorting the $popstem $virusstem virus bam 
		samtools view -h -f 2 -@ 6 $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.dedup.bam | samtools sort -@ 6 - > $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.bam
		
		##for some of the Vesanto bams remapped to non-TIR regions, by this stage no reads are left in the bam as all reads previously in them were stacked in the TIRs
		##removing these bam files from the analysis by updating the sample list for each virus
		bwareadno=$(samtools view -c $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.bam)
		if [ $bwareadno == 0 ]
		then 
			echo sample $popstem has no reads mapping to the $virusstem virus genome without TIRs
		else 
			echo sample $popstem has some reads mapping to the $virusstem virus genome with TIRs
			echo $popstem >> $virusstem.second.sample.list.txt
		fi
	done 
	
	###########################################################################################################
	## final check on each of the filtered bams to see if they have >10 rd in 95% of the genome for analysis ##
	###########################################################################################################
	
	echo Now creating list of samples with sufficient $virusstem virus coverage in the non TIR bam files 
	#fetching the virus sequence of interest with header in single line form, then retreiving the sequence, and measuring its length 
	refseqlength=$(grep -A1 "$virusid" linear.virus.noTIRs.fasta | awk '/^>/{getline; getline; print}' | awk '{print length($1)}')
	#calculating the bp consitituting 95% of the genome length, and rounding value up to whole number using awk
	perc="0.95"
	bptarget=$(awk -v perc="${perc}" -v refseqlength="$refseqlength" 'BEGIN{bpt=(perc*refseqlength); i=int(bpt); print (bpt-i<0.5)?i:i+1}' )
	
	#calculating the read depth across the bam file for each sample bam 
	echo calculating the per position read depth for the $virusstem virus non TIR bams 
	for s in $(cat $virusstem.second.sample.list.txt)
	do
		popstem=$(echo $s)
		#using bedtools to measure per position read depth and count the positions where read depth is > ten 
		depthcount=$(genomeCoverageBed -d -ibam $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.bam | awk '$3 > 9'| grep -w "$virusid" -c)
		echo in the $popstem sample $depthcount bp of the $virusstem virus genome have sufficient read depth, $bptarget bp are needed for analyses
		if [ $depthcount -ge $bptarget ] 
		then
			echo sample $popstem can be used for $virusstem virus diversity analyses
			echo $popstem >> $virusstem.final.sample.list.txt
		else
			echo sample $popstem has insufficient read depth for $virusstem virus diversity analyses 
		fi	
	done
done 
	
######################
## Vesanto segments ##
######################
	
##For each $popstem, checking the $virusstem.final.sample.list.txt files for all the segment variants of Vesanto, to see if its in more than one of them, and if it is checking which haplotype has more reads and selecting that as the only haplotype for that sample
echo Now examining the read depth of the Vesanto virus segments and ensuring only one haplotype is mapped to in each sample
	
for b in $(cat all.sample.list.txt)
do
	popstem=$(echo $b)
	
	for c in $(cat gt5samples.virus.list.txt)
	do
	
	virusid=$(echo $c)
	
	isseg=$(echo $c | grep '^Vesanto')
	
	if [ -n "$isseg" ]
	then
		virusstem=$(echo $c | sed 's/\(_virus\)//' )
		FILE="$virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.bam"
		if [ -f "$FILE" ]
		then
			varreadcount=$(samtools view -b -@ 6 -c $FILE)
			echo -e "$popstem \t $virusstem \t $varreadcount" >> var_read_counts.txt
		else
			echo $popstem has insufficient reads of $virusstem variant 
		fi
	else
		virusstem=$(echo $c | sed -r 's/^[A-Z]{2}[0-9]{6}_//' | sed 's/virus.*/virus/' | sed 's/\(_virus\|_Nudivirus\)//')
		echo no need to do this for $virusstem virus
	fi
	done
done	

#Making txt list of segments to cycle through
cat gt5samples.virus.list.txt | grep '^Vesanto' | sed 's/^.*_S/S/g' | sort | uniq > seg.no.list.txt
	
for d in $(cat seg.no.list.txt)
do
	seg=$(echo $d)
	##making txt file w only the seg of interest
	cat var_read_counts.txt | awk -v seg="$seg" '$2 ~ seg {print}' > var_read_counts_$seg.txt
	for e in $(cat var_read_counts_$seg.txt | cut -f1 | sort | uniq )
	do
		popstem=$(echo $e)
		maxsegpop=$(cat var_read_counts_$seg.txt | awk -v popstem="$popstem" '$1 ~ popstem {print}' | awk -v max=-1 '{if($3>max){want=$0; max=$3}}END{print want} ' | awk '{if($3>0){print}}')
		if [ -n "$maxsegpop" ]
		then
			echo $maxsegpop >> pop_max_vars.txt
		fi
	done
done 	

#######################
## Vesanto remapping ##
#######################
 
for g in $(cat seg.no.list.txt)
do
	seg=$(echo $g)
	
	#checking if there's a directory for this segments analyses, and if not, creating one
	DIR="Vesanto.virus.remapping.$seg"
	if [ -d "$DIR" ] 
	then
		echo Directory "$DIR" already exists
	else
		mkdir "$DIR"
		echo created directory for Vesanto $seg remapping
	fi
	
	##Filtering the competitively mapped Vesanto reads in each population in the segment lists for all reads mapping to each segment
	for h in $(cat var_read_counts_$seg.txt | cut -f1 | sort | uniq )
	do
		popstem=$(echo $h)
		#creating array of virusids for sorting the bam file
		var_array=($(awk -v seg="$seg" '$1 ~ seg {print}' gt5samples.virus.list.txt))
		echo filtering the $popstem virusmap.bam file for reads mapping to all variants of Vesanto $seg
		samtools view -b -@ 6 -f 2 $popstem.virusmap.bam $(echo ${var_array[@]}) > Vesanto.virus.remapping.$seg/$popstem.$seg.bam 
		samtools index Vesanto.virus.remapping.$seg/$popstem.$seg.bam
	
		##Now converting the filtered bam files to fastq, zipping them, splitting them into lanes, and then quality trimming them , before remapping to the focal virus genome
		echo converting $popstem.$seg.bam to fastq
		samtools sort -@8 -n Vesanto.virus.remapping.$seg/$popstem.$seg.bam | samtools bam2fq -@8 -c6 -1 Vesanto.virus.remapping.$seg/$popstem.$seg.1.fastq.gz -2 Vesanto.virus.remapping.$seg/$popstem.$seg.2.fastq.gz -0 /dev/null -s /dev/null -
		
		#for each seg, split the fastq files into lanes + flowcells, and then trim based on quality using trimadapt
		#Splitting the fastq files from the paired end reads into lanes using awk
		awk -v popstem="$popstem" -v seg="$seg" 'BEGIN {FS = ":"} {flowcell=$3 ; lane=$4 ; print > "Vesanto.virus.remapping."seg"/"popstem"."seg"."flowcell"."lane".1.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "Vesanto.virus.remapping."seg"/"popstem"."seg"."flowcell"."lane".1.fastq"}}' <(gzip -dc Vesanto.virus.remapping.$seg/$popstem.$seg.1.fastq.gz)
		awk -v popstem="$popstem" -v seg="$seg" 'BEGIN {FS = ":"} {flowcell=$3 ; lane=$4 ; print > "Vesanto.virus.remapping."seg"/"popstem"."seg"."flowcell"."lane".2.fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "Vesanto.virus.remapping."seg"/"popstem"."seg"."flowcell"."lane".2.fastq"}}' <(gzip -dc Vesanto.virus.remapping.$seg/$popstem.$seg.2.fastq.gz)
		
		#trimming based on quality with cutadapt 
		#first making a list of the unique flowcell and lane combos for this sample in the virus folder 
		ls Vesanto.virus.remapping.$seg/$popstem.$seg.*.*.*.fastq | sed -r "s/Vesanto.virus.remapping.$seg\/$popstem.$seg.//" | sed 's/.1.fastq//' | sed 's/.2.fastq//' | uniq > Vesanto.virus.remapping.$seg/$popstem.$seg.lanelist.txt
		
		##if the samples are from 2014, a Nextseq trim is used to account for two colour chemistry, otherwise a normal -q 18 is used 
		is2014=$(echo $h | grep '^2014_')
		if [ -n "$is2014" ]
		then
			for i in $(cat Vesanto.virus.remapping.$seg/$popstem.$seg.lanelist.txt)
			do
				lane=$(echo $i)
				#zipping the fastq files from each of the lanes 
				gzip Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.1.fastq 
				gzip Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.2.fastq 
				#trimming with cutadapt
				cutadapt -j 4 --nextseq-trim=18 --minimum-length 75 --pair-filter=any -o Vesanto.virus.remapping.$seg/trimmed18_$popstem.$seg.$lane.1.fastq.gz -p Vesanto.virus.remapping.$seg/trimmed18_$popstem.$seg.$lane.2.fastq.gz Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.1.fastq.gz Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.2.fastq.gz
			done
		else
			for i in $(cat Vesanto.virus.remapping.$seg/$popstem.$seg.lanelist.txt)
			do
				lane=$(echo $i)
				#zipping the fastq files from each of the lanes
				gzip Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.1.fastq 
				gzip Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.2.fastq 
				#trimming with cutadapt
				cutadapt -j 4 -q 18 --minimum-length 75 --pair-filter=any -o Vesanto.virus.remapping.$seg/trimmed18_$popstem.$seg.$lane.1.fastq.gz -p Vesanto.virus.remapping.$seg/trimmed18_$popstem.$seg.$lane.2.fastq.gz Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.1.fastq.gz Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.2.fastq.gz
			done 
		fi
	done
	
	#Now remapping the samples to the virus genomes using bwa -mem, and using a fasta of the Vesanto segments with the TIRs removed to re-map to
	#This time, the reads are only mapped to the segment variant fasta with the highest readcount in the competitive mapping for Vesanto
	
	#making a txt file for each segment with the pop_max_vars in it, and making sure its tab separated
	cat pop_max_vars.txt | awk -v seg="$seg" '$2 ~ seg {print}' | sed 's/ /\t/g' > Vesanto.virus.remapping.$seg/pop_max_vars_$seg.txt
	
	while read -r popstem virusstem readcount
	do
		#this virusid is now the segment variant with the greatest read count in the competitive mapping 
		virusid=$(echo $virusstem | sed 's/Vesanto_/Vesanto_virus_/' )
		
		#now mapping these reads to the virus reference genome without TIRs using bwa mem
		#added read group at this point too
		#for library I've just used the $popstem as each sample had its own lib
		#ID and PU now includes the sample for freebayes sample identification - so that multiple samples dont have the same ID if they were on the same lane
		#filtering for only MAPQ >30, no unmapped reads, or reads where the mate is unmapped, and no non-primary alignments
		echo mapping $popstem reads to the $virusstem virus reference genome with no TIRs and adding read groups before sorting and quality filtering the bam using samtools 
		for l in $(cat Vesanto.virus.remapping.$seg/$popstem.$seg.lanelist.txt)
		do
			lane=$(echo $l)
			bwa mem -M -t 12 -v 1 -R $(echo "@RG\tID:$lane.$popstem\tSM:$popstem\tLB:$popstem\tPL:ILLUMINA\tPU:$lane.$popstem") /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.bwa.index/$virusid.noTIRs.bwa.index.fa Vesanto.virus.remapping.$seg/trimmed18_$popstem.$seg.$lane.1.fastq.gz Vesanto.virus.remapping.$seg/trimmed18_$popstem.$seg.$lane.2.fastq.gz | samtools view -b -@ 6 -q 30 -u -F 4 -F 8 -F 0x100 - | samtools sort -n -@ 6 - | samtools fixmate -cm -@ 6 - - | samtools sort -@ 6 -m 6G - > Vesanto.virus.remapping.$seg/$popstem.$virusstem.$lane.bwa.noTIRs.bam 
			#removing the untrimmed fastq.gz files to save space 
			rm Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.1.fastq.gz 
			rm Vesanto.virus.remapping.$seg/$popstem.$seg.$lane.2.fastq.gz 
		done
		
		##Removing Duplicates 
		echo removing duplicates and merging the $virusstem virus bam files from $popstem
		#making an array of the bam files to be merged + dedupped, for each sample
		ls Vesanto.virus.remapping.$seg/$popstem.$virusstem.*.*.bwa.noTIRs.bam > Vesanto.virus.remapping.$seg/$popstem.$virusstem.noTIRs.bam.list.txt
		bam_array=($(sed -e 's/^/INPUT=/' Vesanto.virus.remapping.$seg/$popstem.$virusstem.noTIRs.bam.list.txt))
		#For the 2014 samples, using the default optical duplicate pixel distance of 100, but for 2015 + 2016 using 2500 as they were run on the HiSeq
		is2014=$(echo $popstem | grep '^2014_')
		if [ -n "$is2014" ]
		then
			java -jar /localdisk/home/s1667991/picard-2.22.8/picard.jar MarkDuplicates $(echo ${bam_array[@]}) REMOVE_DUPLICATES=true OUTPUT=Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.dedup.bam METRICS_FILE=Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.dedup.metrics.txt VALIDATION_STRINGENCY=SILENT
		else
			java -jar /localdisk/home/s1667991/picard-2.22.8/picard.jar MarkDuplicates $(echo ${bam_array[@]}) OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 REMOVE_DUPLICATES=true OUTPUT=Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.dedup.bam METRICS_FILE=Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.dedup.metrics.txt VALIDATION_STRINGENCY=SILENT
		fi
		
		##preparing the bam files for final assesment of coverage, and indel realignment by removing any remaining unpaired reads, coordinate sorting and indexing the bam file 
		echo doing a final removal of any unpaired reads and sorting the $popstem $virusstem virus bam which now should contain all Vesanto $seg reads
		samtools view -h -f 2 -@ 6 Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.dedup.bam | samtools sort -@ 6 - > Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.bam
		
		#rename the old second sample list to indicate that it was generated using the comptitive mapping 
		FILE="$virusstem.compmap.second.sample.list.txt"
		if [ -f "$FILE" ]
		then 
			echo $virusstem compmap second sample list already exists
		else
			mv $virusstem.second.sample.list.txt $virusstem.compmap.second.sample.list.txt
		fi
		
		##for some of the Vesanto bams remapped to non-TIR regions, by this stage no reads are left in the bam as all reads previously in them were stacked in the TIRs
		##removing these bam files from the analysis by updating the sample list for each virus
		bwareadno=$(samtools view -c Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.bam)
		if [ $bwareadno == 0 ]
		then 
			echo sample $popstem has no reads mapping to the $virusstem virus genome without TIRs
		else 
			echo sample $popstem has some reads mapping to the $virusstem virus genome without TIRs
			#this sample list 
			echo $popstem >> $virusstem.second.sample.list.txt
		fi
	done < Vesanto.virus.remapping.$seg/pop_max_vars_$seg.txt
		
	###########################################################################################################
	## final check on each of the filtered bams to see if they have >10 rd in 95% of the genome for analysis ##
	###########################################################################################################
	
	for k in $(cat Vesanto.virus.remapping.$seg/pop_max_vars_$seg.txt | cut -f2 | sort | uniq )
	do
		virusstem=$(echo $k)
		virusid=$(echo $k | sed 's/Vesanto_/Vesanto_virus_/' )
		
		#rename the old final sample list to indicate that it was generated using the competitive mapping 
		FILE="$virusstem.compmap.final.sample.list.txt"
		if [ -f "$FILE" ]
		then 
			echo $virusstem compmap final sample list already exists
		else
			mv $virusstem.final.sample.list.txt $virusstem.compmap.final.sample.list.txt
		fi
		
		echo Now creating list of samples with sufficient $virusstem virus coverage in the non TIR bam files 
		
		#fetching the virus sequence of interest with header in single line form, then retreiving the sequence, and measuring its length 
		refseqlength=$(grep -A1 "$virusid" linear.virus.noTIRs.fasta | awk '/^>/{getline; getline; print}' | awk '{print length($1)}')
		#calculating the bp consitituting 95% of the genome length, and rounding value up to whole number using awk
		perc="0.95"
		bptarget=$(awk -v perc="${perc}" -v refseqlength="$refseqlength" 'BEGIN{bpt=(perc*refseqlength); i=int(bpt); print (bpt-i<0.5)?i:i+1}' )
		
		#calculating the read depth across the bam file for each sample bam 
		echo calculating the per position read depth for the $virusstem virus non TIR bams 
		for s in $(cat $virusstem.second.sample.list.txt)
		do
			popstem=$(echo $s)
			#using bedtools to measure per position read depth and count the positions where read depth is > ten 
			depthcount=$(genomeCoverageBed -d -ibam Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.bam | awk '$3 > 9'| grep -w "$virusid" -c)
			echo in the $popstem sample $depthcount bp of the $virusstem virus genome have sufficient read depth, $bptarget bp are needed for analyses
			if [ $depthcount -ge $bptarget ] 
			then
				echo sample $popstem can be used for $virusstem virus diversity analyses
				echo $popstem >> $virusstem.final.sample.list.txt
			else
				echo sample $popstem has insufficient read depth for $virusstem virus diversity analyses 
			fi	
		done
	done
done		

	###############################
	## realignment around indels ## 
	###############################

##On the non-Vesanto viruses 

for a in $(cat gt5samples.virus.list.txt | grep -v 'Vesanto' )
do	
	virusstem=$(echo $a | sed -r 's/^[A-Z]{2}[0-9]{6}_//' | sed 's/virus.*/virus/' | sed 's/\(_virus\|_Nudivirus\)//')
	virusid=$(echo $a)
		
	source activate gatk3_env 
	for s in $(cat $virusstem.final.sample.list.txt)
	do
		popstem=$(echo $s)
		#indexing for realignment around indels 
		samtools index $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.bam $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.bai
		echo realigning sample $popstem mapped to $virusstem virus around indels 
		gatk3 -T RealignerTargetCreator -R /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -I $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.bam -o $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.indeltarget.list
		##Now realigning around these indels
		gatk3 -T IndelRealigner -R /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -I $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.bam -targetIntervals $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs.indeltarget.list -o $virusstem.virus.analyses/$popstem.$virusstem.bwa.noTIRs_InDel.bam
		##this also creates a new bam.bai file 
	done 
	source deactivate 	
done	

##And now on the Vesanto non-competitively mapped bams

for m in $(cat seg.no.list.txt)
do
	source activate gatk3_env
	seg=$(echo $m)
	
	for n in $(cat Vesanto.virus.remapping.$seg/pop_max_vars_$seg.txt | cut -f2 | sort | uniq )
	do
		virusstem=$(echo $n)
		virusid=$(echo $n | sed 's/Vesanto_/Vesanto_virus_/' )
		
		for p in $(cat $virusstem.final.sample.list.txt)
		do
			popstem=$(echo $p)
		
			#indexing for realignment around indels 
			samtools index Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.bam Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.bai	
			
			#echo realigning sample $popstem mapped to $virusstem virus around indels 
			gatk3 -T RealignerTargetCreator -R /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -I Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.bam -o Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.indeltarget.list
			##Now realigning around these indels
			gatk3 -T IndelRealigner -R /mnt/drive3-6tb/Kallithea_Diversity/$virusid.noTIRs.fasta -I Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.bam -targetIntervals Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs.indeltarget.list -o Vesanto.virus.remapping.$seg/$popstem.$virusstem.bwa.noTIRs_InDel.bam
			##this also creates a new bam.bai file 
			
		done
		
	done
	source deactivate
done	
			
