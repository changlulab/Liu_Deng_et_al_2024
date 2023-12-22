#!/bin/bash

start=`date +%s`

QUAL=1

ZIP_TRIM=1

BLACKLIST=1
PREP_INPUT=1
PROC_SAMPLE=1

RUN_CORR=1

RUN_MACS=0
RUN_MACS2=1
RUN_EPIC2=1

REM_OUTPUT=1

GENOME=GENOMEGOESHERE

ODINISTHEBESTDOGEVER

if [ "$ZIP_TRIM" = "1" ]; then
gunzip *.gz
trim_galore *.fastq
fi

cd ..

############################### EXTENSION SETUP #################################

FILES=$PWD/Raw_Data/*.fq
SAM=.sam
BAM=.bam
RBAM=_raw.bam
WIG=.wig
BW=.bw
RBED=_raw.bed
BD=.bed
TRIM=_trimmed.fq
PRO=_promotor.bed
CORPRO=_corr_promotor.bed
WIN=_window.bed
SORT=_sorted.bed
CORSORT=_corr_sorted.bed
NOR=_normalized.bedGraph
CORNOR=_corr_norm.bedGraph
AVG=_Normalized_avg.bed
BED=_Aligned.bed
EXT=_Extended.bed
NORBAM=_Normalized_BAM.bam
NOALBAM=_Normalized_Aligned_BAM.bam
NOWIG=_Normalized.wig
NOBW=_Normalized.bw
NOSHWIG=_Normalized_Short.wig
NOSHBW=_Normalized_Short.bw
FULLSORT=_Full_Sorted.bed
ALBED=_Aligned.bed
FULLNOR=_Full_Normalized.bedGraph
FULLWIN=_Full_Window.bed
FULLPRO=_Full_Promotor.bed
FULLNOBW=_Full_Normalized.bw
FULLNOSHWIG=_Full_Short_Normalized.wig
FULLNOSHBW=_Full_Short_Normalized.bw
FULLNOWIG=_Full_Normalized.wig
FULLNOALBAM=_Full_Normalized_Aligned.bam
FULLNORBAM=_Full_Normalized.bam
FULLAVG=_Full_Normalized_Average.bed
EXTWIN=_Extended_Window.bed
EXTSORT=_Extended_Sort.bed
EXTBW=_Full_Normalized_Extended.bw

################################# FOLDER SETUP ##################################
mkdir Full_Sorted
mkdir Aligned_SAM
mkdir Aligned_BAM
mkdir Aligned_BED
mkdir Normalized
mkdir MACS
mkdir BW
mkdir WIG
mkdir Sort_BAM
mkdir Promotor
mkdir Window_BED
mkdir Sorted
#mkdir Average_Norm
mkdir Normalized_BAM
mkdir Normalized_WIG
mkdir Normalized_BW
#mkdir Normalized_Short_BW
mkdir Full_Normalized
mkdir Full_Window_BED
mkdir Analyzed
mkdir Correlation
mkdir Input
mkdir MACS2
mkdir EPIC2
cd MACS
mkdir P1
mkdir P3
mkdir P5
cd ..

######################### PREPPING THE INPUT FILE ###############################
#Input file is gotten from run script
#Align the input file to the genome
#Isolate unique bam
#Find the 100nt and 4000nt count windows
#Extend the input 100bp on either side of the sequence then 4000nt count windows
#################################################################################
DIR=$HOME/Data/Promotors/$GENOME

input=GET_INPUT

input_base=`basename "$input" .fastq`

echo "Input file: $input_base"

if [ "$PREP_INPUT" = "1" ]; then

bowtie -t -p 16 -S ~/Data/Indexes/$GENOME/$GENOME $PWD/Raw_Data/$input_base$TRIM -S $PWD/Input/Input.sam

if [ "$BLACKLIST" = "1" ]; then
	samtools view -@ 16 -bq $QUAL $PWD/Input/Input.sam > $PWD/Input/Input_raw.bam
	bedtools bamtobed -i $PWD/Input/Input_raw.bam > $PWD/Input/Input_raw.bed
	bedtools intersect -v -b $DIR/blacklist.bed -a $PWD/Input/Input_raw.bed > $PWD/Input/Input.bed
	bedtools bedtobam -i $PWD/Input/Input.bed -g $DIR/Genome.bed > $PWD/Input/Input.bam
fi

if ["$BLACKLIST" = "0" ]; then
	samtools view -@ 16 -bq $QUAL $PWD/Input/Input.sam > $PWD/Input/Input.bam
	bedtools bamtobed -i $PWD/Input/Input.bam > $PWD/Input/Input.bed
fi

bedtools coverage -counts -b $PWD/Input/Input.bam -a $DIR/Promotor_Window.bed > $PWD/Input/Window_Promotor_Input.bed
bedtools coverage -counts -b $PWD/Input/Input.bam -a $DIR/GenomeWindowed.bed > $PWD/Input/Window_Full_Input.bed
bedtools coverage -counts -b $PWD/Input/Input.bam -a $DIR/Promotor.bed > $PWD/Input/Promotor_Input.bed

#bedtools bamtobed -i $PWD/Input/Input.bam > $PWD/Input/Input.bed
bedtools slop -i $PWD/Input/Input.bed -g $DIR/Genome.bed -b 100 > $PWD/Input/Input_Extended.bed
bedtools coverage -counts -b $PWD/Input/Input_Extended.bed -a $DIR/GenomeWindowed.bed > $PWD/Input/Window_Full_Extended_Input.bed

fi

######################## PROCESSING THE SAMPLE FILE #############################
#Get the input length (for normalization)
#Loop through each sample.
#Get the core name of the file
#Align the sample to the genome
#Obtain the unique bam file
#################################################################################
input_length=$(samtools view -c "$PWD/Input/Input.bam")

for fn in $FILES
do

echo `basename "$fn" _trimmed.fq`
f=`basename "$fn" _trimmed.fq`

if [ "$f" = "$input_base" ]; then
continue
fi

if [ "$PROC_SAMPLE" = "1" ]; then
echo "Processing $f"

bowtie -t ~/Data/Indexes/$GENOME/$GENOME -p 16 -S $PWD/Raw_Data/$f$TRIM $PWD/Aligned_SAM/$f$SAM

if [ "$BLACKLIST" = "1" ]; then
	samtools view -@ 16 -bq $QUAL $PWD/Aligned_SAM/$f$SAM > $PWD/Aligned_BAM/$f$RBAM
	bedtools bamtobed -i $PWD/Aligned_BAM/$f$RBAM > $PWD/Aligned_BED/$f$RBED
	bedtools intersect -v -b $DIR/blacklist.bed -a $PWD/Aligned_BED/$f$RBED > $PWD/Aligned_BED/$f$BED
	bedtools bedtobam -i $PWD/Aligned_BED/$f$BED -g $DIR/Genome.bed > $PWD/Aligned_BAM/$f$BAM
fi

if [ "$BLACKLIST" = "0" ]; then
	samtools view -@ 16 -bq $QUAL $PWD/Aligned_SAM/$f$SAM > $PWD/Aligned_BAM/$f$BAM
	bedtools bamtobed -i $PWD/Aligned_BAM/$f$BAM > $PWD/Aligned_BED/$f$BED
fi

fi

############################## CROSS CORRELATION ################################
# This normalization is only for cross-correlation
# Sort the input promotor file (4000nt window)
# normalize the sample promotor file against input promotor file
#################################################################################

if [ "$RUN_CORR" = "1" ]; then

ChIP_length=$(samtools view -c $PWD/Aligned_BAM/$f$BAM)

bedtools coverage -counts -b $PWD/Aligned_BAM/$f$BAM -a $DIR/Promotor.bed > $PWD/Promotor/$f$CORPRO
sort -k1,1 -k2,2g -u -o $PWD/Input/Promotor_Input_Sorted.bed $PWD/Input/Promotor_Input.bed
sort -k1,1 -k2,2g -u -o $PWD/Sorted/$f$CORSORT $PWD/Promotor/$f$CORPRO
paste $PWD/Sorted/$f$CORSORT $PWD/Input/Promotor_Input_Sorted.bed | awk -v OFS="\t" '{print $1,$2,$3,$4,$5/'$ChIP_length'*1000000-$10/'$input_length'*1000000}' > $PWD/Correlation/$f$CORNOR

fi

################################## MACS 1.4 #####################################

if [ "$RUN_MACS" = "1" ]; then

P3=p3
P1=p1
P5=p5

if [ "$GENOME" = "hg19" ] || [ "$GENOME" = "hg38" ]; then
	macs14 -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g hs -n $PWD/MACS/P1/$f -p 1e-1
	macs14 -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g hs -n $PWD/MACS/P3/$f -p 1e-3
	macs14 -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g hs -n $PWD/MACS/P5/$f -p 1e-5
else
	macs14 -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g mm -n $PWD/MACS/P1/$f -p 1e-1
	macs14 -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g mm -n $PWD/MACS/P3/$f -p 1e-3
	macs14 -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g mm -n $PWD/MACS/P5/$f -p 1e-5
fi

fi

#################################### MACS2 ######################################

if [ "$RUN_MACS2" = "1" ]; then

if [ "$GENOME" = "hg19" ] || [ "$GENOME" = "hg38" ]; then
	macs2 callpeak -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g hs -n $f -q 0.05 --outdir $PWD/MACS2
else
	macs2 callpeak -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g mm -n $f -q 0.05 --outdir $PWD/MACS2
fi

fi

################################### EPIC2 #######################################

if [ "$RUN_EPIC2" = "1" ]; then

epic2 -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -gn $GENOME -fs 250 --gaps-allowed 4 -fdr 1e-20 -o $PWD/EPIC2/$f".EPIC2.bed"

fi

################################### BW OUTPUT ###################################

if [ "$REM_OUTPUT" = "1" ]; then

samtools sort -@ 16 -o $PWD/Sort_BAM/$f$BAM $PWD/Aligned_BAM/$f$BAM
samtools index $PWD/Sort_BAM/$f$BAM
#igvtools count -w 25 -f mean -e 200 $PWD/Sort_BAM/$f$BAM $PWD/WIG/$f$WIG $GENOME
#wigToBigWig $PWD/WIG/$f$WIG $DIR/Genome.bed  $PWD/BW/$f$BW


################################ EXTEND NORM ####################################

bedtools slop -i $PWD/Aligned_BED/$f$BED -g ~/Data/Promotors/$GENOME/Genome.bed -b 100 > $PWD/Aligned_BED/$f$EXT
bedtools coverage -counts -b $PWD/Aligned_BED/$f$EXT -a ~/Data/Promotors/$GENOME/GenomeWindowed.bed > $PWD/Aligned_BED/$f$EXTWIN
sort -k1,1 -k2,2g -u  -o $PWD/Aligned_BED/$f$EXTSORT $PWD/Aligned_BED/$f$EXTWIN
paste $PWD/Aligned_BED/$f$EXTSORT $PWD/Input/Window_Full_Extended_Input.bed | awk -v OFS="\t" '{print $1,$2,$3,$4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Normalized_WIG/$f$FULLNOR
bedSort $PWD/Normalized_WIG/$f$FULLNOR $PWD/Normalized_WIG/$f$FULLNOR
bedGraphToBigWig $PWD/Normalized_WIG/$f$FULLNOR ~/Data/Promotors/$GENOME/Genome.bed $PWD/Normalized_BW/$f$EXTBW
bigWigToWig $PWD/Normalized_BW/$f$EXTBW $PWD/Normalized_WIG/$f".wig"

fi

#################################################################################


mv $PWD/Raw_Data/$f".fastq" $PWD/Analyzed
mv $PWD/Raw_Data/$f"_trimmed.fq" $PWD/Analyzed

done


################################## SUMMARIZE ####################################

awk -v OFS="\t" 'BEGIN {print "File","Raw_Reads","Trimmed","Aligned","Aln_Rate","Unique","Unq_Rate","MACS","MACS2","EPIC2"}' > Summary.txt

SAMS=$PWD/Aligned_SAM/*.sam
for f in $SAMS
do
base=`basename "$f" .sam`
samtools view -F 0x904 -c $f >> Aligned.txt
echo $base >> Files.txt
done

wc -l $PWD/Analyzed/*.fastq > RawReads.txt
wc -l $PWD/Analyzed/*_trimmed.fq > TrimReads.txt
wc -l $PWD/Aligned_BED/*Aligned.bed > Unique.txt
wc -l $PWD/MACS/P5/*_peaks.bed > MACS.txt
wc -l $PWD/MACS2/*_peaks.narrowPeak > MACS2.txt
wc -l $PWD/EPIC2/*.EPIC2.bed > EPIC2.txt


#bash ~/Scripts/ROC2.sh $MARK
#~/Scripts/Summary.sh .

paste Files.txt RawReads.txt TrimReads.txt Aligned.txt Unique.txt MACS.txt MACS2.txt EPIC2.txt | awk -v OFS="\t" '{print $1,$2/4,$4/4,$6,$6*4/$4,$7,$7/$6,$9,$11,$13}' >> Summary.txt

rm MACS.txt MACS2.txt EPIC2.txt Aligned.txt TrimReads.txt RawReads.txt Files.txt Unique.txt
#################################################################################

end=`date +%s`
runtime=$((end-start))

echo "Runtime: $runtime"
