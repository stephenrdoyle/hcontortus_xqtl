# XQTL analysis - mapping_variants
- this workbook contains the workflow to map raw sequencing reads and perform variant calling

1. [Preparing the reference](#reference)
2. MAPPING
     1. [Mapping Parental strains](#mapping_parents)  
     2. [Mapping XQTL](#mapping_xqtl)
     3. [Mapping Advanced Intercross](#mapping_ai)
     4. [Mapping Dose response](#mapping_doseresponse)
     5. [Mapping Canadian Field Samples from John Gilleard](#mapping_canada_field)
     6. [Mapping US farm samples from Ray Kaplan](#mapping_us_field)
3. [Analysis](#analysis)




## Project setup
```shell
# working directory
WORKING=/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL

cd ${WORKING}

# project setup
mkdir 00_SCRIPTS
mkdir 01_REFERENCE
mkdir 02_RAW
mkdir 03_MAPPING
mkdir 04_VARIANTS
mkdir 05_ANALYSIS
```


## Preparing the reference <a name="reference"></a>
- reference genome is the latest chromosome-scale assembly for *H. contortus*, available at [WormBse Parasite](https://parasite.wormbase.org/Haemonchus_contortus_prjeb506/Info/Index/) and is described in [Doyle et al 2020 Communications Biology](https://doi.org/10.1038/s42003-020-01377-3).
```shell
cd ${WORKING}/01_REFERENCE

# get reference - note it is the same reference as in WormBase Parasite
# wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS15.genomic.fa.gz

cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa .

# make a index and a dict file - these are used in various steps later on
samtools faidx HAEM_V4_final.chr.fa
samtools dict HAEM_V4_final.chr.fa > HAEM_V4_final.chr.dict
```
[↥ **Back to top**](#top)





## Get and trim the raw reads before mapping
- reads are first trimmed using trimmomatic
- given there is a lot of data, and that pooled sequencing has been performed and therefore downstream analyses perhaps more likely to be affected by poorer quality reads/mapping, trimming is worth doing.

```shell
# make some working spaces
mkdir ${WORKING}/02_RAW/RAW_PARENTS
cd ${WORKING}/02_RAW/RAW_PARENTS

# get the raw reads - note "pathfind" is a Sanger command to retrive raw sequencign data from Sanger archives. The "26582_8" is the lane ID for the sequencing run
pathfind -t lane -i  26582_8 --filetype fastq -l ./ --rename

# make a "samples_lanes.list" tab-delimited file containing the name that will be used for the sample throguhout, and the lane ID. This gets used in the trimming and mapping steps. The file will look like the following, without the hash: 
# MHCO3_P0_L3_n200_01	26582_8_1
# MHCO18_P0_L3_n200_IVM_01	26582_8_2

mkdir TRIM
cd TRIM

rm run_trim_all
while read name lane; do \
echo "bsub.py --queue yesterday --threads 7 20 trim_P ${WORKING}/00_SCRIPTS/run_trimmomatic.sh ${name}_${lane} ../${lane}_1.fastq.gz ../${lane}_2.fastq.gz" >> run_trim_all; \
done < ../samples_lanes.list ; \
chmod a+x run_trim_all; \
./run_trim_all

# get rid of the unpaired fraction
rm *unpaired*

# make a sample list for mapping
mkdir ../../../03_MAPPING/PARENTS

ls -1 *.paired_R1.fastq.gz | sed 's/.paired_R1.fastq.gz//g' > ../../../03_MAPPING/PARENTS/sample.list

```

- where "run_trimmomatic.sh" is:
```shell
#------ Trimmomatic

# tool to trim adapters and poor quality sequence from fastq files

# Requirements: paired end fastq files
# Manual: # http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

# Use: ./run_trimmomatic <out_prefix> <read_1> <read_2>

#sd21@sanger.ac.uk
#V1: April 2016



###########################################################################################################################
# input read1 and read2 files to be processed
out_prefix=$1
read_1=$2
read_2=$3


java -jar /software/pathogen/external/apps/usr/local/Trimmomatic-0.32/trimmomatic-0.32.jar PE \
-threads 7 \
-phred33 \
${read_1} ${read_2} \
${out_prefix}.paired_R1.fastq.gz ${out_prefix}.unpaired_R1.fastq.gz \
${out_prefix}.paired_R2.fastq.gz ${out_prefix}.unpaired_R2.fastq.gz \
ILLUMINACLIP:/nfs/users/nfs_s/sd21/databases/trimmomatic_Illumina-adapters.fa:2:30:10 \
SLIDINGWINDOW:10:20 MINLEN:50
```




###--- map reads to the reference genome
- uses a script called "run_bwamem_splitter.sh" for mapping
- this is a workflow to split fastq reads into chunks, map the chunks separately and in parallel, and then merge everything back together again
- it reads the sample list, and will iteratively work through the list of samples to map each one.
- becasue it checks to see whether a particular job has been started (and if it has, it will skip to the next one), the command can be run multiple times to further speed up the mapping if there are lots of samples to process

```shell
cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/03_MAPPING/PARENTS

screen
while read NAME; do
     ~sd21/bash_scripts/run_bwamem_splitter.sh ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_PARENTS/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_PARENTS/TRIM/${NAME}.paired_R2.fastq.gz;
     done < sample.list &

# once mapping has been completed (dont do this before or while running!), do some clean up to remove unnecessary files

mv *out/*.marked.bam .
mv *out/*.marked.bam.bai .

rm -r *out

```
- where "run_bwamem_splitter.sh" is:
```shell
#!/usr/bin/env bash
########################################################################################
# run_bwamem_splitter.sh
#
# sd21@sanger.ac.uk
# July 2018
# Updated 20200214 - necessary updates to work on farm5, update samtools to 1.6
########################################################################################

sample_name=$1
reference=$2
read1=$3
read2=$4

ID="U$(date +%s)"

if [ "$#" -eq 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	echo ""
    echo "Usage: ~sd21/bash_scripts/run_bwa_splitter <SAMPLE_PREFIX> <REFERENCE> <R1.fastq> <R2.fastq>"
    echo ""
    exit 0
fi

if [ -d "${sample_name}_bwasplitter_out" ]; then
	echo -e "\nThere is already a run started with this sample name. Rename and start again\n"
    exit 0
fi


mkdir ${sample_name}_bwasplitter_out
cd ${sample_name}_bwasplitter_out
mkdir logfiles

# prepare reference and split data
echo -e "# prepare reference and split the raw data
sample_name=$"{1}"
reference=$"{2}"
read1=$"{3}"
read2=$"{4}"
ln -sf $"{reference}" ref.fa
bwa index -b 100000000 ref.fa

if [[ $read1 =~ \.gz$ ]]
then ln -sf $"{read1}" R1.fq.gz; zcat R1.fq.gz | split -d -a 3 -l 4000000 - R1_tmp_
else ln -sf $"{read1}" R1.fq; split -d -a 3 -l 4000000 R1.fq R1_tmp_
fi

if [[ $read2 =~ \.gz$ ]]
then ln -sf $"{read2}" R2.fq.gz; zcat R2.fq.gz | split -d -a 3 -l 4000000 - R2_tmp_
else ln -sf $"{read2}" R2.fq; split -d -a 3 -l 4000000 R2.fq R2_tmp_
fi


#split -d -a 3 -l 4000000 R1.fq R1_tmp_
#split -d -a 3 -l 4000000 R2.fq R2_tmp_
touch step1_FINISHED" > step1_bwasplitter
chmod a+x step1_bwasplitter


# prepare split mapping and set off mapping
echo -e "# prepare the split mapping run files
sample_name=$"{1}"
n=0
for i in \` ls -1 R1_* \` ; do
let \"n+=1\"
echo -e \"bwa mem -t 4 -R '@RG\\\\\\\\\\\tRG:$"{sample_name}"\\\\\\\\\\\tID:$"{sample_name}"\\\\\\\\\\\tSM:$"{sample_name}"' -Y -M ref.fa $"{i}" $"{i/R1/R2}" | samtools view --threads 4 -b - | samtools sort --threads 4 -o $"{i/R1/bwamem}".tmp.sort.bam - \"  > step2.2_bwamem_tmp_$"{n}"; done; chmod a+x step2.2_bwamem_tmp_*
touch step2_FINISHED" > step2.1_bwasplitter

chmod a+x step2.1_bwasplitter



# merge mapping, mark duplicates, generate stats, and finalise

echo -e "# merge mapping, mark duplicates, generate stats, and finalise
sample_name=$"{1}"
ls -1 *.tmp.sort.bam > bam.fofn
samtools merge --threads 4 -cpf -b bam.fofn tmp.merged.sorted.bam
#rm *.tmp.sort.bam
java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=tmp.merged.sorted.bam OUTPUT=merged.sorted.marked.bam METRICS_FILE=tmp.merged.sorted.marked.metrics TMP_DIR=$PWD/tmp
samtools flagstat merged.sorted.marked.bam > $"{sample_name}".merged.sorted.marked.flagstat
#samtools-1.3 stats merged.sorted.marked.bam | grep ^SN | cut -f 2- > $"{sample_name}".merged.sorted.marked.stats
bamtools stats -in merged.sorted.marked.bam > $"{sample_name}".merged.sorted.marked.bamstats
samtools view --threads 4 -F 12 -b merged.sorted.marked.bam -o $"{sample_name}".merged.sorted.marked.bam
samtools view --threads 4 -f 12 merged.sorted.marked.bam -o $"{sample_name}".unmapped.bam
samtools index -b $"{sample_name}".merged.sorted.marked.bam
rm R[12]_*
rm -r *tmp*
mv *.[eo] logfiles/
touch step3_FINISHED
touch bam_splitter_COMPLETE" > step3_bwasplitter
chmod a+x step3_bwasplitter



#----- RELEASE THE KRAKEN!

# run - reference and paired read setup
bsub -q normal -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -n1 -M10000 -J step1_bwasplitter_${ID} -e step1_bwasplitter.e -o step1_bwasplitter.o ./step1_bwasplitter ${sample_name} ${reference} ${read1} ${read2}


# run - prepare mapping scripts
bsub -q normal -w "done(step1_bwasplitter_${ID})" -R'span[hosts=1] select[mem>2500] rusage[mem=2500]' -n1 -M2500 -J step2.1_bwasplitter_${ID} -e step2.1_bwasplitter.e -o step2.1_bwasplitter.o ./step2.1_bwasplitter ${sample_name}

while [ ! -f step1_FINISHED ]
do
  sleep 2
done

jobs=$( ls -1 R1_tmp_* | wc -l )


# run - mapping in array
bsub -q normal -w "done(step2.1_bwasplitter_${ID})"  -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -n6 -M10000 -J step2.2_bwasplitter_${ID}_[1-$jobs] -e step2.2_bwasplitter[1-$jobs].e -o step2.2_bwasplitter[1-$jobs].o ./step2.2_bwamem_tmp_\$LSB_JOBINDEX


# run - merge mapping data into a single bam, mark duplicates, and clean up
bsub -q normal -w "done(step2.2_bwasplitter_${ID}_[1-$jobs])"  -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n6 -M20000 -J step3_bwasplitter_${ID} -e logfiles/step3_bwasplitter.e -o logfiles/step3_bwasplitter.o ./step3_bwasplitter ${sample_name}

```




### Use GATK to realign mapped sequencing reads around indels
- needs to be done given either GATK Unified Genotyper, or Popoolation2 is used for SNPs - note this does not need to be done if GATK Haplotype caller is used
- in this project, both popoolation2 and GATK UG will be used, so worth doing
```
# make a bam list of samples
ls -1 *bam > bam.list

# run indel realigner
~sd21/bash_scripts/run_gatk_indel_realigner.sh HAEM_V4_final.chr.fa bam.list

#--- once indel realignment is completed, remove the original mapping files, keeping only the realigned bam
rm *marked.bam
rm *marked.bam.bai

#--- remove the trimmed fastqs - they are only taking up disk space
rm /lustre/scratch118/infgen/team133/sd21/hc/XQTL/02_RAW/RAW_PARENTS/TRIM/*.gz
```
- where "run_gatk_indel_realigner.sh" is:
```shell
REFERENCE=$1
BAM_LIST=$2

module load gatk/3.7.0

samtools faidx ${REFERENCE}
samtools dict  ${REFERENCE} >  ${REFERENCE%.fa}.dict

echo -e "/software/pathogen/etc/gatk/3.7.0/wrappers/gatk -T RealignerTargetCreator --num_threads 7 -R ${REFERENCE} \\" > run_gatk_RealignerTargetCreator
while read SAMPLE; do
echo -e "--input_file ${SAMPLE} \\" >> run_gatk_RealignerTargetCreator;
done < ${BAM_LIST}
echo -e "-o all_samples.intervals" >> run_gatk_RealignerTargetCreator

chmod a+x run_gatk_RealignerTargetCreator
bsub -q long -n 15 -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J step_01_gatk_indel_realigner -e step_01_gatk_indel_realigner.e -o step_01_gatk_indel_realigner.o ./run_gatk_RealignerTargetCreator

while read SAMPLE; do
echo -e "bsub -q long -w "step_01_gatk_indel_realigner" -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J step_02_gatk_indel_realigner -e step_02_gatk_indel_realigner.e -o  step_02_gatk_indel_realigner.o /software/pathogen/etc/gatk/3.7.0/wrappers/gatk -T IndelRealigner \
-R ${REFERENCE} \\
-I ${SAMPLE} \\
-targetIntervals all_samples.intervals \\
-o ${SAMPLE%.bam}.realigned.bam" >> run_gatk_IndelRealigner; done < ${BAM_LIST}


chmod a+x run_gatk_IndelRealigner
./run_gatk_IndelRealigner
```


## setup for mpileup
```shell
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS

ls $PWD/*bam > /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS/bam.list

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS


#--- run mpileup
~sd21/bash_scripts/run_mpileup.sh XQTL_PARENTS /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list

cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_PARENTS.mpileup &

rm *tmp*
```
- where "run_mpileup.sh" is:
```shell
#!/bin/bash
# run_mpileup
#
# Tool will run an mpileup on all bams in a specified directory, and split the job up by sequence in the reference to parallelise the job. Finally, it shoudl merge all the split jobs into a single file.
# mpileup command based on https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1024-y

(($# == 3)) || { echo -e "\nUsage: $0 <prefix> <reference sequence> <bamlist_file>\n\n"; exit; }

prefix=$1
ref=$2
bamlist=$3

#--- Step 1: prepare input files
cp $ref ref.tmp
samtools faidx ref.tmp
fastaq to_fasta -l 0 ref.tmp ref.tmp2
samtools faidx ref.tmp2
grep ">" ref.tmp | cut -f1  -d" " | sed -e 's/>//g' | cat -n > ref_seq.tmp.list

while read number sequences; do grep -A1 "$sequences" ref.tmp2 > $sequences.tmp.fasta; done < ref_seq.tmp.list


#while read name; do
#	if [ ! -f ${name}.bai ]; then
#	 samtools index -b ${name}
#	fi; done < ${bamlist}

while read number sequences; do
echo -e "samtools mpileup -b $bamlist -r $sequences -f ref.tmp -t DP,SP,AD,ADF,INFO/AD -F0.25 -d500 -E -o $number.$sequences.tmp.mpileup" > run_mpileup.tmp.$number;
done < ref_seq.tmp.list
chmod a+x run_mpileup.tmp*

#while read number sequences; do \
#echo -e "\
#samtools-1.3 mpileup -b $bamlist -r $sequences -f ref.tmp -t DP,SP,AD,ADF,INFO/AD -F0.25 -d500 -E -o $number.$sequences.tmp.mpileup" > run_mp_$sequences.$number; done < ref_seq.tmp.list
jobs=$( wc -l ref_seq.tmp.list | cut -f1 -d" " )
bsub -q long -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J mpileup_array[1-$jobs] -e mpileup_array[1-$jobs].e -o mpileup_array[1-$jobs].o ./run_mpileup.tmp.\$LSB_JOBINDEX



#--- Step 2: run mpileup on split chromosomes/sequences
#chmod a+x run_mpileup.tmp
#./run_mpileup.tmp


#--- Step 3: check to see all mpileup jobs are complete. If so, them merge into a single mpileup file
#a=$( ls -1 *.FINISHED | wc -l )
#b=$( cat ref_seq.tmp.list | wc -l )
#while [[ $a != $b ]]
#	do
#	sleep 100
#	a=$( ls -1 *.FINISHED | wc -l )
#	b=$( cat ref_seq.tmp.list | wc -l )
#	done
#cat $(find ./ -name "*.mpileup" | sort -V) > ${prefix}.mpileup; rm *tmp*



echo -e "prefix=\$1; cat \$(find ./ -name \"*.mpileup\" | sort -V) > \${prefix}.mpileup; touch mpileup_FINISHED" > run_mp_combine
chmod a+x run_mp_combine

bsub -q normal -w "done(mpileup_array)"  -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' -n1 -M1000 -J mp_combine -e mp_combine.e -o mp_combine.o ./run_mp_combine ${prefix}

# if final mpileup throws strange error in popoolation2 for example, can be fixed using the following
#sed 's/\t\t/\t!\t!/g' final.mpileup > final.mpileup2
```


[↥ **Back to top**](#top)




# Mapping - XQTL <a name="mapping_xqtl"></a>
```shell
WORKING=/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL

mkdir ${WORKING}/02_RAW/RAW_XQTL
cd ${WORKING}/02_RAW/RAW_XQTL

# pathfind is a Sanger command to find and download sequencing lanes from Sanger storage
pathfind -t lane -i 23241_1 --filetype fastq -l ./ --rename
pathfind -t lane -i 23241_2 --filetype fastq -l ./ --rename
pathfind -t lane -i 23241_3 --filetype fastq -l ./ --rename
pathfind -t lane -i 23241_5 --filetype fastq -l ./ --rename
pathfind -t lane -i 23241_4 --filetype fastq -l ./ --rename
pathfind -t lane -i 23241_6 --filetype fastq -l ./ --rename
pathfind -t lane -i 23241_8 --filetype fastq -l ./ --rename
pathfind -t lane -i 23241_7 --filetype fastq -l ./ --rename
pathfind -t lane -i 23204_7 --filetype fastq -l ./ --rename
pathfind -t lane -i 23204_8 --filetype fastq -l ./ --rename
pathfind -t lane -i 21529_1 --filetype fastq -l ./ --rename
pathfind -t lane -i 21528_3 --filetype fastq -l ./ --rename
pathfind -t lane -i 21528_2 --filetype fastq -l ./ --rename
pathfind -t lane -i 21528_1 --filetype fastq -l ./ --rename
pathfind -t lane -i 21395_3 --filetype fastq -l ./ --rename
pathfind -t lane -i 21395_2 --filetype fastq -l ./ --rename
pathfind -t lane -i 21395_1 --filetype fastq -l ./ --rename

#---- samples_lanes.list file
/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_XQTL/samples_lanes.list


#--- trim reads
mkdir TRIM
cd TRIM

rm run_trim_all
while read name lane; do \
     echo "bsub.py --threads 7 20 trim_XQTL ${WORKING}/00_SCRIPTS/run_trimmomatic.sh ${name}_${lane} ../${lane}_1.fastq.gz ../${lane}_2.fastq.gz" >> run_trim_all; \
     done < ../samples_lanes.list ; \
chmod a+x run_trim_all; \
./run_trim_all


rm *unpaired*

# make a sample list for mapping
mkdir ../../../03_MAPPING/XQTL

ls -1 *.paired_R1.fastq.gz | sed 's/.paired_R1.fastq.gz//g' > ../../../03_MAPPING/XQTL/sample.list

screen
while read NAME; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_XQTL/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_XQTL/TRIM/${NAME}.paired_R2.fastq.gz; done < sample.list &


#--- cleanup
mv *out/*.marked.bam .
mv *out/*.marked.bam.bai .

rm -r *out

#--- realign indels using GATK - needs to be done given either GATK Unified Genotyper, or Popoolation2 is used for SNPs - note this does not need to be done if GATK Haplotype caller is used
ls -1 *bam > bam.list
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa

~sd21/bash_scripts/run_gatk_indel_realigner HAEM_V4_final.chr.fa bam.list

#--- once indel realignment is completed, remove the original mapping files, keeping only the realigned bam
rm *marked.bam
rm *marked.bam.bai


#--- remove the trimmed fastqs - they are only taking up disk space
rm /lustre/scratch118/infgen/team133/sd21/hc/XQTL/02_RAW/RAW_XQTL/TRIM/*.gz

# setup for mpileup
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/XQTL_ADULT
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/XQTL_BZ
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/XQTL_IVM
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/XQTL_L3_5k
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/XQTL_LEV
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/XQTL_CONTROL

ls $PWD/*adultM*bam > ../../04_VARIANTS/XQTL_ADULT/bam.list
ls $PWD/*control*bam > ../../04_VARIANTS/XQTL_CONTROL/bam.list
ls $PWD/*L3_n5k*bam > ../../04_VARIANTS/XQTL_L3_5k/bam.list
ls $PWD/*n200_IVM*bam > ../../04_VARIANTS/XQTL_IVM/bam.list
ls $PWD/*LEV*bam > ../../04_VARIANTS/XQTL_LEV/bam.list
ls $PWD/*BZ*bam > ../../04_VARIANTS/XQTL_BZ/bam.list

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/XQTL


#--- run mpileup
~sd21/bash_scripts/run_mpileup.sh XQTL_ADULT /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list

cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_PARENTS.mpileup &

rm *tmp*



cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_ADULTS.mpileup &
```
[↥ **Back to top**](#top)






## Advanced Intercross <a name="mapping_ai"></a>

```shell
mkdir 02_RAW/RAW_ADVANCED_INTERCROSS

cd 02_RAW/R
