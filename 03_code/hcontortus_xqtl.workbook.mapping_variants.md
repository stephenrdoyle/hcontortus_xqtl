# XQTL analysis - mapping_variants

Author: Stephen Doyle
Contact: stephen.doyle[at]sanger.ac.uk

- this workbook contains the workflow to map raw sequencing reads and perform variant calling

1. [Preparing the reference](#reference)
2. Parents workflow
     - Get and trim the raw reads before mapping (parent samples)
     - map reads to the reference genome (parent samples)
     - Use GATK to realign mapped sequencing reads around indels (parent samples)
     - Run mpileup in preparation for variant calling (parent samples)
     - Perform SNP calling to generate a VCF file
     - Running NPSTATS
     - Running snpEff
3. XQTL workflow
4. Advanced Intercross workflow
5. Dose response workflow
6. US farm workflow
7. Other
- Canadian Field Samples from John Gilleard

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
     echo -e \"bwa mem -t 4 -R '@RG\\\\\\\\\\\tRG:$"{sample_name}"\\\\\\\\\\\tID:$"{sample_name}"\\\\\\\\\\\tSM:$"{sample_name}"' -Y -M ref.fa $"{i}" $"{i/R1/R2}" | samtools view --threads 4 -b - | samtools sort --threads 4 -o $"{i/R1/bwamem}".tmp.sort.bam - \"  > step2.2_bwamem_tmp_$"{n}";
     done;
     chmod a+x step2.2_bwamem_tmp_*

touch step2_FINISHED" > step2.1_bwasplitter

chmod a+x step2.1_bwasplitter



# merge mapping, mark duplicates, generate stats, and finalise

echo -e "# merge mapping, mark duplicates, generate stats, and finalise

sample_name=$"{1}"

ls -1 *.tmp.sort.bam > bam.fofn

samtools merge --threads 4 -cpf -b bam.fofn tmp.merged.sorted.bam

#rm *.tmp.sort.bam

java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=tmp.merged.sorted.bam

OUTPUT=merged.sorted.marked.bam METRICS_FILE=tmp.merged.sorted.marked.metrics TMP_DIR=$PWD/tmp

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
     -o ${SAMPLE%.bam}.realigned.bam" >> run_gatk_IndelRealigner;
done < ${BAM_LIST}


chmod a+x run_gatk_IndelRealigner
./run_gatk_IndelRealigner
```


## Run mpileup in preparation for variant calling
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

while read number sequences; do
     grep -A1 "$sequences" ref.tmp2 > $sequences.tmp.fasta;
done < ref_seq.tmp.list


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



## Perform SNP calling to generate a VCF file
- in addition to running poolseq based analyses with popoolation2 and NPSTATS, I want to generate a standard VCF file that can be used with VCFTOOLS and snpEff.

```bash
~sd21/bash_scripts/run_mpileup2vcf.sh XQTL_PARENTS HAEM_V4_final.chr.fa bam.list
```
- where "run_mpileup2vcf.sh" is:

```bash
#!/bin/bash
# run_mpileup2vcf
#
# Tool will run an mpileup on all bams in a specified directory, and split the job up by sequence in the reference to parallelise the job. Finally, it shoudl merge all the split jobs into a single file.
# mpileup command based on https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1024-y

# load modules
module load \
samtools/1.6--h244ad75_4 \
fastaq/3.17.0-docker3 \
bcftools/1.9--h68d8f2e_9



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

while read number sequences; do
     grep -A1 "$sequences" ref.tmp2 > $sequences.tmp.fasta;
done < ref_seq.tmp.list


#while read name; do
#	if [ ! -f ${name}.bai ]; then
#	 samtools index -b ${name}
#	fi; done < ${bamlist}

while read number sequences; do
     echo -e "samtools mpileup --ignore-RG -Ou -t DP,SP,AD,ADF,INFO/AD -b $bamlist -r $sequences -f ref.tmp -F0.25 -d500 -E | bcftools call -vm -Oz -o $number.$sequences.tmp.vcf.gz" > run_mpileup2vcf.tmp.$number;
done < ref_seq.tmp.list

chmod a+x run_mpileup2vcf.tmp*


jobs=$( wc -l ref_seq.tmp.list | cut -f1 -d" " )

bsub -q long -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -J mpileup2vcf_array[1-$jobs] -e mpileup2vcf_array[1-$jobs].e -o mpileup2vcf_array[1-$jobs].o ./run_mpileup2vcf.tmp.\$LSB_JOBINDEX



#--- Step 2: bring it together

echo -e "ls -1v *.tmp.vcf.gz > vcffiles.fofn && vcf-concat -f vcffiles.fofn | gzip -c >  ${prefix}.raw.vcf.gz" > run_mp2vcf_combine

chmod a+x run_mp2vcf_combine

bsub -q normal -w "done(mpileup2vcf_array)"  -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' -n1 -M1000 -J mp2vcf_combine -e mp2vcf_combine.e -o mp2vcf_combine.o ./run_mp2vcf_combine ${prefix}

```






## Running NPSTATS
- NPSTATS is a handy tool to generate population genetic metrics from pooled sequnecing data
- Github: https://github.com/lucaferretti/npstat
- Ref: https://doi.org/10.1111/mec.12522

- to use, first need to convert the mpileup into individual pileup files, both by sample and then by chromomsome

```bash
# convert the multisample mpileup into individual pileup files per sample
run_mpileup2pileup XQTL_*.mpileup 3

# once completed, split pileup by sample,  by chromosome
for i in `(grep ">" HAEM_V4_final.chr.fa | sed -e 's/>//g')`; do
     for j in *.final.pileup; do
          grep $i $j > sample_${j%.final.pileup}_chr_${i}.splitpileup;
     done ;
done
```

where "run_mpileup2pileup.sh" is:
```bash
#!/bin/bash

#(($# != 2)) || { echo -e "\nUsage: $0 <file to split> <# columns in each split>\n\n"; exit; }

infile="$1"
columns_per_sample="$2"

cut -f 1-3 $infile > mpileup_positions.tmp
cut -f 4- $infile > mpileup_data.tmp


inc=${columns_per_sample}
ncol=$(awk 'NR==1{print NF}' mpileup_data.tmp)

((inc < ncol)) || { echo -e "\nSplit size >= number of columns\n\n"; exit; }


for((i=1, start=1, end=$inc; i < ncol/inc + 1; i++, start+=inc, end+=inc)); do
     cut -f$start-$end mpileup_data.tmp > "pileup.tmp.$i"
done

for i in pileup.tmp.*; do
     paste mpileup_positions.tmp $i > ${i#pileup.tmp.}.final.pileup;
done

#rm *tmp*

```
- once the data is prepared, run NPSTATS
```bash
for i in *splitpileup; do
     bsub.py --queue long 10 npstats "/nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/npstat/npstat -n 400 -l 10000 -mincov 20 -maxcov 200 -minqual 20 -nolowfreq 2 ${i}" ;
done
```

## Running popoolation2
- popoolation2 is the main tools used for the pariwise analysis of pooled sequencign data
- Code: https://sourceforge.net/p/popoolation2/wiki/Main/
- Paper: https://doi.org/10.1093/bioinformatics/btr589


- to run popoolation2, the following was performed.
```bash
~sd21/bash_scripts/run_mpileup2popoolation2 XQTL_PARENTS ~sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa XQTL_PARENTS.mpileup 200 5000
```
- where :
```bash
PREFIX=$1
FASTA=$2
MPILEUP=$3
POOL_SIZE=$4
WINDOW=$5

ID="U$(date +%s)"

grep ">" ${FASTA} | sed -e 's/>//g' | cat -n > ref.fa.sequence-names.tmp


#--- make sync
echo -e "\
java -jar ~sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/mpileup2sync.jar \
--input ${MPILEUP} \
--output ${PREFIX}.raw.sync \
--min-qual 20 \
--threads 7 \
--fastq-type sanger

while read NUMBER NAME; do
     grep "\${NAME}" ${PREFIX}.raw.sync > \${NUMBER}.\${NAME}.raw.sync.tmp
done < ref.fa.sequence-names.tmp" > run_make_syncronised_file.tmp.${ID}

chmod a+x run_make_syncronised_file.tmp.${ID}


#---- make fst and fet arrays
while read NUMBER NAME; do
     echo -e "\
     perl ~sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/fisher-test.pl --input ${NUMBER}.${NAME}.raw.sync.tmp --output ${NUMBER}.${NAME}.fet.tmp --min-count 4 --min-coverage 30 --max-coverage 2% --suppress-noninformative

     perl ~sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/fst-sliding.pl --pool-size ${POOL_SIZE} --window-size ${WINDOW} --step-size ${WINDOW} --min-count 4 --min-coverage 30 --max-coverage 2% --input ${NUMBER}.${NAME}.raw.sync.tmp --output ${NUMBER}.${NAME}.raw.fst.tmp" > run_pp2_split.tmp.${ID}.${NUMBER};
done < ref.fa.sequence-names.tmp
chmod a+x *run_pp2_split*

echo -e "\
PREFIX=\$1
cat \$(find ./ -name \"*.fet.tmp\" | sort -V) | sed -e \"s/=/\\\t/g\" > \${PREFIX}.merged.fet
cat \$(find ./ -name \"*.fst.tmp\" | sort -V) | sed -e \"s/=/\\\t/g\" > \${PREFIX}.merged.fst
#rm *.tmp*" > run_pp2_combine.tmp.${ID}
chmod a+x run_pp2_combine.tmp.${ID}




jobs=$( wc -l ref.fa.sequence-names.tmp | cut -f1 -d" " )

bsub -q yesterday -R'span[hosts=1] select[mem>10000] rusage[mem=10000]' -M10000 -n 7 -e 01_mp2pp2_makesync.e -o 01_mp2pp2_makesync.o -J 01_mp2pp2_makesync ./run_make_syncronised_file.tmp.${ID}

bsub -q long -w "done(01_mp2pp2_makesync)" -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' -M1000 -J 02_mp2pp2_array[1-$jobs] -e 02_pp2_array[1-$jobs].e -o 02_pp2_array[1-$jobs].o ./run_pp2_split.tmp.${ID}.\$LSB_JOBINDEX

bsub -q normal -w "done(02_mp2pp2_array[1-$jobs])" -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' -n1 -M1000 -J 03_mp2pp2_combine -e 03_pp2_combine.e -o 03_pp2_combine.o ./run_pp2_combine.tmp.${ID} ${PREFIX}

```

## Running snpEff
- snpEff allows functional annotation of a VCF, to enable interpretation fo the effect of a variant on protein coding sequences.
-


- to use snpEff, a database needs to be first made. Note that this only needs to be done once.

```bash
#----- SNPeff analysis of functional variants
# setup for HCON_V4

cd /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff

mkdir data/HCON_V4_20200130

cd  data/HCON_V4_20200130

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa HCON_V4_Dec2019.fa

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/20191212/UPDATED_annotation.gff3 genes.gff

gffread genes.gff -g HCON_V4_20200130.fa -y protein.fa

cp HCON_V4_20200130.fa ../genomes/

# modify config file
cd /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff

echo "

# Haemonchus contortus chromosomes V4
HCON_V4_20200130.genome : HCON_V4_20200130

" > new.genome

cat snpEff.config new.genome > tmp; mv tmp snpEff.config


# build database
bsub.py 10 build  "java -jar snpEff.jar build -v HCON_V4_20200130"

```
- to run snpEff, the following was used. Note that "-no-intergenic -no-downstream -no-upstream" flags are used
```bash
bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130 XQTL_PARENTS.raw.vcf > XQTL_PARENTS.raw.snpeff.vcf"
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

while read NAME; do
     ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_XQTL/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_XQTL/TRIM/${NAME}.paired_R2.fastq.gz;
done < sample.list &


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

cd 02_RAW/RAW_ADVANCED_INTERCROSS

pathfind -t lane -i  27112_2 --filetype fastq -l ./ --rename
pathfind -t lane -i  27112_3 --filetype fastq -l ./ --rename
pathfind -t lane -i  27112_5 --filetype fastq -l ./ --rename
pathfind -t lane -i  27112_6 --filetype fastq -l ./ --rename
pathfind -t lane -i  27112_4 --filetype fastq -l ./ --rename
pathfind -t lane -i  26581_8 --filetype fastq -l ./ --rename

#--- make a sample database called samples_lanes.list
/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_ADVANCED_INTERCROSS/samples_lanes.list

#--- example
head samples_lanes.list
# XQTL_F4_L3_n200_IVM_post_2x_01	26581_8_10
# XQTL_F4_L3_n200_IVM_post_2x_02	26581_8_11
# XQTL_F4_L3_n200_control_pre_01	26581_8_1
# XQTL_F4_L3_n200_IVM_post_2x_03	26581_8_12
# XQTL_F4_L3_n200_control_post_0.5x_01	26581_8_13
# XQTL_F4_L3_n200_control_post_0.5x_02	26581_8_14
# XQTL_F4_L3_n200_control_post_2x_01	26581_8_15
# XQTL_F4_L3_n200_control_post_2x_02	26581_8_16
# XQTL_F4_L3_n200_control_post_2x_03	26581_8_17
# XQTL_F4_L3_n200_control_post_0.5x_03	26581_8_18

#--- trim reads
mkdir TRIM
cd TRIM

rm run_trim_all
while read name lane; do \
     echo "bsub.py --threads 7 20 trim_AI /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/00_SCRIPTS/run_trimmomatic ${name}_${lane} ../${lane}_1.fastq.gz ../${lane}_2.fastq.gz" >> run_trim_all; \
done < ../samples_lanes.list ; \

chmod a+x run_trim_all; \

./run_trim_all



# merge trimmed reads into R1 and R2 files with sensible names
echo -e '"cat ../samples_lanes.list | cut -f1 | sort | uniq | while read -r name; do zcat ${name}_*paired_R1.fastq.gz | gzip > ${name}_R1.fastq.gz; zcat ${name}_*paired_R2.fastq.gz | gzip > ${name}_R2.fastq.gz; done"' > run_cat_reads

chmod a+x run_cat_reads

bsub.py --queue yesterday 1 cat_reads ./run_cat_reads


#--- once reads are concatenated, make a sample list.
#--- note slight tweak to remove the "paired" named samples, to make sure only picking up the merged read samples

mkdir /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/03_MAPPING/ADVANCED_INTERCROSS

ls -1 *R1.fastq.gz | grep -v "paired" | sed 's/_R1.fastq.gz//g' > ../../../03_MAPPING/ADVANCED_INTERCROSS/sample.list


# --- mapping
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/03_MAPPING/ADVANCED_INTERCROSS

screen

while read NAME; do
     ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_ADVANCED_INTERCROSS/TRIM/${NAME}_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_ADVANCED_INTERCROSS/TRIM/${NAME}_R2.fastq.gz;
done < sample.list &


#--- clean up
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
rm /lustre/scratch118/infgen/team133/sd21/hc/XQTL/02_RAW/RAW_ADVANCED_INTERCROSS/TRIM/*.gz



# setup for mpileup
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/ADVANCED_INTERCROSS

ls $PWD/*bam > /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/ADVANCED_INTERCROSS/bam.list

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/ADVANCED_INTERCROSS


#--- run mpileup
~sd21/bash_scripts/run_mpileup.sh XQTL_ADVANCED_INTERCROSS /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list

cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_ADVANCED_INTERCROSS.mpileup &

rm *tmp*
```


```shell
# originally, ran the following - this included both the technical and biological replicate, which is not quite right.

perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,2-6,3-7,4-8 --input XQTL_BZ.raw.sync --output XQTL_BZ.raw.sync.cmh --min-pvalue 0.01
```


```shell
# rerunning without the technical replicate, ie. without 2-6 comparison
# bz
bsub.py 1 --queue long cmh_replicate_x3 "perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,3-7,4-8 --input XQTL_BZ.raw.sync --output XQTL_BZ.rep3.cmh --min-pvalue 0.01"

# lev - note only 2 replicates used
bsub.py 1 --queue long cmh_replicate_x2 "perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,4-8 --input XQTL_LEV.raw.sync --output XQTL_LEV.rep2.cmh --min-pvalue 0.01"

# ivm
bsub.py 1 --queue long cmh_replicate_x3 "perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,3-7,4-8 --input XQTL_IVM.raw.sync --output XQTL_IVM.rep3.cmh --min-pvalue 0.01"
```


```shell
# exploring CMH data with SNPeff results
grep "HIGH\|MODE" XQTL_IVM.raw.snpeff.vcf | awk -F'[\t|]' '{print $1,$2,$14,$9,$10,$17"_"$18}' > high_moderate.pos
```

```R
library(dplyr)

cmh <- read.table("XQTL_IVM.raw.sync.cmh",header=F)
pos <- read.table("high_moderate.pos",header=F)
data <- inner_join(cmh,pos,by=c("V1","V2"))

# grep in R based on gene or position
filter(data,grepl("HCON_00162390", V3.y))
```







## Dose Response workflow
```shell
mkdir 02_RAW/RAW_DOSE_RESPONSE

cd 02_RAW/RAW_DOSE_RESPONSE


# get crams
iget 27606_3#1.cram .
iget 27606_3#2.cram .

# cram to fastq
bsub.py 2 c2fq ~sd21/bash_scripts/run_cram2fastq

ls -1 *gz | while read -r NAME; do
     rename s/#/_/ ${NAME};
done


#--- trim reads
mkdir TRIM
cd TRIM

rm run_trim_all
while read name lane; do \
     echo "bsub.py --threads 7 20 trim_XQTL /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/00_SCRIPTS/run_trimmomatic ${name}_${lane} ../${lane}_1.fastq.gz ../${lane}_2.fastq.gz" >> run_trim_all; \
done < ../samples_lanes.list ; \

chmod a+x run_trim_all; \

./run_trim_all


rm *unpaired*

# make a sample list for mapping
mkdir ../../../03_MAPPING/DOSE_RESPONSE

ls -1 *.paired_R1.fastq.gz | sed 's/.paired_R1.fastq.gz//g' > ../../../03_MAPPING/DOSE_RESPONSE/sample.list

# --- mapping
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/03_MAPPING/DOSE_RESPONSE/


screen

while read NAME; do
     ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_DOSE_RESPONSE/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_DOSE_RESPONSE/TRIM/${NAME}.paired_R2.fastq.gz;
done < sample.list &



#--- clean up
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
rm /lustre/scratch118/infgen/team133/sd21/hc/XQTL/02_RAW/RAW_DOSE_RESPONSE/TRIM/*.gz



# setup for mpileup
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/DOSE_RESPONSE

ls $PWD/*bam > /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/DOSE_RESPONSE/bam.list

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/DOSE_RESPONSE


#--- run mpileup
~sd21/bash_scripts/run_mpileup.sh XQTL_DOSE_RESPONSE /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list &

cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_DOSE_RESPONSE.mpileup &

rm *tmp*


#mpileup to popoolation2
./run_mpileup2popoolation2_lowcov ../../01_REFERENCE/HAEM_V4_final.chr.fa XQTL_DOSE_RESPONSE.mpileup 200 5000 &

```



## Analysis of US farm samples from Ray Kaplan
```bash
mkdir 02_RAW/RAW_US_FIELD
cd 02_RAW/RAW_US_FIELD

pathfind -t lane -i 31192_6 --rename -l ./ --filetype fastq
pathfind -t lane -i 31192_7 --rename -l ./ --filetype fastq
pathfind -t lane -i 31192_8 --rename -l ./ --filetype fastq

cat samples_lanes.list
# Hc_RKUSA_L3_n200_FARM_001	31192_6_1
# Hc_RKUSA_L3_n200_FARM_002	31192_6_2
# Hc_RKUSA_L3_n200_FARM_003	31192_6_3
# Hc_RKUSA_L3_n200_FARM_004	31192_6_4
# Hc_RKUSA_L3_n200_FARM_005	31192_6_5
# Hc_RKUSA_L3_n200_FARM_006	31192_6_6
# Hc_RKUSA_L3_n200_FARM_007	31192_6_7
# Hc_RKUSA_L3_n200_FARM_008	31192_6_8
# Hc_RKUSA_L3_n200_FARM_009	31192_6_9
# Hc_RKUSA_L3_n200_FARM_010	31192_6_10
# Hc_RKUSA_L3_n200_FARM_001	31192_7_1
# Hc_RKUSA_L3_n200_FARM_002	31192_7_2
# Hc_RKUSA_L3_n200_FARM_003	31192_7_3
# Hc_RKUSA_L3_n200_FARM_004	31192_7_4
# Hc_RKUSA_L3_n200_FARM_005	31192_7_5
# Hc_RKUSA_L3_n200_FARM_006	31192_7_6
# Hc_RKUSA_L3_n200_FARM_007	31192_7_7
# Hc_RKUSA_L3_n200_FARM_008	31192_7_8
# Hc_RKUSA_L3_n200_FARM_009	31192_7_9
# Hc_RKUSA_L3_n200_FARM_010	31192_7_10
# Hc_RKUSA_L3_n200_FARM_001	31192_8_1
# Hc_RKUSA_L3_n200_FARM_002	31192_8_2
# Hc_RKUSA_L3_n200_FARM_003	31192_8_3
# Hc_RKUSA_L3_n200_FARM_004	31192_8_4
# Hc_RKUSA_L3_n200_FARM_005	31192_8_5
# Hc_RKUSA_L3_n200_FARM_006	31192_8_6
# Hc_RKUSA_L3_n200_FARM_007	31192_8_7
# Hc_RKUSA_L3_n200_FARM_008	31192_8_8
# Hc_RKUSA_L3_n200_FARM_009	31192_8_9
# Hc_RKUSA_L3_n200_FARM_010	31192_8_10

mkdir TRIM
cd TRIM

rm run_trim_all
while read name lane; do \
     echo "bsub.py --queue yesterday --threads 7 20 trim_P /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/00_SCRIPTS/run_trimmomatic ${name}_${lane} ../${lane}_1.fastq.gz ../${lane}_2.fastq.gz" >> run_trim_all; \
done < ../samples_lanes.list ; \

chmod a+x run_trim_all; \

./run_trim_all



# get rid of the unpaired fraction
rm *unpaired*

# make a sample list for mapping
mkdir ../../../03_MAPPING/US_FIELD

ls -1 *.paired_R1.fastq.gz | sed 's/.paired_R1.fastq.gz//g' > ../../../03_MAPPING/US_FIELD/sample.list




#--- map reads

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/03_MAPPING/US_FIELD

screen

while read NAME; do
     ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_US_FIELD/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_US_FIELD/TRIM/${NAME}.paired_R2.fastq.gz;
done < sample.list &


mv *out/*.marked.bam .
mv *out/*.marked.bam.bai .


# run_merge-bams_multisample
#
# Setup script to collate multiple mapped lanes with same sample ID
# Once all mapping jobs have finished, this script can be run to collate any individual mapped
# samples that shre the same name, e.g. if a sample has been sequenced twice, or, has been
# split into multiple fastq files due to high output

touch run_merger.tmp
for i in $( cat ../../02_RAW/RAW_US_FIELD/samples_lanes.list | cut -f1 | sort | uniq ); do
     echo "bsub.py 10 merge_${i}_tmp ./run_merge.tmp ${i} &" >> run_merger.tmp;
done

echo 'SAMPLE=$1; ls -1 ${SAMPLE}*.bam > ${SAMPLE}.tmp.bamlist ; samtools-1.3 merge -c -b ${SAMPLE}.tmp.bamlist ${SAMPLE}.merge.bam; samtools-1.3 index -b ${SAMPLE}.merge.bam' > run_merge.tmp

chmod a+x run_merge.tmp
chmod a+x run_merger.tmp
./run_merger.tmp



#--- clean up

rm *bamlist
rm *tmp*
rm *marked.bam*



#rm -r *out

#--- realign indels using GATK - needs to be done given either GATK Unified Genotyper, or Popoolation2 is used for SNPs - note this does not need to be done if GATK Haplotype caller is used
ls -1 *bam > bam.list

~sd21/bash_scripts/run_gatk_indel_realigner /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list

#--- once indel realignment is completed, remove the original mapping files, keeping only the realigned bam
rm *marked.bam
rm *marked.bam.bai
rm *.merge.bam
rm *.merge.bam.bai


#--- remove the trimmed fastqs - they are only taking up disk space
rm /lustre/scratch118/infgen/team133/sd21/hc/XQTL/02_RAW/RAW_US_FIELD/TRIM/*.gz



# setup for mpileup
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/US_FIELD

ls $PWD/*bam > /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/US_FIELD/bam.list

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/US_FIELD


#--- run mpileup
~sd21/bash_scripts/run_mpileup XQTL_US_FIELD /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list

rm *tmp*

~sd21/bash_scripts/run_mpileup2popoolation2 XQTL_US_FIELD ~sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa $PWD/XQTL_US_FIELD.mpileup 200 10000
















#XQTL key
#Rep1.1	Rep1.2	Rep2	Rep3
#V13	V27	V39	V49

#Advanced intercross key
#Rep1	Rep2	Rep3
#control - pre v 0.5X	V17		V51		V83
#control - pre v 2X		V11		V135	V161
#IVM pre v 0.5x			V251	V267	V281
#IVM pre v 2x			V287	V297	V305

# run npstats
for i in *splitpileup; do
     bsub.py --queue long 1 npstats "/nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/npstat/npstat -n 400 -l 5000 -mincov 20 -maxcov 200 -minqual 20 -nolowfreq 2 ${i}" ;
done

#--- low coverage
#for i in *splitpileup; do bsub.py --queue long 1 npstats "/nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/npstat/npstat -n 400 -l 5000 -mincov 5 -maxcov 200 -minqual 20 -nolowfreq 2 -annot ../../../GENOME/TRANSCRIPTOME/haemonchus_contortus.PRJEB506.WBPS11.annotations.gff3 ${i}" ; done

for i in $(ls -1 sample_*splitpileup.stats | sed 's/*chr*.*$//g' | sort -V | uniq); do
     echo -e "chr\twindow\tlength\tlength_outgroup\tread_depth\tS\tWatterson\tPi\tTajima_D\tvar_S\tvar_Watterson\tunnorm_FayWu_H\tFayWu_H\tdiv\tnonsyn_pol\tsyn_pol\tnonsyn_div\tsyn_div\talpha" > ${i}.pileup.stats
     for j in $(ls -1 ${i}_*splitpileup.stats | grep -v "mtDNA"); do
          CHR=$( echo ${j} | sed -e 's/sample_[123456789X].*chr_//g' -e 's/.splitpileup.stats//g')
          cat ${j} | grep -v "window" | awk -v CHR=${CHR} '{print CHR,$0}' OFS="\t" >> ${i}.pileup.stats;
          done;
done
```










# run SNPeff
```bash
bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/ls -lev_pretreatment -no-intergenic -no-downstream -no-upstream HCON_V4_20200130  XQTL_AI.raw.vcf >  XQTL_AI.raw.snpeff.vcf"

bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130  XQTL_BZ.raw.vcf >  XQTL_BZ.raw.snpeff.vcf"

bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130 XQTL_IVM.raw.vcf > XQTL_IVM.raw.snpeff.vcf"

bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130 XQTL_LEV.raw.vcf > XQTL_LEV.raw.snpeff.vcf"

bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130 GBX_PARENTS.raw.vcf.gz > GBX_PARENTS.raw.snpeff.vcf.gz"
```
