# XQTL analysis
- lab workbook for the XQTL project

1. [Preparing the reference](#reference)
2. MAPPING
     1. [Mapping Parental strains](#mapping_parents)  
     2. [Mapping XQTL](#mapping_xqtl)
     3. [Mapping Advanced Intercross](#mapping_ai)
     4. [Mapping Dose response](#mapping_doseresponse)
     5. [Mapping Canadian Field Samples from John Gilleard](#mapping_canada_field)
     6. [Mapping US farm samples from Ray Kaplan](#mapping_us_field)
3. [Analysis](#analysis)

working directory: /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL

```shell
# working directory
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL

# project setup
mkdir 00_SCRIPTS
mkdir 01_REFERENCE
mkdir 02_RAW
mkdir 03_MAPPING
mkdir 04_VARIANTS
mkdir 05_ANALYSIS
```


## Step 1 - get and format reference <a name="reference"></a>

```shell
cd 01_REFERENCE

# get reference - note it is the same reference as in WormBase Parasite
# wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS15.genomic.fa.gz

cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa .

# make a index and a dict file
samtools faidx HAEM_V4_final.chr.fa
samtools dict HAEM_V4_final.chr.fa > HAEM_V4_final.chr.dict
```
[↥ **Back to top**](#top)





## Mapping - Parental strains <a name="mapping_parents"></a>
- reads are first trimmed using trimmomatic, followed by mapping with bwa mem

```shell
mkdir 02_RAW/RAW_PARENTS
cd 02_RAW/RAW_PARENTS

pathfind -t lane -i  26582_8 --filetype fastq -l ./ --rename

head samples_lanes.list
# MHCO3_P0_L3_n200_01	26582_8_1
# MHCO18_P0_L3_n200_IVM_01	26582_8_2

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
mkdir ../../../03_MAPPING/PARENTS

ls -1 *.paired_R1.fastq.gz | sed 's/.paired_R1.fastq.gz//g' > ../../../03_MAPPING/PARENTS/sample.list




###--- map reads

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/03_MAPPING/PARENTS
screen
while read NAME; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_PARENTS/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_PARENTS/TRIM/${NAME}.paired_R2.fastq.gz; done < sample.list &


###--- clean up
mv *out/*.marked.bam .
mv *out/*.marked.bam.bai .

rm -r *out

#--- realign indels using GATK - needs to be done given either GATK Unified Genotyper, or Popoolation2 is used for SNPs - note this does not need to be done if GATK Haplotype caller is used
ls -1 *bam > bam.list

~sd21/bash_scripts/run_gatk_indel_realigner HAEM_V4_final.chr.fa bam.list

#--- once indel realignment is completed, remove the original mapping files, keeping only the realigned bam
rm *marked.bam
rm *marked.bam.bai

#--- remove the trimmed fastqs - they are only taking up disk space
rm /lustre/scratch118/infgen/team133/sd21/hc/XQTL/02_RAW/RAW_PARENTS/TRIM/*.gz



## setup for mpileup

mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS

ls $PWD/*bam > /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS/bam.list

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS


#--- run mpileup
~sd21/bash_scripts/run_mpileup XQTL_PARENTS /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list

cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_PARENTS.mpileup &

rm *tmp*
```
[↥ **Back to top**](#top)




# Mapping - XQTL <a name="mapping_xqtl"></a>
```shell  
mkdir 02_RAW/RAW_XQTL
cd 02_RAW/RAW_XQTL


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
echo "bsub.py --threads 7 20 trim_XQTL /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/00_SCRIPTS/run_trimmomatic ${name}_${lane} ../${lane}_1.fastq.gz ../${lane}_2.fastq.gz" >> run_trim_all; \
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

while read NAME; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_ADVANCED_INTERCROSS/TRIM/${NAME}_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_ADVANCED_INTERCROSS/TRIM/${NAME}_R2.fastq.gz; done < sample.list &


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
[↥ **Back to top**](#top)








## XQTL CMH

```bash
# originally, ran the following - this included both the technical and biological replicate, which is not quite right.

perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,2-6,3-7,4-8 --input XQTL_BZ.raw.sync --output XQTL_BZ.raw.sync.cmh --min-pvalue 0.01




# rerunning without the technical replicate, ie. without 2-6 comparison
# bz
bsub.py 1 --queue long cmh_replicate_x3 "perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,3-7,4-8 --input XQTL_BZ.raw.sync --output XQTL_BZ.rep3.cmh --min-pvalue 0.01"

# lev - not only 2 replicates used
bsub.py 1 --queue long cmh_replicate_x2 "perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,4-8 --input XQTL_LEV.raw.sync --output XQTL_LEV.rep2.cmh --min-pvalue 0.01"

# ivm
bsub.py 1 --queue long cmh_replicate_x3 "perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,3-7,4-8 --input XQTL_IVM.raw.sync --output XQTL_IVM.rep3.cmh --min-pvalue 0.01"


# exploring CMH data with SNPeff results

grep "HIGH\|MODE" XQTL_IVM.raw.snpeff.vcf | awk -F'[\t|]' '{print $1,$2,$14,$9,$10,$17"_"$18}' > high_moderate.pos


R
library(dplyr)

cmh <- read.table("XQTL_IVM.raw.sync.cmh",header=F)
pos <- read.table("high_moderate.pos",header=F)
data <- inner_join(cmh,pos,by=c("V1","V2"))

# grep in R based on gene or position
filter(data,grepl("HCON_00162390", V3.y))






```







## Dose Response <a name="mapping_doseresponse"></a>
```shell
mkdir 02_RAW/RAW_DOSE_RESPONSE

cd 02_RAW/RAW_DOSE_RESPONSE


# get crams
iget 27606_3#1.cram .
iget 27606_3#2.cram .

# cram to fastq
bsub.py 2 c2fq ~sd21/bash_scripts/run_cram2fastq

ls -1 *gz | while read -r NAME; do rename s/#/_/ ${NAME}; done


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

while read NAME; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_DOSE_RESPONSE/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_DOSE_RESPONSE/TRIM/${NAME}.paired_R2.fastq.gz; done < sample.list &



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
[↥ **Back to top**](#top)





# ## Canadian Field Samples from John Gilleard <a name="mapping_canada_field"></a>
# ```
# mkdir 02_RAW/RAW_CAN_FIELD
#
# cd 02_RAW/RAW_CAN_FIELD
#
#
# # get crams
# imeta qu -z seq -d id_run = 27606 and lane = 4 | grep "cram" | grep -v phix | awk '{print $2}' >> samples.list
# imeta qu -z seq -d id_run = 27606 and lane = 5 | grep "cram" | grep -v phix | awk '{print $2}' >> samples.list
#
# while read name; do iget /seq/27606/$name ./ ; done < samples.list &
#
# rm *888.cram
# rm *#0.cram
#
#
# # cram to fastq
# bsub.py 2 c2fq ~sd21/bash_scripts/run_cram2fastq
#
#
#
# ls -1 *gz | while read -r NAME; do rename s/#/_/ ${NAME}; done
#
# #--- trim reads
# mkdir TRIM
# cd TRIM
#
# rm run_trim_all
# while read name lane; do \
# echo "bsub.py --threads 7 20 trim_XQTL /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/00_SCRIPTS/run_trimmomatic ${name}_${lane} ../${lane}_1.fastq.gz ../${lane}_2.fastq.gz" >> run_trim_all; \
# done < ../samples_lanes.list ; \
# chmod a+x run_trim_all; \
# ./run_trim_all
#
#
# rm *unpaired*
#
# # make a sample list for mapping
# mkdir ../../../03_MAPPING/CAN_FIELD
#
# ls -1 *_R1.fastq.gz | sed 's/_R1.fastq.gz//g' > ../../../03_MAPPING/CAN_FIELD/sample.list
#
# # --- mapping
# cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/03_MAPPING/CAN_FIELD/
#
#
# screen
#
# while read NAME; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_CAN_FIELD/TRIM/${NAME}_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_CAN_FIELD/TRIM/${NAME}_R2.fastq.gz; done < sample.list &
#
#
#
# #--- clean up
# mv *out/*.marked.bam .
# mv *out/*.marked.bam.bai .
#
# rm -r *out
#
# #--- realign indels using GATK - needs to be done given either GATK Unified Genotyper, or Popoolation2 is used for SNPs - note this does not need to be done if GATK Haplotype caller is used
# ls -1 *bam > bam.list
# ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa
#
# ~sd21/bash_scripts/run_gatk_indel_realigner HAEM_V4_final.chr.fa bam.list
#
#
# #--- once indel realignment is completed, remove the original mapping files, keeping only the realigned bam
# rm *marked.bam
# rm *marked.bam.bai
#
# #--- remove the trimmed fastqs - they are only taking up disk space
# rm /lustre/scratch118/infgen/team133/sd21/hc/XQTL/02_RAW/RAW_DOSE_RESPONSE/TRIM/*.gz
#
#
#
# # setup for mpileup
# mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/DOSE_RESPONSE
#
# ls $PWD/*bam > /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/DOSE_RESPONSE/bam.list
#
# cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/DOSE_RESPONSE
#
#
# #--- run mpileup
# ~sd21/bash_scripts/run_mpileup.sh XQTL_DOSE_RESPONSE /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list &
#
# cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_DOSE_RESPONSE.mpileup &
#
# rm *tmp*
# ```
# [↥ **Back to top**](#top)



# Analysis of US farm samples from Ray Kaplan <a name="mapping_us_field"></a>
```shell
mkdir 02_RAW/RAW_US_FIELD
cd 02_RAW/RAW_US_FIELD

pathfind -t lane -i 31192_6 --rename -l ./ --filetype fastq
pathfind -t lane -i 31192_7 --rename -l ./ --filetype fastq
pathfind -t lane -i 31192_8 --rename -l ./ --filetype fastq

cat samples_lanes.list
Hc_RKUSA_L3_n200_FARM_001	31192_6_1
Hc_RKUSA_L3_n200_FARM_002	31192_6_2
Hc_RKUSA_L3_n200_FARM_003	31192_6_3
Hc_RKUSA_L3_n200_FARM_004	31192_6_4
Hc_RKUSA_L3_n200_FARM_005	31192_6_5
Hc_RKUSA_L3_n200_FARM_006	31192_6_6
Hc_RKUSA_L3_n200_FARM_007	31192_6_7
Hc_RKUSA_L3_n200_FARM_008	31192_6_8
Hc_RKUSA_L3_n200_FARM_009	31192_6_9
Hc_RKUSA_L3_n200_FARM_010	31192_6_10
Hc_RKUSA_L3_n200_FARM_001	31192_7_1
Hc_RKUSA_L3_n200_FARM_002	31192_7_2
Hc_RKUSA_L3_n200_FARM_003	31192_7_3
Hc_RKUSA_L3_n200_FARM_004	31192_7_4
Hc_RKUSA_L3_n200_FARM_005	31192_7_5
Hc_RKUSA_L3_n200_FARM_006	31192_7_6
Hc_RKUSA_L3_n200_FARM_007	31192_7_7
Hc_RKUSA_L3_n200_FARM_008	31192_7_8
Hc_RKUSA_L3_n200_FARM_009	31192_7_9
Hc_RKUSA_L3_n200_FARM_010	31192_7_10
Hc_RKUSA_L3_n200_FARM_001	31192_8_1
Hc_RKUSA_L3_n200_FARM_002	31192_8_2
Hc_RKUSA_L3_n200_FARM_003	31192_8_3
Hc_RKUSA_L3_n200_FARM_004	31192_8_4
Hc_RKUSA_L3_n200_FARM_005	31192_8_5
Hc_RKUSA_L3_n200_FARM_006	31192_8_6
Hc_RKUSA_L3_n200_FARM_007	31192_8_7
Hc_RKUSA_L3_n200_FARM_008	31192_8_8
Hc_RKUSA_L3_n200_FARM_009	31192_8_9
Hc_RKUSA_L3_n200_FARM_010	31192_8_10

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
while read NAME; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_US_FIELD/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_US_FIELD/TRIM/${NAME}.paired_R2.fastq.gz; done < sample.list &


mv *out/*.marked.bam .
mv *out/*.marked.bam.bai .


# run_merge-bams_multisample
#
# Setup script to collate multiple mapped lanes with same sample ID
# Once all mapping jobs have finished, this script can be run to collate any individual mapped
# samples that shre the same name, e.g. if a sample has been sequenced twice, or, has been
# split into multiple fastq files due to high output

touch run_merger.tmp
for i in $( cat ../../02_RAW/RAW_US_FIELD/samples_lanes.list | cut -f1 | sort | uniq ); do echo "bsub.py 10 merge_${i}_tmp ./run_merge.tmp ${i} &" >> run_merger.tmp; done

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

```
[↥ **Back to top**](#top)


###---------------------------------------------------------------------------------------  


XQTL key
#Rep1.1	Rep1.2	Rep2	Rep3
#V13	V27	V39	V49

Advanced intercross key
#						Rep1	Rep2	Rep3
#control - pre v 0.5X	V17		V51		V83
#control - pre v 2X		V11		V135	V161

#IVM pre v 0.5x			V251	V267	V281
#IVM pre v 2x			V287	V297	V305





















































# run npstats

for i in *splitpileup; do bsub.py --queue long 1 npstats "/nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/npstat/npstat -n 400 -l 5000 -mincov 20 -maxcov 200 -minqual 20 -nolowfreq 2 ${i}" ; done


#--- low coverage
#for i in *splitpileup; do bsub.py --queue long 1 npstats "/nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/npstat/npstat -n 400 -l 5000 -mincov 5 -maxcov 200 -minqual 20 -nolowfreq 2 -annot ../../../GENOME/TRANSCRIPTOME/haemonchus_contortus.PRJEB506.WBPS11.annotations.gff3 ${i}" ; done


# bring them all together
for i in $(ls -1 sample_*splitpileup.stats | sed 's/_chr_.*$//g' | sort -V | uniq); do
	echo -e "chr\twindow\tlength\tlength_outgroup\tread_depth\tS\tWatterson\tPi\tTajima_D\tvar_S\tvar_Watterson\tunnorm_FayWu_H\tFayWu_H\tdiv\tnonsyn_pol\tsyn_pol\tnonsyn_div\tsyn_div\talpha" > ${i}.pileup.stats

	for j in $(ls -1 ${i}_*splitpileup.stats | grep -v "mtDNA"); do
	CHR=$( echo ${j} | sed -e 's/sample_[123456789X].*chr_//g' -e 's/.splitpileup.stats//g')
	cat ${j} | grep -v "window" | awk -v CHR=${CHR} '{print CHR,$0}' OFS="\t" >> ${i}.pileup.stats;
	done;
done

rm *splitpileup*


R-3.5.0
library(ggplot2)
library(patchwork)

sample1 <-read.table("sample_1.pileup.stats",header=T,sep="\t")
sample2 <-read.table("sample_2.pileup.stats",header=T,sep="\t")
sample3 <-read.table("sample_3.pileup.stats",header=T,sep="\t")
sample4 <-read.table("sample_4.pileup.stats",header=T,sep="\t")
sample5 <-read.table("sample_5.pileup.stats",header=T,sep="\t")
sample6 <-read.table("sample_6.pileup.stats",header=T,sep="\t")
sample7 <-read.table("sample_7.pileup.stats",header=T,sep="\t")
sample8 <-read.table("sample_8.pileup.stats",header=T,sep="\t")
sample9 <-read.table("sample_9.pileup.stats",header=T,sep="\t")
sample10 <-read.table("sample_10.pileup.stats",header=T,sep="\t")






plot_S1_pi <- ggplot(sample1)+
     geom_point(aes(window*10000,log10(sample1$Pi/sample10$Pi)),size=.1,alpha=0.5,col="blue")+
     geom_point(aes(window*10000,log10(sample1$Pi/sample8$Pi)),size=.1,alpha=0.5,col="red")+
     geom_point(aes(window*10000,log10(sample1$Pi/sample7$Pi)),size=.1,alpha=0.5,col="green")+
     geom_point(aes(window*10000,log10(sample1$Pi/sample6$Pi)),size=.1,alpha=0.5,col="orange")+
     geom_point(aes(window*10000,log10(sample1$Pi/sample5$Pi)),size=.1,alpha=0.5,col="yellow")+
     facet_grid(.~sample1$chr)+theme_bw()+
     labs(title="UGA_S vs R", x="Genomic position", y="log10(S[Pi]/R[Pi])")


plot_S2_pi <- ggplot(sample1)+
          geom_point(aes(window*10000,log10(sample9$Pi/sample10$Pi)),size=.1,alpha=0.5,col="blue")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample8$Pi)),size=.1,alpha=0.5,col="red")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample7$Pi)),size=.1,alpha=0.5,col="green")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample6$Pi)),size=.1,alpha=0.5,col="orange")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample5$Pi)),size=.1,alpha=0.5,col="yellow")+
          facet_grid(.~sample1$chr)+theme_bw()+
          labs(title="Idaho_S vs R", x="Genomic position", y="log10(S[Pi]/R[Pi])")


plot_S3_pi <- ggplot(sample1)+
     geom_point(aes(window*10000,log10(sample1$Pi/sample10$Pi)),size=1,alpha=1,col="blue")+
     geom_point(aes(window*10000,log10(sample1$Pi/sample8$Pi)),size=1,alpha=1,col="red")+
     geom_point(aes(window*10000,log10(sample1$Pi/sample7$Pi)),size=1,alpha=1,col="green")+
     geom_point(aes(window*10000,log10(sample1$Pi/sample6$Pi)),size=1,alpha=1,col="orange")+
     geom_point(aes(window*10000,log10(sample1$Pi/sample5$Pi)),size=1,alpha=1,col="yellow")+
     facet_grid(.~sample1$chr)+theme_bw()+xlim(3.65e7,3.8e7)+
     labs(title="UGA_S vs R", x="Genomic position", y="log10(S[Pi]/R[Pi])")


plot_S4_pi <- ggplot(sample1)+
          geom_point(aes(window*10000,log10(sample9$Pi/sample10$Pi)),size=1,alpha=1,col="blue")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample8$Pi)),size=1,alpha=1,col="red")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample7$Pi)),size=1,alpha=1,col="green")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample6$Pi)),size=1,alpha=1,col="orange")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample5$Pi)),size=1,alpha=1,col="yellow")+
          facet_grid(.~sample1$chr)+theme_bw()+xlim(3.65e7,3.8e7)+
          labs(title="Idaho_S vs R", x="Genomic position", y="log10(S[Pi]/R[Pi])")

plot_S1_pi + plot_S2_pi + plot_S3_pi + plot_S4_pi + plot_layout(ncol=1)





sample1.5 <- sample1[sample1$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample2.5 <- sample2[sample2$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample3.5 <- sample3[sample3$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample4.5 <- sample4[sample4$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample5.5 <- sample5[sample5$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample6.5 <- sample6[sample6$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample7.5 <- sample7[sample7$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample8.5 <- sample8[sample8$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample9.5 <- sample9[sample9$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
sample10.5  <- sample10[sample10$chr=="hcontortus_chr5_Celeg_TT_arrow_pilon",]




plot_S3_pi <- ggplot(sample1.5)+
     geom_point(aes(window*10000,log10(sample1.5$Pi/sample10.5$Pi)),size=1,alpha=1,col="blue")+
     geom_point(aes(window*10000,log10(sample1.5$Pi/sample8.5$Pi)),size=1,alpha=1,col="red")+
     geom_point(aes(window*10000,log10(sample1.5$Pi/sample7.5$Pi)),size=1,alpha=1,col="green")+
     geom_point(aes(window*10000,log10(sample1.5$Pi/sample6.5$Pi)),size=1,alpha=1,col="orange")+
     geom_point(aes(window*10000,log10(sample1.5$Pi/sample5.5$Pi)),size=1,alpha=1,col="yellow")+
     theme_bw()+xlim(3.65e7,3.8e7)+
     labs(title="UGA_S vs R", x="Genomic position", y="log10(S[Pi]/R[Pi])")




#plot labeled windows in chromosome V region
plot1.5 <-     ggplot(sample1.5,aes(window*10000,log10(Pi),label=window*10000))+
                    geom_point(size=1,alpha=1,col="black")+theme_bw()+xlim(3.65e7,3.8e7)+ylim(-5,-1)+labs(title="UGA_S")+
                    geom_text_repel(data=subset(sample1.5,Pi <= 0.001))


plot5.5 <-     ggplot(sample5.5,aes(window*10000,log10(Pi),label=window*10000))+
                    geom_point(size=1,alpha=1,col="black")+theme_bw()+xlim(3.65e7,3.8e7)+ylim(-5,-1)+labs(title="Alday_R")+
                    geom_text_repel(data=subset(sample5.5,Pi <= 0.001))

plot6.5 <-     ggplot(sample6.5,aes(window*10000,log10(Pi),label=window*10000))+
                    geom_point(size=1,alpha=1,col="black")+theme_bw()+xlim(3.65e7,3.8e7)+ylim(-5,-1)+labs(title="Mulligan_R")+
                    geom_text_repel(data=subset(sample6.5,Pi <= 0.001))

plot7.5 <-     ggplot(sample7.5,aes(window*10000,log10(Pi),label=window*10000))+
                    geom_point(size=1,alpha=1,col="black")+theme_bw()+xlim(3.65e7,3.8e7)+ylim(-5,-1)+labs(title="Rau_Shelby_R")+
                    geom_text_repel(data=subset(sample7.5,Pi <= 0.001))

plot8.5 <-     ggplot(sample8.5,aes(window*10000,log10(Pi),label=window*10000))+
                    geom_point(size=1,alpha=1,col="black")+theme_bw()+xlim(3.65e7,3.8e7)+ylim(-5,-1)+labs(title="Pena_R")+
                    geom_text_repel(data=subset(sample8.5,Pi <= 0.001))

plot10.5 <-     ggplot(sample10.5,aes(window*10000,log10(Pi),label=window*10000))+
                    geom_point(size=1,alpha=1,col="black")+theme_bw()+xlim(3.65e7,3.8e7)+ylim(-5,-1)+labs(title="Strickland_R")+
                    geom_text_repel(data=subset(sample10.5,Pi <= 0.001))


plot1.5 + plot5.5 + plot6.5 + plot7.5 + plot8.5 + plot10.5  + plot_layout(ncol=1)


plot_S4_pi <- ggplot(sample1)+
          geom_point(aes(window*10000,log10(sample9$Pi/sample10$Pi)),size=1,alpha=1,col="blue")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample8$Pi)),size=1,alpha=1,col="red")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample7$Pi)),size=1,alpha=1,col="green")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample6$Pi)),size=1,alpha=1,col="orange")+
          geom_point(aes(window*10000,log10(sample9$Pi/sample5$Pi)),size=1,alpha=1,col="yellow")+
          facet_grid(.~sample1$chr)+theme_bw()+xlim(3.65e7,3.8e7)+
          labs(title="Idaho_S vs R", x="Genomic position", y="log10(S[Pi]/R[Pi])")








#----- SNPeff analysis of functional variants
```shell
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
bbsub.py 10 build  "java -jar snpEff.jar build -v HCON_V4_20200130"




# run SNPeff

bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130  XQTL_AI.raw.vcf >  XQTL_AI.raw.snpeff.vcf"

bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130  XQTL_BZ.raw.vcf >  XQTL_BZ.raw.snpeff.vcf"

bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130 XQTL_IVM.raw.vcf > XQTL_IVM.raw.snpeff.vcf"

bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130 XQTL_LEV.raw.vcf > XQTL_LEV.raw.snpeff.vcf"
```
[↥ **Back to top**](#top)










#L3_5000
```
L3_5000 <- read.table("XQTL_L3_5k/XQTL_L3_5k.merged.fst",header=F)
data<-L3_5000

ggplot(data)+geom_point(aes(1:nrow(data)*5000,V7,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8))+
  theme_bw()+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

L3_5000_c2<-L3_5000[L3_5000$V1=="hcontortus_chr2_Celeg_TT_arrow_pilon",]
data<-L3_5000_c2
ggplot(data)+geom_point(aes(V2,V7,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.035)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous()+
  xlim(2.4e6,3.5e6)+
  theme_bw()+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

L3_5000_c5<-L3_5000[L3_5000$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon",]
data<-L3_5000_c5
ggplot(data)+geom_point(aes(V2,V7,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous()+
  xlim(35e6,40e6)+
  theme_bw()+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
```






#snpeff
bsub -q long -R "select[mem>5000] rusage[mem=5000]" -M5000 -o snpeff_new.o -e snpeff_new.e "java -Xmx4g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/COMPARATIVE_GENOMICS/snpEff/snpEff.jar -no-intergenic -no-downstream -no-upstream HCON_V4_20200130 GBX_PARENTS.raw.vcf.gz > GBX_PARENTS.raw.snpeff.vcf.gz"


```
cd ~/lustre118_link/hc/XQTL/04_VARIANTS/PARENTS

R
library(ggplot2,dplyr,ggrepel,data.table, patchwork)

parents <-read.table("XQTL_PARENTS.merged.fst",header=F)
xqtl <- read.table("../XQTL_IVM/XQTL_IVM.merged.fst",header=F)
data <- inner_join(parents,xqtl,c("V1","V2"))

data_chr5 <- data[data$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon",]

ggplot(data,aes(V7.x,V13,col=V1,label=V2))+geom_point() +
geom_text_repel(data = subset(data, V7.x > 0.4 & V13 > 0.03  ))


ggplot(data,aes(V2,V7.x*V13,label=V2))+geom_point(size=0.1)+facet_grid(V1~.)+geom_text_repel(data = subset(data, V7.x > 0.4 & V13 > 0.03  ))




xqtl_fet <- fread("xqtl_chr5.2.fet",header=F)
parents_fet <- fread("parents_chr5.fet",header=F)
data_fet <- inner_join(parents_fet,xqtl_fet,c("V1","V2"))

ggplot(data_fet,aes(V2,V7*V4.y,label=V2))+geom_point(size=0.1)+facet_grid(V1~.)



p1 <- ggplot(data_chr5,aes(V2,V7.x*V13,label=V2))+geom_point(size=0.1)+geom_text_repel(data = subset(data_chr5, V7.x > 0.4 & V13 > 0.03  ))+xlim(25e6,55e6)

p2 <- ggplot(data_fet2,aes(V2,V7*V4.y,label=V2))+geom_point(size=0.1)+geom_text_repel(data = subset(data_fet, V7*V4.y>150  ))+xlim(30e6,45e6)
p1 + p2 + plot_layout(ncol=1)















working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_BZ

# get some testdata

 cut -f1,2,13 XQTL_BZ.merged.fst > peakfinding.testdata

 cut -f1,2,13 ../XQTL_CONTROL/XQTL_CONTROL.merged.fst > peakfinding.control.testdata

1. all data - genome wide average
     - pre-treatment vs post =
          cat peakfinding.testdata | datamash mean 3
          0.015070584074789

     - time matched control pre/post =
          cat peakfinding.control.testdata | datamash mean 3
          0.011916732789999


2. list of manually curated "peaks"



hcontortus_chr1_Celeg_TT_arrow_pilon	6992500	0.09456757


3. extract peak coordinates

#!/bin/bash

PREFIX=$1
PEAKFILE=$2
FST_DATA=$3
WINDOW=$4

#eg.
#PREFIX=TEST
#PEAKFILE=peak.data
#FST_DATA=peakfinding.testdata
#WINDOW=500000


printf CHR"\t"PEAK_COORD"\t"PEAK_FST"\t"PEAK_START_COORD"\t"PEAK_END_COORD"\t"PEAK_WINDOW_SIZE"\t"PEAK_START_FST"\t"PEAK_END_FST"\n" > ${PREFIX}.peak_windows

while read COL CHR PEAK FST; do
     # ignore comment lines
     [[ "$COL" =~ ^# ]] && continue

     # extract data from master file
     cut -f1,2,${COL} ${FST_DATA} > data.tmp

     # extract a window of data around the peak
     grep "${CHR}" data.tmp | awk -v PEAK=$PEAK -v WINDOW=$WINDOW '{if($2>(PEAK-WINDOW) && $2< (PEAK+WINDOW)) print}' > data.tmp2;

     # extract from the windowed data Fst values that are half the peak height, taking into account the genome wide average
     GENOME_AVERAGE=$(cat data.tmp | datamash mean 3 sstdev 3 | awk '{print $1+(3*$2)}')

     awk -v GENOME_AVERAGE=$GENOME_AVERAGE -v FST=$FST '{if($3>((FST+GENOME_AVERAGE)/2)) print}' data.tmp2 > data.tmp3;

     # print the boundaries of the peak
     PEAK_START=$(head -n 1 data.tmp3 | cut -f2);
     PEAK_START_FST=$(head -n 1 data.tmp3 | cut -f3);
     PEAK_END=$(tail -n 1 data.tmp3 | cut -f2);
     PEAK_END_FST=$(tail -n 1 data.tmp3 | cut -f3);
     PEAK_WINDOW_SIZE=$(($PEAK_END-PEAK_START));
     printf ${CHR}"\t"${PEAK}"\t"${FST}"\t"$PEAK_START"\t"$PEAK_END"\t"${PEAK_WINDOW_SIZE}"\t"${PEAK_START_FST}"\t"${PEAK_END_FST}"\n" >> ${PREFIX}.peak_windows;
     rm *tmp*
     done < ${PEAKFILE}


XQTL
V13 1:5
V27 2:6
V39 3:7
V49 4:8

AI
#						Rep1	Rep2	Rep3
#control - pre v 0.5X	V17		V51		V83
#control - pre v 2X		V11		V135	V161

#IVM pre v 0.5x			V251	V267	V281
#IVM pre v 2x			V287	V297	V305


Dose response
V7 1:2

cut -f1,2,13 XQTL_*.merged.fst | sort -k3 | tail -n 100

#bz
hcontortus_chr1_Celeg_TT_arrow_pilon	6992500	0.09456757
hcontortus_chr1_Celeg_TT_arrow_pilon	8807500	0.08951118
hcontortus_chr1_Celeg_TT_arrow_pilon	12722500	0.07817395
#ivm-XQTL
hcontortus_chr5_Celeg_TT_arrow_pilon	34007500	0.04864623
hcontortus_chr5_Celeg_TT_arrow_pilon	36252500	0.08151205
hcontortus_chr5_Celeg_TT_arrow_pilon	37467500	0.07105946
hcontortus_chr5_Celeg_TT_arrow_pilon	45397500	0.04814783
hcontortus_chr2_Celeg_TT_arrow_pilon	2947500	0.05130916
#ivm-ai
hcontortus_chr1_Celeg_TT_arrow_pilon	11047500	0.07026732
hcontortus_chr3_Celeg_TT_arrow_pilon	38012500	0.10238732
#doseresponse
hcontortus_chr5_Celeg_TT_arrow_pilon	34397500	0.50000000




# XQTL_IVM
# Rep1.1 - V13
13   hcontortus_chr5_Celeg_TT_arrow_pilon	37467500	0.07105946
13   hcontortus_chr2_Celeg_TT_arrow_pilon	2947500	0.05130916
#Rep1.2 - V27
27   hcontortus_chr5_Celeg_TT_arrow_pilon	33842500	0.06678443
27   hcontortus_chr5_Celeg_TT_arrow_pilon	46717500	0.04836741
#Rep2 - V39
#Rep3 - V49
49   hcontortus_chr5_Celeg_TT_arrow_pilon	36317500	0.10849596
49   hcontortus_chr5_Celeg_TT_arrow_pilon	35822500	0.10612988

./run_find_peak_windows.sh XQTL_IVM peak.data XQTL_IVM.merged.fst 500000


# XQTL_LEV  
# Rep1.1 - V13
13   hcontortus_chr4_Celeg_TT_arrow_pilon	14817500	0.07934698
13   hcontortus_chr5_Celeg_TT_arrow_pilon	31467500	0.10024437
#Rep1.2 - V27
27   hcontortus_chr5_Celeg_TT_arrow_pilon	22157500	0.06093688
27   hcontortus_chr5_Celeg_TT_arrow_pilon	28482500	0.08273666
#Rep2 - V39
#Rep3 - V49



# XQTL_BZ
# Rep1.1 - V13
13   hcontortus_chr1_Celeg_TT_arrow_pilon	6992500	0.09456757
13   hcontortus_chr1_Celeg_TT_arrow_pilon	8807500	0.08951118
13   hcontortus_chr1_Celeg_TT_arrow_pilon	12722500	0.07817395
# Rep1.2 - V27
# Rep2 - V39
# Rep3 - V49
49	hcontortus_chr1_Celeg_TT_arrow_pilon	8337500	0.23776054




# AI
# Rep1 - V287	 
287  hcontortus_chr1_Celeg_TT_arrow_pilon	11047500	0.07026732
# Rep 2 - V297
297  hcontortus_chr2_Celeg_TT_arrow_pilon	7147500	0.21312042
297  hcontortus_chr4_Celeg_TT_arrow_pilon	29077500	0.06586891
# Rep 3 - V305
305  hcontortus_chr3_Celeg_TT_arrow_pilon	38012500	0.10238732


# DOSE_RESPONSE
# rep 1 - V7
7    hcontortus_chr5_Celeg_TT_arrow_pilon	34052500	0.35451348


# XQTL_ADULTS
# rep 1 - V7
7    hcontortus_chr5_Celeg_TT_arrow_pilon	37377500	0.19438392
7    hcontortus_chr2_Celeg_TT_arrow_pilon	27157500	0.12806700




COL=7
DATA=XQTL_PARENTS.merged.fst

CUTOFF=$(cut -f1,2,${COL} *.merged.fst | sort -k3 | datamash mean 3 sstdev 3 | awk '{print $1+(3*$2)}')

cut -f1,2,${COL} *.merged.fst | sort -k3 | awk -v CUTOFF=$CUTOFF -v COL=$COL '{if($3>CUTOFF) print COL,$0}' OFS="\t"  > all.high.peaks

./run_find_peak_windows.sh XQTL_allhightest all.high.peaks ${DATA} 500000

sort -k1,1 -k2,2n  XQTL_allhightest.peak_windows | awk '{if($6!=0 && NF==8) print}' > sorted


R ${DATA} ${COL}
args <- commandArgs()

library(ggplot2)
library(viridis)

name=args[2]
columns=args[3]

a<-read.table("sorted",header=T)
ggplot(a)+
     geom_rect(aes(xmin=a$PEAK_START_COORD,ymin=1:nrow(a)-0.5,xmax=a$PEAK_END_COORD,ymax=1:nrow(a)+0.5,fill=PEAK_FST))+
     facet_grid(a$CHR~.)+
     xlim(0,50e6)+
     labs(title=paste0("file=",name,", ","column=",columns))+
     scale_fill_viridis(direction=-1,limits=c(0,0.5))+
     theme_bw()

ggsave(paste0("predictedpeaks_",name,"_","column_",columns,".pdf"))


******
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
