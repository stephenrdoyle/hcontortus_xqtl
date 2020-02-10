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
```


## Step 1 - get and format reference <a name="reference"></a>

```shell
cd 01_REFERENCE

# get reference - note it is the same reference as in WBP13
cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa .

# make a index and a dict file
samtools faidx HAEM_V4_final.chr.fa
samtools dict HAEM_V4_final.chr.fa > HAEM_V4_final.chr.dict
```
[↥ **Back to top**](#top)





## Mapping - Parental strains <a name="mapping_parents"></a>

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
```



###--- map reads
```
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

```
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



## Dose Response <a name="mapping_doseresponse"></a>

```
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








## Canadian Field Samples from John Gilleard <a name="mapping_canada_field"></a>


```
mkdir 02_RAW/RAW_CAN_FIELD

cd 02_RAW/RAW_CAN_FIELD


# get crams
imeta qu -z seq -d id_run = 27606 and lane = 4 | grep "cram" | grep -v phix | awk '{print $2}' >> samples.list
imeta qu -z seq -d id_run = 27606 and lane = 5 | grep "cram" | grep -v phix | awk '{print $2}' >> samples.list

while read name; do iget /seq/27606/$name ./ ; done < samples.list &

rm *888.cram
rm *#0.cram


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
mkdir ../../../03_MAPPING/CAN_FIELD

ls -1 *_R1.fastq.gz | sed 's/_R1.fastq.gz//g' > ../../../03_MAPPING/CAN_FIELD/sample.list

# --- mapping
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/03_MAPPING/CAN_FIELD/


screen

while read NAME; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_CAN_FIELD/TRIM/${NAME}_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_CAN_FIELD/TRIM/${NAME}_R2.fastq.gz; done < sample.list &



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
```
[↥ **Back to top**](#top)


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













# Analyses <a name="analysis"></a>

```bash
# load R - using 3.6.0
R
```

```R
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(scales)
library(patchwork)

pal<-palette(c("cornflowerblue","blue"))


###---------------------------------------------------------------------------------------        
# Parents


parents<-read.table("PARENTS/XQTL_PARENTS.merged.fst",header=F)
data<-parents


# genome wide levels of significance
#--- mean of replicates Fst
gwide_sig3_fst_p0 <- mean(c(parents$V7))+(3*sd(c(parents$V7)))
gwide_sig5_fst_p0 <- mean(c(parents$V7))+(5*sd(c(parents$V7)))


ggplot(data)+geom_point(aes(1:nrow(data)*5000,V7,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept = gwide_sig3_fst_p0, linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept = gwide_sig5_fst_p0, linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



###---------------------------------------------------------------------------------------        
#L3_5000

L3_5000 <- read.table("XQTL_L3_5k/XQTL_L3_5k.merged.fst",header=F)
data<-L3_5000

ggplot(data)+geom_point(aes(1:nrow(data)*5000,V7,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
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
  scale_x_continuous(,labels=comma)+
  #xlim(45000000,55000000)+
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
  scale_x_continuous(,labels=comma)+
  #xlim(45000000,55000000)+
  theme_bw()+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


###---------------------------------------------------------------------------------------        

# XQTL

xqtl_control <- read.table("XQTL_CONTROL/XQTL_CONTROL.merged.fst",header=F)
data<-xqtl_control

# genome wide levels of significance
#--- mean of replicates Fst
gwide_sig3_fst <- mean(c(xqtl_control$V21,xqtl_control$V13,xqtl_control$V29))+(3*sd(c(xqtl_control$V21,xqtl_control$V13,xqtl_control$V29)))
gwide_sig5_fst <- mean(c(xqtl_control$V21,xqtl_control$V13,xqtl_control$V29))+(5*sd(c(xqtl_control$V21,xqtl_control$V13,xqtl_control$V29)))


#Rep1	Rep2	Rep3
#V11	V21	V29

control_r1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V11,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_r2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V21,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_r3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V29,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_r1 + control_r2 + control_r3 + plot_layout(ncol=1)

###---------------------------------------------------------------------------------------        
#BZ

xqtl_bz <- read.table("XQTL_BZ/XQTL_BZ.merged.fst",header=F)
data<-xqtl_bz

#Rep1.1	Rep1.2	Rep2	Rep3
#V13	V27	V39	V49

bz_r1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V13,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

bz_r2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V27,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

bz_r3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V39,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

bz_r4 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V49,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

bz_r1 + bz_r2 + bz_r3 + bz_r4 +plot_layout(ncol=1)


###---------------------------------------------------------------------------------------        
# XQTL_LEV


xqtl_lev <- read.table("XQTL_LEV/XQTL_LEV.merged.fst",header=F)
data<-xqtl_lev

#Rep1.1	Rep1.2	Rep2	Rep3
#V13	V27	V39	V49

lev_r1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V13,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

lev_r2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V27,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

lev_r3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V39,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.3)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

lev_r4 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V49,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.3)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

lev_r1 + lev_r2 + lev_r3 + lev_r4 +plot_layout(ncol=1)


###---------------------------------------------------------------------------------------        
# XQTL_IVM

xqtl_ivm <- read.table("XQTL_IVM/XQTL_IVM.merged.fst",header=F)
data<-xqtl_ivm

#Rep1.1	Rep1.2	Rep2	Rep3
#V13	V27	V39	V49

ivm_r1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V13,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.075)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_r2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V27,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.075)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_r3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V39,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.075)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_r4 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V49,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.075)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_r1 + ivm_r2 + ivm_r3 + ivm_r4 + plot_layout(ncol=1)


###---------------------------------------------------------------------------------------        
# ADVANCED_INTERCROSS

xqtl_ai <- read.table("ADVANCED_INTERCROSS/XQTL_ADVANCED_INTERCROSS.merged.fst",header=F)
data<-xqtl_ai

# genome wide levels of significance
#--- mean of replicates Fst
gwide_sig3_fst_0.5X <- mean(c(xqtl_ai$V17,xqtl_ai$V51,xqtl_ai$V83))+(3*sd(c(xqtl_ai$V21,xqtl_ai$V51,xqtl_ai$V83)))
gwide_sig5_fst_0.5X <- mean(c(xqtl_ai$V21,xqtl_ai$V51,xqtl_ai$V83))+(5*sd(c(xqtl_ai$V21,xqtl_ai$V51,xqtl_ai$V83)))

gwide_sig3_fst_2X <- mean(c(xqtl_ai$V11,xqtl_ai$V135,xqtl_ai$V161))+(3*sd(c(xqtl_ai$V11,xqtl_ai$V135,xqtl_ai$V161)))
gwide_sig5_fst_2X <- mean(c(xqtl_ai$V11,xqtl_ai$V135,xqtl_ai$V161))+(5*sd(c(xqtl_ai$V11,xqtl_ai$V135,xqtl_ai$V161)))

#						Rep1	Rep2	Rep3
#control - pre v 0.5X	V17		V51		V83
#control - pre v 2X		V11		V135	V161

#IVM pre v 0.5x			V251	V267	V281
#IVM pre v 2x			V287	V297	V305


#rep1
control_0.5x.1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V17,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_0.5X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_0.5X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_2X.1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V11,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_2X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_2X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_0.5X.1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V251,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_0.5X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_0.5X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_2X.1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V287,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_2X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_2X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_0.5x.1 + ivm_0.5X.1 + control_2X.1 + ivm_2X.1 + plot_layout(ncol=2)


#rep2
control_0.5x.2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V51,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_0.5X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_0.5X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_2X.2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V135,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_2X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_2X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_0.5X.2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V267,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_0.5X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_0.5X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_2X.2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V297,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_2X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_2X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_0.5x.2 + ivm_0.5X.2 + control_2X.2 + ivm_2X.2 + plot_layout(ncol=2)


#rep2
control_0.5x.3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V83,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_0.5X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_0.5X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_2X.3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V161,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_2X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_2X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_0.5X.3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V281,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_0.5X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_0.5X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ivm_2X.3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V305,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.1)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  geom_hline(yintercept=gwide_sig3_fst_2X,linetype = "dashed", colour="orange",size=0.5)+
  geom_hline(yintercept=gwide_sig5_fst_2X,linetype = "dashed", colour="red",size=0.5)+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

control_0.5x.3 + ivm_0.5X.3 + control_2X.3 + ivm_2X.3 + plot_layout(ncol=2)




# summary of XQTL and AI endpoints
ivm_r1 + ivm_2X.1 + ivm_2X.2 + ivm_2X.3 + plot_layout(ncol=1)

###---------------------------------------------------------------------------------------        
#XQTL F2 adults

xqtl_adultM <- read.table("XQTL_ADULT/XQTL_ADULTS.merged.fst",header=F)
data<-xqtl_adultM

adultM <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V7,alpha=V4,colour = ifelse(as.numeric(V1) %% 2 ==1, "0", "1")),size=0.1)+
  ylim(0,0.2)+
  xlab("Relative window position in genome")+
  ylab("Fst")+
  scale_color_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,3e8,0.5e8),labels=comma)+
  theme_bw()+
  theme(legend.position="none",
        panel.background=element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
adultM








# XQTL_AI
#IVM pre v 2x
xqtl_ai <- read.table("ADVANCED_INTERCROSS/XQTL_ADVANCED_INTERCROSS.merged.fst",header=F)
pairwise_ai_R1vR2 <-	ggplot(xqtl_ai,aes(V287,V297,col=V1,label=V2))+
 						geom_point(alpha=0.2,size=0.1)+
 						geom_text_repel(data = subset(xqtl_ai, V287 > quantile(xqtl_ai$V287,0.999) & V297 > quantile(xqtl_ai$V297,0.999) ))+
 						theme_bw()+ theme(legend.position="none")

pairwise_ai_R1vR3 <-	ggplot(xqtl_ai,aes(V287,V305,col=V1,label=V2))+
 						geom_point(alpha=0.2,size=0.1)+
 						geom_text_repel(data = subset(xqtl_ai, V287 > quantile(xqtl_ai$V287,0.999) & V305 > quantile(xqtl_ai$V297,0.999)  ))+
 						theme_bw()+ theme(legend.position="none")

pairwise_ai_R1vR2 + pairwise_ai_R1vR3 + plot_layout(ncol=2)



# XQTL_BZ
#Rep1.1	Rep1.2	Rep2	Rep3
#V13	V27	V39	V49

xqtl_bz <- read.table("XQTL_BZ/XQTL_BZ.merged.fst",header=F)
pairwise_bz_R1vR2 <-	ggplot(xqtl_bz,aes(V13,V39,col=V1,label=V2))+
 						geom_point(alpha=0.2,size=0.1)+
 						geom_text_repel(data = subset(xqtl_bz, V13 > quantile(xqtl_bz$V13,0.999) & V39 > quantile(xqtl_bz$V39,0.999) ))+
 						theme_bw()+ theme(legend.position="none")

pairwise_bz_R1vR3 <-	ggplot(xqtl_bz,aes(V13,V49,col=V1,label=V2))+
 						geom_point(alpha=0.2,size=0.1)+
 						geom_text_repel(data = subset(xqtl_bz, V13 > quantile(xqtl_bz$V13,0.999) & V49 > quantile(xqtl_bz$V49,0.999) ))+
 						theme_bw()+ theme(legend.position="none")
pairwise_bz_R1vR2 + pairwise_bz_R1vR3 + plot_layout(ncol=2)






# XQTL_LEV

#Rep1.1	Rep1.2	Rep2	Rep3
#V13	V27	V39	V49
xqtl_lev <- read.table("XQTL_LEV/XQTL_LEV.merged.fst",header=F)
pairwise_lev_R1vR2 <-	ggplot(xqtl_lev,aes(V13,V39,col=V1,label=V2))+
 						geom_point(alpha=0.2,size=0.1)+
 						geom_text_repel(data = subset(xqtl_lev, V13 > quantile(xqtl_lev$V13,0.995) & V39 > quantile(xqtl_lev$V39,0.995) ))+
 						theme_bw()+ theme(legend.position="none")

pairwise_lev_R1vR3 <-	ggplot(xqtl_lev,aes(V13,V49,col=V1,label=V2),stroke = NA)+
 						geom_point(alpha=0.2,size=0.1)+
 						geom_text_repel(data = subset(xqtl_lev, V13 > quantile(xqtl_lev$V13,0.995) & V49 > quantile(xqtl_lev$V49,0.995) ))+
 						theme_bw()+ theme(legend.position="none")
pairwise_lev_R1vR2 + pairwise_lev_R1vR3 + plot_layout(ncol=2)




# XQTL_IVM

#Rep1.1	Rep1.2	Rep2	Rep3
#V13	V27	V39	V49
xqtl_ivm <- read.table("XQTL_IVM/XQTL_IVM.merged.fst",header=F)
pairwise_ivm_R1vR2 <-	ggplot(xqtl_ivm,aes(V13,V39,col=V1,label=V2))+
 						geom_point(alpha=0.2,size=0.1)+
 						geom_text_repel(
 						 & V39 > quantile(xqtl_ivm$V39,0.995) ))+
 						theme_bw()+ theme(legend.position="none")

pairwise_ivm_R1vR3 <-	ggplot(xqtl_ivm,aes(V13,V49,col=V1,label=V2),stroke = NA)+
 						geom_point(alpha=0.2,size=0.1)+
 						geom_text_repel(data = subset(xqtl_ivm, V13 > quantile(xqtl_ivm$V13,0.995) & V49 > quantile(xqtl_ivm$V49,0.995) ))+
 						theme_bw()+ theme(legend.position="none")
pairwise_ivm_R1vR2 + pairwise_ivm_R1vR3 + plot_layout(ncol=2)






###########################################################################################
# allele frequencies for Jamie


# BZ
cd ~/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_BZ

awk '{if($1=="hcontortus_chr1_Celeg_TT_arrow_pilon" && $2>=5000000 && $2<=10000000) print $0}' OFS="\t" XQTL_BZ_pwc > XQTL_BZ_pwc_chr1_5-10Mb
awk '{if($1=="hcontortus_chr1_Celeg_TT_arrow_pilon" && $2>=5000000 && $2<=10000000) print $0}' OFS="\t" XQTL_BZ.merged.fst > XQTL_BZ.merged_25-35Mb.fst
R-3.4.0
library(ggplot2)
library(patchwork)


a<-read.table("XQTL_BZ_pwc_chr1_5-10Mb",header=F,na.strings="na")
a_fst<-read.table("XQTL_BZ.merged_25-35Mb.fst",header=F,na.strings="na")
b<-a[a$V12>0.4 & a$V19>0.4 & a$V25>0.4 & a$V30>0.4,]
b<-na.omit(b)
c<-a[a$V2==7029790,]

ggplot()+geom_point(aes(a$V2,a$V12),alpha=0.2,cex=0.1)+geom_point(aes(b$V2,b$V12),col="red")+theme_bw()+geom_vline(xintercept=7029790,col="blue")+ylim(0,1)

plot<-ggplot()+geom_vline(xintercept=7029790,col="blue")+geom_point(aes(a$V2,a$V12),alpha=0.2,cex=0.1)+geom_point(aes(b$V2,b$V12),col="red")+theme_bw()+ylim(0,1)+ylab("Allele freq diff. - pre vs post treatment")
fst_plot<-ggplot()+geom_vline(xintercept=7029790,col="blue")+geom_point(aes(a_fst$V2,a_fst$V13),alpha=0.5)+theme_bw()+ylab("Fst - pre vs post treatment")

fst_plot + plot + plot_layout(ncol=1)


write.table(b,"chr1.allelefreqdiff.gt0.4.txt",sep="\t",quote=F,row.names=F)


while read chr pos rest; do \
	awk -v chr="$chr" -v pos="$pos" -F'[\t/]' '{if($1==chr && $2==pos && $4=="2") print $1,$2,$3,$5"/"$6,$9,$11/$12,$13/$14,$15/$16,$17/$18,$19/$20,$21/$22,$23/$24,$25/$26}' OFS="\t" XQTL_BZ_rc >> XQTL_BZ_peakSNPs.txt; \
done < chr1.allelefreqdiff.gt0.4.txt &




# LEV - chr 5
cd ~/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV

awk '{if($1=="hcontortus_chr5_Celeg_TT_arrow_pilon" && $2>=25000000 && $2<=35000000) print $0}' OFS="\t" XQTL_LEV_pwc > XQTL_BZ_pwc_chr5_25-35Mb
awk '{if($1=="hcontortus_chr5_Celeg_TT_arrow_pilon" && $2>=25000000 && $2<=35000000) print $0}' OFS="\t" XQTL_LEV.merged.fst > XQTL_LEV.merged_25-35Mb.fst

R-3.4.0
library(ggplot2)
library(patchwork)
a<-read.table("XQTL_BZ_pwc_chr5_25-35Mb",header=F,na.strings="na")
a_fst<-read.table("XQTL_LEV.merged_25-35Mb.fst",header=F,na.strings="na")
b<-a[a$V12>0.4 & a$V19>0.4 & a$V30>0.4,]
b<-na.omit(b)



plot<-ggplot()+geom_point(aes(a$V2,a$V12),alpha=0.2,cex=0.1)+geom_point(aes(b$V2,b$V12),col="red")+theme_bw()+ylim(0,1)+ylab("Allele freq diff. - pre vs post treatment")
fst_plot<-ggplot()+geom_point(aes(a_fst$V2,a_fst$V13),alpha=0.5)+theme_bw()+ylab("Fst - pre vs post treatment")

fst_plot + plot + plot_layout(ncol=1)



# LEV - chr 4
cd ~/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV

awk '{if($1=="hcontortus_chr4_Celeg_TT_arrow_pilon" && $2>=10000000 && $2<=20000000) print $0}' OFS="\t" XQTL_LEV_pwc > XQTL_LEV_pwc_chr4_10-20Mb
awk '{if($1=="hcontortus_chr4_Celeg_TT_arrow_pilon" && $2>=10000000 && $2<=20000000) print $0}' OFS="\t" XQTL_LEV.merged.fst > XQTL_LEV.merged_10-20Mb.fst

R-3.4.0
library(ggplot2)
library(patchwork)
a<-read.table("XQTL_LEV_pwc_chr4_10-20Mb",header=F,na.strings="na")
a_fst<-read.table("XQTL_LEV.merged_10-20Mb.fst",header=F,na.strings="na")
b<-a[a$V12>0.4 & a$V19>0.4 & a$V30>0.4,]
b<-na.omit(b)

plot<-ggplot()+geom_point(aes(a$V2,a$V12),alpha=0.2,cex=0.1)+geom_point(aes(b$V2,b$V12),col="red")+theme_bw()+ylim(0,1)+ylab("Allele freq diff. - pre vs post treatment")
fst_plot<-ggplot()+geom_point(aes(a_fst$V2,a_fst$V13),alpha=0.5)+theme_bw()+ylab("Fst - pre vs post treatment")

fst_plot + plot + plot_layout(ncol=1)






# IVM - chr 2
cd ~/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_VM

awk '{if($1=="hcontortus_chr2_Celeg_TT_arrow_pilon" && $2>=0 && $2<=10000000) print $0}' OFS="\t" XQTL_IVM_pwc > XQTL_IVM_pwc_chr1_0-10Mb
awk '{if($1=="hcontortus_chr2_Celeg_TT_arrow_pilon" && $2>=0 && $2<=10000000) print $0}' OFS="\t" XQTL_IVM.merged.fst > XQTL_IVM.merged_0-10Mb.fst

R-3.4.0
library(ggplot2)
library(patchwork)
a<-read.table("XQTL_IVM_pwc_chr1_0-10Mb",header=F,na.strings="na")
a_fst<-read.table("XQTL_IVM.merged_0-10Mb.fst",header=F,na.strings="na")
b<-a[a$V12>0.3 & a$V19>0.3,]
b<-na.omit(b)

plot<-ggplot()+geom_point(aes(a$V2,a$V12),alpha=0.2,cex=0.1)+geom_point(aes(b$V2,b$V12),col="red")+theme_bw()+ylim(0,1)+ylab("Allele freq diff. - pre vs post treatment")
fst_plot<-ggplot()+geom_point(aes(a_fst$V2,a_fst$V13),alpha=0.5)+theme_bw()+ylab("Fst - pre vs post treatment")

fst_plot + plot + plot_layout(ncol=1)



# IVM - chr 5
cd ~/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_VM

awk '{if($1=="hcontortus_chr5_Celeg_TT_arrow_pilon" && $2>=30000000 && $2<=50000000) print $0}' OFS="\t" XQTL_IVM_pwc > XQTL_IVM_pwc_chr5_30-50Mb
awk '{if($1=="hcontortus_chr5_Celeg_TT_arrow_pilon" && $2>=30000000 && $2<=50000000) print $0}' OFS="\t" XQTL_IVM.merged.fst > XQTL_IVM.merged_30-50Mb.fst

R-3.4.0
library(ggplot2)
library(patchwork)
a<-read.table("XQTL_IVM_pwc_chr5_30-50Mb",header=F,na.strings="na")
a_fst<-read.table("XQTL_IVM.merged_30-50Mb.fst",header=F,na.strings="na")
b<-a[a$V12>0.4 & a$V19>0.4 & a$V25>0.4 & a$V30>0.4,]
b<-na.omit(b)

plot<-ggplot()+geom_point(aes(a$V2,a$V12),alpha=0.2,cex=0.1)+geom_point(aes(b$V2,b$V12),col="red")+theme_bw()+ylim(0,1)+ylab("Allele freq diff. - pre vs post treatment")
fst_plot<-ggplot()+geom_point(aes(a_fst$V2,a_fst$V13),alpha=0.5)+theme_bw()+ylab("Fst - pre vs post treatment")

fst_plot + plot + plot_layout(ncol=1)
```











################################################################################
## Dose response data - quick look
```shell
cd ~/lustre118_link/hc/XQTL/04_VARIANTS/DOSE_RESPONSE

R
```

```R
library(ggplot2)
library(ggrepel)


data <- read.table("DOSE_RESPONSE.merged.fst",header=F)


ggplot(data)+geom_point(aes(data$V2,data$V7,col=data$V1),alpha=0.5,size=0.5)+facet_grid(.~data$V1)+theme_bw()

chr5 <- data[data$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon",]

ggplot(chr5,aes(V2,V7,label=V2))+geom_point()+xlim(2.5e7,4.5e7)+geom_text_repel(data = subset(chr5,V7 > 0.25))+theme_bw()

```









################################################################################


chromosome.labels <- c("I","II","III","IV","V", "X" )
names(chromosome.labels) <- c("hcontortus_chr1_Celeg_TT_arrow_pilon","hcontortus_chr2_Celeg_TT_arrow_pilon","hcontortus_chr3_Celeg_TT_arrow_pilon","hcontortus_chr4_Celeg_TT_arrow_pilon","hcontortus_chr5_Celeg_TT_arrow_pilon","hcontortus_chrX_Celeg_TT_arrow_pilon")



library(ggplot2)
library(patchwork)


data<-read.table("XQTL_IVM.merged.fst",header=F)
data5 <- data[data$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon",]

top1 <- data[data$V13 > quantile(data$V13,prob=1-0.5/100),]
top1_5 <- top1[top1$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon",]

ggplot(data,aes(V2,V13,group,V1))+
	facet_grid(V1~.,labeller = labeller(V1 = chromosome.labels))+
	geom_point(data=top1,aes(top1$V2,top1$V13,alpha=top1$V4),col="red")+
	geom_point(aes(alpha=V4),size=0.5)+
	theme_bw()+
	labs(x="Genome position (bp)", y="Fst (XQTL_IVM 1:5 (V13))", alpha="window coverage")






plot1<-ggplot(data=data5,aes(V2,V13))+
	theme_bw()+
	 geom_rect(aes(xmin=27640082,ymin=0,xmax=27661088,ymax=Inf),color="grey",fill="grey")+  # glc-3 HCON_00148840
   	 geom_rect(aes(xmin=28367524,ymin=0,xmax=28394349,ymax=Inf),color="grey",fill="grey")+  # osm-6 HCON_00149350
	geom_rect(aes(xmin=36159873,ymin=0,xmax=36160068,ymax=Inf),color="grey",fill="grey")+  # microsat 8a20 - V3, approx in V4
   	geom_rect(aes(xmin=37487982,ymin=0,xmax=37497398,ymax=Inf),color="red",fill="red")+  # cky-1 NPAS4 HCON_00155390
  	geom_rect(aes(xmin=37591595,ymin=0,xmax=37598257,ymax=Inf),color="grey",fill="grey")+  # HCON_00155440  GPCR serpentin "chemoreceptor"
	geom_rect(aes(xmin=37538652,ymin=0,xmax=37556986,ymax=Inf),color="grey",fill="grey")+  # lin-12 HCON_00155420
	geom_rect(aes(xmin=39459752,ymin=0,xmax=39485324,ymax=Inf),color="grey",fill="grey")+  # egl2
   	geom_rect(aes(xmin=42448291,ymin=0,xmax=42464088,ymax=Inf),color="grey",fill="grey")+  # lgc-55 HCON_00158990
   	geom_rect(aes(xmin=45761945,ymin=0,xmax=45765591,ymax=Inf),color="grey",fill="grey")+  # HCON_00161700 - neurotransmitter gated channel ligand domain containing protein
  	geom_rect(aes(xmin=45965820,ymin=0,xmax=45985698,ymax=Inf),color="grey",fill="grey")+  # dop-5 HCON_00161830
  	geom_rect(aes(xmin=45986412,ymin=0,xmax=46020767,ymax=Inf),color="grey",fill="grey")+  # hen-1 HCON_00161834  	
  	geom_rect(aes(xmin=47243601,ymin=0,xmax=47259265,ymax=Inf),color="grey",fill="grey")+  # pgp-11 HCON_0016278
  	geom_point(data=top1_5,aes(top1_5$V2,top1_5$V13,alpha=top1_5$V4),col="red")+
  	geom_point(aes(alpha=V4),size=0.5)+
  	labs(x="Genome position (bp)", y="Fst (XQTL_IVM 1:5 (V13))", alpha="window coverage")

plot1<-plot1+xlim(3.7e7,3.8e7)
plot1<-plot1+xlim(4.5e7,4.9e7)

plot2 <- ggplot(genes)+ geom_rect(aes(xmin=genes$V2,ymin=0,xmax=genes$V3,ymax=1),color="grey",fill="grey")+geom_text_repel(aes(label=genes$V1,x=genes$V2+(genes$V3-genes$V2),y=0.5),angle = 90,size=2)

plot2<-plot2+xlim(3.7e7,3.8e7)

plot1+plot2+plot_layout(ncol=1)

genes <-read.table("chr5_geneIDs.txt",header=F)

plot + geom_rect(aes(xmin=genes$V2,ymin=0,xmax=genes$V3,ymax=1),color="grey",fill="grey")
plot+geom_vline(aes(x=genes$V2))











library(ggplot2)
library(patchwork)


fst <-read.table("XQTL_US_FIELD.merged.fst",header=F)

# susceptibles
plot_1 <- ggplot(fst,aes(V2,V21))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Idaho_S")+theme_bw()

# UGA_S vs moderates
plot_2 <- ggplot(fst,aes(V2,V7))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Ober_NY_Mod")+theme_bw()
plot_3 <- ggplot(fst,aes(V2,V9))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Histon_Maryland_Mod")+theme_bw()
plot_4 <- ggplot(fst,aes(V2,V11))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Krushing_Ohio_Mod")+theme_bw()


# Idaho_S vs moderates
plot_5 <- ggplot(fst,aes(V2,V37))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="Idaho_S vs Ober_NY_goat_Mod")+theme_bw()
plot_6 <- ggplot(fst,aes(V2,V51))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="Idaho_S vs Histon_Maryland_Mod")+theme_bw()
plot_7 <- ggplot(fst,aes(V2,V63))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="Idaho_S vs Krushing_Ohio_Mod")+theme_bw()


# UGA_S vs high resistance
plot_8 <- ggplot(fst,aes(V2,V13))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Adlay_Georgia_R")+theme_bw()
plot_9 <- ggplot(fst,aes(V2,V15))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Mulligan_Georgia_R")+theme_bw()
plot_10 <- ggplot(fst,aes(V2,V17))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Rau-Shelby_Georgia_R")+theme_bw()
plot_11 <- ggplot(fst,aes(V2,V19))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Pena_NY_alpaca_R")+theme_bw()
plot_12 <- ggplot(fst,aes(V2,V23))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="UGA_S vs Strickland_Georgia_R")+theme_bw()


# Idaho_S vs high resistance
plot_13 <- ggplot(fst,aes(V2,V73))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="Idaho_S vs Adlay_Georgia_R")+theme_bw()
plot_14 <- ggplot(fst,aes(V2,V81))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="Idaho_S vs Mulligan_Georgia_R")+theme_bw()
plot_15 <- ggplot(fst,aes(V2,V87))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="Idaho_S vs Rau-Shelby_Georgia_R")+theme_bw()
plot_16 <- ggplot(fst,aes(V2,V91))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="Idaho_S vs Pena_NY_alpaca_R")+theme_bw()
plot_17 <- ggplot(fst,aes(V2,V95))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+ylim(0,1)+labs(title="Idaho_S vs Strickland_Georgia_R")+theme_bw()




plot_1 + plot_2 + plot_3 + plot_4 + plot_layout(ncol=1)
plot_1 + plot_5 + plot_6 + plot_7 + plot_layout(ncol=1)


plot_8 + plot_9 + plot_10 + plot_11 + plot_12 + plot_layout(ncol=1)
plot_13 + plot_14 + plot_15 + plot_16 + plot_17 + plot_layout(ncol=1)

plot_1 + plot_2 + plot_3 + plot_4 + plot_8 + plot_9 + plot_10 + plot_11 + plot_12 + plot_layout(ncol=1)






















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












#------------ PLAY
plot7.5 <-     ggplot(sample7.5,aes(window*10000,log10(Pi),label=window*10000))+
                              geom_point(size=1,alpha=1,col="black")+theme_bw()+ylim(-5,-1)+labs(title="Rau_Shelby_R")+
                              geom_text_repel(data=subset(sample7.5,Pi <= 0.001))
acr8 31508821	31530834

plot_10 + plot7.5 +plot_layout(ncol=1)




plot10in <- plot_10 + xlim(3e7,3.3e7)

lev_plo1v5 <- ggplot(lev,aes(V2,V13))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+labs(title="XQTL_LEV_1v5")+theme_bw()+xlim(3e7,3.3e7)

lev_plot2v6 <- ggplot(lev,aes(V2,V27))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+labs(title="XQTL_LEV_1v5")+theme_bw()+xlim(3e7,3.3e7)

lev_plot3v7 <- ggplot(lev,aes(V2,V39))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+labs(title="XQTL_LEV_1v5")+theme_bw()+xlim(3e7,3.3e7)

lev_plot4v8 <- ggplot(lev,aes(V2,V49))+geom_point(size=0.5,alpha=0.5)+facet_grid(.~V1)+labs(title="XQTL_LEV_1v5")+theme_bw()+xlim(3e7,3.3e7)


plot10in + lev_plo1v5 + lev_plot2v6 + lev_plot3v7 + lev_plot4v8 +  plot_layout(ncol=1)

















# mapping gilleard ISE, CAVR, WRS RNAseq to compare with XQTL data



#--- map reads

for i in `ls -1 | grep "_1.fastq.gz" | sed 's/_1.fastq.gz//g'  `; do \
bsub.py 20 --threads 8 starmap_${i} /nfs/users/nfs_s/sd21/lustre118_link/software/TRANSCRIPTOME/STAR/bin/Linux_x86_64/STAR \
--runThreadN 8 \
--genomeDir /nfs/users/nfs_s/sd21/lustre118_link/REFERENCE_SEQUENCES/haemonchus_contortus/hc_v4chr_rnaseq_star_index \
--readFilesIn ${i}_1.fastq.gz ${i}_2.fastq.gz \
--readFilesCommand zcat \
--alignIntronMin 10 \
--outTmpDir starmap_${i}_tmp \
--outFileNamePrefix starmap_${i}_other_out \
--outSAMtype BAM SortedByCoordinate \
--outWigType bedGraph \
--twopassMode Basic \
; done





#----- SNPeff
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
