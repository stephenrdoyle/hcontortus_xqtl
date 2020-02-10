# XQTL analysis

working directory: /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL

```shell
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL


mkdir 00_SCRIPTS
mkdir 01_REFERENCE
mkdir 02_RAW
mkdir 03_MAPPING
mkdir 04_VARIANTS
```

################################################################################################
# Step 1 -reference
################################################################################################
```shell
cd 01_REFERENCE

cp /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/REF/HAEM_V4_final.chr.fa .
samtools faidx HAEM_V4_final.chr.fa
samtools dict HAEM_V4_final.chr.fa > HAEM_V4_final.chr.dict
```






################################################################################################
# Parents
################################################################################################

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



#--- map reads

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/03_MAPPING/PARENTS
screen
while read NAME; do ~sd21/bash_scripts/run_bwamem_splitter ${NAME} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_PARENTS/TRIM/${NAME}.paired_R1.fastq.gz /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/02_RAW/RAW_PARENTS/TRIM/${NAME}.paired_R2.fastq.gz; done < sample.list &


#--- clean up
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



# setup for mpileup
mkdir /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS

ls $PWD/*bam > /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS/bam.list

cd /lustre/scratch118/infgen/team133/sd21/hc/XQTL/04_VARIANTS/PARENTS


#--- run mpileup
~sd21/bash_scripts/run_mpileup XQTL_PARENTS /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list

cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_PARENTS.mpileup &

rm *tmp*





################################################################################################
# XQTL
################################################################################################


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
pathfind -t lane -i 23204_8 --filet