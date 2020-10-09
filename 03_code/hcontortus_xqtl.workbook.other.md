hcontortus_xqtl.workbook.other



###########################################################################################
# allele frequencies for Jamie


# BZ
cd ~/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_BZ

awk '{if($1=="hcontortus_chr1_Celeg_TT_arrow_pilon" && $2>=5000000 && $2<=10000000) print $0}' OFS="\t" XQTL_BZ_pwc > XQTL_BZ_pwc_chr1_5-10Mb
awk '{if($1=="hcontortus_chr1_Celeg_TT_arrow_pilon" && $2>=5000000 && $2<=10000000) print $0}' OFS="\t" XQTL_BZ.merged.fst > XQTL_BZ.merged_25-35Mb.fst


```R
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
```


```shell
while read chr pos rest; do \
	awk -v chr="$chr" -v pos="$pos" -F'[\t/]' '{if($1==chr && $2==pos && $4=="2") print $1,$2,$3,$5"/"$6,$9,$11/$12,$13/$14,$15/$16,$17/$18,$19/$20,$21/$22,$23/$24,$25/$26}' OFS="\t" XQTL_BZ_rc >> XQTL_BZ_peakSNPs.txt; \
done < chr1.allelefreqdiff.gt0.4.txt &
```



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




# mapping gilleard ISE, CAVR, WRS RNAseq to compare with XQTL data
```
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
```




### mapping Gilleard parental strains to compare against XQTL data
- gillard BC was mapped against the V3 version of the genome, so needs to be done again.





```shell
cd ~/lustre118_link/hc/XQTL/03_MAPPING/GILLEARD_BC/DNASEQ

parental_ISE_susceptible_5764_1	5764_1
parental_ISE_susceptible_6514_1	6514_1
parental_ISE_susceptible_7756_1	7756_1
parental_WRS_resistant_5764_2	5764_2
parental_WRS_resistant_6474_1	6474_1
parental_WRS_resistant_7756_2	7756_2
parental_CAVR_resistant_5764_3	5764_3
parental_CAVR_resistant_6474_2	6474_2
parental_CAVR_resistant_7756_3	7756_3



while read name id; do pf data --rename --filetype fastq -l ./ --type lane --id ${id}; done < sample_lanes.list

while read name lane; do ~sd21/bash_scripts/run_bwamem_splitter ${name} /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa $PWD/${lane}_1.fastq.gz $PWD/${lane}_2.fastq.gz; done < sample_lanes.list



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




#--- run mpileup
ls -1 *bam > bam.list

~sd21/bash_scripts/run_mpileup.sh GBX_PARENTS /lustre/scratch118/infgen/team133/sd21/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa bam.list &

cat $(find . -name "*tmp.mpileup" | sort -V) | sed 's/\t\t/\t!\t!/g' > XQTL_DOSE_RESPONSE.mpileup &

rm *tmp*


#mpileup to popoolation2
~sd21/bash_scripts/run_mpileup2popoolation2 GBX /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa GBX.mpileup 50 5000 &
```
