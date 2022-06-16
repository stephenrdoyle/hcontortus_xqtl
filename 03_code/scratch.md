# scratch space for random code than needs a home




## playing around with correlations between EC50 and US farm Fst data

```bash
/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/IVM

ln -s ../../04_VARIANTS/US_FIELD/XQTL_US_FIELD.merged.fst
```


```R
library(tidyverse)

# manually input the EC50 data from Table SX
# farm order = F9,F2,F3,F4,F5,F6,F7,F8,F10
ivm_conc <- c(0.977,12.04,9.15,11.27,297.9,619,312.1,NA,259.4)
bz_conc <- c(9.67,19.93,15,7.71,15,61.57,29.6,NA,29.6)
lev_conc <- c(1.5,1.5,0.58,0.57,0.91,0.96,9.36,NA,0.88)
mox_conc <- c(26.6,221.4,166.1,205.3,6038,10000,3987,NA,5263)


data <- read.table("XQTL_US_FIELD.merged.fst",header=F)

# farm order = F9,F2,F3,F4,F5,F6,F7,F8,F10
fst_data <- data[,c(1,2,21,7,9,11,13,15,17,19,23)]


ivm_slopes <-
     apply(as.matrix(fst_data[,-c(1,2)]),1,
          function(y) {
               fit <- lm(y~ivm_conc,data=data.frame(y,ivm_conc))
               coef(fit)[2]
          })

bz_slopes <-
     apply(as.matrix(fst_data[,-c(1,2)]),1,
          function(y) {
               fit <- lm(y~bz_conc,data=data.frame(y,bz_conc))
               coef(fit)[2]
          })

lev_slopes <-
     apply(as.matrix(fst_data[,-c(1,2)]),1,
          function(y) {
               fit <- lm(y~lev_conc,data=data.frame(y,lev_conc))
               coef(fit)[2]
          })

mox_slopes <-
     apply(as.matrix(fst_data[,-c(1,2)]),1,
          function(y) {
               fit <- lm(y~mox_conc,data=data.frame(y,mox_conc))
               coef(fit)[2]
          })


data_out <- data.frame(fst_data,ivm_slopes,bz_slopes,lev_slopes,mox_slopes,check.names=FALSE)


ggplot(data_out,aes(V2,ivm_slopes)) + geom_point() + facet_grid(V1~.)
```





# generate CMH data
perl /nfs/users/nfs_s/sd21/lustre118_link/software/POOLSEQ/popoolation2_1201/cmh-test.pl --min-count 4 --min-coverage 30 --max-coverage 2% --population 1-5,4-8 --input XQTL_LEV.raw.sync --output XQTL_LEV.rep2.cmh --min-pvalue 0.01


# SNPeff positional data
#grep "HIGH\|MODE" <(zcat XQTL_LEV.raw.snpeff.vcf.gz) | awk -F'[\t|]' '{print $1,$2,$14,$9,$10}' > snpeff_high_moderate.pos

#grep "HIGH\|MODE" <(zcat XQTL_LEV.raw.snpeff.vcf.gz) | awk -F'[\t;|]' '{ print $1,$2,$4,$5,$6,$29,$24,$25,$32,$33}' OFS="\t" > snpeff_filtered_data.txt

zcat XQTL_LEV.raw.snpeff.vcf.gz | grep -v "#" | grep -v "prime_UTR" | grep -v "intron/|intra" | awk -F'[\t;|]' '{ print $1,$2,$4,$5,$6,$29,$24,$25,$32,$33}' OFS="\t" > snpeff_coding_variants.filtered.txt


zcat XQTL_LEV.raw.snpeff.vcf.gz | awk -F'[\t;|]' '{ print $1,$2,$4,$5,$6,$29,$24,$25,$32,$33}' OFS="\t" | grep "HCON" | grep -v "#" > genic.variants.data


cat XQTL_IVM.raw.snpeff.vcf | awk -F'[\t;|]' '{ print $1,$2,$4,$5,$6,$29,$24,$25,$32,$33}' OFS="\t" | grep "HCON" | grep -v "#" > genic.variants.data

```R
merge_cmh_snpeff <- function(snpeff_data,cmh_data){
library(tidyverse)
library(data.table)

cmh_data <- fread(cmh_data)
snpeff_data <- fread(snpeff_data)


# join datasets together, keeping only positions in both
data <- inner_join(snpeff_data, cmh_data, by=c("V1","V2"))


colnames(data) <- c("CHR", "SNP_POS","REF_ALLELE", "VARINAT_ALLELE", "SNP_QUAL", "TRANSCRIPT", "VARIANT_TYPE", "VARIANT_EFFECT", "DNA_CHANGE", "AA_CHANGE", "REF_BASE", "FREQ_POP1", "FREQ_POP2", "FREQ_POP3", "FREQ_POP4", "FREQ_POP5", "FREQ_POP6", "FREQ_POP7", "FREQ_POP8", "CMH")

nrow(data)
#[1] 619751

data_coding <- filter(data, VARIANT_TYPE!="intron_variant" & VARIANT_TYPE!="5_prime_UTR_variant" & VARIANT_TYPE!="3_prime_UTR_variant")
write.table(data_coding, "snpeff_coding_variants_per_gene.cmh.table", row.names=FALSE, sep="\t", quote = FALSE)


data_noncoding <- filter(data, VARIANT_TYPE=="intron_variant" | VARIANT_TYPE=="5_prime_UTR_variant" | VARIANT_TYPE=="3_prime_UTR_variant")
write.table(data_noncoding, "snpeff_noncoding_variants_per_gene.cmh.table", row.names=FALSE, sep="\t", quote = FALSE)
}

# ivm
merge_cmh_snpeff("genic.variants.data","XQTL_IVM.rep3.cmh")



# remove poorer quality variants
data_Q900 <- filter(data,SNP_QUAL>900)
nrow(data_Q900)
#[1] 306960

write.table(data_Q900, "snpeff_coding_variants.qual900.cmh.table", row.names=FALSE, sep="\t", quote = FALSE)
```

grep "HCON_00107690" snpeff_coding_variants.cmh.table
grep "HCON_00107690" snpeff_coding_variants.qual900.cmh.table





snpeff_missense_pos <- read.table("lev_missense.pos")

data_missense <- inner_join(cmh, snpeff_missense_pos, by=c("V1","V2"))
colnames(data_missense) <- c("CHR", "SNP", "REF_BASE", "FREQ_POP1", "FREQ_POP2", "FREQ_POP3", "FREQ_POP4", "FREQ_POP5", "FREQ_POP6", "FREQ_POP7", "FREQ_POP8", "CMH")







data_high_mod %>% grepl("HCON_00107690",TRANSCRIPT)
#> nothing found
data_high_mod %>% grepl("HCON_00107700",TRANSCRIPT)
#> nothing found

# HCON_00107690: 14977822	14984960
data_high_mod %>% filter(CHR=="hcontortus_chr4_Celeg_TT_arrow_pilon" & SNP>=14977822 & SNP<=14984960)

#HCON_00107700: 14991116	15003109
data_high_mod %>% filter(CHR=="hcontortus_chr4_Celeg_TT_arrow_pilon" & SNP>=14991116 & SNP<=15003109)
data_missense %>% filter(SNP>14991116 & SNP<=15003109)
}



get_sig_gene_variants  <- function(snpff_data, cmh_data, gene_ID){

library(tidyverse)

cmh <- read.table(cmh_data)
snpeff_pos <- read.table(snpff_data)

# join datasets together, keeping only positions in both
data_high_mod <- inner_join(cmh, snpeff_pos, by=c("V1","V2"))
colnames(data_high_mod) <- c("CHR", "SNP", "REF_BASE", "FREQ_POP1", "FREQ_POP2", "FREQ_POP3", "FREQ_POP4", "FREQ_POP5", "FREQ_POP6", "FREQ_POP7", "FREQ_POP8", "CMH", "TRANSCRIPT")

data_high_mod %>% grepl(x=data_high_mod, pattern=gene_ID, TRANSCRIPT)
}



get_sig_gene_variants("high_moderate.pos", "XQTL_LEV.raw.sync.cmh", "HCON_00107690")
get_sig_gene_variants("high_moderate.pos", "XQTL_LEV.raw.sync.cmh", "HCON_00107700")




# extracting top SNPs in protein kinase genes within the chr5 QTL
for i in HCON_00155230 \
HCON_00155340 \
HCON_00155280 \
HCON_00155380 \
HCON_00155295 \
HCON_00155270 \
HCON_00155190 \
HCON_00155240; do
grep "${i}" snpeff_coding_variants_per_gene.cmh.table | sort -k20,20g | head -n1;
grep "${i}" snpeff_noncoding_variants_per_gene.cmh.table | sort -k19,19g | head -n1;
done | sort -k19,19g




for i in  HCON_00155230 \
HCON_00155377 \
HCON_00155340 \
HCON_00155280 \
HCON_00155380 \
HCON_00155375 \
HCON_00155295 \
HCON_00155270 \
HCON_00155190 \
HCON_00155220 \
HCON_00155215 \
HCON_00155240 \
HCON_00155290 \
HCON_00155390 \
HCON_00155400 \
HCON_00155200 \
HCON_00155210 \
HCON_00155300 \
HCON_00155310 \
HCON_00155320 \
HCON_00155330 \
HCON_00155360 \
HCON_00155402 \
HCON_00155405 \
HCON_00155410; do
grep "${i}" snpeff_coding_variants_per_gene.cmh.table | sort -k20,20g | head -n1;
grep "${i}" snpeff_noncoding_variants_per_gene.cmh.table | sort -k19,19g | head -n1;
done | sort -k19,19g
