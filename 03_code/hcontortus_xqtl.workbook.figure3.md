# Figure 3 - Benzimidazole

1. [Figure 3 A](#figure3a)
2. [Figure 3 B](#figure3b)
3. [Figure 3 C](#figure3c)

---

## Figure 3 A  <a name="figure3a"></a>

Aim is to show genetic differentiation between pre/post benzimidazole treatment on chromosome 1 where the major peak lies

```shell
# working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/BZ
```

### R to plot
```R
# load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# import fst data
xqtl_bz_fst <- read.table("XQTL_BZ.merged.fst", header = F)

# reformat
xqtl_bz_fst <- dplyr::select(xqtl_bz_fst,  V1,  V2, V13, V39, V49)
xqtl_bz_fst <- xqtl_bz_fst %>% mutate(mean_FST = rowMeans(select(.,V13,V39,V49)))
colnames(xqtl_bz_fst) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN")



# calculate a genome wide significance cutoff
gw_r1 <- mean(xqtl_bz_fst$FST_R1)+3*sd(xqtl_bz_fst$FST_R1)
gw_r2 <- mean(xqtl_bz_fst$FST_R2)+3*sd(xqtl_bz_fst$FST_R2)
gw_r3 <- mean(xqtl_bz_fst$FST_R3)+3*sd(xqtl_bz_fst$FST_R3)
gw_mean <- mean(xqtl_bz_fst$FST_MEAN)+3*sd(xqtl_bz_fst$FST_MEAN)


# extract chromosome 1 data
xqtl_bz_fst_chr1 <- xqtl_bz_fst[xqtl_bz_fst$CHR == "hcontortus_chr1_Celeg_TT_arrow_pilon", ]

#xqtl_bz_fst_chr1_peaks <- xqtl_bz_fst_chr1[(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD),]


# # get predicted peak windows
#
# peaks <- read.table("peak.windows.bed", header = T)
#
# peak_subset <- xqtl_bz_fst_chr1[(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD), ]

# set colours for data, based on thresholds and replicates
xqtl_bz_fst_chr1 <- mutate(xqtl_bz_fst_chr1,
       point_colour = case_when(
       (FST_R1 < gw_r1 & FST_R2 < gw_r2 & FST_R3 < gw_r3) ~ "#6699FFFF",
       (FST_R1 > gw_r1 & FST_R2 < gw_r2 & FST_R3 < gw_r3) ~ "#FFCC00FF",
       (FST_R1 < gw_r1 & FST_R2 > gw_r2 & FST_R3 < gw_r3) ~ "#FFCC00FF",
       (FST_R1 < gw_r1 & FST_R2 < gw_r2 & FST_R3 > gw_r3) ~ "#FFCC00FF",
       (FST_R1 > gw_r1 & FST_R2 > gw_r2 & FST_R3 < gw_r3) ~ "#FF9900FF",
       (FST_R1 > gw_r1 & FST_R2 < gw_r2 & FST_R3 > gw_r3) ~ "#FF9900FF",
       (FST_R1 > gw_r1 & FST_R2 > gw_r2 & FST_R3 > gw_r3) ~ "#FF0000FF",))

plot_a <- ggplot(xqtl_bz_fst_chr1, aes(POS, FST_MEAN, colour=point_colour, size=ifelse(FST_MEAN>gw_mean,1,0.3))) +
               geom_hline(yintercept = gw_mean, linetype = "dashed", col = "black") +
               geom_vline(xintercept = 7029790, linetype = "dashed", col = "darkgrey", size=1) +                  
               geom_point() + facet_grid(CHR~.) + scale_color_identity() + scale_size_identity() +
               ylim(0, 0.1) + xlim(0, 50e6) +
               theme_bw() + theme(legend.position = "none", text = element_text(size = 10)) +
               labs(title = "A", x = "Genomic position (bp)", y = expression(paste("Genetic differentiation between pre- and post-treatment", " (",~italic(F)[ST],")"))) +
               facet_grid(CHR ~ .)




# make the plot
plot_a <- ggplot(xqtl_bz_fst_chr1) +
     geom_hline(yintercept = gw_mean, linetype = "dashed", col = "black") +
     geom_vline(xintercept = 7029790, linetype = "dashed", col = "grey")+
     geom_point(aes(POS, FST_MEAN, group = CHR), col = "cornflowerblue", size = 0.5) +
     geom_point(data = subset(xqtl_bz_fst_chr1,(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD[1]) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD[1]) & (V13 > genomewide_sig)), aes(V2, V13), col = "red", size = 1) +
     geom_point(data = subset(xqtl_bz_fst_chr1,(xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD[2]) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD[2]) & (V13 > genomewide_sig)), aes(V2, V13), col = "red", size = 1) +
     geom_point(data = subset(xqtl_bz_fst_chr1, (xqtl_bz_fst_chr1$V2 >= peaks$PEAK_START_COORD[3]) & (xqtl_bz_fst_chr1$V2 <= peaks$PEAK_END_COORD[3]) & (V13 > genomewide_sig)), aes(V2, V13), col = "red", size = 1) +
     ylim(0, 0.1) + xlim(0, 50e6) +
     theme_bw() + theme(legend.position = "none", text = element_text(size = 10)) +
     labs(title = "A", x = "Genomic position (bp)", y = expression(paste("Genetic differentiation between pre- and post-treatment", " (",~italic(F)[ST],")")))
     facet_grid( V1 ~ .)
```
---



## Figure 3 B - beta tubulin isotype 1 data <a name="figure3b"></a>

Aim is to show allele frequencies of Phe167Tyr and Phe200Tyr pre and post treatment, including non-treated time matched control

### prepare the data
```shell
#working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_BZ

# get sample IDs from bam.lists
grep "pre" bam.list  > bz_pretreatment_samples.list
grep "post" bam.list  > bz_posttreatment_samples.list

grep "pre" ../XQTL_CONTROL/bam.list  > control_pretreatment_samples.list
grep "post" ../XQTL_CONTROL/bam.list  > control_posttreatment_samples.list


# extract allele counts
vcftools --vcf XQTL_BZ.raw.snpeff.vcf --keep bz_pretreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out bz_pretreatment
vcftools --vcf XQTL_BZ.raw.snpeff.vcf --keep bz_posttreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out bz_posttreatment

vcftools --gzvcf ../XQTL_CONTROL/XQTL_CONTROL.raw.vcf.gz --keep control_pretreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out control_pretreatment
vcftools --gzvcf ../XQTL_CONTROL/XQTL_CONTROL.raw.vcf.gz --keep control_posttreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out control_posttreatment

#where "btub.positions" contains:
#P167
#hcontortus_chr1_Celeg_TT_arrow_pilon 7029569
#P198
#hcontortus_chr1_Celeg_TT_arrow_pilon 7029788
# P200
#hcontortus_chr1_Celeg_TT_arrow_pilon 7029790

# calculate variant frequencies
for i in `ls bz*AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done

for i in `ls control*AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done

# didnt use this, but keeping for another day
#
# hcontortus_chr1_Celeg_TT_arrow_pilon	7029569	0.125	0.197368	0.223301	0.264463	0.15	0.128571	0.131148	0.155172
# hcontortus_chr1_Celeg_TT_arrow_pilon	7029790	0.721311	0.805556	0.866071	0.926316	0.534884	0.477612	0.446281	0.347368
#
#
#
# awk '{print $1,$2,$7,$8,$9,$10,"bz_pre"}' OFS="\t" bz.freq > bz.freq2
# awk '{print $1,$2,$3,$4,$5,$6,"bz_post"}' OFS="\t" bz.freq >> bz.freq2
#
#
# awk '{print $1,$2,$7,$3,"BZ_treated","R1"}' OFS="\t" bz.freq > bz.freq2
# awk '{print $1,$2,$8,$4,"BZ_treated","R2"}' OFS="\t" bz.freq >> bz.freq2
# awk '{print $1,$2,$9,$5,"BZ_treated","R3"}' OFS="\t" bz.freq >> bz.freq2
# awk '{print $1,$2,$10,$6,"BZ_treated","R4"}' OFS="\t" bz.freq >> bz.freq2
#
#
#
# ggplot(a) + geom_segment(aes(x="1", xend="2", y=V3, yend=V4,col=factor(V2),group=V2), size=.75)+facet_grid(V5~V2)
#
# grep "^hcon" out.AD.FORMAT | awk -F '[\t,]' '{for(i=3;i<=NF;i+=2) print $1,$2,(i+=1)/($i+(i+=1))}'
#
# for ((i=3,j=4;i<=j;i+=2,j+=2)); do \
#      grep "^hcon" out.AD.FORMAT | \
#      awk -v i=$i -v j=$j -F '[\t,]' '{if(i<NF) print $1,$2,$j/($i+$j)}'
# done

```

### R to plot
```R
# load required libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)

# load and reformat the data
bz_pre <- read.table("bz_pretreatment.ADfreq")
bz_pre <- filter(bz_pre,V1=="hcontortus_chr1_Celeg_TT_arrow_pilon")
bz_pre <- dplyr::select(bz_pre,  V1,  V2,  V3 , V5, V6)
colnames(bz_pre) <- c("CHR", "POS", "R1", "R2", "R3")
bz_pre <- melt(bz_pre,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")

bz_post <- read.table("bz_posttreatment.ADfreq")
bz_post <- filter(bz_post,V1=="hcontortus_chr1_Celeg_TT_arrow_pilon")
bz_post <- dplyr::select(bz_post,  V1,  V2,  V3 , V5, V6)
colnames(bz_post) <- c("CHR", "POS", "R1", "R2", "R3")
bz_post <- melt(bz_post,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")

bz_data <- dplyr::full_join(bz_pre,  bz_post,  by = c("CHR", "POS", "SAMPLE_ID"))
bz_data$TREATMENT <- "2. Benzimidazole"
colnames(bz_data) <- c("CHR", "POS", "SAMPLE_ID", "PRE_TREATMENT", "POST_TREATMENT", "TREATMENT")


control_pre <- read.table("control_pretreatment.ADfreq")
control_pre <- filter(control_pre,V1=="hcontortus_chr1_Celeg_TT_arrow_pilon")
colnames(control_pre) <- c("CHR", "POS", "R1", "R2", "R3")
control_pre <- melt(control_pre,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")

control_post <- read.table("control_posttreatment.ADfreq")
control_post <- filter(control_post,V1=="hcontortus_chr1_Celeg_TT_arrow_pilon")
colnames(control_post) <- c("CHR", "POS", "R1", "R2", "R3")
control_post <- melt(control_post,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")

control_data <- dplyr::full_join(control_pre,  control_post,  by = c("CHR", "POS", "SAMPLE_ID"))
control_data$TREATMENT <- "1. Untreated"
colnames(control_data) <- c("CHR", "POS", "SAMPLE_ID", "PRE_TREATMENT", "POST_TREATMENT", "TREATMENT")


# bring datasets together
data <- dplyr::bind_rows(control_data,  bz_data)

# change the labels
data <- data %>%
  mutate(POS = str_replace_all(POS,  c("7029569", "7029790"),  c("Phe167Tyr", "Phe200Tyr")))



# make the plot
plot <- ggplot(data) +
     geom_segment(aes(x = "1.PRE",  xend = "2.POST",  y = PRE_TREATMENT,  yend = POST_TREATMENT, col = factor(SAMPLE_ID), group = POS),  size = 1) +
     labs(title = "B", x = "Sampling time-point", y = "Resistant allele frequency", col = "Replicate") +
     ylim(0, 1) +
     facet_grid(POS ~ TREATMENT) +
     theme_bw() + theme(text = element_text(size = 10))

# perform pairwise t tests between pre/post for each SNP on BZ treated samples
bz_data_stats <- bz_data %>%
  gather(key = "TREATMENT",  value = "FREQ",  PRE_TREATMENT,  POST_TREATMENT)

bz_stat.test <- bz_data_stats %>%
    group_by(POS) %>%
    pairwise_t_test(
      FREQ ~ TREATMENT,  paired = TRUE,
      p.adjust.method = "bonferroni"
      ) %>%
    select(-df,  -statistic,  -p) # Remove details

bz_stat.test$TREATMENT <- "2. Benzimidazole"
bz_stat.test <- bz_stat.test %>%
  mutate(POS = str_replace(POS,  c("7029569", "7029790"),  c("Phe167Tyr", "Phe200Tyr")))

# perform pairwise t tests between pre/post for each SNP on control samples
control_data_stats <- control_data %>%
  gather(key = "TREATMENT",  value = "FREQ",  PRE_TREATMENT,  POST_TREATMENT)

control_stat.test <- control_data_stats %>%
    group_by(POS) %>%
    pairwise_t_test(
      FREQ ~ TREATMENT,  paired = TRUE,
      p.adjust.method = "bonferroni"
      ) %>%
    select(-df,  -statistic,  -p) # Remove details


control_stat.test$TREATMENT <- "1. Untreated"
control_stat.test <- control_stat.test %>%
  mutate(POS = str_replace(POS,  c("7029569", "7029790"),  c("Phe167Tyr", "Phe200Tyr")))

p.data <- dplyr::bind_rows(control_stat.test,  bz_stat.test)


# make new plot with p values annotated on it
plot_b <- plot +
     geom_text(data = p.data,  aes(x = 1.5,  y = 0.95,  group = POS,  label = paste('P = ', p.adj)), size = 2.5)
```

---

## Figure 3 C <a name="figure3c"></a>

Aim is the show correlation between Glu198Val variant in beta tubulin isotype 2 and Benzimidazole EC50

### Prepare the data
```shell
working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/US_FIELD/VCF

echo "hcontortus_chr2_Celeg_TT_arrow_pilon 13435823" > btub.positions

cat bam.list > samples.list

vcftools --vcf 2.hcontortus_chr2_Celeg_TT_arrow_pilon.snpeff.vcf --keep samples.list --positions btub.positions --extract-FORMAT-info AD --out us_farms_btub2

for i in `ls *AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10),$12/($11+$12),$14/($13+$14),$16/($15+$16),$18/($17+$18),$20/($19+$20),$22/($21+$22)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
done
```

### R to plot
```R
# load required libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)
library(ggrepel)

us_btub2 <- read.table("us_farms_btub2.ADfreq")
colnames(us_btub2) <- c("CHR", "POS", "Farm 1", "Farm 2", "Farm 3", "Farm 4", "Farm 5", "Farm 6", "Farm 7", "Farm 8", "Farm 9", "Farm 10")
us_btub2 <- melt(us_btub2,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")

bz_conc <- c(1, 19.93, 15, 7.71, 15, 61.57, 29.6, NA, 9.67, 29.6)

us_btub2$BZ_CONCENTRATION <- bz_conc
colnames(us_btub2) <- c("CHR", "POS", "SAMPLE_ID", "ALLELE_FREQ", "BZ_CONCENTRATION")

us_btub2 <- us_btub2 %>%
  mutate(POS = str_replace(POS,  c("13435823"),  c("Glu198Val")))

# calculate correlation coefficient between alllele frequency and concentration
af_bz_cor <- cor.test(us_btub2$ALLELE_FREQ,  us_btub2$BZ_CONCENTRATION,  method = "pearson",  use = "complete.obs")

# make the plot
plot_c <- ggplot(us_btub2)+
     geom_smooth(aes(BZ_CONCENTRATION, ALLELE_FREQ), method = 'lm', col = 'grey')+
     geom_jitter(aes(BZ_CONCENTRATION, ALLELE_FREQ, col = SAMPLE_ID), size = 2)+
     geom_text(aes(15, 0.95, label = paste('r = ', signif(af_bz_cor$estimate, 3), '\n', 'P = ', signif(af_bz_cor$p.value, 3))), size = 2.5)+
     geom_text_repel(aes(BZ_CONCENTRATION, ALLELE_FREQ, label = SAMPLE_ID, col = SAMPLE_ID), size = 3.5)+
     labs(title = "C", y = "Variant Allele Frequency", x = "Benzimidazole EC50 (uM)", col = "US farm ID") +
     ylim(-0.05, 1) +
     #facet_grid(. ~ POS) +
     facet_grid(POS ~ "US farm") +
     theme_bw()+
     theme(legend.position = "none", text = element_text(size = 10))
```

```R
# bring the panels together

library(patchwork)

# plot it
plot_a / (plot_b | plot_c)


# save it
ggsave("Figure_benzimidazole.pdf",  useDingbats = FALSE, width = 170, height = 130, units = "mm")
ggsave("Figure_benzimidazole.png")
```
![](04_analysis/Figure_benzimidazole.png)


******
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
