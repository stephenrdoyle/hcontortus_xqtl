# Supplementary Figures

1. [Figure S1](#Figure_S1)
2. [Figure S2](#Figure_S2)
3. [Figure S3](#Figure_S3)
4. [Figure S4](#Figure_S4)
5. [Figure S5](#Figure_S5)



#-----------------------------------------------------------------------------------------

#-----


library(ggplot2)
library(patchwork)

# chromosome colours
chr_colours <- c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

fst <-read.table("/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/US_FIELD/XQTL_US_FIELD.merged.fst",header=F)
fst <- fst[fst$V1!="hcontortus_chr_mtDNA_arrow_pilon",]

# susceptibles
plot_1 <- ggplot(fst,aes(1:nrow(fst)*10000,V21,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 9 (ivermectin susceptible)",  y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

# UGA_S vs moderates
plot_2 <- ggplot(fst,aes(1:nrow(fst)*10000,V7,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 2 (moderate ivermectin resistance)", y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_3 <- ggplot(fst,aes(1:nrow(fst)*10000,V9,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 3 (moderate ivermectin resistance)",   y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_4 <- ggplot(fst,aes(1:nrow(fst)*10000,V11,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 4 (moderate ivermectin resistance)",  y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())


# UGA_S vs high resistance
plot_5 <- ggplot(fst,aes(1:nrow(fst)*10000,V13,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 5 (high ivermectin resistance)",  y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_6 <- ggplot(fst,aes(1:nrow(fst)*10000,V15,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 6 (high ivermectin resistance)",  y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_7 <- ggplot(fst,aes(1:nrow(fst)*10000,V17,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 7 (high ivermectin resistance)",  y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_8 <- ggplot(fst,aes(1:nrow(fst)*10000,V19,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 8 (high ivermectin resistance)",   y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_9 <- ggplot(fst,aes(1:nrow(fst)*10000,V23,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 10 (high ivermectin resistance)", x = "Chromosomal position (bp)",  y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8))

plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 + plot_8 + plot_9 + plot_layout(ncol=1)

ggsave("XQTL_SupplementaryFigure_USfarm1.pdf", useDingbats = FALSE, width = 170, height = 250, units = "mm")
ggsave("XQTL_SupplementaryFigure_USfarm1.png")



# susceptibles
plot_1 <- ggplot(fst,aes(1:nrow(fst)*10000,V21,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 1 (ivermectin susceptible) vs Farm 9 (ivermectin susceptible)",  y = "")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

# UGA_S vs moderates
plot_2 <- ggplot(fst,aes(1:nrow(fst)*10000,V37,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 9 (ivermectin susceptible) vs Farm 2 (moderate ivermectin resistance)", y = "")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_3 <- ggplot(fst,aes(1:nrow(fst)*10000,V51,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 9 (ivermectin susceptible) vs Farm 3 (moderate ivermectin resistance)",   y = "")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_4 <- ggplot(fst,aes(1:nrow(fst)*10000,V63,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 9 (ivermectin susceptible) vs Farm 4 (moderate ivermectin resistance)",  y = "")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())


# UGA_S vs high resistance
plot_5 <- ggplot(fst,aes(1:nrow(fst)*10000,V73,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 9 (ivermectin susceptible) vs Farm 5 (high ivermectin resistance)",  y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_6 <- ggplot(fst,aes(1:nrow(fst)*10000,V81,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 9 (ivermectin susceptible) vs Farm 6 (high ivermectin resistance)",  y = "")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_7 <- ggplot(fst,aes(1:nrow(fst)*10000,V87,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 9 (ivermectin susceptible) vs Farm 7 (high ivermectin resistance)",  y = "")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_8 <- ggplot(fst,aes(1:nrow(fst)*10000,V91,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 9 (ivermectin susceptible) vs Farm 8 (high ivermectin resistance)",   y = "")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

plot_9 <- ggplot(fst,aes(1:nrow(fst)*10000,V95,colour = V1))+
     geom_point(size=0.5,alpha=0.5)+
     ylim(0,1)+
     labs(title = "Farm 9 (ivermectin susceptible) vs Farm 10 (high ivermectin resistance)", x = "Chromosomal position (bp)",  y = "Genetic differentiation (Fst)")+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+theme(legend.position = "none", text = element_text(size = 8))

plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7 + plot_8 + plot_9 + plot_layout(ncol=1)

ggsave("XQTL_SupplementaryFigure_USfarm2.pdf", useDingbats = FALSE, width = 170, height = 250, units = "mm")
ggsave("XQTL_SupplementaryFigure_USfarm2.png")




#---------------------------------------------------------------------------------------

## Figure S1 - isotype 1 data from US farms <a name="Figure_S1"></a>

Aim is to determine the frequency of Phe167Tyr & Phe200Tyr variants in the US farm dataset

## Prepare the data
```shell
# working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/US_FIELD/VCF

# file containing position of variants
cp ../../XQTL_BZ/btub.positions btub1.positions

# extract allele count data for each farm
vcftools --vcf 1.hcontortus_chr1_Celeg_TT_arrow_pilon.snpeff.vcf --keep samples.list --positions btub1.positions --extract-FORMAT-info AD --out us_farms_btub1

# convert allele count data to variant frequency
grep "^hcon" us_farms_btub1.AD.FORMAT | awk -F '[\t, ]' '{print $1, $2, $4/($3+$4), $6/($5+$6), $8/($7+$8), $10/($9+$10), $12/($11+$12), $14/($13+$14), $16/($15+$16), $18/($17+$18), $20/($19+$20), $22/($21+$22)}' OFS="\t" > us_farms_btub1.ADfreq
```

### R to plot
```R
# load libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(rstatix)

# reformat data
us_btub1 <- read.table("us_farms_btub1.ADfreq")
colnames(us_btub1) <- c("CHR", "POS", "Farm 1", "Farm 2", "Farm 3", "Farm 4", "Farm 5", "Farm 6", "Farm 7", "Farm 8", "Farm 9", "Farm 10")
us_btub1 <- melt(us_btub1,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")
colnames(us_btub1) <- c("CHR", "POS", "SAMPLE_ID", "ALLELE_FREQ")
us_btub1 <- us_btub1 %>%
  mutate(POS = str_replace(POS,  c("7029569", "7029790"),  c("Phe167Tyr", "Phe200Tyr")))

# make the figure
ggplot(us_btub1, aes(x = SAMPLE_ID, y = ALLELE_FREQ, fill = factor(POS))) +
     geom_bar(position = "dodge",  stat = "identity") +
     labs(title = "A",  x="Sampling location",  y="Resistant allele frequency",  fill = "Variant") +
     theme_bw() + theme(text = element_text(size = 10))

# save it
ggsave("FigureSX_USfarm_btub1.pdf",  useDingbats=FALSE, width=170, height=100, units="mm")
ggsave("FigureSX_USfarm_btub1.png")
```
![](04_analysis/FigureSX_USfarm_btub1.png)



#-----------------------------------------------------------------------------------------

## Figure S2 - panel a - evidence of ivermectin selection on btub variants <a name="Figure_S2"></a>

Aim is to show evidence (or lack thereof) of Phe167Tyr & Phe200Tyr variants in the ivermectin XQTL data.

### Prepare the data
```shell
#working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_IVM


grep "pre" bam.list  > ivm_pretreatment_samples.list
grep "post" bam.list  > ivm_posttreatment_samples.list

cp ../XQTL_BZ/btub.positions .

# extract allele count data for for beta tubulin variants in ivermectin XQTL experiment
vcftools --vcf XQTL_IVM.raw.snpeff.vcf --keep ivm_pretreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out ivm_btub_pretreatment
vcftools --vcf XQTL_IVM.raw.snpeff.vcf --keep ivm_posttreatment_samples.list --positions btub.positions --extract-FORMAT-info AD --out ivm_btub_posttreatment

# convert allele count data to variant frequency
for i in `ls ivm*AD.FORMAT`; do
      grep "^hcon" ${i} | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10)}' OFS="\t" > ${i%.AD.FORMAT}.ADfreq;
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

# reform the data
pre <- read.table("ivm_btub_pretreatment.ADfreq")
colnames(pre) <- c("CHR", "POS", "R1", "R1.2", "R2", "R3")
pre <- melt(pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

post <- read.table("ivm_btub_posttreatment.ADfreq")
colnames(post) <- c("CHR", "POS", "R1", "R1.2", "R2", "R3")
post <- melt(post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

data <- dplyr::full_join(pre, post, by = c("CHR", "POS", "SAMPLE_ID"))
data$TREATMENT <- "Ivermectin"
colnames(data) <- c("CHR", "POS", "SAMPLE_ID", "PRE_TREATMENT", "POST_TREATMENT", "TREATMENT")

# change the labels
data <- data %>%
  mutate(POS = str_replace(POS, c("7029569", "7029790"), c("Phe167Tyr", "Phe200Tyr")))


# make the plot
plot_a <- ggplot(data) +
     geom_segment(aes(x = "1.PRE", xend = "2.POST", y = PRE_TREATMENT, yend = POST_TREATMENT, col = factor(SAMPLE_ID), group = POS), size = 1) +
     labs(title = "A", x = "Sampling time-point", y = "Resistant allele frequency", col = "Replicate") +
     ylim(-0.05, 1.05) +
     facet_grid(TREATMENT ~ POS) +
     theme_bw() + theme(text = element_text(size = 10))



# perform pairwise t tests between pre/post for each SNP on BZ treated samples
data_stats <- data %>%
  gather(key = "TREATMENT", value = "FREQ", PRE_TREATMENT, POST_TREATMENT)

stat.test <- data_stats %>%
    group_by(POS) %>%
    pairwise_t_test(
      FREQ ~ TREATMENT, paired = TRUE,
      p.adjust.method = "bonferroni"
      ) %>%
    select(-df, -statistic, -p) # Remove details

stat.test$TREATMENT <- "Ivermectin"

p.data <- stat.test


# make new plot with p values annotated on it
plot_a <- plot_a +
     geom_text(data = p.data, aes(x = 1.5, y = 0.95, group = POS, label = paste('P = ',p.adj)), size = 3)

plot_a
```

## Figure SX - panel b - correlation between btubulin isotype 1 and ivermectin concentration

Aim is to show frequency of Phe167Tyr Phe200Tyr variants against EC50 for ivermectin in US farms

### Prepare the data
```shell
#working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/US_FIELD/VCF

cp ../../XQTL_BZ/btub.positions btub1.positions

cat bam.list > samples.list

vcftools --vcf 1.hcontortus_chr1_Celeg_TT_arrow_pilon.snpeff.vcf --keep samples.list --positions btub1.positions --extract-FORMAT-info AD --out us_farms_btub1

grep "^hcon" us_farms_btub1.AD.FORMAT | awk -F '[\t,]' '{print $1,$2,$4/($3+$4),$6/($5+$6),$8/($7+$8),$10/($9+$10),$12/($11+$12),$14/($13+$14),$16/($15+$16),$18/($17+$18),$20/($19+$20),$22/($21+$22)}' OFS="\t" > us_farms_btub1.ADfreq
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
library(patchwork)

# reformat the data
us_btub1 <- read.table("us_farms_btub1.ADfreq")
colnames(us_btub1) <- c("CHR", "POS", "Farm 1", "Farm 2", "Farm 3", "Farm 4", "Farm 5", "Farm 6", "Farm 7", "Farm 8", "Farm 9", "Farm 10")

us_btub1 <- melt(us_btub1,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")

# manually input the EC50 data from Table SX
ivm_conc <- c(1.51, 1.51, 12.04, 12.04, 9.15, 9.15, 11.27, 11.27, 297.9, 297.9, 619, 619, 312.1, 312.1, NA, NA, 0.977, 0.977, 259.4, 259.4)

us_btub1$IVM_CONCENTRATION <- ivm_conc
colnames(us_btub1) <- c("CHR", "POS", "SAMPLE_ID", "ALLELE_FREQ", "IVM_CONCENTRATION")

us_btub1 <- us_btub1 %>%
  mutate(POS = str_replace(POS,  c("7029569", "7029790"),  c("Phe167Tyr", "Phe200Tyr")))

# calculate correlation coefficient between alllele frequency and concentration
af_ivm_cor <- cor.test(us_btub1$ALLELE_FREQ,  us_btub1$IVM_CONCENTRATION,  method = "pearson",  use = "complete.obs")

P167 <- us_btub1[us_btub1$POS != "Phe167Tyr", ]
P200 <- us_btub1[us_btub1$POS != "Phe200Tyr", ]

P167cor <- cor.test(P167$ALLELE_FREQ,  P167$IVM_CONCENTRATION,  method = "pearson",  use = "complete.obs")
P200cor <- cor.test(P200$ALLELE_FREQ,  P200$IVM_CONCENTRATION,  method = "pearson",  use = "complete.obs")

P167cor.data <- as.data.frame(P167cor$estimate)
P167cor.data$pvalue <- P167cor$p.value
colnames(P167cor.data) <- c("COR", "PVALUE")

P200cor.data <- as.data.frame(P200cor$estimate)
P200cor.data$pvalue <- P200cor$p.value
colnames(P200cor.data) <- c("COR", "PVALUE")

cor.data <- dplyr::bind_rows(P167cor.data,  P200cor.data)
cor.data$CHR <- "hcontortus_chr1_Celeg_TT_arrow_pilon"
cor.data$POS <-  c("Phe167Tyr", "Phe200Tyr")

colnames(cor.data) <- c("COR", "PVALUE", "CHR", "POS")


# make the plot
plot_b <- ggplot(us_btub1) +
     geom_smooth(aes(IVM_CONCENTRATION, ALLELE_FREQ), method = 'lm', col = 'grey', se = FALSE) +
     geom_jitter(aes(IVM_CONCENTRATION, ALLELE_FREQ, col = SAMPLE_ID), size = 2) +
     geom_text_repel(aes(IVM_CONCENTRATION, ALLELE_FREQ, label = SAMPLE_ID, col = SAMPLE_ID), size = 2) +
     labs(title = "B", y = "Resistant Allele Frequency", x = "Ivermectin EC50 (nM)", col = "US farm ID") +
     ylim(-0.05, 1.05) +
     facet_grid(. ~ POS) +
     theme_bw() + theme(legend.position = "none", text = element_text(size = 10))

plot_b <- plot_b + geom_text(data = cor.data,  aes(x = 500,  y = 1,  group = POS,  label = paste('r = ', signif(COR, 3), '\n', 'P = ', signif(PVALUE, 3))), size = 3)

# combine panels a and b using patchwork
plot_a + plot_b + plot_layout(ncol = 2)

# save it
ggsave("FigureSX_USfarm_btub1vsIVM.pdf",  useDingbats = FALSE,  width = 170,  height = 100,  units = "mm")
ggsave("FigureSX_USfarm_btub1vsIVM.png")

```
![](04_analysis/FigureSX_USfarm_btub1vsIVM.png)



#-------------------------------------------------------------------------------
## Supplementary Figure - multiple seqeunce alignment of acr-8 showing Ser168Thr conservation <a name="Figure_S3"></a>

Aim is to show amino acid alignments of acr8 from all available clade V nematodes, and highlight conservation of the Serine at the Ser168Th variant position putatively linked to levamisole resistance

### Prepare the data  
```shell
# downloaded protein sequences from WBP frmo ortholog set of C. elegans acr-8

WBGene00036931|cabrigprjna10731  Caenorhabditis_briggsae
Cni-acr-8|canigoprjna384657 Caenorhabditis_nigoni
WBGene00058468|caremaprjna53967  Caenorhabditis_remanei
Csp11.Scaffold629.g13743|catropprjna53597  Caenorhabditis_tropicalis
ACAC_0000398801|ancantprjeb493  Angiostrongylus_cantonensis
ANCDUO_01121|anduodprjna72581   Ancylostoma_duodenale
ACOC_0000247401|ancostprjeb494  Angiostrongylus_costaricensis
Cang_2012_03_13_00065.g3190|caangaprjna51225   Caenorhabditis_angaria
arcanus-mkr-S_1590-0.0-mRNA-1|prarcaprjeb27334  Pristionchus_arcanus
Acey_s0072.g653|anceylprjna231479 Ancylostoma_ceylanicum
HPOL_0000787701|hepolyprjeb15396 Heligmosomoides_polygyrus
DCO_008622|dicoroprjdb3143  Diploscapter_coronatus
FL83_01328|calateprjna248912 Caenorhabditis_latens
HPLM_0000547301|haplacprjeb509  Haemonchus_placei
fissidentatus-mkr-S40-3.15-mRNA-1|prfissprjeb27334  Pristionchus_fissidentatus
maxplancki-mkr-S100-2.104-mRNA-1|prmaxpprjeb27334   Pristionchus_maxplancki
HCON_00151270|hacontprjeb506 Haemonchus_contortus
japonicus-mkr-S267-0.26-mRNA-1|prjapoprjeb27334 Pristionchus_japonicus
OESDEN_01101|oedentprjna72579   Oesophagostomum_dentatum
MicoRS5524-ag_msk-S79-2.36-mRNA-1|mijapoprjeb27334  Micoletzkya_japonica
WBGene00000047|caelegprjna13758  Caenorhabditis_elegans
WBGene00092568|prpaciprjna12644  Pristionchus_pacificus
SVUK_0001037501|stvulgprjeb531  Strongylus_vulgaris
WR25_23289|dipachprjna280107 Diploscapter_pachys
Sp34_X0071520|casp34prjdb5687   Caenorhabditis_inopinata
NECAME_16202|neamerprjna72135   Necator_americanus
OTIPU.nOt.2.0.1.g01642|ostipuprjeb15512   Oscheius_tipulae
nAv.1.0.1.g09901|acviteprjeb1697 Acanthocheilonema_viteae
mayeri-mkr-S55-0.24-mRNA-1|prmayeprjeb27334 Pristionchus_mayeri
NBR_0000568801|nibrasprjeb511   Nippostrongylus_brasiliensis

# removed - poor / missing alignment in region of interest
entomophagus-mkr-S570-0.22-mRNA-1|prentoprjeb27334  Pristionchus_entomophagus
exspectatus-mkr-S_1791-0.14-mRNA-1|prexspprjeb24288  Pristionchus_exspectatus
CGOC_0000275601|cygoldprjeb498  Cylicostephanus_goldi
DICVIV_07636|diviviprjna72587   Dictyocaulus_viviparus
ANCCAN_17789|ancaniprjna72585   Ancylostoma_caninum
Csp5_scaffold_00005.g324|casiniprjna194557 Caenorhabditis_sinica
WBGene00139702|cabrenprjna20035  Caenorhabditis_brenneri
mbelari.g26242|mebelaprjeb30104  Mesorhabditis_belari


# substitute names from transcript ID to species name
while read old new; do sed -i "s/${old}/${new}/" wbp_acr8.fa; done < rename

# extract fasta sequences based on the new names to a new fasta
while read old new; do samtools faidx wbp_acr8.fa ${new}; done < rename > wb_cladeV_acr8.fa

# make an alignment using mafft
module load mafft/7.407=1

mafft --localpair --maxiterate 16 --reorder "wb_cladeV_acr8.fa" > "wb_cladeV_acr8.aln"
```

### R to plot
```R
# load required libraries
library(ggmsa)
library(ggplot2)

# define the alignment window - chosen to flank the Ser168Thr position in H. contortus
ALIGNMENT_START = 225
ALIGNMENT_END = 275

# load data and clean up for plotting
data <- tidy_msa("wb_cladeV_acr8.aln", ALIGNMENT_START, ALIGNMENT_END)


# colour scheme - based on amino acid properties
<!-- Aliphatic #CD2127  A
Aliphatic #CD2127  G
Aliphatic #CD2127  I
Aliphatic #CD2127  L
Aliphatic #CD2127  P
Aliphatic #CD2127  V
Aromatic #AAC64B  F
Aromatic #AAC64B  W
Aromatic #AAC64B  Y
Acidic  #ED6823  D
Acidic  #ED6823  E
Basic   #759CD1  R
Basic   #759CD1  H
Basic   #759CD1  K
Hydroxylic   #F09AC1  S
Hydroxylic   #F09AC1  T
Sulfur-containing  #EDAB20  C
Sulfur-containing  #EDAB20  M
Amidic  #283983  N
Amidic  #283983  Q -->

colours <- c("A" = "#CD2127", "G" = "#CD2127", "I" = "#CD2127", "L" = "#CD2127", "P" = "#CD2127", "V" = "#CD2127", "F" = "#AAC64B", "W" = "#AAC64B", "Y" = "#AAC64B", "D" = "#ED6823", "E" = "#ED6823", "R" = "#759CD1", "H" = "#759CD1", "K" = "#759CD1", "S" = "#F09AC1", "T" = "#F09AC1", "C" = "#EDAB20", "M" = "#EDAB20", "N" = "#283983", "Q" = "#283983", "-" = "white")

# make the plot
ggplot(data, aes(y = name, x = position)) +
   geom_tile(aes(fill = character), col = "grey", na.rm = TRUE) +
   geom_text(aes(label = character), col = "black", size = 2) +
   geom_text(aes(x = 253, y = 31), label = "*") +
   geom_text(aes(x = 253, y = 29.5), label = "H. contortus\nSer168Thr", size = 3) +
   theme_minimal() + theme(legend.position = "none")+
   labs(y = "Species", x = "Alignment position", text = element_text(size = 10)) +
   scale_fill_manual(values = colours)

# save it
ggsave("acr8_multiple_sequence_alignment.pdf", width = 170, height = 180, units = "mm")
ggsave("acr8_multiple_sequence_alignment.png")
```
![](04_analysis/acr8_multiple_sequence_alignment.png)
#-------------------------------------------------------------------------------



5000 L3 per pool example
#-----
library(ggplot2)
library(patchwork)



data <- read.table("/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_L3_5k/XQTL_L3_5k.merged.fst",header=F)
data <- L3_5000[L3_5000$V1!="hcontortus_chr_mtDNA_arrow_pilon",]

# chromosome colours
chr_colours <- c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")


plot_a <- ggplot(data)+
     geom_point(aes(1:nrow(data) * 5000, V7, alpha = V4, colour = V1),size = 0.1)+
     ylim(0, 0.1) +
     labs(title = "A", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+ theme(legend.position = "none", text = element_text(size = 10))


data_chr5 <- data[data$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon",]

plot_b <- ggplot(data_chr5)+
     geom_point(aes(V2,V7),size = 0.5, colour = "cornflowerblue")+
  ylim(0,0.1)+
  labs(title = "B",  x = "Genomic position (bp)",  y = "Genetic differentiation (Fst)") +
  scale_x_continuous() +
  xlim(36.5e6, 38.5e6) +
  theme_bw() + theme(legend.position = "none", text = element_text(size = 10))

plot_a + plot_b + plot_layout(ncol = 1)

ggsave("XQTL_SupplementaryFigure_5000L3.pdf", useDingbats = FALSE, width = 170, height = 100, units = "mm")
ggsave("XQTL_SupplementaryFigure_5000L3.png")






#-------------------------------------------------------------------------------

#-----
```R
#XQTL F2 adults

library(ggplot2)
library(patchwork)
library(dplyr)

# chromosome colours
chr_colours <- c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

# load data1
male_data <- read.table("/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_ADULT/XQTL_ADULTS.merged.fst",header=F)
male_data <- male_data[male_data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]

plot_a <- ggplot(male_data)+
     geom_point(aes(1:nrow(male_data)*5000, V7, colour = V1), size=0.1)+
     #ylim(0,0.2)+
     labs(title = "A", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \nsurviving adult males and untreated L3 (Fst)") +
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+
     theme(legend.position = "none", text = element_text(size = 8))


chr5_maledata <- filter(male_data,male_data$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon")
plot_b <- ggplot(chr5_maledata)+
     geom_point(aes(V2, V7, colour = V1), size=1)+
     #ylim(0,0.2)+
     xlim(36.5e6,38.5e6)+
     labs(title = "B", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \nsurviving adult males and untreated L3 (Fst)") +
     scale_color_manual(values = chr_colours) +
     theme_bw()+
     theme(legend.position = "none", text = element_text(size = 8))

# dose response curve
dose_response_data <- read.table("/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/IVM/DrenchRiteF5_IVM_pooled.txt", header=T)
ivm_dose_response_data <- dplyr::filter(dose_response_data,dose_response_data$Treatment=="IVM_aglycone" & dose_response_data$Total>0)


#probit model
#logr_probit <- glm(IVM_logged ~ (X._developed_to_L3/100), data=ivm_dose_response_data, family=binomial(link="probit"))
#preddat <- data.frame(IVM_logged = seq(0, 3.5, .01) )
#preddat$pred <- predict(logr_probit, newdata = preddat, type = "response")







plot_c <- ggplot(ivm_dose_response_data,aes(IVM_logged,X._developed_to_L3/100)) +
     geom_point(aes(color=as.factor(Plate)))+
     #geom_line(data = preddat, aes(y = pred), color = "red")+
     geom_smooth(method="glm",method.args = list(family = quasibinomial(link = "probit")),colour="blue")+
     labs(title = "C", y = "Proportion developing to L3", x = "Ivermectin concentration [log10(nM)]", colour = "Replicate") +
     theme_bw()+
     theme(text = element_text(size = 8))


# genome wide dose response

dose_data <- read.table("/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/DOSE_RESPONSE/DOSE_RESPONSE.merged.fst",header=F)
dose_data <- dose_data[dose_data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]

plot_d <- ggplot(dose_data)+
     geom_point(aes(1:nrow(dose_data)*5000, V7, colour = V1), size=0.1)+
     #ylim(0,0.2)+
     labs(title = "D", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \n<EC25 and >EC75 L3 pools (Fst)") +
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw()+
     theme(legend.position = "none", text = element_text(size = 8))


     # chr5_dose_data <- filter(dose_data,dose_data$V1=="hcontortus_chr5_Celeg_TT_arrow_pilon")
     # ggplot(chr5_dose_data)+
     #      geom_point(aes(V2, V7, colour = V1), size=1)+
     #      #ylim(0,0.2)+
     #      xlim(36.5e6,38.5e6)+
     #      labs(title = "B", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
     #      scale_color_manual(values = chr_colours) +
     #      theme_bw()+
     #      theme(legend.position = "none", text = element_text(size = 8))


plot_a + (plot_b | plot_c) + plot_d + plot_layout(ncol=1)

ggsave("XQTL_SupplementaryFigureX_male_doseresponse.pdf", useDingbats = FALSE, width = 170, height = 200, units = "mm")
ggsave("XQTL_SupplementaryFigureX_male_doseresponse.png")
```
******
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
