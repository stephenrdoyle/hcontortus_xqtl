# Levamisole analyses


## Figure 4 A <a name="figure4a"></a>

Aim is show genetic differentiation between rpe and post treatment levamisole samples on chromosomes 4 and 5

```shell
# working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/LEV
```

### R to plot
```R
# load required libraries
library(ggplot2)
library(dplyr)
library(stringr)

# load and reformat the data
xqtl_lev_fst <- read.table("XQTL_LEV.merged.fst", header = F)
xqtl_lev_fst <- xqtl_lev_fst[xqtl_lev_fst$V1! = "hcontortus_chr_mtDNA_arrow_pilon", ]
xqtl_lev_fst <- dplyr::select(xqtl_lev_fst,  V1,  V2,  V13 , V39, V49)
xqtl_lev_fst <- xqtl_lev_fst %>% mutate(mean_FST = rowMeans(select(.,V13)))
colnames(xqtl_lev_fst) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN")

# calculate a genome wide significance cutoff
gw_r1 <- mean(xqtl_lev_fst$FST_R1) + 3*sd(xqtl_lev_fst$FST_R1)
gw_r2 <- mean(xqtl_lev_fst$FST_R2) + 3*sd(xqtl_lev_fst$FST_R2)
gw_r3 <- mean(xqtl_lev_fst$FST_R3) + 3*sd(xqtl_lev_fst$FST_R3)
gw_mean <- mean(xqtl_lev_fst$FST_MEAN) + 3*sd(xqtl_lev_fst$FST_MEAN)


# extract chromosome 4 and 5 data
lev_chr4 <- xqtl_lev_fst[xqtl_lev_fst$CHR == "hcontortus_chr4_Celeg_TT_arrow_pilon", ]
lev_chr5 <- xqtl_lev_fst[xqtl_lev_fst$CHR == "hcontortus_chr5_Celeg_TT_arrow_pilon", ]

# fix names
lev_chr45_data <- dplyr::bind_rows(lev_chr4, lev_chr5)
lev_chr45_data <- lev_chr45_data %>%
 mutate(CHR = str_replace_all(CHR, c("hcontortus_chr4_Celeg_TT_arrow_pilon" = "Chromosome 4", "hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))


lev_chr45_data <- mutate(lev_chr45_data,
        point_colour = case_when(
        (FST_R1 < gw_r1 & FST_R2 < gw_r2 & FST_R3 < gw_r3) ~ "#6699FFFF",
        (FST_R1 > gw_r1 & FST_R2 < gw_r2 & FST_R3 < gw_r3) ~ "#FFCC00FF",
        (FST_R1 < gw_r1 & FST_R2 > gw_r2 & FST_R3 < gw_r3) ~ "#FFCC00FF",
        (FST_R1 < gw_r1 & FST_R2 < gw_r2 & FST_R3 > gw_r3) ~ "#FFCC00FF",
        (FST_R1 > gw_r1 & FST_R2 > gw_r2 & FST_R3 < gw_r3) ~ "#FF9900FF",
        (FST_R1 > gw_r1 & FST_R2 < gw_r2 & FST_R3 > gw_r3) ~ "#FF9900FF",
        (FST_R1 > gw_r1 & FST_R2 > gw_r2 & FST_R3 > gw_r3) ~ "#FF0000FF",))


# make the plot
plot_a <- ggplot(lev_chr45_data, aes(POS, FST_MEAN, colour=point_colour, size=ifelse(FST_MEAN>gw_mean,0.6,0.3))) +
               geom_hline(yintercept = gw_mean, linetype = "dashed", col = "black") +
               geom_point() + facet_grid(CHR~.) + scale_color_identity() + scale_size_identity() +
               xlim(0, 50e6) +
               theme_bw() + theme(legend.position = "none", text = element_text(size = 10)) +
               labs(title = "A", x = "Genomic position (bp)", y = expression(paste("Genetic differentiation between pre- and post-treatment", " (",~italic(F)[ST],")"))) +
               facet_grid(CHR ~ .)


#plot_a
```

## Figure 4 B <a name = "figure4b"></a>

Aim is to show a gene model of acr-8, highlighing the position of the lev-resistance associated indel and the new Ser168Thr variant

```shell
# working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/LEV

#gff
ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/20200407/UPDATED_annotation.gff3 ANNOTATION.gff
```


### R to plot
```R
# load required R libraries
library(data.table)
# note - "data.table" has a function called "fread" which is great for quickly loading really large datasets
library(ggplot2)
library(dplyr)

# load data
# --- genome annotation
gff <- fread(cmd = "grep ^hcontortus ANNOTATION.gff")
colnames(gff) <- c("chr", "source", "feature", "start", "end", "point1", "strand", "frame", "info")


# select gene ID - ACR-8 and the cuticle collegen gene that sits in the intron of acr-8
gene1 = 'HCON_00151270'
gene2 = 'HCON_00151260'


# filter data to select chromosome and mRNA
mrna1_data <- gff[grep(gene1, gff$info), ]
mrna1_data <- mrna1_data[mrna1_data$feature == 'mRNA', ]
mrna1_data <- cbind(mrna1_data, read.table(text = as.character(mrna1_data$info), sep = ";"))
mrna1_data <- cbind(mrna1_data, read.table(text = as.character(mrna1_data$V3), sep = "=", col.names = c("ID", "unique_ID")))
mrna1_id <- head(data.frame(mrna1_data$unique_ID), 1) # gives 1st isoform if multiple
chromosome <- mrna1_data[1, 1]
colnames(chromosome)<-c("chromosome_ID")

data1 <- gff[grep(mrna1_id$mrna1_data.unique_ID, gff$info), ]

# filter by feature type
cds1 <- data1[data1$feature == "CDS", ]
mrna1 <- data1[data1$feature == "mRNA", ]

# generate intron and utr data
intron1 <- data.frame(head(cds1$end, -1), tail(cds1$start, -1), (tail(cds1$start, -1)-head(cds1$end, -1))/2)
colnames(intron1) <- c("start", "end", "midpoint")

utr5 <- data.frame(head(mrna1$start, 1), head(sort(cds1$start), 1))
colnames(utr5) <- c("start", "end")
utr3 <- data.frame(head(mrna1$end, 1), tail(sort(cds1$end), 1))
colnames(utr3) <- c("start", "end")
utr1 <- rbind(utr5, utr3)

# repeat above now for gene2
mrna2_data <- gff[grep(gene2, gff$info), ]
mrna2_data <- mrna2_data[mrna2_data$feature == 'mRNA', ]
mrna2_data <- cbind(mrna2_data, read.table(text = as.character(mrna2_data$info), sep = ";"))
mrna2_data <- cbind(mrna2_data, read.table(text = as.character(mrna2_data$V3), sep = "=", col.names = c("ID", "unique_ID")))
mrna2_id <- head(data.frame(mrna2_data$unique_ID), 1) # gives 1st isoform if multiple
chromosome <- mrna2_data[1, 1]
colnames(chromosome) <- c("chromosome_ID")

data2 <- gff[grep(mrna2_id$mrna2_data.unique_ID, gff$info), ]

cds2 <- data2[data2$feature == "CDS", ]
mrna2 <- data2[data2$feature == "mRNA", ]

intron2 <- data.frame(head(cds2$end, -1), tail(cds2$start, -1), (tail(cds2$start, -1)-head(cds2$end, -1))/2)
colnames(intron2) <- c("start", "end", "midpoint")

utr52 <- data.frame(head(mrna2$start, 1), head(sort(cds2$start), 1))
colnames(utr52) <- c("start", "end")
utr32 <- data.frame(head(mrna2$end, 1), tail(sort(cds2$end), 1))
colnames(utr32) <- c("start", "end")
utr2 <- rbind(utr52, utr32)


# indel coordinates (size)
# 31527022 to 31527119 (97 bp) / 31527121 (99 bp)


# make plot
plot_b <- ggplot()+
#gene1
 geom_rect(data = utr1, aes(xmin = utr1$start, ymin = 0.5, xmax = utr1$end, ymax = 1.5), fill = NA, col = "grey", size = 0.4) +
 geom_segment(data = intron1, aes(x = intron1$start, xend = intron1$start + intron1$midpoint, y = 1, yend = 1.5), size = 0.5) +
 geom_segment(data = intron1, aes(x = intron1$start + intron1$midpoint, xend = intron1$end, y = 1.5, yend = 1), size = 0.5) +
 geom_rect(data = cds1, aes(xmin = cds1$start, ymin = 0.5, xmax = cds1$end, ymax = 1.5), fill = "grey", col = NA) +
 geom_text(aes(x = mrna1$end+(0.15*(mrna1$end-mrna1$start)), y = 1, label = gene1),size = 2) +
 geom_segment(aes(x = 31521884, xend = 31521884, y = 0.5, yend = 1.5), size = 1, col = "orange") + # Ser168Thr
 geom_text(aes(x = 31521884, y = 0.35), label = "Ser168Thr", size = 2) +
#gene2
 geom_rect(data = utr2, aes(xmin = utr2$start, ymin = 2, xmax = utr2$end, ymax = 3), fill = NA, col = "grey", size = 0.4) +
 geom_segment(data = intron2, aes(x = intron2$start, xend = intron2$start+intron2$midpoint, y = 2.5, yend = 3), size = 0.5) +
 geom_segment(data = intron2, aes(x = intron2$start+intron2$midpoint, xend = intron2$end, y = 3, yend = 2.5), size = 0.5) +
 geom_rect(data = cds2, aes(xmin = cds2$start, ymin = 2, xmax = cds2$end, ymax = 3), fill = "grey", col = NA) +
 geom_text(aes(x = mrna1$end+(0.15*(mrna1$end-mrna1$start)), y = 2.5, label = gene2), size = 2) +
 # acr-8 indel
 geom_segment(aes(x = 31527022, xend = 31527022, y = 0.5, yend = 3), size = 1, col = "red") +
 geom_text(aes(x = 31527022, y = 0.35), label = "indel", size = 2) +
 # plot layout
 theme_classic()+
 #xlab("Genome position (bp)")+
 labs(title = "B", x = paste("Chromosome 5 position (bp)"))+
 xlim(mrna1$start - (0.1 * (mrna1$end - mrna1$start)), mrna1$end + (0.25 * (mrna1$end-mrna1$start))) +
 scale_y_reverse(lim = c(3, 0.35)) +
 scale_fill_discrete(guide = FALSE) +
 theme(axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = 10))

#plot_b
```

## Figure 4 C <a name = "figure4c"></a>

Aim is to show a zoomed in view of the acr-8 indel, showing sequnecing reads spanning the deletion and that they have been soft clipped when mapped because of the indel. Two slightly different deletion coordinates are found, from: 31527022 to 31527119 (97 bp) and 31527022 to 31527121 (99 bp).

### R to plot
```R
# load required R libraries
library(GenomicAlignments)
library(viridis)
library(stringr)

# load in the required data and reformat
pre_bam <- "pre.bam"
#pre_bam <- "XQTL_F3_L3_n200_LEV_pre_03_23241_8_1.merged.sorted.marked.realigned.bam"

data_pre <- as.data.frame(readGAlignmentPairs(pre_bam, use.names = TRUE, param = ScanBamParam(which = GRanges("hcontortus_chr5_Celeg_TT_arrow_pilon", IRanges(31525841, 31529149)))))
data_pre$number <- rep(1:100, length.out = nrow(data_pre))
data_pre$treatment <- "1. Pre-treatment"

post_bam <- "post.bam"
#post_bam <- "XQTL_F3_L3_n200_LEV_post_03_23241_8_2.merged.sorted.marked.realigned.bam"

data_post <- as.data.frame(readGAlignmentPairs(post_bam, use.names = TRUE, param = ScanBamParam(which = GRanges("hcontortus_chr5_Celeg_TT_arrow_pilon", IRanges(31525841, 31529149)))))
data_post$number <- rep(1:100, length.out = nrow(data_post))
data_post$treatment <- "2. Levamisole treated"

#data_ise <- as.data.frame(readGAlignmentPairs("ISE.bam", use.names = TRUE, param = ScanBamParam(which = GRanges("hcontortus_chr5_Celeg_TT_arrow_pilon", IRanges(31525841, 31529149)))))
#data_ise$number <- rep(1:100, length.out = nrow(data_ise))
#data_ise$treatment <- "1.ISE"

#data_uga <- as.data.frame(readGAlignmentPairs("UGA.bam", use.names = TRUE, param = ScanBamParam(which = GRanges("hcontortus_chr5_Celeg_TT_arrow_pilon", IRanges(31525841, 31529149)))))
#data_uga$number <- rep(1:100, length.out = nrow(data_uga))
#data_uga$treatment <- "2.UGA"

#data <- dplyr::bind_rows(data_ise, data_uga, data_pre, data_post)
data <- dplyr::bind_rows(data_pre, data_post)


# find reads that have been softclipped, and determine the length of the softclip
softclip_match <- "([:digit:]{1,})S"

# for read 1
cigars1 <- data$cigar.first
softclip1_data <- str_extract(cigars1, softclip_match)
softclip1_data <- gsub("S", "", softclip1_data)
softclip1_data[is.na(softclip1_data)] <- 0
softclip1_data <- as.numeric(unlist(softclip1_data))

data$softclip1_length <- softclip1_data

# for read 2
cigars2 <- data$cigar.last
softclip2_data <- str_extract(cigars2, softclip_match)
softclip2_data <- gsub("S", "", softclip2_data)
softclip2_data[is.na(softclip2_data)] <- 0
softclip2_data <- as.numeric(unlist(softclip2_data))

data$softclip2_length <- softclip2_data


# make the plot
plot_c <- ggplot(data) +
   geom_rect(aes(xmin = start.first, ymin = number-0.4, xmax = end.first, ymax = number + 0.4, fill = softclip1_length)) +
   geom_rect(aes(xmin = start.last, ymin = number-0.4, xmax = end.last, ymax = number + 0.4, fill = softclip2_length)) +
   # if paired end links are required, could add the following. Gets a bit messy at this scale so have removed
   #geom_curve(aes(x = start.first, xend = end.last+1, y = number, yend = number+0.1), curvature = -0.05, col = "grey")+
   #xlim(31526800, 31527250)+
   coord_cartesian(xlim = c(31526850, 31527250)) +
   labs(title = "C", x = "Chromosomal position (bp)") +
   theme_bw() + theme(legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6),
   axis.title.y = element_blank(),
   axis.text.y = element_blank(),
   axis.ticks.y = element_blank()) +
   facet_grid(treatment ~ .) +
   scale_fill_viridis(option = "plasma")

#plot_c

```


## Figure 4 D

Aim is to show the allele frequency change at position Ser168Thr pre and post treatment in levamisole treated relative to control

### Generate the required data
```shell
#working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV


grep "pre" bam.list > lev_pretreatment_samples.list
grep "post" bam.list > lev_posttreatment_samples.list

grep "pre" ../XQTL_CONTROL/bam.list > control_pretreatment_samples.list
grep "post" ../XQTL_CONTROL/bam.list > control_posttreatment_samples.list

vcftools --vcf XQTL_LEV.raw.snpeff.vcf --keep lev_pretreatment_samples.list --positions acr8.positions --extract-FORMAT-info AD --out lev_pretreatment
vcftools --vcf XQTL_LEV.raw.snpeff.vcf --keep lev_posttreatment_samples.list --positions acr8.positions --extract-FORMAT-info AD --out lev_posttreatment

vcftools --gzvcf ../XQTL_CONTROL/5.hcontortus_chr5_Celeg_TT_arrow_pilon.tmp.vcf.gz --keep control_pretreatment_samples.list --positions acr8.positions --extract-FORMAT-info AD --out control_pretreatment
vcftools --gzvcf ../XQTL_CONTROL/5.hcontortus_chr5_Celeg_TT_arrow_pilon.tmp.vcf.gz --keep control_posttreatment_samples.list --positions acr8.positions --extract-FORMAT-info AD --out control_posttreatment

#where "acr8.positions" contains:
#P 31521884
#hcontortus_chr5_Celeg_TT_arrow_pilon 31521884

# convert allele counts to varinat frequency
for i in `ls lev*AD.FORMAT`; do
   grep "^hcon" ${i} | awk -F '[\t, ]' '{print $1, $2, $4/($3+$4), $6/($5+$6), $8/($7+$8), $10/($9+$10)}' OFS = "\t" > ${i%.AD.FORMAT}.ADfreq;
done

for i in `ls control*AD.FORMAT`; do
   grep "^hcon" ${i} | awk -F '[\t, ]' '{print $1, $2, $4/($3+$4), $6/($5+$6), $8/($7+$8)}' OFS = "\t" > ${i%.AD.FORMAT}.ADfreq;
done

cd ~/lustre118_link/hc/XQTL/05_ANALYSIS/LEV

ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV/lev_pretreatment.ADfreq
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV/lev_posttreatment.ADfreq
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV/control_pretreatment.ADfreq
ln -s /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_LEV/control_posttreatment.ADfreq
```

```R
# load the required R libraries
library(tidyverse)
library(rstatix)
library(reshape2)

# read data and reformat
lev_pre <- read.table("lev_pretreatment.ADfreq")
colnames(lev_pre) <- c("CHR", "POS", "R1", "R1.2", "R2", "R3")
lev_pre <- melt(lev_pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

lev_post <- read.table("lev_posttreatment.ADfreq")
colnames(lev_post) <- c("CHR", "POS", "R1", "R1.2", "R2", "R3")
lev_post <- melt(lev_post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

lev_data <- dplyr::full_join(lev_pre, lev_post, by = c("CHR", "POS", "SAMPLE_ID"))
lev_data$TREATMENT <- "Levamisole treated"
colnames(lev_data) <- c("CHR", "POS", "SAMPLE_ID", "PRE_TREATMENT", "POST_TREATMENT", "TREATMENT")

control_pre <- read.table("control_pretreatment.ADfreq")
colnames(control_pre) <- c("CHR", "POS", "R1", "R2", "R3")
control_pre <- melt(control_pre, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

control_post <- read.table("control_posttreatment.ADfreq")
colnames(control_post) <- c("CHR", "POS", "R1", "R2", "R3")
control_post <- melt(control_post, id = c("CHR", "POS"), variable.name = "SAMPLE_ID")

control_data <- dplyr::full_join(control_pre, control_post, by = c("CHR", "POS", "SAMPLE_ID"))
control_data$TREATMENT <- "Untreated control"
colnames(control_data) <- c("CHR", "POS", "SAMPLE_ID", "PRE_TREATMENT", "POST_TREATMENT", "TREATMENT")


# bring datasets together and fix the labels
data <- dplyr::bind_rows(control_data, lev_data)

data <- data %>%
 mutate(POS = str_replace_all(POS, c("31521884" = "Ser168Thr")))



# make the plot
plot_d <- ggplot(data) +
   geom_segment(aes(x = "1.Pre-treatment", xend = "2.Post-treatment", y = PRE_TREATMENT, yend = POST_TREATMENT, col = factor(TREATMENT), group = POS), size = 1) +
   labs(title = "D", x = "Sampling time-point", y = "Variant allele frequency", col = "Replicate") +
   ylim(0, 1) +
   theme_bw() + theme(text = element_text(size = 10), legend.position = "none")
   facet_grid(. ~ "Ser168Thr")

# perform pairwise t tests between pre/post for each SNP on BZ treated samples
lev_data_stats <- lev_data %>%
 gather(key = "TREATMENT", value = "FREQ", PRE_TREATMENT, POST_TREATMENT)

lev_stat.test <- lev_data_stats %>%
  group_by(POS) %>%
  pairwise_t_test(
   FREQ ~ TREATMENT, paired = TRUE,
   p.adjust.method = "bonferroni"
   ) %>%
  select(-df, -statistic, -p) # Remove details

lev_stat.test$TREATMENT <- "Levamisole treated"
lev_stat.test <- lev_stat.test %>%
 mutate(POS = str_replace_all(POS, c("31521884" = "Ser168Thr")))

# perform pairwise t tests between pre/post for each SNP on control samples
control_data_stats <- control_data %>%
 gather(key = "TREATMENT", value = "FREQ", PRE_TREATMENT, POST_TREATMENT)

control_stat.test <- control_data_stats %>%
  group_by(POS) %>%
  pairwise_t_test(
   FREQ ~ TREATMENT, paired = TRUE,
   p.adjust.method = "bonferroni"
   ) %>%
  select(-df, -statistic, -p) # Remove details


control_stat.test$TREATMENT <- "Untreated control"
control_stat.test <- control_stat.test %>%
 mutate(POS = str_replace_all(POS, c("31521884" = "Ser168Thr")))

p.data <- dplyr::bind_rows(control_stat.test, lev_stat.test)


# make new plot with p values annotated on it
plot_d <- plot_d +
   geom_text(data = p.data, aes(x = 1.5, y = 0.90, label = paste('P = ', p.adj[2])), size = 3, col = "red")+
   geom_text(data = p.data, aes(x = 1.5, y = 0.15, label = paste('P = ', p.adj[1])), size = 3, col = "blue")

#plot_d
```

```R
# bring the panels together in a single figure

library(patchwork)

# temporary blank plot
plot_e <- ggplot() + geom_blank()


plot_a / (plot_b | plot_d) / (plot_c | plot_e)  + plot_layout(ncol = 1, height = c(4, 1.5, 3))


ggsave("XQTL_Figure_4.pdf", useDingbats = FALSE, width = 170, height = 200, units = "mm")
ggsave("XQTL_Figure_4.png")
```
![](04_analysis/XQTL_Figure_4.png)








### INCOMPLETE - get frequency of acr8 variant from global dataset
```shell
#-----
grep "pre" bam.list > lev_pretreatment_samples.list
grep "post" bam.list > lev_posttreatment_samples.list

grep "pre" ../XQTL_CONTROL/bam.list > control_pretreatment_samples.list
grep "post" ../XQTL_CONTROL/bam.list > control_posttreatment_samples.list

for i in $(ls *.list); do
vcftools --gzvcf 5.hcontortus_chr5_Celeg_TT_arrow_pilon.cohort.vcf.gz --keep ${i} --positions acr8.positions --freq --out global_${i}; done
```



# Analysis of acr8 protein structure

Analysis is focused on the Ser168Thr variant that has a strong correlation with levamisole resistance




>Haemonchus_contortus
MRAFGIVLVTIVVYIVESNRYAEQLYEDLLYYYNKNVRPVKNASESIKVKFGASLIRIID
VDEVNQVLTTNLWLEMQWFDYKLAWEPTRWGGIKKLHIPSDQIWIPDILLYNNADGEPHI
TIMSDALVYHNGLIVWKPPSIYKSFCPINIEYFPYDSQSCSLKFGGWSYNGFLLDVRQLP
SADSPIINREDESGKEYQYLEKGMDLSGFYQSLEWDLMSLTSQRHEQLYPGCCGQDFYID
VTFDIALRRKTLFYTVNLVIPCMLIAILTTFVFYIPPIEHKMTFSISILVSLTVFYLVLI
ELIPPTSLVIPLIGRYLLFTMFLVAMSILLSVLSLNYYRRDGSAHPMPHWMQAVFIETLP
KYIGIKKAEEDEISDRGSSTSEIGQLLDSSRRPSPYFLTVPTSIQDDGQKRTRGEINKMR
LSQLAQLRGMHPDLIRRMIDNVSFIADHFRAMKKEDKISSDWSYVANVIDRLVLIIFTFV
NLVGTVLIIQNSPMLFDHTEPMKIGSASRPLSGDTFESTFDANFTQESWWKNSEGL*

# sequence search at Protein Data Bank (PDB)
https://www.rcsb.org/search/advanced/sequence



- top hit : 4AQ9 (https://www.rcsb.org/3d-view/4AQ9)

Gating movement in acetylcholine receptor analysed by time- resolved electron cryo-microscopy (open class)
Unwin, N., Fujiyoshi, Y.

(2012) J Mol Biol 422: 617

Released	2012-08-01
Method
ELECTRON MICROSCOPY 6.2 Å
Organisms
Torpedo marmorata
Macromolecule
ACETYLCHOLINE RECEPTOR BETA SUBUNIT (protein)
ACETYLCHOLINE RECEPTOR DELTA SUBUNIT (protein)
ACETYLCHOLINE RECEPTOR GAMMA SUBUNIT (protein)
ACETYLCHOLINE RECEPTOR SUBUNIT ALPHA (protein)




- using the "3D view", it is possible to find the homology position of Ser168Thr at position 174 which is a Thr residue in the alpha subunit of the pentamer
- The W[S/T]Y is clearly highly conserved - in parasitic nematodes, is WSY, in Caenorhabditis it is WTY, and in this model (marbled electric ray), it is WTY
- C. elegans is a bit weird in this sense in that it is WTY, but the Caenorhabditis species all look a bit more different relative to the other parasitic species.
- Loss of C.elegans acr-8 doesnt result in levamisole sensitivity, however, Haem acr-8 does substitute for C. elegans lev-8, and lev-9 has the WSY motif.

- The shared Trp(W) residue is conserved immediately downstream, ie pos 173
     - the Trp is a really important residue
     - Notably, the levamisole sensitivity of UNC-29 and UNC-38 nAChR isoforms seems to depend on a glutamate side chain at the position four downstream from TrpB (Rayes et al., 2004), perhaps indicating that, like other agonists, levamisole recognition is sensitive to Loop B-Loop C interactions. (https://www.frontiersin.org/articles/10.3389/fphys.2014.00160/full)

- these residues lie at the interface between the alpha subunit and the gamma subunit



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

******
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
