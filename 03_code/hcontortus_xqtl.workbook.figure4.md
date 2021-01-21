# Figure 4 - levamisole analyses

1. [Figure 4 A](#figure4a)
2. [Figure 4 B](#figure4b)
3. [Figure 4 C](#figure4c)
4. [Figure 4 D](#figure4d)
5. [Figure 4 E](#figure4e)


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
lev <- read.table("XQTL_LEV.merged.fst", header = F)
lev <- lev[lev$V1! = "hcontortus_chr_mtDNA_arrow_pilon", ]
lev <- dplyr::select(lev, V1, V2, V13)
colnames(lev) <- c("CHR", "POS", "FST")
lev$LABEL <- "Levamisole"
lev$ROW_ID <- 1:nrow(lev)

lev_chr4 <- lev[lev$CHR == "hcontortus_chr4_Celeg_TT_arrow_pilon", ]
lev_chr5 <- lev[lev$CHR == "hcontortus_chr5_Celeg_TT_arrow_pilon", ]

data <- dplyr::bind_rows(lev_chr4, lev_chr5)
data <- data %>%
 mutate(CHR = str_replace_all(CHR, c("hcontortus_chr4_Celeg_TT_arrow_pilon" = "Chromosome 4", "hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))

peaks <- read.table("XQTL_LEV.peak_windows", header = T)
peaks <- peaks %>%
 mutate(CHR = str_replace_all(CHR, c("hcontortus_chr4_Celeg_TT_arrow_pilon" = "Chromosome 4", "hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))

genes <- read.table("candidategenes.data", header = T)
genes <- genes %>%
 mutate(CHR = str_replace_all(CHR, c("hcontortus_chr4_Celeg_TT_arrow_pilon" = "Chromosome 4", "hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))


# set colour for chromosomes
chr_colours<-c("cornflowerblue")

# genome wide signficance per sample
data_gws <- data %>%
  group_by(CHR) %>%
  summarise(GWS = mean(FST) + 3*sd(FST))

# make the plot
plot_a <- ggplot(data)+
   geom_hline(data = data_gws, aes(yintercept = GWS), linetype = "dashed", col = "black")+
   geom_vline(data = genes, aes(xintercept = POS), col = "darkgrey", size = 1)+
   geom_point(aes(POS, FST, group = CHR), size = 0.5, col = "cornflowerblue")+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[1]) & (data$POS <= peaks$PEAK_END_COORD[1]) & (FST > data_gws$GWS[1])), aes(POS, FST), col = "red", size = 1)+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[2]) & (data$POS <= peaks$PEAK_END_COORD[2]) & (FST > data_gws$GWS[2])), aes(POS, FST), col = "red", size = 1)+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[3]) & (data$POS <= peaks$PEAK_END_COORD[3]) & (FST > data_gws$GWS[2])), aes(POS, FST), col = "red", size = 1)+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[4]) & (data$POS <= peaks$PEAK_END_COORD[4]) & (FST > data_gws$GWS[2])), aes(POS, FST), col = "red", size = 1)+
   ylim(0, 0.1)+xlim(0, 50e6)+
   labs(title = "A", x = "Chromosomal position (bp)", y = "Genetic differentiation between \npre- and post-treatment (Fst)")+
   theme_bw()+theme(legend.position = "none", text = element_text(size = 10))+
   facet_grid(CHR~.)

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


plot_a + plot_b + (plot_c | (plot_d / plot_e)) + plot_layout(ncol = 1, height = c(3, 1, 3))
plot_a / ((plot_b / plot_c) | (plot_d / plot_e))  + plot_layout(ncol = 1, height = c(3, 3))
plot_a / (plot_b | plot_d) / (plot_c | plot_e)  + plot_layout(ncol = 1, height = c(3, 1.5, 3))


ggsave("XQTL_Figure_4.pdf", useDingbats = FALSE, width = 250, height = 300, units = "mm")
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


******
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
