# Figure 5 - ivermectin analyses

1. [Figure 5 A](#figure5a)
2. [Figure 5 B](#figure5b)
3. [Figure 5 C](#figure5c)
4. [Figure 5 D](#figure5d)
5. [Figure 5 E](#figure5e)


## Figure 5 A <a name="figure5a"></a>

Aim is show genetic differentiation between pre and post treatment ivermectin samples on chromosome 5

```shell
# working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/IVM
```

ln -s ../../04_VARIANTS/XQTL_IVM/XQTL_IVM.merged.fst



### R to plot
```R
# load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

# load and reformat the data
ivm <- read.table("XQTL_IVM.merged.fst", header = F)
ivm <- ivm[ivm$V1! = "hcontortus_chr_mtDNA_arrow_pilon", ]
ivm <- dplyr::select(ivm, V1, V2, V13)
colnames(ivm) <- c("CHR", "POS", "FST")
ivm$LABEL <- "Ivermectin"
ivm$ROW_ID <- 1:nrow(ivm)

ivm_chr5 <- ivm[ivm$CHR == "hcontortus_chr5_Celeg_TT_arrow_pilon", ]

data <- ivm_chr5
data <- data %>%
 mutate(CHR = str_replace_all(CHR, c("hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))

peaks <- read.table("XQTL_IVM.peak_windows", header = T)
peaks <- peaks %>%
 mutate(CHR = str_replace_all(CHR, c("hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))

#genes <- read.table("candidategenes.data", header = T)
#genes <- genes %>%
# mutate(CHR = str_replace_all(CHR, c("hcontortus_chr4_Celeg_TT_arrow_pilon" = "Chromosome 4", "hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))


# set colour for chromosomes
chr_colours<-c("cornflowerblue")

# genome wide signficance per sample
data_gws <- data %>%
  group_by(CHR) %>%
  summarise(GWS = mean(FST) + 3*sd(FST))

# make the plot
plot_a <- ggplot(data)+
   geom_hline(data = data_gws, aes(yintercept = GWS), linetype = "dashed", col = "black")+
   #geom_vline(data = genes, aes(xintercept = POS), col = "darkgrey", size = 1)+
   geom_point(aes(POS, FST, group = CHR), size = 0.5, col = "cornflowerblue")+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[1]) & (data$POS <= peaks$PEAK_END_COORD[1]) & (FST > data_gws$GWS[1])), aes(POS, FST), col = "red", size = 1)+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[3]) & (data$POS <= peaks$PEAK_END_COORD[3]) & (FST > data_gws$GWS[1])), aes(POS, FST), col = "red", size = 1)+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[4]) & (data$POS <= peaks$PEAK_END_COORD[4]) & (FST > data_gws$GWS[1])), aes(POS, FST), col = "red", size = 1)+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[5]) & (data$POS <= peaks$PEAK_END_COORD[5]) & (FST > data_gws$GWS[1])), aes(POS, FST), col = "red", size = 1)+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[6]) & (data$POS <= peaks$PEAK_END_COORD[6]) & (FST > data_gws$GWS[1])), aes(POS, FST), col = "red", size = 1)+
   geom_point(data = subset(data, (data$POS >= peaks$PEAK_START_COORD[7]) & (data$POS <= peaks$PEAK_END_COORD[7]) & (FST > data_gws$GWS[1])), aes(POS, FST), col = "red", size = 1)+
   ylim(0, 0.1)+xlim(0, 50e6)+
   labs(title = "A", x = "Chromosomal position (bp)", y = "Genetic differentiation between \npre- and post-treatment (Fst)")+
   theme_bw()+theme(legend.position = "none", text = element_text(size = 10))+
   facet_grid(CHR~.)

#plot_a
```

## Figure 4 B <a name = "figure4b"></a>

Aim is to show zoomed in region around chromosome 5 main peak, highlighting farm data, genes present, whether genes are differentially expressed

```shell
# working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/IVM

#gff
ln -fs /nfs/users/nfs_s/sd21/lustre118_link/hc/GENOME/TRANSCRIPTOME/TRANSCRIPTOME_CURATION/20200407/UPDATED_annotation.gff3 ANNOTATION.gff
```

grep "gene" ANNOTATION.gff | sed -e 's/owner=irisadmin@local.host;//g' -e 's/date_last_modified=[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9];//g' | awk -F '[\t=;]' '{print $12,$1,$4,$5,$7}' OFS="\t" > genes.positions
```


```
genes <- read.table("genes.positions",header=F)
colnames(genes) <- c("name_PUGAxPISE","chr","start","end","dir")

rnaseq <- read.table("summary_DESEQ_WBPS15gff_alpha0.01_males_all3IVMcomps.tabular", header=T)

rnaseq_data <- dplyr::full_join(genes,rnaseq,by="name_PUGAxPISE")

rnaseq_data_chr5 <- rnaseq_data[rnaseq_data$chr == "hcontortus_chr5_Celeg_TT_arrow_pilon", ]
rnaseq_data_chr5 <- rnaseq_data_chr5 %>%
 mutate(chr = str_replace_all(chr, c("hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))
```


```
# parental ISE vs parent UGA
rnaseq_plot1 <- ggplot(rnaseq_data_chr5) +
       geom_point(aes(start,log2FoldChange),size=0.5,col="lightgrey")+
       geom_point(aes(start,log2FoldChange, col=-log10(padj), size=-log10(padj)))+
       geom_text_repel(data=subset(rnaseq_data_chr5,log2FoldChange > 2 | log2FoldChange < -2),aes(start,log2FoldChange, label=name_PUGAxPISE)) +
       scale_colour_viridis(direction=-1, limits = c(0, 30))+
       scale_size_continuous(limits = c(0, 30)) +
       theme_bw()+
       xlim(3.65e7,3.80e7)+
       labs(Colour="-log10(adjusted p-value)", Size="-log10(adjusted p-value)", x="Genomics position", y= "log2(fold change): Pre vs Post treatment")


# parental ISE vs post treatmetn F3
rnaseq_plot2 <- ggplot(rnaseq_data_chr5) +
      geom_point(aes(start,log2FoldChange.2),size=0.5,col="lightgrey")+
      geom_point(aes(start,log2FoldChange.2,col=-log10(padj.2),size=-log10(padj.2)))+
      geom_text_repel(data=subset(rnaseq_data_chr5,log2FoldChange.2 > 2 | log2FoldChange.2 < -2),aes(start,log2FoldChange.2, label=name_PUGAxPISE)) +
      scale_colour_viridis(direction=-1,limits = c(0, 30))+
      scale_size_continuous(limits = c(0, 30)) +
      theme_bw()+
      xlim(3.65e7,3.80e7)+
      labs(Colour="-log10(adjusted p-value)", Size="-log10(adjusted p-value)", x="Genomics position", y= "log2(fold change): Pre vs Post treatment")


# pre vs post treatment
rnaseq_plot3 <-      ggplot(rnaseq_data_chr5) +
            geom_point(aes(start,log2FoldChange.1),size=0.5,col="lightgrey")+
            geom_point(aes(start,log2FoldChange.1,col=-log10(padj.1),size=-log10(padj.1)))+
            geom_text_repel(data=subset(rnaseq_data_chr5,log2FoldChange.1 > 2 | log2FoldChange.1 < -2),aes(start,log2FoldChange.1, label=name_PUGAxPISE)) +
            scale_colour_viridis(direction=-1,limits = c(0, 30))+
            scale_size_continuous(limits = c(0, 30)) +
            theme_bw()+
            xlim(3.65e7,3.80e7)+
            labs(Colour="-log10(adjusted p-value)", Size="-log10(adjusted p-value)", x="Genomics position", y= "log2(fold change): Pre vs Post treatment")

rnaseq_plot1 + rnaseq_plot2 + rnaseq_plot3 + plot_layout(ncol=1)


### R to plot
```R
# load required R libraries
library(data.table)
# note - "data.table" has a function called "fread" which is great for quickly loading really large datasets
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(viridis)

# load data
# --- genome annotation
gff <- fread(cmd = "grep ^hcontortus ANNOTATION.gff")
colnames(gff) <- c("chr", "source", "feature", "start", "end", "point1", "strand", "frame", "info")




us_farm1 <-read.table("sample_1.pileup.stats",header=T,sep="\t")
us_farm1$sample <- "us_farm1"
us_farm1$ivm_EC50 <- "1.51"
us_farm2 <-read.table("sample_2.pileup.stats",header=T,sep="\t")
us_farm2$sample <- "us_farm2"
us_farm2$ivm_EC50 <- "12.04"
us_farm3 <-read.table("sample_3.pileup.stats",header=T,sep="\t")
us_farm3$sample <- "us_farm3"
us_farm3$ivm_EC50 <- "9.15"
us_farm4 <-read.table("sample_4.pileup.stats",header=T,sep="\t")
us_farm4$sample <- "us_farm4"
us_farm4$ivm_EC50 <- "11.27"
us_farm5 <-read.table("sample_5.pileup.stats",header=T,sep="\t")
us_farm5$sample <- "us_farm5"
us_farm5$ivm_EC50 <- "297.9"
us_farm6 <-read.table("sample_6.pileup.stats",header=T,sep="\t")
us_farm6$sample <- "us_farm6"
us_farm6$ivm_EC50 <- "619"
us_farm7 <-read.table("sample_7.pileup.stats",header=T,sep="\t")
us_farm7$sample <- "us_farm7"
us_farm7$ivm_EC50 <- "312.1"
us_farm8 <-read.table("sample_8.pileup.stats",header=T,sep="\t")
us_farm8$sample <- "us_farm8"
us_farm8$ivm_EC50 <- "NA"
us_farm9 <-read.table("sample_9.pileup.stats",header=T,sep="\t")
us_farm9$sample <- "us_farm9"
us_farm9$ivm_EC50 <- "0.977"
us_farm10 <-read.table("sample_10.pileup.stats",header=T,sep="\t")
us_farm10$sample <- "us_farm10"
us_farm10$ivm_EC50 <- "259.4"


us_farm_data <-  dplyr::bind_rows(us_farm1, us_farm2, us_farm3, us_farm4, us_farm5, us_farm6, us_farm7, us_farm8, us_farm9, us_farm10)

us_farm_data_chr5 <- us_farm_data[us_farm_data$chr == "hcontortus_chr5_Celeg_TT_arrow_pilon", ]
us_farm_data_chr5 <- us_farm_data_chr5 %>%
 mutate(chr = str_replace_all(chr, c("hcontortus_chr5_Celeg_TT_arrow_pilon" = "Chromosome 5")))



plot_fst <- ggplot(data)+
   geom_hline(data = data_gws, aes(yintercept = GWS), linetype = "dashed", col = "black")+
   #geom_vline(data = genes, aes(xintercept = POS), col = "darkgrey", size = 1)+
   geom_point(aes(POS, FST, group = CHR), size = 1, col = "cornflowerblue", alpha=0.5)+
   xlim(3.65e7,3.85e7)+
   labs(title = "B", x="Genomic position", y="Genetic differentiation between \npre- and post-treatment (Fst)")+
   theme_bw()+theme(text = element_text(size = 10))

plot_pi <- ggplot(us_farm_data_chr5)+
     geom_point(aes((window*10000)-5000,(Pi),col=as.numeric(ivm_EC50)), size=1, alpha=0.5)+
     xlim(3.65e7,3.85e7)+
     scale_colour_viridis(direction=1,limits = c(0, 800))+
     labs(title = "C", x="Genomic position", y="US Farm \nNucleotide diversity (Pi)", colour="Ivermectin EC50 per farm")+
     theme_bw()+theme(text = element_text(size = 10))

#plot_tajD <- ggplot(us_farm_data_chr5)+
#         geom_point(aes((window*10000)-5000,(Tajima_D),col=as.numeric(ivm_EC50)),size=1.5,alpha=0.5)+
#         theme_bw()+
#         xlim(3.65e7,3.85e7)+
#         scale_colour_viridis(direction=1, guide = FALSE,limits = c(0, 800))+
#         labs(x="Genomic position", y="Tajima's D")+
#         theme()

rnaseq_plot2 <- ggplot(rnaseq_data_chr5) +
               geom_point(aes(start,log2FoldChange.2),size=0.5,col="lightgrey")+
               geom_point(aes(start,log2FoldChange.2,col=-log10(padj.2),size=-log10(padj.2)*0.5))+
               #geom_text_repel(data=subset(rnaseq_data_chr5,log2FoldChange.2 > 2 | log2FoldChange.2 < -2),aes(start,log2FoldChange.2, label=name_PUGAxPISE),size=3) +
               scale_colour_viridis(direction=1,limits = c(0, 30), option="magma")+
               scale_size_continuous(limits = c(0, 30)) +
               xlim(3.65e7,3.85e7)+ylim(-7,7)+
               labs(title = "D", colour="-log2(adjusted p-value)", size="-log2(adjusted p-value)", x="Genomic position", y= "RNAseq log2(fold change): \nPost treatment vs ISE parent")+
               theme_bw()+ theme(text = element_text(size = 10))



# bring it all together

plot_a / ((plot_fst / plot_pi / rnaseq_plot2) |  (plot_spacer() / guide_area())) + plot_layout(ncol=1, guides='collect',heights = c(1, 3)) & theme(legend.position = 'bottom',legend.direction = "horizontal", legend.box="vertical")


ggsave("XQTL_Figure_5.pdf", useDingbats = FALSE, width = 250, height = 300, units = "mm")
ggsave("XQTL_Figure_5.png")


```

## Figure 5 C <a name = "figure4c"></a>



******
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
