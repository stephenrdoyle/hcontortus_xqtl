# Figure 2 - Fst plot for parent strains

1. [Figure 2 A](#figure2a)
2. [Figure 2 B](#figure2b)


---
## Figure 2 A <a name="figure2a"></a>

Aim is to show genetic differentiation between the two parental strains,  MHco3(ISE) and MHco18(UGA)

```shell
# my working dir:
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/GENOMEWIDE
```

```
# download data needed to make the figures
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/sd21/hcontortus_xqtl/figure2/XQTL_BZ.merged.fst
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/sd21/hcontortus_xqtl/figure2/XQTL_CONTROL.merged.fst
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/sd21/hcontortus_xqtl/figure2/XQTL_IVM.merged.fst
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/sd21/hcontortus_xqtl/figure2/XQTL_LEV.merged.fst
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/sd21/hcontortus_xqtl/figure2/XQTL_PARENTS.merged.fst
```

### R to plot
```R
# load required libraries
library(tidyverse)
library(patchwork)


# reformat the data
parents <- read.table("XQTL_PARENTS.merged.fst",  header = F)
parents <- parents[parents$V1 != "hcontortus_chr_mtDNA_arrow_pilon",  ]
parents <- dplyr::select(parents,  V1,  V2,  V7)
colnames(parents) <- c("CHR",  "POS",  "FST")
parents$LABEL <- "MHco3(ISE) vs MHco18(UGA)"
parents$ROW_ID <- 1:nrow(parents)
colnames(parents) <- c("CHR",  "POS",  "FST",  "LABEL",  "ROW_ID")

data <- parents

# genome wide significance per sample
data_gws <- data %>%
    group_by(LABEL) %>%
    summarise(GWS = mean(FST) + 3*sd(FST))

# chromosome colours
chr_colours <- c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

# make the plot
plot_a <- ggplot(data)+
     geom_hline(data = data_gws,  aes(yintercept = GWS),  linetype = "dashed",  col = "black") +
     geom_point(aes(ROW_ID * 5000,  FST,  colour = CHR,  group = LABEL),  size = 0.1) +
     ylim(0, 1) +
     labs(title = "A",  x = "Genomic position (bp)",  y = expression(paste("Genetic differentiation"," (",~italic(F)[ST],")"))) +
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
     theme_bw() + theme(legend.position = "none",  text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     facet_grid(LABEL ~ .)

```

## Figure 2 B <a name="figure2b"></a>

Aim is to show genetic differentiation between pre and post treatment F3 generation for each of the control,  benzimidazole,  levamisole,  and ivermectin XQTL datasets

```R
# load required libraries
library(tidyverse)
library(patchwork)



# control
control <- read.table("XQTL_CONTROL.merged.fst", header=F)
control <- control[control$V1!="hcontortus_chr_mtDNA_arrow_pilon", ]
control <- dplyr::select(control,  V1,  V2,  V11 , V21, V29)
control <- control %>% mutate(mean_FST = rowMeans(select(.,V11,V21,V29)))
control$LABEL <- "1. Control"
control$ROW_ID <- 1:nrow(control)
colnames(control) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN",  "LABEL",  "ROW_ID")

# benzimidazole
bz <- read.table("XQTL_BZ.merged.fst", header = F)
bz <- bz[bz$V1 != "hcontortus_chr_mtDNA_arrow_pilon", ]
bz <- dplyr::select(bz,  V1,  V2,  V13 , V39, V49)
bz <- bz %>% mutate(mean_FST = rowMeans(select(.,V13,V39,V49)))
bz$LABEL <- "2. Benzimidazole"
bz$ROW_ID <- 1:nrow(bz)
colnames(bz) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN",  "LABEL",  "ROW_ID")

data <- dplyr::bind_rows(control, bz)


# levamisole
lev <- read.table("XQTL_LEV.merged.fst", header = F)
lev <- lev[lev$V1 != "hcontortus_chr_mtDNA_arrow_pilon", ]
lev <- dplyr::select(lev,  V1,  V2,  V13 , V39, V49)
lev <- lev %>% mutate(mean_FST = rowMeans(select(.,V13)))
lev$LABEL <- "3. Levamisole"
lev$ROW_ID <- 1:nrow(lev)
colnames(lev) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN",  "LABEL",  "ROW_ID")


data <- dplyr::bind_rows(data, lev)

#ivermectin
ivm <- read.table("XQTL_IVM.merged.fst", header = F)
ivm <- ivm[ivm$V1 != "hcontortus_chr_mtDNA_arrow_pilon", ]
ivm <- dplyr::select(ivm,  V1,  V2,  V13 , V39, V49)
ivm <- ivm %>% mutate(mean_FST = rowMeans(select(.,V13,V39,V49)))
ivm$LABEL <- "4. Ivermectin"
ivm$ROW_ID <- 1:nrow(ivm)
colnames(ivm) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN",  "LABEL",  "ROW_ID")


data <- dplyr::bind_rows(data, ivm)


# genome wide signficance per sample
data_gws <- data %>%
    group_by(LABEL) %>%
    summarise(GWS = mean(FST_MEAN) + 3*sd(FST_MEAN))



# set the colours
chr_colours<-c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

# make the plot
plot_b <- ggplot(data) +
     geom_hline(data = data_gws,  aes(yintercept = GWS),  linetype = "dashed", col = "black") +
     geom_point(aes(ROW_ID * 5000,  FST_MEAN,  colour = CHR,  group = LABEL), size = 0.1) +
     #ylim(0, 0.12) +
     labs(title = "B", x = "Genomic position (bp)",  y = expression(paste("Genetic differentiation between pre- and post-treatment", " (",~italic(F)[ST],")")))+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0, 3e8, 0.5e8), limits = c(0, 300e6)) + ylim(0,0.15) +
     theme_bw() + theme(legend.position = "none", text = element_text(size=10)) +
     facet_grid(LABEL ~ .)


# combine the parental and XQTL plots into a single figure
plot_a + plot_b + plot_layout(ncol = 1,  heights = c(1, 4))

# save it
ggsave("genomewide_fst_plots.pdf", useDingbats = FALSE, width = 170, height = 200, units = "mm")
ggsave("genomewide_fst_plots.png")
```
![](04_analysis/genomewide_fst_plots.png)


# replicate 1

data_gws <- data %>%
    group_by(LABEL) %>%
    summarise(GWS = mean(FST_R1) + 3*sd(FST_R1))

ggplot(data) +
     geom_hline(data = data_gws,  aes(yintercept = GWS),  linetype = "dashed", col = "black") +
     geom_point(aes(ROW_ID * 5000,  FST_R1,  colour = CHR,  group = LABEL), size = 0.1) +
     #ylim(0, 0.12) +
     labs(title = "XQTL replicate set 1", x = "Genomic position (bp)",  y = expression(paste("Genetic differentiation between pre- and post-treatment", " (",~italic(F)[ST],")")))+
     scale_color_manual(values = chr_colours) +
     scale_x_continuous(breaks = seq(0, 3e8, 0.5e8), limits = c(0, 300e6)) +
     theme_bw() + theme(legend.position = "none", text = element_text(size=10)) +
     facet_grid(LABEL ~ ., scales="free_y")

ggsave("genomewide_fst_plots_R1_supplement.pdf", useDingbats = FALSE, width = 170, height = 200, units = "mm")
ggsave("genomewide_fst_plots_R1_supplement.png")

# replicate 2
data_gws <- data %>%
    group_by(LABEL) %>%
    summarise(GWS = mean(FST_R2) + 3*sd(FST_R2))

ggplot(data) +
          geom_hline(data = data_gws,  aes(yintercept = GWS),  linetype = "dashed", col = "black") +
          geom_point(aes(ROW_ID * 5000,  FST_R2,  colour = CHR,  group = LABEL), size = 0.1) +
          #ylim(0, 0.12) +
          labs(title = "XQTL replicate set 2", x = "Genomic position (bp)",  y = expression(paste("Genetic differentiation between pre- and post-treatment", " (",~italic(F)[ST],")")))+
          scale_color_manual(values = chr_colours) +
          scale_x_continuous(breaks = seq(0, 3e8, 0.5e8), limits = c(0, 300e6)) +
          theme_bw() + theme(legend.position = "none", text = element_text(size=10)) +
          facet_grid(LABEL ~ ., scales="free_y")

ggsave("genomewide_fst_plots_R2_supplement.pdf", useDingbats = FALSE, width = 170, height = 200, units = "mm")
ggsave("genomewide_fst_plots_R2_supplement.png")

# replicate 3
data_gws <- data %>%
    group_by(LABEL) %>%
    summarise(GWS = mean(FST_R3) + 3*sd(FST_R3))

          ggplot(data) +
               geom_hline(data = data_gws,  aes(yintercept = GWS),  linetype = "dashed", col = "black") +
               geom_point(aes(ROW_ID * 5000,  FST_R3,  colour = CHR,  group = LABEL), size = 0.1) +
               #ylim(0, 0.12) +
               labs(title = "XQTL replicate set 3", x = "Genomic position (bp)",  y = expression(paste("Genetic differentiation between pre- and post-treatment", " (",~italic(F)[ST],")")))+
               scale_color_manual(values = chr_colours) +
               scale_x_continuous(breaks = seq(0, 3e8, 0.5e8), limits = c(0, 300e6)) +
               theme_bw() + theme(legend.position = "none", text = element_text(size=10)) +
               facet_grid(LABEL ~ ., scales="free_y")

ggsave("genomewide_fst_plots_R3_supplement.pdf", useDingbats = FALSE, width = 170, height = 200, units = "mm")
ggsave("genomewide_fst_plots_R3_supplement.png")






******
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
<!--

### other stuff that wasn't used but worth keeping


```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_PARENTS.merged.fst", header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon", ]

data$LABEL <- "Parents: MHco3(ISE) vs MHco18(UGA)"
# genome wide levels of significance


chr_colours<-c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

fst_column <- data$V7

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_a <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig,  linetype="dashed", col="black")+
     geom_point(aes(1:nrow(data)*5000,  fst_column,  colour = V1), size=0.1)+
     ylim(0, 1)+
     labs(title="A", x="Chromosome position (5 kbp window)",  y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0, 3e8, 0.5e8))+
     theme_bw()+theme(legend.position="none", text = element_text(size=10))+
     facet_grid(LABEL~.)


#ggsave("XQTL_parents_fst.pdf", useDingbats=FALSE, width=170, height=50, units="mm")

```

```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_CONTROL.merged.fst", header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon", ]
data$LABEL <- "Control"
# genome wide levels of significance


chr_colours<-c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

# Replicates
# 1:4 = V11
# 2:5 = V21
# 3:6 = V29

fst_column <- data$V11

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_b <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig,  linetype="dashed", col="black")+
     geom_point(aes(1:nrow(data)*5000,  fst_column,  colour = V1), size=0.1)+
     ylim(0, 0.1)+
     labs(title="B", x="Chromosome position (5 kbp window)",  y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0, 3e8, 0.5e8))+
     theme_bw()+theme(legend.position="none", text = element_text(size=10))+
     facet_grid(LABEL~.)

```
```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_BZ.merged.fst", header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon", ]
data$LABEL <- "Benzimidazole"
# genome wide levels of significance


chr_colours<-c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

# Replicates
# 1:5 = V13
# 2:6 = V27
# 3:7 = V39
# 4:8 = V49


fst_column <- data$V13

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_c <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig,  linetype="dashed", col="black")+
     geom_point(aes(1:nrow(data)*5000,  fst_column,  colour = V1), size=0.1)+
     ylim(0, 0.1)+
     labs(title="C", x="Chromosome position (5 kbp window)",  y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0, 3e8, 0.5e8))+
     theme_bw()+theme(legend.position="none", text = element_text(size=10))+
     facet_grid(LABEL~.)
```

```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_LEV.merged.fst", header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon", ]
data$LABEL <- "Levamisole"
# genome wide levels of significance


chr_colours<-c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

# Replicates
# 1:5 = V13
# 2:6 = V27
# 3:7 = V39
# 4:8 = V49


fst_column <- data$V13

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_d <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig,  linetype="dashed", col="black")+
     geom_point(aes(1:nrow(data)*5000,  fst_column,  colour = V1), size=0.1)+
     ylim(0, 0.1)+
     labs(title="D", x="Chromosome position (5 kbp window)",  y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0, 3e8, 0.5e8))+
     theme_bw()+theme(legend.position="none", text = element_text(size=10))+
     facet_grid(LABEL~.)
```

## Ivermectin
```R
R
library(ggplot2)

# get data
data <- read.table("XQTL_IVM.merged.fst", header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon", ]
data$LABEL <- "Ivermectin"
# genome wide levels of significance


chr_colours<-c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")

# Replicates
# 1:5 = V13
# 2:6 = V27
# 3:7 = V39
# 4:8 = V49


fst_column <- data$V13

genomewide_sig <- mean(fst_column)+(3*sd(fst_column))

plot_e <- ggplot(data)+
     geom_hline(yintercept=genomewide_sig,  linetype="dashed", col="black")+
     geom_point(aes(1:nrow(data)*5000,  fst_column,  colour = V1), size=0.1)+
     ylim(0, 0.1)+
     labs(title="E", x="Chromosome position (5 kbp window)",  y="Genetic differentiation (Fst)")+
     scale_color_manual(values=chr_colours)+
     scale_x_continuous(breaks=seq(0, 3e8, 0.5e8))+
     theme_bw()+theme(legend.position="none", text = element_text(size=10))+
     facet_grid(LABEL~.)
```

# bring it together
```R
library(patchwork)
plot_a + plot_b + plot_c + plot_d + plot_e + plot_layout(ncol=1)


``` -->
