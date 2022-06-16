
library(tidyverse)
library(zoo)
library(ggridges)

# control
control <- read.table("XQTL_CONTROL.merged.fst", header=F)
control <- control %>% filter(V1 != "hcontortus_chr_mtDNA_arrow_pilon")
control <- control %>% filter(V1 != "hcontortus_chrX_Celeg_TT_arrow_pilon")
control <- dplyr::select(control,  V1,  V2,  V11 , V21, V29)
control <- control %>% mutate(mean_FST = rowMeans(select(.,V11,V21,V29)))
control <- control %>% mutate(roll_mean = rollmean(mean_FST, k = 200, fill = NA))
control$LABEL <- "1. Control"
control$LABEL2 <- "0.175"
control$ROW_ID <- 1:nrow(control)
colnames(control) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN", "FST_ROLLMEAN", "LABEL", "LABEL2", "ROW_ID")


# benzimidazole
bz <- read.table("XQTL_BZ.merged.fst", header = F)
bz <- bz %>% filter(bz$V1 != "hcontortus_chr_mtDNA_arrow_pilon")
bz <- bz %>% filter(bz$V1 != "hcontortus_chrX_Celeg_TT_arrow_pilon")
bz <- dplyr::select(bz,  V1,  V2,  V13 , V39, V49)
bz <- bz %>% mutate(mean_FST = rowMeans(select(.,V13,V39,V49)))
bz <- bz %>% mutate(roll_mean = rollmean(mean_FST, k = 200, fill = NA))
bz$LABEL <- "2. Fenbendazole"
bz$LABEL2 <- "0.15"
bz$ROW_ID <- 1:nrow(bz)
colnames(bz) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN", "FST_ROLLMEAN", "LABEL","LABEL2",  "ROW_ID")

data <- dplyr::bind_rows(control, bz)


# levamisole
lev <- read.table("XQTL_LEV.merged.fst", header = F)
lev <- lev %>% filter(lev$V1 != "hcontortus_chr_mtDNA_arrow_pilon")
lev <- lev %>% filter(lev$V1 != "hcontortus_chrX_Celeg_TT_arrow_pilon")
lev <- dplyr::select(lev,  V1,  V2,  V13 , V39, V49)
lev <- lev %>% mutate(mean_FST = rowMeans(select(.,V13)))
lev <- lev %>% mutate(roll_mean = rollmean(mean_FST, k = 200, fill = NA))
lev$LABEL <- "3. Levamisole"
lev$LABEL2 <- "0.125"
lev$ROW_ID <- 1:nrow(lev)
colnames(lev) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN", "FST_ROLLMEAN", "LABEL", "LABEL2", "ROW_ID")


data <- dplyr::bind_rows(data, lev)

#ivermectin
ivm <- read.table("XQTL_IVM.merged.fst", header = F)
ivm <- ivm %>% filter(ivm$V1 != "hcontortus_chr_mtDNA_arrow_pilon")
ivm <- ivm %>% filter(ivm$V1 != "hcontortus_chrX_Celeg_TT_arrow_pilon")
ivm <- dplyr::select(ivm,  V1,  V2,  V13 , V39, V49)
ivm <- ivm %>% mutate(mean_FST = rowMeans(select(.,V13,V39,V49)))
ivm <- ivm %>% mutate(roll_mean = rollmean(mean_FST, k = 200, fill = NA))
ivm$LABEL <- "4. Ivermectin"
ivm$LABEL2 <- "0.1"
ivm$ROW_ID <- 1:nrow(ivm)
colnames(ivm) <- c("CHR",  "POS",  "FST_R1", "FST_R2", "FST_R3", "FST_MEAN", "FST_ROLLMEAN", "LABEL", "LABEL2", "ROW_ID")

data <- dplyr::bind_rows(data, ivm)


ggplot(data) + geom_point(aes(ROW_ID,as.numeric(LABEL2)+FST_ROLLMEAN,col=FST_ROLLMEAN, group=LABEL)) +
     scale_colour_gradient(low = "royalblue", high = "red", na.value = NA, guide="none") +
     labs(x="", y="") +
     theme_classic() +
     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

ggsave()








# load required libraries
library(tidyverse)
library(rstatix)
library(ggrepel)
library(reshape2)
library(patchwork)

us_btub2 <- read.table("/nfs/users/nfs_s/sd21/lustre118_link/haemonchus_contortus/XQTL/04_VARIANTS/US_FIELD/VCF/us_farms_btub2.ADfreq")
colnames(us_btub2) <- c("CHR", "POS", "Farm 1", "Farm 2", "Farm 3", "Farm 4", "Farm 5", "Farm 6", "Farm 7", "Farm 8", "Farm 9", "Farm 10")
us_btub2 <- melt(us_btub2,  id = c("CHR",  "POS"),  variable.name = "SAMPLE_ID")

bz_conc <- c(0.05, 19.93, 15, 7.71, 15, 61.57, 29.6, 18.47, 9.67, 29.6)

us_btub2$BZ_CONCENTRATION <- bz_conc
colnames(us_btub2) <- c("CHR", "POS", "SAMPLE_ID", "ALLELE_FREQ", "BZ_CONCENTRATION")

us_btub2 <-
     us_btub2 %>%
     mutate(POS = str_replace(POS,  c("13435823"),  c("Glu198Val")))

# calculate correlation coefficient between alllele frequency and concentration
af_bz_cor <- cor.test(us_btub2$ALLELE_FREQ,  us_btub2$BZ_CONCENTRATION,  method = "pearson",  use = "complete.obs")

# make the plot
plot_btub2_EC50 <-
     ggplot(us_btub2)+
     geom_smooth(aes(BZ_CONCENTRATION, ALLELE_FREQ), method = 'lm', col = 'grey')+
     geom_jitter(aes(BZ_CONCENTRATION, ALLELE_FREQ, col = BZ_CONCENTRATION), size = 4)+
     ylim(-0.05, 1) +
     scale_colour_gradient(low = "royalblue", high = "red", na.value = NA, guide="none") +
     theme_bw() + theme(legend.position = "none", text = element_text(size = 10))

plot_btub2_EC50
