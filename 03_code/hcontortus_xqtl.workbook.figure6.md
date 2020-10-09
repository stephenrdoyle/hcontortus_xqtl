
#-----
library(ggplot2)
library(patchwork)

data <- read.table("/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/ADVANCED_INTERCROSS/XQTL_ADVANCED_INTERCROSS.merged.fst",header=F)
data <- data[data$V1!="hcontortus_chr_mtDNA_arrow_pilon",]

# chromosome colours
chr_colours <- c("blue", "cornflowerblue", "blue", "cornflowerblue", "blue", "cornflowerblue")




#						Rep1	Rep2	Rep3
#control - pre v 0.5X	V17		V51		V83
#control - pre v 2X		V11		V135	V161

#IVM pre v 0.5x			V251	V267	V281
#IVM pre v 2x			V287	V297	V305


#rep1
# control_0.5x.1 <- ggplot(data)+
#      geom_point(aes(1:nrow(data)*5000, V17,  colour = V1),size=0.05)+
#      ylim(0,0.1) +
#      labs(x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
#      scale_color_manual(values = chr_colours) +
#      scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
#      theme_bw()+
#      theme(legend.position = "none", text = element_text(size = 8))
#
# control_2X.1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V11,colour = V1),size=0.05)+
#   ylim(0,0.1)+
#   labs(x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
#   scale_color_manual(values = chr_colours) +
#   scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
#   theme_bw()+
#   theme(legend.position = "none", text = element_text(size = 8))
#
# ivm_0.5X.1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V251,colour = V1),size=0.05)+
#   ylim(0,0.1)+
#   labs(x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
#   scale_color_manual(values = chr_colours) +
#   scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
#   theme_bw()+
#   theme(legend.position = "none", text = element_text(size = 8))

ivm_2X.1 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V287,colour = V1),size=0.05)+
  ylim(0,0.1)+
  labs(title = "F", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
  scale_color_manual(values = chr_colours) +
  scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 8))

#control_0.5x.1 + ivm_0.5X.1 + control_2X.1 + ivm_2X.1 + plot_layout(ncol=2)


#rep2
control_0.5x.2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V51,colour = V1),size=0.05)+
  ylim(0,0.1)+
  labs(title = "A", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
  scale_color_manual(values = chr_colours) +
  scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 8))

control_2X.2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V135,colour = V1),size=0.05)+
  ylim(0,0.1)+
  labs(title = "B", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
  scale_color_manual(values = chr_colours) +
  scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 8))

ivm_0.5X.2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V267,colour = V1),size=0.05)+
  ylim(0,0.1)+
  labs(title = "C", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
  scale_color_manual(values = chr_colours) +
  scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 8))

ivm_2X.2 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V297,colour = V1),size=0.05)+
  ylim(0,0.1)+
  labs(title = "D", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
  scale_color_manual(values = chr_colours) +
  scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 8))

 ivm_2X.22 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V297,colour = V1),size=0.05)+
    ylim(0,0.1)+
    labs(title = "E", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
    scale_color_manual(values = chr_colours) +
    scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
    theme_bw()+
    theme(legend.position = "none", text = element_text(size = 8))



#control_0.5x.2 + ivm_0.5X.2 + control_2X.2 + ivm_2X.2 + plot_layout(ncol=2)


#rep3
# control_0.5x.3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V83,colour = V1),size=0.1)+
#   ylim(0,0.1)+
#   labs(x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
#   scale_color_manual(values = chr_colours) +
#   scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
#   theme_bw()+
#   theme(legend.position = "none", text = element_text(size = 10))
#
# control_2X.3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V161,colour = V1),size=0.1)+
#   ylim(0,0.1)+
#   labs(x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
#   scale_color_manual(values = chr_colours) +
#   scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
#   theme_bw()+
#   theme(legend.position = "none", text = element_text(size = 10))
#
# ivm_0.5X.3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V281,colour = V1),size=0.1)+
#   ylim(0,0.1)+
#   labs(x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
#   scale_color_manual(values = chr_colours) +
#   scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
#   theme_bw()+
#  theme(legend.position = "none", text = element_text(size = 10))

ivm_2X.3 <- ggplot(data)+geom_point(aes(1:nrow(data)*5000,V305,colour = V1),size=0.05)+
  ylim(0,0.1)+
  labs(title = "G", x = "Chromosomal position (bp)",  y = "Genetic differentiation between \npre- and post-treatment (Fst)") +
  scale_color_manual(values = chr_colours) +
  scale_x_continuous(breaks = seq(0,  3e8,  0.5e8), limits = c(0,  300e6)) +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 8))

#control_0.5x.3 + ivm_0.5X.3 + control_2X.3 + ivm_2X.3 + plot_layout(ncol=2)




# summary of XQTL and AI endpoints
(control_0.5x.2 | ivm_0.5X.2 ) / (
control_2X.2 | ivm_2X.2) / ivm_2X.22 + ivm_2X.1 + ivm_2X.3 + plot_layout(ncol=1)

ggsave("XQTL_Figure6.pdf", useDingbats = FALSE, width = 170, height = 250, units = "mm")
ggsave("XQTL_Figure6.png")
