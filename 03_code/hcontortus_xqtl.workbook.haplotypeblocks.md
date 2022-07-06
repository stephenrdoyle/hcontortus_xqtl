# Estimation of haplotype blocks in the XQTL experiment
- Aim: estimate average haplotype block length that would be maintained from the parental population to the sampled F3
- Rationale: want to know what the minimum QTL size in genome, below which is likely to be noise, above which provides greater confidence of a selective sweep

Overview (see cross design schematic of XQTL paper)
- 100 P0 susceptible females X 100 P0 resistant males
- 5000 F1 progeny used for passage (F2 progeny will be 1st generation of recombinants)
- X F2 adults producing 2x200 F3 L3 samples which were sequenced (F3 progeny will be 2nd generation of recombinants)

```R
# 5000 F1 parasites were using in the passage
f1_autosomes <- 2 * 5000 * 5

# from genetic map, expect 0.69 crossovers per chromosome on average
crossovers <- 0.69
f2_recombinants <- crossovers * f1_autosomes
f2_recombinants
#> 34,500 F2 recombination events

# 5000 F3 parasites sampled from treated F2 adults
f2_autosomes <- 2 * 400 * 5
f3_recombinants <- crossovers * f2_autosomes
f3_recombinants
#> 2760 F3 recombination events

# total recombination events

total_recombinants <- f2_recombinants + f3_recombinants
total_recombinants


# PER POPULATION of 400 sampled
autosome_size <- 237412613

pop_haplotype_block <- (autosome_size * 2)/total_recombinants
pop_haplotype_block
#> 12743.56 bp

# PER INDIVIDUAL of the 400 sampled
co_per_indiv <- total_recombinants / 400
co_per_indiv
#> 93.15 recombination events per individual

indv_haplotype_block <- (autosome_size * 2) / co_per_indiv
indv_haplotype_block
#> 5.097 Mb haplotype blocks from parent in each individual
```




***
# poisson distribution calculation of groups of windows above threshold
```R


vline.data <- data %>%
              group_by(id) %>%
              summarize(mean_fst_3sd = mean(V13)+3*sd(V13))


sample <- bz

#
# frequency of a point about the threshold if randomly distributed in genome, based on number of total Fst windows
f_windows <- nrow(filter(sample, V13 > (mean(V13) + 3*sd(V13)))) / nrow(sample)

# frequency of a point about the threshold if randomly distributed in genome, based on expected haplotype block size
f_haplotypes <- nrow(filter(sample, V13 > (mean(V13) + 3*sd(V13)))) / pop_haplotype_block

# calculate pvalue from poisson
# ppois(1, lambda=bz_freq, lower=F)

# loop over values from 1 to 10
data  <- NULL;
for (i in 1:10)
 {
      p_window <- NULL
      p_haploblock <- NULL
  p_window[1] <- ppois(i, lambda=f_windows, lower=F)
  p_window[2] <- i
  p_haploblock[1] <- ppois(i, lambda=f_haplotypes, lower=F)
  p_haploblock[2] <- i
  data <- rbind(data, p_window, p_haploblock)
 }

myDF <- data.frame(name = row.names(data), data)
colnames(myDF) <- c("test","p_value","count")


ggplot(myDF,aes(count,-log10(p_value),col=test)) +
     geom_point() +
     #facet_grid(test~.) +
     theme_bw() +
     geom_hline(yintercept=-log10((0.01/nrow(sample))),col="blue") +
     geom_hline(yintercept=-log10((0.01/pop_haplotype_block)),col="red")+
     labs(title="Possion calculation of significance based on number of windows found in a group above a threshold", x="Number of windows in a group", col="Region tested")

# save it
ggsave("xqtl_poisson_windows_sig.png")
```
![](../04_analysis/xqtl_poisson_windows_sig.png)
- note: threshold is bonferroni threshold at 0.01 based on number of windows (red) or haploblocks (blue)
- suggests three windows above threshold sufficiently significant
- also suggests could colour runs of windows by significance value.
