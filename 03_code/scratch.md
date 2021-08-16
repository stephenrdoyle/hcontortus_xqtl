# scratch space for random code than needs a home




## playing around with correlations between EC50 and US farm Fst data

/nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/05_ANALYSIS/IVM

ln -s ../../04_VARIANTS/US_FIELD/XQTL_US_FIELD.merged.fst

# manually input the EC50 data from Table SX

```R
library(tidyverse)

# farm order = F9,F2,F3,F4,F5,F6,F7,F8,F10
ivm_conc <- c(1.51,11.27,9.15,11.27,297.9,619,312.1,300,259.4)
bz_conc <- c(9.67,19.93,15,7.71,15,61.57,29.6,NA,29.6)
lev_conc <- c(1.5,1.5,0.58,0.57,0.91,0.96,9.36,NA,0.88)
mox_conc <- c(26.6,221.4,166.1,205.3,6038,10000,3987,NA,5263)


data <- read.table("XQTL_US_FIELD.merged.fst",header=F)
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


ggplot(data_out,aes(V2,mox_slopes)) + geom_point() + facet_grid(V1~.)
```
