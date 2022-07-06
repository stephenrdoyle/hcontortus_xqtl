## Peak finding

Author: Stephen Doyle
Contact: stephen.doyle[at]sanger.ac.uk

- didnt use this in the end, but keeping for future reference
```shell

working dir: /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/04_VARIANTS/XQTL_BZ

# get some testdata

 cut -f1,2,13 XQTL_BZ.merged.fst > peakfinding.testdata

 cut -f1,2,13 ../XQTL_CONTROL/XQTL_CONTROL.merged.fst > peakfinding.control.testdata

1. all data - genome wide average
     - pre-treatment vs post =
          cat peakfinding.testdata | datamash mean 3
          0.015070584074789

     - time matched control pre/post =
          cat peakfinding.control.testdata | datamash mean 3
          0.011916732789999


2. list of manually curated "peaks"



hcontortus_chr1_Celeg_TT_arrow_pilon	6992500	0.09456757


3. extract peak coordinates


XQTL
V13 1:5
V27 2:6
V39 3:7
V49 4:8

AI
#						Rep1	Rep2	Rep3
#control - pre v 0.5X	V17		V51		V83
#control - pre v 2X		V11		V135	V161

#IVM pre v 0.5x			V251	V267	V281
#IVM pre v 2x			V287	V297	V305


Dose response
V7 1:2

cut -f1,2,13 XQTL_*.merged.fst | sort -k3 | tail -n 100

#bz
hcontortus_chr1_Celeg_TT_arrow_pilon	6992500	0.09456757
hcontortus_chr1_Celeg_TT_arrow_pilon	8807500	0.08951118
hcontortus_chr1_Celeg_TT_arrow_pilon	12722500	0.07817395
#ivm-XQTL
hcontortus_chr5_Celeg_TT_arrow_pilon	34007500	0.04864623
hcontortus_chr5_Celeg_TT_arrow_pilon	36252500	0.08151205
hcontortus_chr5_Celeg_TT_arrow_pilon	37467500	0.07105946
hcontortus_chr5_Celeg_TT_arrow_pilon	45397500	0.04814783
hcontortus_chr2_Celeg_TT_arrow_pilon	2947500	0.05130916
#ivm-ai
hcontortus_chr1_Celeg_TT_arrow_pilon	11047500	0.07026732
hcontortus_chr3_Celeg_TT_arrow_pilon	38012500	0.10238732
#doseresponse
hcontortus_chr5_Celeg_TT_arrow_pilon	34397500	0.50000000




# XQTL_IVM
# Rep1.1 - V13
13   hcontortus_chr5_Celeg_TT_arrow_pilon	37467500	0.07105946
13   hcontortus_chr2_Celeg_TT_arrow_pilon	2947500	0.05130916
#Rep1.2 - V27
27   hcontortus_chr5_Celeg_TT_arrow_pilon	33842500	0.06678443
27   hcontortus_chr5_Celeg_TT_arrow_pilon	46717500	0.04836741
#Rep2 - V39
#Rep3 - V49
49   hcontortus_chr5_Celeg_TT_arrow_pilon	36317500	0.10849596
49   hcontortus_chr5_Celeg_TT_arrow_pilon	35822500	0.10612988

./run_find_peak_windows.sh XQTL_IVM peak.data XQTL_IVM.merged.fst 500000


# XQTL_LEV
# Rep1.1 - V13
13   hcontortus_chr4_Celeg_TT_arrow_pilon	14817500	0.07934698
13   hcontortus_chr5_Celeg_TT_arrow_pilon	31467500	0.10024437
#Rep1.2 - V27
27   hcontortus_chr5_Celeg_TT_arrow_pilon	22157500	0.06093688
27   hcontortus_chr5_Celeg_TT_arrow_pilon	28482500	0.08273666
#Rep2 - V39
#Rep3 - V49



# XQTL_BZ
# Rep1.1 - V13
13   hcontortus_chr1_Celeg_TT_arrow_pilon	6992500	0.09456757
13   hcontortus_chr1_Celeg_TT_arrow_pilon	8807500	0.08951118
13   hcontortus_chr1_Celeg_TT_arrow_pilon	12722500	0.07817395
# Rep1.2 - V27
# Rep2 - V39
# Rep3 - V49
49	hcontortus_chr1_Celeg_TT_arrow_pilon	8337500	0.23776054




# AI
# Rep1 - V287
287  hcontortus_chr1_Celeg_TT_arrow_pilon	11047500	0.07026732
# Rep 2 - V297
297  hcontortus_chr2_Celeg_TT_arrow_pilon	7147500	0.21312042
297  hcontortus_chr4_Celeg_TT_arrow_pilon	29077500	0.06586891
# Rep 3 - V305
305  hcontortus_chr3_Celeg_TT_arrow_pilon	38012500	0.10238732


# DOSE_RESPONSE
# rep 1 - V7
7    hcontortus_chr5_Celeg_TT_arrow_pilon	34052500	0.35451348


# XQTL_ADULTS
# rep 1 - V7
7    hcontortus_chr5_Celeg_TT_arrow_pilon	37377500	0.19438392
7    hcontortus_chr2_Celeg_TT_arrow_pilon	27157500	0.12806700




COL=7
DATA=XQTL_PARENTS.merged.fst

CUTOFF=$(cut -f1,2,${COL} *.merged.fst | sort -k3 | datamash mean 3 sstdev 3 | awk '{print $1+(3*$2)}')

cut -f1,2,${COL} *.merged.fst | sort -k3 | awk -v CUTOFF=$CUTOFF -v COL=$COL '{if($3>CUTOFF) print COL,$0}' OFS="\t"  > all.high.peaks

./run_find_peak_windows.sh XQTL_allhightest all.high.peaks ${DATA} 500000

sort -k1,1 -k2,2n  XQTL_allhightest.peak_windows | awk '{if($6!=0 && NF==8) print}' > sorted


R ${DATA} ${COL}
```
args <- commandArgs()

library(ggplot2)
library(viridis)

name=args[2]
columns=args[3]

a<-read.table("sorted",header=T)
ggplot(a)+
     geom_rect(aes(xmin=a$PEAK_START_COORD,ymin=1:nrow(a)-0.5,xmax=a$PEAK_END_COORD,ymax=1:nrow(a)+0.5,fill=PEAK_FST))+
     facet_grid(a$CHR~.)+
     xlim(0,50e6)+
     labs(title=paste0("file=",name,", ","column=",columns))+
     scale_fill_viridis(direction=-1,limits=c(0,0.5))+
     theme_bw()

ggsave(paste0("predictedpeaks_",name,"_","column_",columns,".pdf"))
```



- where "run_find_peak_windows.sh"
```shell
#!/bin/bash

PREFIX=$1
PEAKFILE=$2
FST_DATA=$3
WINDOW=$4

#eg.
#PREFIX=TEST
#PEAKFILE=peak.data
#FST_DATA=peakfinding.testdata
#WINDOW=500000


printf CHR"\t"PEAK_COORD"\t"PEAK_FST"\t"PEAK_START_COORD"\t"PEAK_END_COORD"\t"PEAK_WINDOW_SIZE"\t"PEAK_START_FST"\t"PEAK_END_FST"\n" > ${PREFIX}.peak_windows

while read COL CHR PEAK FST; do
     # ignore comment lines
     [[ "$COL" =~ ^# ]] && continue

     # extract data from master file
     cut -f1,2,${COL} ${FST_DATA} > data.tmp

     # extract a window of data around the peak
     grep "${CHR}" data.tmp | awk -v PEAK=$PEAK -v WINDOW=$WINDOW '{if($2>(PEAK-WINDOW) && $2< (PEAK+WINDOW)) print}' > data.tmp2;

     # extract from the windowed data Fst values that are half the peak height, taking into account the genome wide average
     GENOME_AVERAGE=$(cat data.tmp | datamash mean 3 sstdev 3 | awk '{print $1+(3*$2)}')

     awk -v GENOME_AVERAGE=$GENOME_AVERAGE -v FST=$FST '{if($3>((FST+GENOME_AVERAGE)/2)) print}' data.tmp2 > data.tmp3;

     # print the boundaries of the peak
     PEAK_START=$(head -n 1 data.tmp3 | cut -f2);
     PEAK_START_FST=$(head -n 1 data.tmp3 | cut -f3);
     PEAK_END=$(tail -n 1 data.tmp3 | cut -f2);
     PEAK_END_FST=$(tail -n 1 data.tmp3 | cut -f3);
     PEAK_WINDOW_SIZE=$(($PEAK_END-PEAK_START));
     printf ${CHR}"\t"${PEAK}"\t"${FST}"\t"$PEAK_START"\t"$PEAK_END"\t"${PEAK_WINDOW_SIZE}"\t"${PEAK_START_FST}"\t"${PEAK_END_FST}"\n" >> ${PREFIX}.peak_windows;
     rm *tmp*
     done < ${PEAKFILE}
```
