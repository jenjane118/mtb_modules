# run apero
# https://github.com/Simon-Leonard/APERO

install.packages("reshape2")
install.packages("snowfall")
devtools::install_github("Simon-Leonard/APERO")

library(Rsamtools)
library(reshape2)
library(snowfall)

library(APERO)

# Load annotation file (needs to be in .ptt form?)
#ptt=read.csv("E:/APERO/NC_016810.ptt",sep="\t",skip=2,header=T,stringsAsFactors = F)

# 5'end detection
res=APERO_start_detection(work_dir = ".", bam_name = "~/bam_files/PRJNA327080_15/SRR3725585_sorted.bam",
                          wmax = 10, min_dist = 10, enrichment = 0.1, min_read_number = 0, genome_size = 4411532
                          )


#3'end detection
res2=APERO_end_detection(work_dir = "./APERO/", start_table = res, mTEX_bam = "~/bam_files/PRJNA327080_15/SRR3725585_sorted.bam",
                         readthrough_proportion = 0.01, Fmin=NA, thread_number = 8, genome_size = 4411532
                         )