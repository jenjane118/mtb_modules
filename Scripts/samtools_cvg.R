# make barplots from coverage data for individual strands 

# coverage files generated using samtools depth



files <- list.files("PRJNA390669_12/pos_cov", pattern="*.coverage", full.names=TRUE)
files

nchar("SRR5689235")

# create new dataframe to store all coverage info
#using readr::read_tsv as it is faster than read.csv
temp <- do.call(cbind,lapply(files,readr::read_tsv, 
                col_names=c("Chromosome", "pos", "cov")))
#extract coverage column from each sample for pos cov df
ncol(temp)
pos_cov_df <- temp[,c(1:3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36)]
head(pos_cov_df)
rm(temp)

colnames(pos_cov_df) <- c("Chromosome", "pos", substr(files, 24, 33))
head(pos_cov_df)
ncol(pos_cov_df)

#show that the lower percentiles are always within a given small range regardless of the remaining values
#produces sample quantiles corresponding to the given probabilities. 
# quantile returns estimates of underlying distribution quantiles based on one or two order statistics from
# the supplied elements in x at probabilities in probs.
p_hi <- c(0.05, 0.1, 0.15, 0.20, 0.25, 0.50, 0.75, 1.0)
p_lo <- c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)


lapply(pos_cov_df[,3:14], quantile, probs=c(0.05, 0.1, 0.15, 0.20, 0.25, 0.50, 0.75, 1.0))
lapply(pos_cov_df[,3:14], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20))

#save the percentiles in a df
pos_percent.df <- as.data.frame(lapply(pos_cov_df[,3:14], 
                        quantile, probs=p_lo))
head(pos_percent.df)
write.table(pos_percent.df, "pos_percentiles_PRJNA390669.tsv", quote=F, sep="\t")

# boxplot
#each boxplot is the distribution of number of reads at given percentile for the 12 samples in the dataset

#par(yaxp  = c(0, 150, 15)) 
boxplot(t(pos_percent.df), show.names=TRUE, ylim=c(0,100),
        ylab="Number of reads in given percentile (raw)", 
        xlab="Percentiles of number of reads PRJNA390669_12 samples",
        yaxp  = c(0, 100, 10),
        main = "positive strand reads"
        )

png(file="PRJNA390669_pos_cvg.png", width=480, height = 480)
boxplot(t(pos_percent.df), show.names=TRUE, ylim=c(0,100),
        ylab="Number of reads in given percentile (raw)", 
        xlab="Percentiles of number of reads PRJNA390669_12 samples",
        yaxp  = c(0, 100, 10),
        main = "positive strand reads"
)
dev.off()

## for higher percentiles
hi_pos_percent.df <- as.data.frame(lapply(pos_cov_df[,3:14], 
                                       quantile, probs=p_hi))

png(file="hi_PRJNA390669_pos_cvg.png")
boxplot(t(hi_pos_percent.df), show.names=TRUE,
        ylab="Number of reads in given percentile (raw)", 
        xlab="Percentiles of number of reads PRJNA390669_12 samples",
        main = "positive strand reads"
)
dev.off()

boxplot(t(hi_pos_percent.df), show.names=TRUE,
        ylab="Number of reads in given percentile (raw)", 
        xlab="Percentiles of number of reads PRJNA390669_12 samples",
        main = "positive strand reads"
)

## repeat for negative strand

files <- list.files("PRJNA390669_12/neg_cov", pattern="*.coverage", full.names=TRUE)
files
# create new dataframe to store all coverage info
#using readr::read_tsv as it is faster than read.csv
temp <- do.call(cbind,lapply(files,readr::read_tsv, 
                             col_names=c("Chromosome", "pos", "cov")))
#extract coverage column from each sample for pos cov df
ncol(temp)
neg_cov_df <- temp[,c(1:3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36)]
head(neg_cov_df)
rm(temp)

colnames(neg_cov_df) <- c("Chromosome", "pos", substr(files, 24, 33))
head(neg_cov_df)
ncol(neg_cov_df)

#show that the lower percentiles are always within a given small range regardless of the remaining values
#produces sample quantiles corresponding to the given probabilities. 
# quantile returns estimates of underlying distribution quantiles based on one or two order statistics from
# the supplied elements in x at probabilities in probs.
p_hi <- c(0.1, 0.15, 0.20, 0.25, 0.50, 0.75, 0.95)
p_lo <- c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)


lapply(neg_cov_df[,3:14], quantile, probs=c(0.05, 0.1, 0.15, 0.20, 0.25, 0.50, 0.75, 1.0))
lapply(neg_cov_df[,3:14], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20))

#save the percentiles in a df
neg_percent.df <- as.data.frame(lapply(neg_cov_df[,3:14], 
                                       quantile, probs=p_lo))
write.table(neg_percent.df, "neg_percentiles_PRJNA390669.tsv", quote=F, sep="\t")

# boxplot
#each boxplot is the distribution of number of reads at given percentile for the 12 samples in the dataset

boxplot(t(neg_percent.df), show.names=TRUE, ylim=c(0,100),
        main ="negative strand reads",
        ylab="Number of reads in given percentile (raw)", 
        xlab="Percentiles of number of reads PRJNA390669_12 samples",
        yaxp  = c(0, 100, 10)
)

png(file="PRJNA390669_neg_cvg.png")
boxplot(t(neg_percent.df), show.names=TRUE, ylim=c(0,100),
        ylab="Number of reads in given percentile (raw)", 
        xlab="Percentiles of number of reads PRJNA390669_12 samples",
        yaxp  = c(0, 100, 10),
        main = "negative strand reads"
)
dev.off()

# for higher percentage

hi_neg_percent.df <- as.data.frame(lapply(neg_cov_df[,3:14], 
                                       quantile, probs=p_hi))
boxplot(t(hi_neg_percent.df), show.names=TRUE, ylab="Number of reads in given percentile (raw)", 
        xlab="Percentiles of number of reads PRJNA390669_12 samples",
        main = "negative strand reads"
        )

png(file="hi_PRJNA390669_neg_cvg.png")
boxplot(t(hi_neg_percent.df), show.names=TRUE, ylab="Number of reads in given percentile (raw)", 
        xlab="Percentiles of number of reads PRJNA390669_12 samples",
        main = "negative strand reads"
)
dev.off()


####

# do for single sample
## make plot of positive/neg strands using samtools covg for same dataset
#positive strand
all_probs<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
# import cvg files into dataframe
temp<-readr::read_tsv("PRJNA390669_12/pos_cov/SRR5689224_sorted_pos.coverage",
                      col_names=c("Chromosome", "pos", "cov"))
pos_strand.df<-temp[,c(2:3)]
per_strand.df<-as.data.frame(quantile(pos_strand.df$cov, probs=all_probs))
head(per_strand.df)

rm(temp)

temp<-readr::read_tsv("PRJNA390669_12/neg_cov/SRR5689224_sorted_neg.coverage",
                      col_names = c("Chromosome", "pos", "cov"))
neg_strand.df<-temp[,c(2:3)]
per_strand.df<-cbind(per_strand.df, quantile(neg_strand.df$cov, probs=all_probs))
colnames(per_strand.df)<-c("pos", "neg")
per_strand.df
rm(temp)



plot(all_probs, per_strand.df$pos,
     cex = 0.6, col="blue", 
     xlab="Percentiles of reads per bp",
     ylab = "Number of reads in a given percentile (Raw)",
     main = "Strand specific distribution of bp coverage",
     type="p"
)
points(all_probs, per_strand.df$neg, col="red", cex = 0.6)

