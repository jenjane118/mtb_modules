## Jennifer J Stiens
## 19/10/20

## work with bedgraph coverage files
## make plot with each sample from PRJNA390669_12


files <- list.files("PRJNA390669_12", pattern="*.bg", full.names=TRUE)
files

# create new dataframe to store all bedgraph info
#using readr::read_tsv as it is faster than read.csv
temp <- do.call(cbind,lapply(files,readr::read_tsv, col_names=c("Chromosome", "pos", "cov")))
#extract coverage column from each sample for df
ncol(temp)
cov_df <- temp[,c(1:3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36)]
head(cov_df)
rm(temp)
#colnames(cov.df) <- c("Chromosome", "pos", gsub(pattern='_sorted_cov.bg$', replacement="", files))
colnames(cov_df) <- c("Chromosome", "pos", substr(files, 16, 25))
head(cov_df)
ncol(cov_df)

#show that the lower percentiles are always within a given small range regardless of the remaining values
#produces sample quantiles corresponding to the given probabilities. 
# quantile returns estimates of underlying distribution quantiles based on one or two order statistics from
# the supplied elements in x at probabilities in probs.
lapply(cov_df[,3:14], quantile, probs=c(0.05, 0.1, 0.15, 0.20, 0.25, 0.50, 0.75, 1.0))
lapply(cov_df[,3:14], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20))

#save the percentiles in a df
percent.df <- as.data.frame(lapply(cov_df[,3:14], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)))
View(percent.df)
write.table(percent.df, "percentiles_PRJNA390669.tsv", quote=F, sep="\t")

#keep the number of reads mapped for each sample (from samtools flagstat)
depths <- matrix(nrow=12, ncol=1)
rownames(depths) <- colnames(cov_df)[3:14]
head(depths)

# match mapped read numbers from flagstat files
library(stringr)

files <- list.files("PRJNA390669_12/PRJNA390669_flagstats", pattern="*.txt", full.names=TRUE)

temp_line=NULL
for (i in 1:length(files)){
  temp_line<-readLines(files[i], n=5)[5]
  depths[i,1]<-str_match(temp_line, "^([0-9]+) \\+ .*")[, 2]
}

depths
depths_df<-as.data.frame(depths, stringsAsFactors = F)
depths_df$V1<-as.numeric(depths_df$V1)
class(depths_df$V1)

boxplot(t(percent.df), show.names=TRUE, ylim=c(0,150),ylab="Number of reads in given percentile (raw)", xlab="Percentiles of number of reads across 12 M.tb samples" )
t(percent.df)

library(beeswarm)
beeswarm(as.data.frame(t(percent.df)), cex=1, add=T, vertical=T)

#each boxplot is the distribution of number of reads at given percentile for the 12 samples in the dataset
png(file="boxplots_low_percentiles.png")
boxplot(t(as.data.frame(lapply(cov_df[,3:14], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)))), show.names=TRUE, ylim=c(0,150),  ylab="Number of reads in given percentile (raw)", xlab="Percentiles of number of reads across 12 M.tb samples")
dev.off()

boxplot(t(as.data.frame(lapply(cov_df[,3:14], quantile, probs=c(0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65)))), show.names=TRUE, ylim=c(0,800))


##################################################################
### repeat with another dataset for comparison
files <- list.files("PRJNA327080_15", pattern="*.bg", full.names=TRUE)
files

# create new dataframe to store all bedgraph info
#using readr::read_tsv as it is faster than read.csv
temp <- do.call(cbind,lapply(files,readr::read_tsv, col_names=c("Chromosome", "pos", "cov")))
#extract coverage column from each sample for df
ncol(temp)
cov_df_15 <- temp[,c(1:3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45)]
head(cov_df_15)
rm(temp)
#colnames(cov.df) <- c("Chromosome", "pos", gsub(pattern='_sorted_cov.bg$', replacement="", files))
colnames(cov_df_15) <- c("Chromosome", "pos", substr(files, 16, 25))
head(cov_df_15)
ncol(cov_df_15)

#show that the lower percentiles are always within a given small range regardless of the remaining values
#produces sample quantiles corresponding to the given probabilities. 
# quantile returns estimates of underlying distribution quantiles based on one or two order statistics from
# the supplied elements in x at probabilities in probs.
lapply(cov_df_15[,3:17], quantile, probs=c(0.05, 0.1, 0.15, 0.20, 0.25, 0.50, 0.75, 1.0))
lapply(cov_df_15[,3:17], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20))

#save the percentiles in a df
percent_df_15 <- as.data.frame(lapply(cov_df_15[,3:17], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)))
View(percent_df_15)
write.table(percent_df_15, "percentiles_PRJNA327080.tsv", quote=F, sep="\t")

#keep the number of reads mapped for each sample (from samtools flagstat)
depths <- matrix(nrow=15, ncol=1)
rownames(depths) <- colnames(cov_df_15)[3:17]
head(depths)

# match mapped read numbers from flagstat files
library(stringr)

files <- list.files("PRJNA327080_15/PRJNA327080_flagstats", pattern="*.txt", full.names=TRUE)
files
for (i in 1:length(files)){
  temp_line<-readLines(files[i], n=5)[5]
  depths[i,1]<-str_match(temp_line, "^([0-9]+) \\+ .*")[, 2]
}
depths
depths_df<-as.data.frame(depths, stringsAsFactors = F)
depths_df$V1<-as.numeric(depths_df$V1)
class(depths_df$V1)

boxplot(t(percent_df_15), show.names=TRUE, ylim=c(0,60),ylab="Number of reads in given percentile (raw)", xlab="Percentiles of number of reads across 15 M.tb samples (PRJNA327080)",  )
t(percent_df_15)

#####################################################

#do coverage graph of one sample with cds only

temp <- do.call(cbind,lapply("SRR5689224_cds_cov.bg",readr::read_tsv))
#extract coverage column for df
ncol(temp)
rm(temp)
# i want second to last column for reads, but use type to filter for CDS
cov_12_cds<-NULL
cov_12_cds <- temp[,c(1,3, 12)]
colnames(cov_12_cds) <- c("chrom", "type", "reads")
head(cov_12_cds)
cov_12_cds <- cov_12_cds[which(cov_12_cds$type=="CDS"),]
head(cov_12_cds)

# make percentiles of one sample
quantile(cov_12_cds[,3], probs=c(0.05, 0.1, 0.15, 0.20, 0.25, 0.50, 0.75, 1.0))
#save the percentiles in a df
cds_percent_df <- as.data.frame(quantile(cov_12_cds[,3], probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)))
#make boxplot
boxplot(t(cds_percent_df), show.names=TRUE, ylim=c(0,600),ylab="Number of reads in given percentile (raw)", xlab="Percentiles of number of reads in one sample",  )


##################################################################

# make positive and negative strand barplots

### repeat with another dataset for comparison
files_pos <- list.files("PRJNA390669_12/pos_cov/", pattern="*.bg", full.names=TRUE)

files_neg <- list.files("PRJNA390669_12/neg_cov/", pattern="*.bg", full.names=TRUE)


# create new dataframe to store all bedgraph info
#using readr::read_tsv as it is faster than read.csv
temp <- do.call(cbind,lapply(files_neg,readr::read_tsv, col_names=c("Chromosome", "pos", "cov")))
#extract coverage column from each sample for df
ncol(temp)
poscov_df <- temp[,c(1:3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39)]
negcov_df <- temp[,c(1:3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39)]
head(negcov_df)
rm(temp)

#pos strand df
colnames(poscov_df) <- c("Chromosome", "pos", substr(files, 25, 34))
head(poscov_df)
ncol(poscov_df)
#save the percentiles in a df
pos_percent_df <- as.data.frame(lapply(poscov_df[,3:14], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)))
View(pos_percent_df)
#write.table(pos_percent_df, "percentiles_PRJNA390669_pos.tsv", quote=F, sep="\t")
boxplot(t(pos_percent_df), show.names=TRUE, ylim=c(0,150),ylab="Number of reads in given percentile (raw)", xlab="Percentiles of number of reads on + strand, PRJNA390669")
t(pos_percent_df)

#neg strand df
colnames(negcov_df) <- c("Chromosome", "pos", substr(files, 25, 34))
head(negcov_df)
neg_percent_df <- as.data.frame(lapply(negcov_df[,3:14], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)))
View(neg_percent_df)
#write.table(neg_percent_df, "percentiles_PRJNA390669_neg.tsv", quote=F, sep="\t")
boxplot(t(neg_percent_df), show.names=TRUE, ylim=c(0,500),ylab="Number of reads in given percentile (raw)", xlab="Percentiles of number of reads on - strand, PRJNA390669")

##################################################################

# analyse cds coverage only
files_cds <- list.files("PRJNA390669_12/cds_cov/", pattern="*.bg", full.names=TRUE)
files_cds

temp <- do.call(cbind,lapply(files_cds,readr::read_tsv))
#extract coverage column from each sample for df
ncol(temp)
head(temp)
cds_12_df <- temp[,c(1, 3, 10, 23, 36, 49, 62, 75, 88, 101, 114, 127, 140, 153)]
head(cds_12_df)
rm(temp)

colnames(cds_12_df) <- c("Chromosome", "type", substr(files_cds, 25, 34))
head(cds_12_df)
# select only CDS entries
cds_12_df <- cds_12_df[which(cds_12_df$type=="CDS"),]
View(cds_12_df)
percents_cds_df<-NULL
percents_cds_df <- as.data.frame(lapply(cds_12_df[,3:14], quantile, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)))
View(percents_cds_df)
#write.table(neg_percent_df, "percentiles_PRJNA390669_neg.tsv", quote=F, sep="\t")
boxplot(t(percents_cds_df), show.names=TRUE, ylim=c(0,1500),ylab="Number of reads in given percentile (raw)", xlab="Percentiles of number of reads from coding only, PRJNA390669")

## all reads are zero

