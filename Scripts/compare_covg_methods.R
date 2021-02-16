## Jennifer J. Stiens
## 28/10/20


## script to cycle through all samples (.bam files) and make boxplots from entire
## dataset
## uses genomic alignments method used in baerhunter (see quantile_calc.R)

#cycle through all bam files in a dataset and make boxplot
# 1) make alignment of .bam
# 2) make coverage vector for each strand for each .bam
# 3) store covg info in pos and neg strand dataframes
# 3) can't unlist all covg info--too big. must use compressed cvg vector and apply 
# percentile function and store in dataset of percentiles (can run again for hi percentiles)
# 4) use pos and neg strand percentile info for entire dataset (12 rows) 
#     to make boxplot of percentiles for each strand

library(IRanges)
library(Rsamtools)
library(GenomicAlignments)


lo_probs <- c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)
all_probs <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.8, 0.9, 0.95)

dataset <- "PRJEB65014_3"

files <- list.files("test_bams", pattern="*.bam$", full.names=TRUE, ignore.case = TRUE)
files


pos_per.df <- as.data.frame(matrix(0, nrow=length(all_probs), ncol=length(files)), 
                            row.names = all_probs)
head(pos_per.df)
neg_per.df <- as.data.frame(matrix(0, nrow=length(all_probs), ncol=length(files)),
                            row.names = all_probs)

for (i in 1:length(files)){
  f<-files[i]
  # paired end reads
  # reversely stranded
  file_alignment <- readGAlignmentPairs(f, strandMode = 2)

  # do this for each strand
  strands <- c("+", "-")
  for (j in 1:2){
    target_strand = strands[j]
    strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]
    strand_cvg <- coverage(strand_alignment)
    list_components <- names(strand_cvg)
    target <- c()
    if (length(list_components)==1) {
      target <- list_components
    } else {
      return(paste("Invalid BAM file:",f, sep = " "))
    }
    vals <- runValue(strand_cvg)

  # add to relevant dataframe
    if (target_strand == "+"){
      per_row <- lapply(vals, quantile, probs=all_probs)
      #add to pos_cvg.df
      pos_per.df[,i]<-per_row
    }else{
      per_row <- lapply(vals, quantile, probs=all_probs)
      #add to neg_cvg.df
      neg_per.df[,i]<-per_row
    }
  }
}
head(pos_per.df)
head(neg_per.df)

write.table(pos_per.df, file=paste(dataset, "_pos.tsv", sep=""), quote=F, sep=" ")
write.table(neg_per.df, file=paste(dataset, "_neg.tsv", sep=""), quote=F, sep=" ")

# make boxplots of each strand
# PRJEB65014_3
# positive strand boxplot

pos_per.df<-read.delim("percentiles/PRJEB65014_3_pos.tsv", sep = " ",
                       stringsAsFactors = F)
neg_per.df<-read.delim("percentiles/PRJEB65014_3_neg.tsv", sep = " ",
                       stringsAsFactors = F)

# t(pos_per.df)
# 
# par(mfrow=c(1,1))
boxplot(t(pos_per.df), show.names=TRUE,
        main = paste("Percentile distrbution in Dataset",dataset),
        ylab="Number of reads in given percentile (raw)",
        ylim = c(1,500),
        col = "light blue",
        xlab= "Percentiles of number of reads"
)
par(new=T)
boxplot(t(neg_per.df), show.names=TRUE,
        col = "pink",
        ylim = c(1,500)
       # ylab="Number of reads in given percentile (raw)"
       # xlab= paste("Percentiles of number of reads in dataset ", dataset, sep="")
)
legend(0.3, 300, legend=c("positive", "negative"), col=c("light blue", "pink"), pch=20, cex = 1.0)
# 
# use percentile dataframes from running on server to make boxplots:

dataset="PRJNA327080_15"

pos_per.df<-read.delim("percentiles/PRJNA327080_15_pos.tsv", sep = " ",
                       stringsAsFactors = F)
neg_per.df<-read.delim("percentiles/PRJNA327080_15_pos.tsv", sep = " ",
                       stringsAsFactors = F)


boxplot(t(pos_per.df), show.names=TRUE,
        main ="BH_GAlignment method, positive strand",
        ylab="Number of reads in given percentile (raw)",
        xlab= paste("Percentiles of number of reads in dataset ", dataset, sep="")
)

boxplot(t(neg_per.df), show.names=TRUE,
        main ="BH_GAlignment method, negative strand",
        ylab="Number of reads in given percentile (raw)",
        xlab= paste("Percentiles of number of reads in dataset ", dataset, sep="")
)

dataset = "PRJNA390669_12"

pos_per.df<-read.delim("percentiles/PRJNA390669_12_pos.tsv", sep = " ",
                       stringsAsFactors = F)
neg_per.df<-read.delim("percentiles/PRJNA390669_12_neg.tsv", sep = " ",
                       stringsAsFactors = F)

boxplot(t(pos_per.df), show.names=TRUE,
        main = paste("Percentile distrbution in Dataset",dataset),
        ylab="Number of reads in given percentile (raw)",
        ylim = c(1,4000),
        col = "light blue",
        xlab= "Percentiles of number of reads"
)
par(new=T)
boxplot(t(neg_per.df), show.names=TRUE,
        col = "pink",
        ylim = c(1,4000)
        # ylab="Number of reads in given percentile (raw)"
        # xlab= paste("Percentiles of number of reads in dataset ", dataset, sep="")
)
legend(0.3, 3800, legend=c("positive", "negative"), col=c("light blue", "pink"), pch=20, cex = 1.0)




## make plot of positive strands using samtools covg and bh-ga for same dataset
# calculate sam_per_strand.df (per_strand.df) from samtools_cvg.R
# calculate bh_per_strand.df (hi_df) from quantile_calc.R


sam_per_strand.df <- per_strand.df
bh_per_strand.df <- hi_df

plot(sam_per_strand.df$pos, bh_per_strand.df$pos,
     cex = 0.6, col="blue", 
     xlab="Percentiles of reads per bp samtools-depth",
     xlim = c(0, 4200),
     xaxp = c(0, 4000, 4),
     ylab = "Percentiles of reads per bp BH-GAlignments",
     ylim = c(0, 4200),
     main = "Comparison of percentiles for BH-GA and samtools-depth for sample SRR5689224",
     type="p"
)

plot(all_probs, sam_per_strand.df$pos,
     cex = 0.6, col="blue",
     ylim = c(0, 4200),
     xlab = "Percentiles of reads per bp",
     ylab = "Number of reads in given percentile (raw)",
     main = "Comparison of percentiles for BH-GA and samtools-depth for sample SRR5689224",
     type = 'p')
points(all_probs, bh_per_strand.df$pos, col="red", cex = 0.6)
legend(0.2, 3500, col = c("blue", "red"), 
       pch = 1,
       c("samtools-depth", "BH-GAlignments"), cex=0.8)
