if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("IRanges")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicAlignments")
install.packages("ape")
BiocManager::install("GenomicRanges")


library(IRanges)
library(Rsamtools)
library(GenomicAlignments)
library(ape)
library(GenomicRanges)


# function to find coverage quantiles of individual bam files

# x is a class of GAlignments, Each row is a read
f <- "PRJNA327080_15/SRR3725585_sorted.bam"

#x <- readGAlignments(f)
#Coverage of the reads : this will generate a RleList
#Rle is a run-length encoding
#xcov <- coverage(x)
#head(xcov)
#xnum <- as.numeric(xcov$AL123456.3)   #Uncompress the coverage will freeze



# paired end reads
# reversely stranded
file_alignment <- readGAlignmentPairs(f, strandMode = 2)

## do this for - strand
target_strand = "-"
strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]
strand_alignment
#ycov_pos <- coverage(y)
#head(ycov_pos)
#ynum_pos <- as.numeric(ycov_pos$AL123456.3)
# this causes computer to freeze

## Create a strand coverage vector and extract coverage values for each strand
## IRanges function 'coverage' counts the number of ranges over each bp.
strand_cvg <- coverage(strand_alignment)
list_components <- names(strand_cvg) ## "AL123456.3"
target <- c()
if (length(list_components)==1) {
  target <- list_components
} else {
  return(paste("Invalid BAM file:",f, sep = " "))
}

# makes an atomic list from covg values
vals <- runValue(strand_cvg)

max(vals)
#2898199 
range(vals)
median(vals)
mean(vals)


lo_probs <- c(0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
all_probs <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.8, 0.9, 0.95)

lo_percent.df <- as.data.frame(quantile(vals, probs=lo_probs))
hi_percent.df <- as.data.frame(quantile(vals, probs=all_probs))

lo_df <- as.data.frame(matrix(0, 
                              nrow=nrow(lo_percent.df), ncol = 2))
colnames(lo_df)<-c("pos", "neg")
rownames(lo_df)<-rownames(lo_percent.df)
lo_df$neg<-lo_percent.df$AL123456.3

hi_df <- as.data.frame(matrix(0,
                              nrow=nrow(hi_percent.df), ncol = 2))
colnames(hi_df)<-c("pos", "neg")
rownames(hi_df)<-rownames(hi_percent.df)
hi_df$neg<-hi_percent.df$AL123456.3

lo_df
hi_df

## repeat for pos strand
target_strand = "+"
strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]
strand_cvg <- coverage(strand_alignment)
list_components <- names(strand_cvg) ## "AL123456.3"
target <- c()
if (length(list_components)==1) {
  target <- list_components
} else {
  return(paste("Invalid BAM file:",f, sep = " "))
}
vals <- runValue(strand_cvg)
head(vals)
cvg<-vals$AL123456.3
head(vals[[target]])

head(cvg)
lo_percent.df <- as.data.frame(quantile(vals, probs=lo_probs))
hi_percent.df <- as.data.frame(quantile(vals, probs=all_probs))
head(lo_percent.df)
head(hi_percent.df)


# add to rel dataframe
lo_df$pos<-lo_percent.df$AL123456.3
hi_df$pos<-hi_percent.df$AL123456.3
lo_df
hi_df


plot(lo_probs, lo_df$pos,
     cex = 0.6, col="blue", 
     xlab="Percentiles of reads per bp",
     ylab = "Number of reads in a given percentile (Raw)",
     main = "Strand specific distribution of bp coverage",
)
points(lo_probs, lo_df$neg, col="red", cex = 0.6)
abline(a=-12, b=576, col="black")





plot(all_probs, hi_df$pos,
     cex = 0.6, col="blue", 
     xlab="Percentiles of reads per bp",
     ylab = "Number of reads in a given percentile (Raw)",
     main = "Strand distr of covg sample SRR5689224, BH-GAligment method",
     cex.main = 0.9,
     type="p"
     )
points(all_probs, hi_df$neg, col="red", cex = 0.6)
legend(0.3, 3500, legend=c("pos", "neg"), pch=1, col=c("blue","red"))

# find slope between 5-10%
pl<-(lo_df$pos[2]-lo_df$pos[1])/0.05
pl  #420
# for sample from PRJNA327080
#120

# find slope between 5-30%
sl<-(lo_df$pos[6]-lo_df$pos[1])/0.25
sl    #576
# slope between 25-30
ml<-(lo_df$pos[6]-lo_df$pos[5])/0.05
ml   #740
# slope between 30-35
sl<-(lo_df$pos[7]-lo_df$pos[6])/0.05
sl   #900
# for sample from PRJNA327080 #200
# slope between 30-50
ns<-(lo_df$pos[10]-lo_df$pos[6])/0.20
ns    #1110
# for sample from PRJNA327080 270

# plot differences between percentiles of 5%
x<-seq(10,50, 5)
x
diff_10<-diff(lo_df$pos)
diff_10_neg<-diff(lo_df$neg)
diff_10  #21 24 29 33 37 45 50 59 68
plot(x, diff_10,
     cex = 0.6, col="blue", 
     main = "Plot of differences between percentiles",
     xlab="Percentiles of reads per bp",
     ylab = "Differences in number of reads between percentiles",
     cex.main = 0.9,
     type="p"
     )
points(x, diff_10_neg, cex=0.6, col="red")


# when difference between percentiles more than double over 5%, make threshold?
# different threshold for each strand
# make diff vector diff_5p<-diff(lo_df$pos)
# if i > (i-1)*2, threshold=lo_df[i]

diff_5p<-diff(lo_df$pos)
diff_5p
low_coverage_cutoff <- 5  #default value
found<-FALSE
# difference between 5-10% as baseline difference
baseline<-lo_df$pos[2]-lo_df$pos[1]
for (i in 1:length(diff_5p)){
  if (found != TRUE & diff_5p[i] > baseline*2){
    low_coverage_cutoff <- lo_df$pos[i]
    found<-TRUE
    }
}
low_coverage_cutoff 
#257
## for sample from PRJNA327080 81
rownames(lo_df)[lo_df$pos == low_coverage_cutoff]
#"40%"
# # for sample from PRJNA327080 45%

# do again with neg strand
diff_5n<-diff(lo_df$neg)
diff_5n
low_coverage_cutoff <- 5  #default value
found<-FALSE
# difference between 5-10% as baseline difference
baseline<-lo_df$neg[2]-lo_df$neg[1]
for (i in 1:length(diff_5n)){
  if (found != TRUE & diff_5n[i] > baseline*2){
    low_coverage_cutoff <- lo_df$neg[i]
    found == TRUE
  }
}
low_coverage_cutoff 
#302
# for sample from PRJNA327080 #81
rownames(lo_df)[lo_df$neg == low_coverage_cutoff]
#45%
# for sample from PRJNA327080 45%




hi_df
p<-c(0.95,0.96,0.97,0.98,0.99)
hi<-quantile(vals, probs=p)
hi
plot(p,hi)


# plot pos vs neg strand
plot(hi_df$neg, hi_df$pos,
     cex = 0.6,
     col = "red",
     cex.main = 0.8,
     main = "Strand coverage percentile comparison
              0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.8, 0.9, 0.95",
     xlab = "Number of reads per percentile, Negative strand",
     xlim = c(0, 4200),
     xaxp = c(0, 4200, 40),
     ylab = "Number of reads per percentile, Positive strand",
     yaxp = c(0, 4200, 40),
     type="p"
)

par(mfrow=c(2,2))
png("strand_perc_compare", width=480, height = 500)
plot(hi_df$neg, hi_df$pos,
     cex = 0.6,
     col = "red",
     cex.main = 0.8,
     main = "Strand coverage percentile comparison
              0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.8, 0.9, 0.95",
     xlab = "Number of reads per percentile, Negative strand",
     xlim = c(0, 4200),
     ylim = c(0, 4200),
     xaxp = c(0, 4200, 20),
     ylab = "Number of reads per percentile, Positive strand",
     yaxp = c(0, 4200, 20),
     type="p"
)
dev.off()

# to make custom axes
#axis(side, at=, labels=, pos=, lty=, col=, las=, tck=, ...)

hi_df$pos[5]
hi_df$neg[5]
