## baerhunter implementation and testing
## uses functions from baerhunter_code.R

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

# 1) write function that finds percentage of genome not likely to be transcribed using .gffs

# write function that uses gff as parameter
# output: % of genome likely to be non-transcribed and therefore can be used to estimate noise

mtb_gff<-read.gff("MtbH37RvNC_000962.3.gff")
head(mtb_gff)

#test data
#mtb_gff <- read.gff("short_mtb.gff")

#bovis_gff<-read.gff("LT708304_updated_aug19.gff")
#head(bovis_gff)
#rm(bovis_gff)

# make df selecting only coding regions
# this is likely to be problematic using different .gff formats
# original script uses 'parent' to exclude the redundant features
# but this designation not in all gffs
# 'type' column might be more useful designation of 'major features'?

# using other elements of attribute column will be problematic. The ID=gene is 
# also not present in bovis gff, whereas 'type' is

cd_mtb_gff <- mtb_gff[which(mtb_gff$type=="CDS"),]
View(cd_mtb_gff)

# subset non cds elements
nc_mtb_gff <- mtb_gff[which(mtb_gff$type %in% c("tRNA", "rRNA", "sRNA")),]
View(nc_mtb_gff)

# get total bp from first row of gff file 
# this can be type:'remark', 'source' or source: 'annotation' or 'feature'
# or get from header: ##sequence-region LT708304.1 1 4349904?
total_genome_len <- mtb_gff[1,5]
total_genome_len
#test with 13,000bp
#total_genome_len <- 13000

#need to use genomic ranges to determine length of non-overlapped regions

#create grange object with subsetted gff for only coding genes
mtb_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(cd_mtb_gff$start, end = cd_mtb_gff$end),
  strand   = Rle(strand(cd_mtb_gff$strand))
)
mtb_gr

#create grange object with nc rnas
nc_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(nc_mtb_gff$start, end = nc_mtb_gff$end),
  strand   = Rle(strand(nc_mtb_gff$strand))
  )

nc_gr

# need to add bp on each end of coding genes
five_p <-flank(mtb_gr, 50)
three_p <- flank(mtb_gr, 30, start=FALSE)
#take union of these 5' and 3' utrs and coding genes
m1 <- union(mtb_gr, five_p)
utrs_mtb_gr <- union(m1, three_p)
# now take union with ncRNAs
com_mtb_gr <- union(utrs_mtb_gr, nc_gr)

# to find regions that are not transcribed on either strand:

#gaps(g) should represent non-transcribed regions. need to specify to ignore strand
# have to use reduce method to do this first
# "the reduce method will align the ranges and merge overlapping ranges to produce 
# a simplified set."
red_mtb_gr <- reduce(com_mtb_gr, ignore.strand=T)
non_trx_mtb <- gaps(red_mtb_gr)
non_trx_mtb
# get width of each gap
# doesn't use gap from last gene to end of genome, need to add this in manually?
w<-width(non_trx_mtb)
# add widths together to get length of non-transcribed region
sum(w)

# add in distance from last entry in red_mtb_gr to end of genome
red_df<-as.data.frame(red_mtb_gr)
#last gap is end of last row to total genome length
last_gap <- total_genome_len - red_df[nrow(red_df), 3]

(sum(w) + last_gap)/total_genome_len
# [1] 0.06431417
# looks like about 6.4% non-transcribed at all?

###################################################

# to identify amount transcribed on each strand:
# use gr object that has all nc and coding and utrs, with strand info
# find gaps
str_spec_nontrx <- gaps(com_mtb_gr)
str_spec_nontrx

# need to calculate widths for each strand
gaps_neg_strand <- str_spec_nontrx[strand(str_spec_nontrx) == "-"]
w<-width(gaps_neg_strand)
# add final gap
com_mtb_df <- as.data.frame(com_mtb_gr)
neg_ends <- com_mtb_df[which(com_mtb_df$strand == "-"),3]
last_gap_neg<- total_genome_len - max(neg_ends) -1
(sum(w) + last_gap_neg)/total_genome_len
# 0.5186652
gaps_pos_str <- str_spec_nontrx[strand(str_spec_nontrx) == "+"]
W<-width(gaps_pos_str)
# add final gap
pos_df <- as.data.frame(gaps_pos_str)
pos_ends <- com_mtb_df[which(com_mtb_df$strand == "+"), 3]
last_gap_pos <- total_genome_len - max(pos_ends) -1
(sum(W) + last_gap_pos)/total_genome_len
# 0.5377225

# overall non-transcribed percentage
(sum(w) + sum(W) + last_gap_neg + last_gap_pos) / (2*total_genome_len)
# 0.5281938

#######################################################

# function to find coverage quantiles of individual bam files

#BiocManager::install("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
#library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

# x is a class of GAlignments, Each row is a read
f <- "SRR5689224_sorted.bam"

# single end reads (need to get an example file)
x <- readGAlignments(f)
#Coverage of the reads : this will generate a RleList
#Rle is a run-length encoding
xcov <- coverage(x)
head(xcov)
xnum <- as.numeric(xcov$AL123456.3)   #Uncompress the coverage
head(xnum)


# paired end reads
# reversely stranded
# need to run this for each strand 'target_strand'
target_strand = "-"
file_alignment <- readGAlignmentPairs(f, strandMode = 2)

## do this for each strand
strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]

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
vals <- runValue(strand_cvg)
percent.df <- as.data.frame(quantile(vals, probs=c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20)))
head(percent.df)
## choose cutoff based on determination of untranscribed region of genome above
## for mtb(above) choose 10-20% for cuttoff?
#10%
low_cutoff <- percent.df[6,1]
# 20%
hi_cutoff <- percent.df[16,1]
low_cutoff; hi_cutoff

#repeat for both strands and put in dataframe
cutoff_df <- as.data.frame(matrix(0, nrow=nrow(percent.df), ncol = 2))
rownames(cutoff_df)<-rownames(percent.df)
colnames(cutoff_df)<-c("+", "-")
if (target_strand == "+"){
  cutoff_df$`+`<-percent.df$AL123456.3
}else{
    cutoff_df$`-`<-percent.df$AL123456.3
  }

cutoff_df
t(cutoff_df)
target_strand

boxplot(t(cutoff_df), show.names=TRUE,ylab="Number of reads in given percentile (raw)", xlab=c("Percentiles of number of reads in sample"), main=c("read distr for sense/antisense strands of one sample"))

