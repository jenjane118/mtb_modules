## find proportion of Mtb genome that is transcribed

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

