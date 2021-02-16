#Comparing GFF3 files. The created GFF3 files were compared for similarity, the code for subsetting ncRNAs and the Jaccard Index code were obtained from baerhunter scripts.

#r comparing .GFF3 files which were created using feature_file_editor of baerhunter package

library("Rsamtools")
library("rtracklayer")
#import .gff3s using rtracklayer. import 4 .gff3 files at a time to compare high expression filtered .gff3s
#for the dataset at the same time
annot_file1 <- file.path("E:/mtb_5_10_67035.GFF3")
annot_file2 <- file.path("E:/filtered_high_expression_ncRNA_mtb_5_10_67035.GFF3")
annot_file3 <- file.path("E:/mtb_5_10_47863.GFF3")
annot_file4 <- file.path("E:/filtered_high_expression_ncRNA_mtb_5_10_47863.GFF3")
annot <- import.GFF3(annot_file1)
annot2 <- import.GFF3(annot_file2)
annot3 <- import.GFF3(annot_file3)
annot4 <- import.GFF3(annot_file4)
# create subsets of the putative sRNAs and UTRs first
pred.sRNA1 <- subset(annot, type == "putative_sRNA")
pred.UTR1 <- subset(annot, type == "putative_UTR")
num.pred.sRNA.plus1 <- length(ranges(subset(annot, (type == "putative_sRNA") & (strand == "+")) ) )
num.pred.sRNA.minus1 <- length(ranges(subset(annot, (type == "putative_sRNA") & (strand == "-")) ) )
num.pred.UTRs.plus1 <- length(ranges(subset(annot, (type == "putative_UTR") & (strand == "+")) ) )
num.pred.UTRs.minus1 <- length(ranges(subset(annot, (type == "putative_UTR") & (strand == "-")) ) )
# create a set with both sRNAs and UTRs
pred.all1 <- subset(annot, (type == "putative_sRNA" | type == "putative_UTR"))


#2nd annot file
pred.sRNA2 <- subset(annot2, type == "putative_sRNA")
pred.UTR2 <- subset(annot2, type == "putative_UTR")
num.pred.sRNA.plus2 <- length(ranges(subset(annot2, (type == "putative_sRNA") & (strand == "+")) ) )
num.pred.sRNA.minus2 <- length(ranges(subset(annot2, (type == "putative_sRNA") & (strand == "-")) ) )
num.pred.UTRs.plus2 <- length(ranges(subset(annot2, (type == "putative_UTR") & (strand == "+")) ) )
num.pred.UTRs.minus2 <- length(ranges(subset(annot2, (type == "putative_UTR") & (strand == "-")) ) )
# create a set with both sRNAs and UTRs
pred.all2 <- subset(annot2, (type == "putative_sRNA" | type == "putative_UTR"))

# 3rd annot file
pred.sRNA3 <- subset(annot3, type == "putative_sRNA")
pred.UTR3 <- subset(annot3, type == "putative_UTR")
num.pred.sRNA.plus3 <- length(ranges(subset(annot3, (type == "putative_sRNA") & (strand == "+")) ) )
num.pred.sRNA.minus3 <- length(ranges(subset(annot3, (type == "putative_sRNA") & (strand == "-")) ) )
num.pred.UTRs.plus3 <- length(ranges(subset(annot3, (type == "putative_UTR") & (strand == "+")) ) )
num.pred.UTRs.minus3 <- length(ranges(subset(annot3, (type == "putative_UTR") & (strand == "-")) ) )
# create a set with both sRNAs and UTRs
pred.all3 <- subset(annot3, (type == "putative_sRNA" | type == "putative_UTR"))

#4th annot file
pred.sRNA4 <- subset(annot4, type == "putative_sRNA")
pred.UTR4 <- subset(annot4, type == "putative_UTR")
num.pred.sRNA.plus4 <- length(ranges(subset(annot4, (type == "putative_sRNA") & (strand == "+")) ) )
num.pred.sRNA.minus4 <- length(ranges(subset(annot4, (type == "putative_sRNA") & (strand == "-")) ) )
num.pred.UTRs.plus4 <- length(ranges(subset(annot4, (type == "putative_UTR") & (strand == "+")) ) )
num.pred.UTRs.minus4 <- length(ranges(subset(annot4, (type == "putative_UTR") & (strand == "-")) ) )
# create a set with both sRNAs and UTRs
pred.all4 <- subset(annot4, (type == "putative_sRNA" | type == "putative_UTR"))



#creating a hits object that contains the locations of overlaps between query and subject.
# query, subject. hits produced gives index of overlapping ranges, not actual ranges. can get ranges from
#actual objects using the index.
sRNA.hits <- findOverlaps(pred.sRNA1, pred.sRNA3, type="any",
                          minoverlap=1L, ignore.strand=FALSE)
head(sRNA.hits)
UTR.hits <- findOverlaps(pred.UTR3, pred.UTR1, type="any",
                         minoverlap=1L, ignore.strand=FALSE)

all.hits <- findOverlaps(pred.all3, pred.all1, type="any",
                         minoverlap=1L, ignore.strand=FALSE)
head(all.hits)
num_UTR_confirmed <- length(unique(subjectHits(UTR.hits)))
num_sRNA_confirmed <- length(unique(subjectHits(sRNA.hits)))

#how can I get the widths from overlapping ranges, and the widths of the ranges that the overlaps come from?
#one attempt was with intersect but this did not work so is hashed.
#gr3 <- intersect(pred.sRNA2, pred.sRNA4)
#width_overlaps <- width(gr3)

#can look at which ranges from predicted.sRNAs formed the hits object, s.RNA.hits. This gives the whole
#range to look at, to compare with the overlapping range
pred.sRNA4[subjectHits(sRNA.hits)]
pred.sRNA2[queryHits(sRNA.hits)]

#find overlapping ranges between ranges of predicted sRNAs, and hold the intersection of these regions.
orhits <- overlapsRanges(ranges(pred.sRNA1), ranges(pred.sRNA3), hits=sRNA.hits)
#check that the ranges produced are within the ranges of the query and subject
orhits[1:50,]
#create lists of widths of these ranges.
ranges_width_67035 <- width(ranges(pred.sRNA1[queryHits(sRNA.hits)]))
ranges_width_47863 <- width(pred.sRNA3[subjectHits(sRNA.hits)])
#create and export a dataframe of the widths of the ranges and the widths of the intersections of those
#ranges.
ranges_width_47863[1:50]
df <- data.frame(width(orhits), ranges_width_67035, ranges_width_47863)
write.csv(df, "./67035_47863_overlaps.csv")

#Experimental code
#I experimented with making line graphs of the overlapping ranges to initially visualise the
# data before exporting it
#a line graph is not the appropriate way to display the data but I wanted
#a quick way to look at the differences.
x <- seq(1,47,1)
df <- data.frame(x, width_overlaps, ranges_witdth_he1616, ranges_width_he67035)
ggplot(data=df, aes(x=x), xlab="width n") +
  geom_line(aes(y = width_overlaps, colour = "overlapping ranges")) +
  geom_line(aes(y = ranges_witdth_he1616, colour = "E-GEOD-47863")) +
  geom_line(aes(y= ranges_width_he67035, colour="E-GEOD-67035"))

#Merging GFF3 files to obtain a representative file. The following lines were written where BEDOPS was used to convert GFF3 files to BED files so that BEDtools intersect could be used to merge putative feature ranges in the 3 GFF3 files.
#BEDOPS was used to convert the GFF3 file to a BED file
#Gff2bed < gff1 > gff1.bed
#grep desired features from each GFF3 file, eg:
#Grep “putative_sRNA” gff1.bed > sRNAs_gff1.bed
#BEDtools intersect was used with parameter choices:
# -s to force strandedness
#bedtools intersect -s -a high_exp_47863.bed -b high_exp_67035.bed high_exp_1616.bed > combined_u.bed
#obtain putative features that are not overlapping with intersections:
#Bedtools intersect -s -v -a bed_file.bed -b bed_file2.bed bed_file3.bed > file1_exclusive_ncRNAs.bed

#cat the file containing intersections to the original annotation file
#cat the files with exclusively expressed features to this file, sort based on genomic ranges, and then order the columns to convert back to GFF3 file

#Sort -nk2 combined_features.bed > sorted_combined_features.bed

#awk '{ print $1 "\t" $7 "\t" $8 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $9 "\t" $10 }' sorted_combined_features.bed > combined_features.gff3


##### R was used to change base-0 format of BED start coordinates of ranges to base-1 format used in GFF3
gf <- import.gff3("E:/sorted_combined_features.gff3")
start(gf) <- start(gf) + 1
export.gff3(gf, "combined_mtb_5_10.gff3")

#################

#baerhunter was run here in the same way as above, now with the representative GFF3 file. This produced a final count matrix.

#the following script was written in R to normalise the data and visualise it. 
#Normalising the count matrix.
#‘’’{r counts}
#script to carry out DEseq2's VST normalisation on count matrix
#load libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
