#combining .gff3 files with predicted sRNA/UTR annotations to get one master list

library(GenomicRanges)
library(rtracklayer)
library(baerhunter)

# 1. Subset elements by ncRNAs (sRNA and UTR)
# 2. make genomic ranges object for each dataset
# 3. combine all ncRNAs together in one grange object
# 4. reduce to simplified set of ncRNA coordinates
# 5. map back to get names of ncRNAs? if combine sRNAs and UTRs, to identify which? 

# import each .gff and subset by ncRNAs
ref<-import.gff3("ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3")

g1<-import.gff3("bh_gffs/output_BH_28_02/output_BH_28_02_PRJNA327080_15.gff3")
g2<-import.gff3("bh_gffs/output_BH_28_02/output_BH_28_02_PRJNA390669_12.gff3")
g3<-import.gff3("bh_gffs/output_BH_28_02/output_BH_28_02_PRJEB65014.gff3")
g4<-import.gff3("bh_gffs/output_BH_01_03_PRJNA278760_22.gff3")

# add all sRNAs and utrs to common list (do these separately)
srnas_gr1<-g1[elementMetadata(g1)[,"type"] %in% "putative_sRNA"]
# filter for sRNAs/UTRs > 1000 nt
srnas_gr1<-srnas_gr1[which(width(srnas_gr1)<=1000),]
utrs_gr1<-g1[elementMetadata(g1)[,"type"] %in% "putative_UTR"]
utrs_gr1<-utrs_gr1[which(width(utrs_gr1)<=1000),]
srnas_gr2<-g2[elementMetadata(g2)[,"type"] %in% "putative_sRNA"]
srnas_gr2<-srnas_gr2[which(width(srnas_gr2)<=1000),]
utrs_gr2<-g2[elementMetadata(g2)[,"type"] %in% "putative_UTR"]
utrs_gr2<-utrs_gr2[which(width(utrs_gr2)<=1000),]
srnas_gr3<-g3[elementMetadata(g3)[,"type"] %in% "putative_sRNA"]
srnas_gr3<-srnas_gr3[which(width(srnas_gr3)<=1000),]
utrs_gr3<-g3[elementMetadata(g3)[,"type"] %in% "putative_UTR"]
utrs_gr3<-utrs_gr3[which(width(utrs_gr3)<=1000),]
srnas_gr4<-g4[elementMetadata(g4)[,"type"] %in% "putative_sRNA"]
srnas_gr4<-srnas_gr4[which(width(srnas_gr4)<=1000),]
utrs_gr4<-g4[elementMetadata(g4)[,"type"] %in% "putative_UTR"]
utrs_gr4<-utrs_gr4[which(width(utrs_gr4)<=1000),]


total_srnas_gr<-c(srnas_gr1, srnas_gr2, srnas_gr3, srnas_gr4)
length(total_srnas_gr)
#1552 v 785

# reduce to align ranges and merge overlapping ranges, those with specified gap between 
# will not be merged, revmap returns list of overlapping ranges
red_srnas_gr<-reduce(total_srnas_gr, min.gapwidth=1L)
## with 1L min.gapwidth==2879 ranges, with 3L==2874 ranges, with 5L==2866 ranges


# do the same with utrs
total_utrs_gr<-c(utrs_gr1, utrs_gr2, utrs_gr3, utrs_gr4)
red_utrs_gr<-reduce(total_utrs_gr, min.gapwidth=1L)


## how to get names back? 
#extract ID info as a dataframe:
# total_srna_ids<-mcols(total_srnas_gr)$ID
# total_srna_ids[1:10]
# #get list of all overlapping ranges for each reduced ncrna range, can see which 
# # srnas are in how many datasets
# total_ncrnas_gr[5221]$ID
# total_ncrnas_gr[5221]$upstream_feature
# ## use findoverlaps 'within' ?
# ncrnas_gr<-NULL
# ov_gr1<-findOverlaps(ncrnas_gr1, red_ncrnas_gr, type="any")
# ncrnas_gr<-ncrnas_gr1[queryHits(ov_gr1)]
# length(ncrnas_gr) 
# length(ncrnas_gr1)
# # maybe easier to use union?
# union_ncrnas<-NULL
# union_ncrnas<-union(ncrnas_gr1, ncrnas_gr2)
# union_ncrnas<-union(union_ncrnas, ncrnas_gr3)
# union_ncrnas<-union(union_ncrnas, ncrnas_gr4)
# union_ncrnas
# #2879 ranges, also no metadata columns!
# make into data frame and add back at least id, like in bh? still will need a naming
# step, but maybe don't need to be named unless verified?

length(red_srnas_gr)
length(red_utrs_gr)

# subset by strand
pos_srnas_gr<-subset(red_srnas_gr, strand=="+")
neg_srnas_gr<-subset(red_srnas_gr, strand=="-")
# use code from srna_calc and utr_calc to re-name the new elements 
# (needed to add 'as.integer' to remove whitespace before coordinates)
names(pos_srnas_gr) <- apply(as.data.frame(pos_srnas_gr), 1, function(x) paste("ID=putative_sRNA:p", as.integer(x[2]), "_", as.integer(x[3]), ";", sep = ""))
names(neg_srnas_gr) <- apply(as.data.frame(neg_srnas_gr), 1, function(x) paste("ID=putative_sRNA:m", as.integer(x[2]), "_", as.integer(x[3]), ";", sep = ""))

pos_utrs_gr<-subset(red_utrs_gr, strand=="+")
neg_utrs_gr<-subset(red_utrs_gr, strand=="-")
names(pos_utrs_gr) <- apply(as.data.frame(pos_utrs_gr),1, function(x) paste("ID=putative_UTR:p", as.integer(x[2]), "_", as.integer(x[3]), ";", sep = ""))
names(neg_utrs_gr) <- apply(as.data.frame(neg_utrs_gr),1, function(x) paste("ID=putative_UTR:m", as.integer(x[2]), "_", as.integer(x[3]), ";", sep = ""))

# add to major features with strand feature editor

# make 'major features file' for each strand (use major_features func from BH)
pos_features<-major_features("ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3", ".", "+", "ncRNA")
neg_features<-major_features("ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3", ".", "-", "ncRNA")

# strand_feature_editor(target_strand, sRNA_IRanges, UTR_IRanges, major_strand_features) 
pos_strand<-strand_feature_editor(target_strand = "+", pos_srnas_gr, pos_utrs_gr, pos_features)
neg_strand<-strand_feature_editor(target_strand = "-", neg_srnas_gr, neg_utrs_gr, neg_features)

# use last step of ffe to create new file

## Creating the final annotation dataframe by combining both strand dataframe and adding missing information like child features from the original GFF3 file.

gff <- read.delim("ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3", header = FALSE, comment.char = "#")
annotation_dataframe <- rbind(gff, pos_strand, neg_strand)
## Remove all the repeating information.
annotation_dataframe <- unique(annotation_dataframe)
## Order the dataframe by feature start coordinates.
annotation_dataframe <- annotation_dataframe[order(annotation_dataframe[,4]),]

## Restore the original header.
f <- readLines("ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3")
header <- c()
i <- 1
while (grepl("#",f[i])==TRUE) {
  f_line <- f[i]
  header <- c(header,f_line)
  i <- i+1
}
# add a line to indicate the origin of the file (single # commas should be ignored by programs)
header <- c(header, "# produced by baerhunter")

## Create the final GFF3 file.
output_file<-"combined_gffs_01_03.gff3"
write.table(header, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(annotation_dataframe, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

#how many in combined list are confirmed?
# function from verifying_BHgffs.R
compare_known_with_predicted <- function(annotation_file, known_RNA_file, minoverlap=minoverlap) {
  # import the homemade ncRNA file using rtracklayer's import function
  library(rtracklayer)
  # creates a GRanges object
  annot <- import.gff3(annotation_file)
  head(annot)
  
  # create subsets of the putative sRNAs and UTRs first
  pred.sRNA <- subset(annot, type == "putative_sRNA")
  head(pred.sRNA)
  #min(width(pred.sRNA))
  pred.UTR <- subset(annot, type == "putative_UTR")
  #num.pred.sRNA.plus <- length(ranges(subset(annot, (type == "putative_sRNA") & (strand == "+")) ) ) 
  #num.pred.sRNA.minus <- length(ranges(subset(annot, (type == "putative_sRNA") & (strand == "-")) ) ) 
  #num.pred.UTRs.plus <- length(ranges(subset(annot, (type == "putative_UTR") & (strand == "+")) ) ) 
  #num.pred.UTRs.minus <- length(ranges(subset(annot, (type == "putative_UTR") & (strand == "-")) ) ) 
  
  # read in the data from ncRNAs_and_TSS.txt file made with Gerrick and TSS annotations
  known.ori <- read.table(file=known_RNA_file, header=FALSE,   
                          sep = " ", skip=1) 
  # select necessary columns
  known <- known.ori[,1:7]
  colnames(known)<-c("Chromosome", "sRNA_ID", "start", "end", "strand", "distance", "TSS")
  head(known)
  
  #seqnames need to match what's in predicted file (AL123456.3)
  # used different annotation file first time ran BH (NC)
  known.all <- GRanges(
    seqnames ="AL123456.3", 
    #seqnames = "NC_000962.3",
    ranges   =IRanges(start=known$start, end=known$end), 
    strand   =Rle(strand(known$strand)),
    sRNA_ID  =known$sRNA_ID
  )
  head(known.all)
  # should we filter this for sRNAs with TSS sites?
  
  
  # need to make a new UTR list from TSS site list
  # read in data from "comb_cortesTSS_srna.txt" for TSS sites
  known_TSS_file = "comb_cortesTSS_srna.txt"
  tss <- read.table(file=known_TSS_file, header=T, sep = " ")
  #head(tss)
  # create 5'UTR list from TSS sites?
  utr.tss <- GRanges(
    seqnames = "AL123456.3",
    #seqnames = "NC_000962.3",
    ranges   = IRanges(start=tss$Genome.position, end=tss$Genome.position),
    strand   = Rle(strand(tss$Strand))
  )
  
  # we ignore the strand and count features on both strands (ignore.strand=T).
  # In the original GRanges for known.all, the starts and ends are strand specific, 
  # ie. forward start is < end and for reverse, start > end. In my list, use strand 
  # and all starts are < ends).
  
  sRNA.hits <- findOverlaps(known.all, pred.sRNA, type="any", minoverlap=minoverlap, ignore.strand=F) 
  #UTR.hits <- findOverlaps(utr.tss, pred.UTR, type="any", minoverlap=1, ignore.strand=FALSE ) 
  #filter for tss within first 20nts of start (should we be doing this with sRNAs, too?)
  UTR.hits <- findOverlaps(utr.tss, pred.UTR, type="start", maxgap=20, minoverlap=1, ignore.strand=F)
  
  #how many predicted hits  are confirmed?
  num_UTR_confirmed <- length(unique(subjectHits(UTR.hits)))
  num_sRNA_confirmed <- length(unique(subjectHits(sRNA.hits)))
  num_all_confirmed <- length(unique(c(subjectHits(UTR.hits), subjectHits(sRNA.hits))))
  
  num_UTR_confirmed
  num_sRNA_confirmed
  num_all_confirmed
  
  num_pred_UTR <- length(ranges(pred.UTR))
  num_pred_sRNA <- length(ranges(pred.sRNA))
  num_pred_UTR
  num_pred_sRNA 
  res <- c(num_pred_sRNA, num_sRNA_confirmed,
           num_pred_UTR, num_UTR_confirmed)
  names(res) <- c('num_pred_sRNA', 'num_sRNA_confirmed',
                  'num_pred_UTR', 'num_UTR_confirmed')
  return(res)
}

res_combined_01_03<-compare_known_with_predicted(annotation_file = "combined_gffs_01_03.gff3", 
                                                 minoverlap = 5L,
                                                 known_RNA_file = "ncRNA_verified.txt"
)
res_combined_01_03
#num_pred_sRNA num_sRNA_confirmed       num_pred_UTR  num_UTR_confirmed 
#  510                 17               1958                136 

annot <- import.gff3("combined_gffs_01_03.gff3")
pred.sRNA <- subset(annot, type == "putative_sRNA")
pred.UTR <- subset(annot, type == "putative_UTR")
w<-width(pred.sRNA)
mean(w)
min(w)
max(w)
long_srnas<-as.data.frame(pred.sRNA[which(width(pred.sRNA)>1000),2:3])
long_srnas
#> long_srnas
#seqnames   start     end width strand          type score
#1 AL123456.3  767366  768899  1534      + putative_sRNA    NA
#2 AL123456.3 1271316 1272321  1006      - putative_sRNA    NA
#3 AL123456.3 2058035 2059284  1250      - putative_sRNA    NA
> 
#write.table(long_srnas, "long_srnas.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

# 56 are longer than 1000 nts (with earlier parameters)
#this is a really long one:
#AL123456.3	.	putative_sRNA	1065730	1066752	.	-	.	ID=putative_sRNA:m1065730_1066752;upstream_feature=m1067194_1073291;downstream_feature=m1064963_1065474

# we can filter out longest srnas?

# I wonder if there should be a bigger gap so it merges some utrs that have only
# 1 nt between them? Changed to 5L for UTRs and only reduced predicted by 2 and 
# confirmed by 1.

# need to make a naming program to implement nomenclature from lamichlane
