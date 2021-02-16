## compare presrat mtb non-genic predictions with 'stable RNAs' from mycobrowser

p_table<-read.table("Data/non-genic_presrat_mtb.txt", sep="\t", header=T)
head(p_table)
p_table$Strand<-gsub('ve', '',p_table$Strand)
head(p_table)
summary(p_table)

#divide by strands
f_p_table<-p_table[which(p_table$Strand == '+'),]
r_p_table<-p_table[which(p_table$Strand == '-'),]
# create genomic ranges object
library(GenomicRanges)
f_presrat_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(f_p_table$Start, end = f_p_table$End),
  strand   = Rle(strand(f_p_table$Strand)),
  rank     = f_p_table$Rank
)
r_presrat_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(r_p_table$End, end = r_p_table$Start),
  strand   = Rle(strand(r_p_table$Strand)),
  rank     = r_p_table$Rank
)
# combine granges objects
presrat_gr<-c(f_presrat_gr, r_presrat_gr)



#check overlap with stable_ncrna_myco.txt list
library(rtracklayer)
myco_table<-readGFF("Data/stable_rna_myco.txt")
nrow(myco_table)
myco_gr<-GRanges(
  seqnames = 'H37Rv',
  ranges   = IRanges(myco_table$start, end=myco_table$end),
  strand   = Rle(strand(myco_table$strand)),
  name     = myco_table$Name
)
ov_pres<-findOverlaps(myco_gr, presrat_gr, minoverlap = 10, ignore.strand=FALSE)
length(ov_pres)

subjectHits(ov_pres)
pres_hits<-as.data.frame(presrat_gr[subjectHits(ov_pres)])
pres_hits
myco_hits<-as.data.frame(myco_gr[queryHits(ov_pres)])
myco_hits

# overlaps with 7 from list
