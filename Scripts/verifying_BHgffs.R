## Overlap of predicted and known ncRNAs
#Here, we check whether the features predicted by baerhunter have any 
#overlap with confirmed small RNAs from list made for
#Laura with verified TSS (Gerrick,2018; Cortes, 2013; Shell, 2015)
# Gerrick, E. R. et al. Small RNA profiling in mycobacterium tuberculosis 
#identifies mrsi as necessary for an anticipatory iron sparing response. 
#Proc. Natl. Acad. Sci. U. S. A. 115, 6464–6469 (2018).

library(rtracklayer)


# the first time I used BH, I used refseq file NC_000962.3, so need
# to change this when making predicted grange objects


# We will repeat this process several times so worth packaging in a function
compare_known_with_predicted <- function(annotation_file, known_RNA_file, minoverlap=minoverlap) {
  # import the homemade ncRNA file using rtracklayer's import function
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
  #res
  return(res)
}


# run function with BH run with default parameters (low=5,hi=10,min_srna=40,min_utr=50)
res_6.11_1 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA327080_15.gff3", 
  minoverlap=1L, 
  known_RNA_file ="ncRNA_verified.txt"
  )
res_6.11_1

res_6.11_5 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA327080_15.gff3", 
  minoverlap=5L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_6.11_5
#num_pred_sRNA num_sRNA_confirmed       num_pred_UTR  num_UTR_confirmed 
#  1674                  7               3066                140 


res_6.11_10 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA327080_15.gff3", 
  minoverlap=10L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_6.11_20 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA327080_15.gff3", 
  minoverlap=20L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_6.11_50 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA327080_15.gff3", 
  minoverlap=50L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_6.11_5
# Plot the results as barplots
install.packages("ggplot2")
library(ggplot2)
v.y  <- c(res_6.11_5[1], res_6.11_5[2], res_6.11_10[2], res_6.11_20[2])
v.labels <- c("predicted", "verified.5", "verified.10", "verified.20")
v.type <- c(rep("sRNA", 4))
df <- data.frame(cbind(v.labels, v.type, v.y))
colnames(df) <- c("labels", "type", "num_predictions")
df
#put the v.types factor in the right order in the plot
df$labels <- factor(df$labels, levels=c('predicted', 'verified.5', 'verified.10', 'verified.20'))
# df$v.y is factor and needs to be numeric
p <- ggplot(data=df, aes(x=type, y=as.numeric(as.character(num_predictions)), fill=labels)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_y_continuous(limits=c(0, 2000)) +
  #xlab("sRNA") +
  ylab("Number of predictions") +
  labs(fill='') +
  scale_fill_grey() +
  labs(title="Baerhunter (6-11) predictions", caption="verified.x = predictions with min x nt overlap with ncRNAs from verified list") +
  theme(plot.caption = element_text(hjust=0.5))
p

head(df)


res_PRJNA278760_5 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA278760_22.gff3", 
  minoverlap=5L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_PRJNA278760_5
#  num_pred_sRNA num_sRNA_confirmed       num_pred_UTR  num_UTR_confirmed 
#   2851                  8               3032                197 
res_PRJNA278760_10 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA278760_22.gff3", 
  minoverlap=10L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_PRJNA278760_20 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA278760_22.gff3", 
  minoverlap=20L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_PRJNA278760_50 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA278760_22.gff3", 
  minoverlap=50L, 
  known_RNA_file ="ncRNA_verified.txt"
)


v.y  <- c(res_PRJNA278760_5[1], res_PRJNA278760_5[2], res_PRJNA278760_10[2], res_PRJNA278760_20[2])
v.labels <- c("predicted", "verified.5", "verified.10", "verified.20")
v.type <- c(rep("sRNA", 4))
PRJNA278760.df <- data.frame(cbind(v.labels, v.type, v.y))
colnames(PRJNA278760.df) <- c("labels", "type", "num_predictions")
PRJNA278760.df


## pretty much all finding same low numbers of sRNAs, despite predicting loads
## more with larger datasets. But my list isn't that long, 205 ncRNAs.

# how many sRNAs are predicted by really high depth dataset, PRJNA390669?
res_PRJNA390669_5 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA390669.gff3", 
  minoverlap=5L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_PRJNA390669_5
#num_pred_sRNA num_sRNA_confirmed       num_pred_UTR  num_UTR_confirmed 
#232                  3               2699                134 


res_PRJNA390669_10 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA390669.gff3", 
  minoverlap=10L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_PRJNA390669_20 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA390669.gff3", 
  minoverlap=20L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_PRJNA390669_50 <- compare_known_with_predicted(
  annotation_file="BH_gffs_default/PRJNA390669.gff3", 
  minoverlap=50L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_PRJNA390669_5
## this dataset predicted way fewer sRNAs at these parameters
v.y  <- c(res_PRJNA390669_5[1], res_PRJNA390669_5[2], res_PRJNA390669_10[2], res_PRJNA390669_20[2])
v.labels <- c("predicted", "verified.5", "verified.10", "verified.20")
v.type <- c(rep("sRNA", 4))
PRJNA390669.df <- data.frame(cbind(v.labels, v.type, v.y))
colnames(PRJNA390669.df) <- c("labels", "type", "num_predictions")
PRJNA390669.df

# run small dataset, PRJEB65014_3 run at 50 and 100 (low and hi)
res_50.100_PRJEB65014_3 <- compare_known_with_predicted(
  annotation_file="output_BH_13_11/PRJEB56014_3.gff3", 
  minoverlap=1L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_50.100_PRJEB65014_3
# 112 predicted, 6 confirmed sRNAs


# high depth dataset with 300/450 parameters
res_300_450_PRJNA390669_12 <- compare_known_with_predicted(
  annotation_file = "PRJNA390669_12_RC.gff3",
  minoverlap = 1L,
  known_RNA_file = "ncRNA_verified.txt"
)
res_300_450_PRJNA390669_12
# sRNAs 676 predicted, 10 confirmed
# UTRs  1764 predicted, 58 confirmed

## what was original results using default parameters for this dataset?
res_default_PRJNA390669_12 <-compare_known_with_predicted(
  annotation_file = "BH_gffs_default/PRJNA390669.gff3",
  minoverlap = 1L,
  known_RNA_file = "ncRNA_verified.txt"
)
res_default_PRJNA390669_12
#> res_default_PRJNA390669_12
#num_pred_sRNA num_sRNA_confirmed       num_pred_UTR  num_UTR_confirmed 
#232                  4               2699                134 
#so clearly much better with higher parameters

# ran one bam to look at on artemis: PRJNA327080_15: SRR3725585_sorted.bam
res_27_11_med<-compare_known_with_predicted(
  annotation_file = "output_BH_27_11/PRJNA327080_585.gff3",
  minoverlap = 5L,
  known_RNA_file = "ncRNA_verified.txt"
)
res_27_11_med

#test on single bam from high density
res_test.3.12<-compare_known_with_predicted(annotation_file = "test2t.3.12.gff3",
                                            minoverlap = 5L,
                                            known_RNA_file = "ncRNA_verified.txt"
                                            )
res_test.3.12

res_07_12_327080<-compare_known_with_predicted(annotation_file = "output_BH_07_12/PRJNA327080_15.gff3",
                                               minoverlap = 5L,
                                               known_RNA_file = "ncRNA_verified.txt"
                                               )
res_07_12_327080
#num_pred_sRNA num_sRNA_confirmed       num_pred_UTR  num_UTR_confirmed 
#574                 10               1623                 52 
res_07_12_390669<-compare_known_with_predicted(annotation_file = "output_BH_07_12/PRJNA390669_12.gff3",
                                               minoverlap = 5L,
                                               known_RNA_file = "ncRNA_verified.txt"
                                               )
res_07_12_390669
#num_pred_sRNA num_sRNA_confirmed       num_pred_UTR  num_UTR_confirmed 
#495                  9               1843                 61 

#how many in combined list are confirmed?
res_combined_08_12<-compare_known_with_predicted(annotation_file = "combined_gffs_08_12", 
                                                 minoverlap = 5L,
                                                 known_RNA_file = "ncRNA_verified.txt"
                                                 )
res_combined_08_12
#num_pred_sRNA num_sRNA_confirmed       num_pred_UTR  num_UTR_confirmed 
# 972                 25               2170                150 
