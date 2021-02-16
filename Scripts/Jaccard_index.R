## find way to determine length of overlap (intersection) and total bp of sRNAs 
## in order to calculate Jaccard index for each sRNA 'hit'

library(rtracklayer)
library(GenomicRanges)
#annotation_file<-"output_BH_25_11/PRJNA327080_15.gff3"
#known_RNA_file<-"ncRNA_verified.txt"
#minoverlap<-1
# We will repeat this process several times so worth packaging in a function
jaccard_known_with_predicted <- function(annotation_file, known_RNA_file, minoverlap=1L) {
  # import the homemade ncRNA file using rtracklayer's import function
  annot <- import.gff3(annotation_file)
  
  # create subsets of the putative sRNAs and UTRs first
  pred.sRNA <- subset(annot, type == "putative_sRNA")
  pred.sRNA
  pred.UTRs.plus <- subset(annot, (type == "putative_UTR") & (strand == "+"))
  num_pred.UTRs.plus <- length(ranges(subset(annot, (type == "putative_UTR") & (strand == "+")) ))
  num_pred.UTRs.plus
  pred.UTRs.minus <- subset(annot, (type == "putative_UTR") & (strand == "-"))
  num_pred.UTRs.minus <- length(ranges(subset(annot, (type == "putative_UTR") & (strand == "-")) ))
  num_pred.UTRs.minus
  
  # read in the data from ncRNAs_and_TSS.txt file made with Gerrick and TSS annotations
  known.ori <- read.table(file=known_RNA_file, header=FALSE,   
                          sep = " ", skip=1) 
  # select necessary columns
  known <- known.ori[,1:7]
  
  colnames(known)<-c("Chromosome", "sRNA_ID", "start", "end", "strand", "distance", "TSS")
  known
  # read in data from "comb_cortesTSS_srna.txt" for TSS sites
  known_TSS_file = "comb_cortesTSS_srna.txt"
  tss <- read.table(file=known_TSS_file, header=T, sep = " ")
  
  known.sRNA <- GRanges(
    seqnames ="AL123456.3", 
    ranges   =IRanges(start=known$start, end=known$end), 
    strand   =Rle(strand(known$strand)),
    sRNA_ID  =known$sRNA_ID
  )
  
  # create 5'UTR list from TSS sites?
  utr.tss <- GRanges(
    seqnames = "AL123456.3",
    ranges   = IRanges(start=tss$Genome.position, end=tss$Genome.position),
    strand   = Rle(strand(tss$Strand))
  )
  
  # In the original GRanges for known.all, the starts and ends are strand specific, 
  # ie. forward start is < end and for reverse, start > end. In my list, use strand 
  # and all starts are < ends). 
  sRNA.hits <- findOverlaps(known.sRNA, pred.sRNA, type="any", minoverlap=minoverlap, ignore.strand=FALSE) 

  #UTR.hits <- findOverlaps(utr.tss, pred.UTR, type="any", minoverlap=minoverlap, ignore.strand=FALSE) 
  # filter for position of TSS within UTR (like in first 20nt?) 
  # this needs to be strand specific--for + strands need to use 'start' 
  # for - strand, need to use 'end'
  plus_UTR.hits <- findOverlaps(utr.tss, pred.UTRs.plus, type="start", maxgap=20, 
                                minoverlap=1, ignore.strand=F)
  minus_UTR.hits <- findOverlaps(utr.tss, pred.UTRs.minus, type="end", maxgap=20, 
                                 minoverlap=1, ignore.strand=F)
  
  
  
  #how is this different than findOverlaps? finds same hits
  subov<-subsetByOverlaps(known.sRNA, pred.sRNA, type="any", minoverlap = 1, ignore.strand=FALSE)
  subov         ##known sRNAs (coordinates) that overlap predicted sRNAs
  
  
  
  s<-subjectHits(sRNA.hits)
  pred.sRNA[s]   ## pred sRNAs that overlap at all with known sRNAs (some are overlapping more than one) #195
  unique(pred.sRNA[s])     ## pred sRNAs that overlap at all with known sRNAs (no repeats) #147
  s.df<-as.data.frame(pred.sRNA[s])
  
  # # these subject hits in test set are way too long for sRNAs. Hitting too many known sRNAs 
  # # just because too long (too much noise with default params)
  q<- queryHits(sRNA.hits) ## known that overlap any pred sRNAS, some overlap more than one #195
  unique(queryHits(sRNA.hits))  ## known sRNAs that overlap any predicted sRNAs (no repeats) #191
  
  q.df<-as.data.frame(known.sRNA[q])
  q.df<-data.frame(lapply(q.df, as.character), stringsAsFactors = F)

  # need to find parallel intersection/union using pintersect and punion for sRNA 
  # hits (or all predicted sRNAs?)
  pi<-as.data.frame(pintersect(known.sRNA[q], pred.sRNA[s], ignore.strand=F))
  pu<-as.data.frame(punion(known.sRNA[q], pred.sRNA[s], ignore.strand=F))
  
  #calculate jaccard index which is intersection of ranges/union of ranges
  #these will include non-unique subject hits--predicted sRNAs that overlap with
  # more than one known sRNA (and may have more entries than total number of unique
  # sRNA hits. Jaccard index will vary depending on which sRNA it is matching with. 
  JIndex<-round(pi$width/pu$width, 3)
  
  #new dataframe for hits and jaccard index
  Jindex_df<-NULL
  Jindex_df<- as.data.frame(cbind(q.df$sRNA_ID, q.df$start, q.df$end, q.df$strand))
  # bind index column separately as it's numeric
  Jindex_df<-cbind(Jindex_df, JIndex)
  colnames(Jindex_df)<-c("sRNA_ID", "start", "end", "strand", "index")
  Jindex_df
  
  # for each 
  #how many predicted hits  are confirmed?
  num_UTR_confirmed <- length(unique(subjectHits(plus_UTR.hits)))+ 
    length(unique(subjectHits(minus_UTR.hits)))
  num_sRNA_confirmed <- length(unique(subjectHits(sRNA.hits)))
  
  num_pred.sRNA <- length(ranges(pred.sRNA))
  num_pred.UTRs <- length(ranges(subset(annot, (type=="putative_UTR"))))
  
  res <- c(num_pred.sRNA, num_sRNA_confirmed, num_pred.UTRs, num_UTR_confirmed)
  names(res) <- c('num_pred_sRNA', 'num_conf_sRNA', 'num_pred_UTR', 
                  'num_conf_UTR')
  ret_list<-list(res, Jindex_df)
  return(ret_list)
}

# run on server after fixed bug and put in automatic low/high threshold:
res_25.11_15 <- jaccard_known_with_predicted(
  annotation_file="output_BH_25_11/PRJNA327080_15.gff3", 
  minoverlap=5L, 
  known_RNA_file ="ncRNA_verified.txt"
)
res_25.11_15[1]
#num_pred_sRNA num_conf_sRNA  num_pred_UTR  num_conf_UTR 
#3152           147          8086           387 

#ridiculous number of pred_UTRs!

PRJNA327080_newBH<-jaccard_known_with_predicted(
  annotation_file = "PRJNA327080_newBH.gff3",
  minoverlap = 5L,
  known_RNA_file = "ncRNA_verified.txt"
)
PRJNA327080_newBH
## num_pred_sRNA num_conf_sRNA  num_pred_UTR  num_conf_UTR 
## 959             8          2005            71 

# very low number of overlap with confirmed sRNA. HOw does it decrease with the 
# full dataset? also, jaccard index for those overlapping is really low.



res_prjna327080_new <- jaccard_known_with_predicted(
  annotation_file = "bh_gffs/PRJNA327080_newBH.gff3",
  minoverlap = 1L,
  known_RNA_file = "ncRNA_verified.txt"
)
# extract summary of hits
res_prjna327080_new[1]
# extract dataframe with hits and jaccard index
j<-as.data.frame(res_prjna327080_new[2])
# extract vector of jaccard index
ji<-as.numeric(j[,5])
mean(ji)


# look at PRJNA327080_15 with 80 and 160?
res_80_160_PRJNA327080_15<-jaccard_known_with_predicted(
  annotation_file = "output_BH_16_11/PRJNA327080_15.gff3",
  minoverlap = 5L,
  known_RNA_file = "ncRNA_verified.txt"
)
res_80_160_PRJNA327080_15
## num_pred_sRNA num_conf_sRNA  num_pred_UTR  num_conf_UTR 
## 1158             9          1947            75 


## use two sample test case of edited bh with moving parameters
test_PRJNA327080_2<-jaccard_known_with_predicted(
  annotation_file = "BH_test_4.gff3",
  minoverlap = 5L,
  known_RNA_file = "ncRNA_verified.txt"
)
test_PRJNA327080_2
## num_pred_sRNA num_conf_sRNA  num_pred_UTR  num_conf_UTR 
## 700            12          1348            49 

PRJNA327080_newBH<-jaccard_known_with_predicted(
  annotation_file = "PRJNA327080_newBH.gff3",
  minoverlap = 5L,
  known_RNA_file = "ncRNA_verified.txt"
  )
PRJNA327080_newBH
## num_pred_sRNA num_conf_sRNA  num_pred_UTR  num_conf_UTR 
## 959             8          2005            71 

# very low number of overlap with confirmed sRNA. HOw does it decrease with the 
# full dataset? also, jaccard index for those overlapping is really low.

prjna327080_07_12<-jaccard_known_with_predicted(
  annotation_file = "output_BH_07_12/PRJNA327080_15.gff3",
  minoverlap = 5L,
  known_RNA_file = "ncRNA_verified.txt"
)
prjna327080_07_12
