## function for verifying srnas with presence of tss


tss_srnas <- function(gene_list, refseq_name, tss_file) {
  # verify srna by searching for tss in sequence
  #gene_list <- mod_genes
  library(GenomicRanges)
  library(stringr)
  # create grange object of srna coordinates
  # parse list for sRNAs predicted from baerhunter
  og_SRNAs <- gene_list[grepl("sRNA", gene_list)]
  if (length(og_SRNAs)==0){
    return(NULL)
  }else{
    srnas    <- str_split_fixed(og_SRNAs, "sRNA:", 2)[,2]
    SRNA_df           <- data.frame(matrix(0, nrow=length(srnas), ncol=4),stringsAsFactors = F)
    colnames(SRNA_df) <- c("pred_srna", "start", "stop", "strand")
    SRNA_df$pred_srna  <- og_SRNAs
  }
  # parse coordinates and strand from snra name
  for (i in 1:length(srnas)){
    strand <- substr(srnas[i], 1, 1)
    start  <- str_split_fixed(srnas[i], "_", 2)[,1]
    start  <- sub(".", "", start)
    stop   <- str_split_fixed(srnas[i], "_", 2)[,2]
    SRNA_df$strand[i] <- strand
    SRNA_df$start[i]  <- as.numeric(start)
    SRNA_df$stop[i]   <- as.numeric(stop)
  }
  SRNA_df$strand <- sub('p', '+', SRNA_df$strand)
  SRNA_df$strand <- sub('m', '-', SRNA_df$strand)
  
  # create grange object of srna coordinates
  srna_Granges<-GRanges(seqnames = refseq_name,
                        ranges   = IRanges(SRNA_df$start, 
                                           end = SRNA_df$stop),
                        strand   = Rle(strand(SRNA_df$strand)))
  
  # find TSS in srnas
  # subset srna_granges by strand 
  srna_fwd_Granges<-srna_Granges[strand(srna_Granges) == "+"]
  fwd_df<-as.data.frame(srna_fwd_Granges)
  srna_rev_Granges<-srna_Granges[strand(srna_Granges) =="-"]
  rev_df<-as.data.frame(srna_rev_Granges)
  # read in tss and make into granges obj
  tss<-read.delim(tss_file, sep=" ")
  tss_gr<-GRanges(seqnames = refseq_name,
                  ranges   = IRanges(start = tss$Genome.position, 
                                     end = tss$Genome.position),
                  strand   = Rle(strand(tss$Strand))
  )
  ## find overlap of each to tss within 10 bp of start 
  # forward strand tss
  fwd_tss_overlaps <- findOverlaps(
    srna_fwd_Granges,
    tss_gr,
    type = "start",
    maxgap = 20,
    select = "first",
    ignore.strand = FALSE
  )
  qh<-as.matrix(fwd_tss_overlaps)
  fwd_df$tss<-ifelse(is.na(qh), FALSE, tss[qh])
  
  # reverse strand tss
  rev_tss_overlaps<- findOverlaps(
    srna_rev_Granges,
    tss_gr,
    type = "end",
    maxgap = 20,
    select = "first",
    ignore.strand = FALSE
  )
  qh<-as.matrix(rev_tss_overlaps)
  rev_df$tss<-ifelse(is.na(qh), FALSE, tss[qh])
  #put back together
  df_tss<-rbind(fwd_df, rev_df)
  # order by start
  df_tss_o<-df_tss[order(df_tss$start),]
  # add to dataframe
  SRNA_df$TSS <- df_tss_o$tss
  SRNA_df$pred_srna <- ifelse(SRNA_df$TSS!=FALSE, paste(SRNA_df$pred_srna, "*", sep=""), SRNA_df$pred_srna)
  return(SRNA_df)
}

