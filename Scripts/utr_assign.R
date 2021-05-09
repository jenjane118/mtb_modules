## function for assigning categories to utrs

utr_assign <- function(gene_list, refseq_name, annot_file, tss_file) {
  # returns a dataframe with each predicted utr given an assigment
  # 5', 3' and related gene or 'NA'
  # gene_df: dataframe of a list of genes including predicted utrs from baerhunter
  # refseq_name: name of reference sequence that will be on granges, ex: "AL123456.3"
  # ref_seq: reference annotation file for making granges object
  
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  
  # read in gff for ref genome (use same one as for BH prediction) as genomic ranges object
  ref_granges <- import(annot_file)
  #isCircular(ref_granges)<-TRUE
  # filter for type='gene'
  ref_granges <- ref_granges[mcols(ref_granges)$type=="gene"]
  ref_genes   <- as.data.frame(ref_granges)
  # col'ID' has 'gene:Rv0001' and 'gene_id' has just ORF name 'Rv0001'
  gene_list <- mod_genes
  # parse list for UTRs predicted from baerhunter
  UTRs <- gene_list[grepl("UTR", gene_list)]
  if (length(UTRs)==0){
    return(NULL)
  }else{
    utrs <- str_split_fixed(UTRs, "UTR.", 2)[,2]
    UTR_df           <- data.frame(matrix(0, nrow=length(utrs), ncol=4),stringsAsFactors = F)
    colnames(UTR_df) <- c("pred_utr", "start", "stop", "strand")
    UTR_df$pred_utr  <- UTRs
  }
  # parse coordinates and strand from utr name
  for (i in 1:length(utrs)){
    strand <- substr(utrs[i], 1, 1)
    start  <- str_split_fixed(utrs[i], "_", 2)[,1]
    start  <- sub(".", "", start)
    stop   <- str_split_fixed(utrs[i], "_", 2)[,2]
    UTR_df$strand[i] <- strand
    UTR_df$start[i]  <- as.numeric(start)
    UTR_df$stop[i]   <- as.numeric(stop)
  }
  UTR_df$strand <- sub('p', '+', UTR_df$strand)
  UTR_df$strand <- sub('m', '-', UTR_df$strand)
  
  # create grange object of utr coordinates
  utr_Granges<-GRanges(seqnames = refseq_name,
                         ranges   = IRanges(UTR_df$start, 
                                            end = UTR_df$stop),
                         strand   = Rle(strand(UTR_df$strand))
                         )

  # searching whether utr co-locates with gene?
  # precede is gene that "is preceded" by the UTR (comes after UTR==downstream)
  pg<-GenomicRanges::precede(utr_Granges, ref_granges,ignore.strand=F)
  utr_precede_genes<-ref_genes$gene_id[pg]
  # precede doesn't pay attention to circular genome, use first gene in gff
  utr_precede_genes[is.na(utr_precede_genes)] <- ref_genes$gene_id[1]
  # follow is gene that 'is followed' by UTR (comes before UTR==upstream)
  ug<-GenomicRanges::follow(utr_Granges, ref_granges, ignore.strand=F)
  utr_follow_genes<-ref_genes$gene_id[ug]
  utr_follow_genes[is.na(utr_follow_genes)] <- ref_genes$gene_id[nrow(ref_genes)]
  ## nearest(utrs, refseq, select="all", ignore.strand = F) will give me nearest gene feature.
  # select=all means returns both nearest features in case of a tie, arbitrary chooses one.
  # see how many ties
  a<-nearest(utr_Granges, ref_granges, select="all", ignore.strand=F)
  b<-as.matrix(a)
  #show frequency of hits for each utr
  c<-as.data.frame(table(b[,1]))
  utr_tied<-c$Freq
  # see nearest gene (arbitrarily choose in ties)
  ng<-nearest(utr_Granges, ref_granges, select="arbitrary", ignore.strand=F)
  utr_nearest_genes<-ref_genes$gene_id[ng]

  #include preceding (downstream), following (upstream), and nearest genes in df 
  # (and indicate if there is a tie)
  UTR_df$downstream <-utr_precede_genes
  UTR_df$upstream   <-utr_follow_genes
  UTR_df$nearest    <-utr_nearest_genes
  UTR_df$tied       <-utr_tied
  
  # if nearest == precede, 3' UTR; if nearest ==follow, 5'UTR; if tied==look for tss
  # find TSS in utrs
  # subset utr_granges by strand 
  utr_fwd_Granges<-utr_Granges[strand(utr_Granges) == "+"]
  fwd_df<-as.data.frame(utr_fwd_Granges)
  utr_rev_Granges<-utr_Granges[strand(utr_Granges) =="-"]
  rev_df<-as.data.frame(utr_rev_Granges)
  # read in tss and make into granges obj
  tss<-read.delim(tss_file, sep=" ")
  tss_gr<-GRanges(seqnames = refseq_name,
                  ranges   = IRanges(tss$Genome.position, 
                                     end = tss$Genome.position),
                  strand   = Rle(strand(tss$Strand))
  )
  
  ## find overlap of each to tss within 10 bp of start 
  # forward strand tss
  fwd_tss_overlaps <-
    findOverlaps(
      utr_fwd_Granges,
      tss_gr,
      type = "start",
      maxgap = 10,
      select = "first",
      ignore.strand = FALSE
    )
  qh<-as.matrix(fwd_tss_overlaps)
  fwd_df$tss<-ifelse(is.na(qh), FALSE, tss[qh])
 
  # reverse strand tss
  rev_tss_overlaps<-
    findOverlaps(
      utr_rev_Granges,
      tss_gr,
      type = "end",
      maxgap = 10,
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
  UTR_df$TSS <- df_tss_o$tss

  # limit identified 5' UTRs to those with TSS (5' UTRs only) and 3' UTRs 
  # that are upstream of a gene (finish 20-40 nt before downstream gene)
  
  # only need start of downstream gene to see if 3' UTR (end of UTR > 20nts from downstream start)
  for (i in 1:nrow(UTR_df)){
    if (UTR_df$strand[i] == "+"){
      UTR_df$downstream_start[i]  <-
        ref_genes$start[match(UTR_df$downstream[i], ref_genes$gene_id)]
      UTR_df$dist_to_start[i] <-
        UTR_df$downstream_start[i] - UTR_df$stop[i]
    }else{
      # for minus strand, need end coordinate of downstream gene
      UTR_df$downstream_start[i]  <-
        ref_genes$end[match(UTR_df$downstream[i], ref_genes$gene_id)]
      UTR_df$dist_to_start[i] <-
        UTR_df$start[i] - UTR_df$downstream_start[i]}
  }
  # assign 5', 3' or NA to each utr
  for (i in 1:nrow(UTR_df)){
    # if tied (same dist between two genes, located between two genes)
    if (UTR_df$tied[i] == 2){
      # check for TSS
      if (UTR_df$TSS[i] != FALSE){
        UTR_df$utr[i] <- paste("5UTR", UTR_df$downstream[i], sep="_")
      }else{
        UTR_df$utr[i] <- "NA"}
    #not tied (closer to either upstream or downstream gene)
    }else{
      # check for TSS (if closer to downstream gene, needs tss to be 5')
      if (UTR_df$TSS[i] != FALSE){
        UTR_df$utr[i] <- paste("5UTR", UTR_df$downstream[i], sep="_")
      }else{
        if (UTR_df$nearest[i] == UTR_df$upstream[i]) {
          # check distance to downstream gene
          if (UTR_df$dist_to_start[i] > 40){
            UTR_df$utr[i] <- paste("3UTR", UTR_df$upstream[i], sep="_")
          }else{
            UTR_df$utr[i] <- "NA"}
        }else{
          UTR_df$utr[i] <- "NA"}
      }
    }
   }
  return(UTR_df)
}


