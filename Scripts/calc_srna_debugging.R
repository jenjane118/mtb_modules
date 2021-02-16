## figure out how this function works and how to fix


msf<-m_strand_features[1:10,]
upr<-union_peak_ranges[1:10]
ts<-"+"


#sRNA_calc <- function(major_strand_features, target_strand, union_peak_ranges) {
  ## This function predicts sRNAs.
  ## Convert strand feature coordinates into IRanges.
strand_IRange <- IRanges(start = msf[,4], end = msf[,5])
strand_IRange
# what in upr overlaps with msf
s<-subsetByOverlaps(upr, strand_IRange, maxgap=1L)
s
# what of upr is in subset? logical vector of entire upr (TTFTFFFFFF)
n<-upr %in% s
n
# uprs in subset== same as subset (s) why do this?
upr[n]

# this just reverses the vector(T to F, etc), so uprs that are NOT in subset?
!n           #FFTFTTTTT
srnas<-upr[!n,]
srnas
## alternatively, upr that is NOT in subset (original code)
m<-!upr %in% s  #FFTFTTTTTT  (like !n)
m
srnas2<-upr[!m,]
srnas2

## so original code is exactly the same, with or without brackets. just reverses
## logical vector

  ## Select only the ranges that do not overlap the annotated features. Also, disregard the ranges that finish/start 1 position before the genomic feature, because they should be considered as UTRs.
IGR_sRNAs <- upr[! (upr %in% subsetByOverlaps(upr, strand_IRange, maxgap = 1L)),]
# original line of code, does this change results?
#IGR_sRNAs <- union_peak_ranges[! union_peak_ranges %in% subsetByOverlaps(union_peak_ranges, strand_IRange, maxgap = 1L),]  
  ## Construct the IDs for the new sRNAs to be added into the attribute colmn of the annotation.
if (target_strand=="+") {
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste("ID=putative_sRNA:p", x[1], "_", x[2], ";", sep = ''))
  } else if (target_strand== "-") {
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste("ID=putative_sRNA:m", x[1], "_", x[2], ";", sep = ''))
  } else {
    return("Select strand")
}
  

# what about utr function? could problem be there?


s<-c(2,10,12,25)
e<-c(7,15,20,29)
genes_ir<-IRanges(start=s, end=e)
genes_ir

p1<-c(1,8,22,28,45)
p2<-c(3,12,24,40,50)
peaks_ir<-IRanges(start=p1, end=p2)
peaks_ir

subs_over<-IRanges::subsetByOverlaps(peaks_ir, genes_ir, maxgap = 1L)

subs_over
peaks_ir %in% subs_over
peaks_ir[!(peaks_ir %in% subs_over)]

peaks_in_subs<-findOverlaps(peaks_ir, subs_over, type="equal")
peaks_in_subs

peaks_not_in_subs


#UTR_calc <- function(major_strand_features, target_strand, union_peak_ranges, min_UTR_length) {
  ## This function predicts UTRs.
  ## Convert strand feature coordinates into IRanges.
#short_IRange <- IRanges(start = msf[,4], end = msf[,5])
#short_IRange
  ## Find the peak union ranges that overlap with genomic features. Also, include the ranges that do not overlap the features but start/finish 1 position away from it.
o_features <- subsetByOverlaps(peaks_ir, genes_ir, maxgap = 1L)
o_features   ## same as overlaps for sRNAs
  ## Join the overlapping features with the genomic ones for further cutting.
o_features <- c(o_features, genes_ir)
o_features  ## all of ranges in one irange object. overlapping + all of features
  ## Cut the overlapping features on teh border where they overlap with the genomic features.
split_f <- disjoin(o_features)
split_f[6:7]

subsetByOverlaps(split_f, genes_ir)
  ## Now select only the UTR "overhangs" that are created by cutting overlapping features on the border.
  #added brackets after the ! to correct error--not sure that solved anything
UTRs <- split_f[! (split_f %in% subsetByOverlaps(split_f, genes_ir))]
UTRs
  ## Select only UTRs that satisfy the minimum length condition.
min_utr<-2
UTRs <- UTRs[width(UTRs)>=min_utr,]
UTRs  
## Construct the IDs for the new UTRs to be added into the attribute colmn of the annotation.
  if (target_strand=="+") {
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste("ID=putative_UTR:p", x[1], "_", x[2],";", sep = ''))
  } else if (target_strand== "-") {
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste("ID=putative_UTR:m", x[1], "_", x[2],";", sep = ''))
  } else {
    return("Select strand")
  }
  return(UTRs)
  
}


##Error in h(simpleError(msg, call)) : 
#error in evaluating the argument 'i' in selecting a method for function '[': 'match' requires vector arguments
