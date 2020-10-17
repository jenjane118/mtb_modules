library(dplyr)

# Dejesus predicted RNAs 
dejesus_table <-
  read.csv(
    "dejesus_predicted_sRNA_list.csv",
    sep = ",",
    header = T,
    skip = 1,
    stringsAsFactors = F
  )
head(dejesus_table)

# make df from table
dejesus_df<-NULL
dejesus_df <-
    as.data.frame(matrix(0, nrow = nrow(dejesus_table), ncol = 4))
dejesus_df <- select(dejesus_table, 2, 3, 4, 5)
ref <- rep(c("DeJesus, 2017"))
dejesus_df <- cbind(dejesus_df, ref)
colnames(dejesus_df) <-
  c("sRNA_ID", "start", "stop", "strand", "reference")
head(dejesus_df)

# make fwd and rev df's
f_dejesus_df<-dejesus_df[dejesus_df$strand == "+", ]
r_dejesus_df<-dejesus_df[dejesus_df$strand == "-", ]

## see which ncRNAs aren't in Gerrick
unique_dejesus_df<-NULL

for (i in 1:nrow(dejesus_df)){
  if (!(dejesus_df$sRNA_ID[i] %in% gerrick_df$sRNA_ID)){
    unique_dejesus_df<-rbind(unique_dejesus_df, dejesus_df[i,])
  }
}
unique_dejesus_df
nrow(unique_dejesus_df)
#20

library(GenomicRanges)
# create a genomic range out of dejesus data (use stop==start for distance)
#but need to create different one for each strand

F_dejesus_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(f_dejesus_df$start, end = f_dejesus_df$start),
  strand   = Rle(strand(f_dejesus_df$strand)),
  sRNA_ID  = f_dejesus_df$sRNA_ID
)

R_dejesus_gr <- GRanges(
    seqnames = "H37Rv",
    ranges   = IRanges(r_dejesus_df$stop, end = r_dejesus_df$stop),
    strand   = Rle(strand(r_dejesus_df$strand)),
    sRNA_ID  = r_dejesus_df$sRNA_ID
)

# use gr of all tss sites 
head(total_tss_df)
class(total_tss_gr)

## for forward strand
## find distance to TSS (hits object) (need to use different gr with just start nucleotide)
tss_nearest_fwd<-NULL
tss_nearest_fwd<-
  distanceToNearest(F_dejesus_gr, total_tss_gr, select=c("arbitrary"))

# creating distance column from metadata in hits object
dist_to_tss<-NULL
dist_to_tss<-as.data.frame(mcols(tss_nearest_fwd))
head(dist_to_tss)
nrow(dist_to_tss)
#39

# dataframe of hits from subjectHits of hits object 
# matched to actual tss gr object to get tss sites
tss_hits<-NULL
tss_hits<-as.data.frame(total_tss_gr[subjectHits(tss_nearest_fwd)])
nrow(tss_hits)
head(tss_hits)
#39

# extract column of actual tss sites (one column) from dataframe above
TSS_fwd<-as.data.frame(tss_hits$start)
colnames(TSS_fwd)<-c("TSS")
head(TSS_fwd)


#create dataframe with genomic range object 
f_dist_dejesus<-NULL
f_dist_dejesus <-
  as.data.frame(F_dejesus_gr[queryHits(tss_nearest_fwd)])
head(f_dist_dejesus)

f_dist_dejesus<- cbind(f_dist_dejesus, dist_to_tss, TSS_fwd)
f_dist_dejesus$end<-f_dejesus_df$stop
head(f_dist_dejesus)
nrow(f_dist_dejesus)

#remove column 'width'
f_dist_dejesus$width<-NULL
# reorder columns
f_dist_dejesus<-f_dist_dejesus[,c(1,5,2,3,4,6,7)]

## reverse columns

## find distance to TSS (hits object) (need to use different gr with just start nucleotide)
tss_nearest_rev<-NULL
tss_nearest_rev<-
  distanceToNearest(R_dejesus_gr, total_tss_gr, select=c("arbitrary"))

# creating distance column from metadata in hits object
dist_to_tss<-NULL
dist_to_tss<-as.data.frame(mcols(tss_nearest_rev))
head(dist_to_tss)
nrow(dist_to_tss)
#39, 23

# dataframe of hits from subjectHits of hits object 
# matched to actual tss gr object to get tss sites
tss_hits<-NULL
tss_hits<-as.data.frame(total_tss_gr[subjectHits(tss_nearest_rev)])
nrow(tss_hits)
head(tss_hits)
#39, 23

# extract column of actual tss sites (one column) from dataframe above
TSS_rev<-NULL
TSS_rev<-as.data.frame(tss_hits$start)
colnames(TSS_rev)<-c("TSS")
nrow(TSS_rev)
#39, 23

#create dataframe with genomic range object 
r_dist_dejesus<-NULL
r_dist_dejesus <-
  as.data.frame(R_dejesus_gr[queryHits(tss_nearest_rev)])
head(r_dist_dejesus)

r_dist_dejesus<- cbind(r_dist_dejesus, dist_to_tss, TSS_rev)
r_dist_dejesus$end<-r_dejesus_df$stop
head(r_dist_dejesus)
nrow(r_dist_dejesus)

#remove column 'width'
r_dist_dejesus$width<-NULL
# reorder columns
r_dist_dejesus<-r_dist_dejesus[,c(1,5,2,3,4,6,7)]

## bind fwd and rev 
dejesus_tss_df<-NULL
dejesus_tss_df<-rbind(f_dist_dejesus, r_dist_dejesus)
nrow(dejesus_tss_df)
#62

#add source of TSS
# add source col
dejesus_tss_df$tss_source<- with(dejesus_tss_df,
                                  ifelse (TSS %in% cortes_df3$Genome.position &
                                            TSS %in% shell_df$Genome.position, c("Cortes, Shell"),
                                          ifelse (TSS %in% shell_df$Genome.position, 
                                                  c("Shell"), c("Cortes"))))

head(dejesus_tss_df)

# add column for citation
dejesus_tss_df$citation<-rep(c("Dejesus, 2017"))
head(dejesus_tss_df)

## add to larger list of ncRNAs

for (i in 1:nrow(dejesus_tss_df)){
  if (!(dejesus_tss_df$sRNA_ID[i] %in% ncRNA_TSS_df$sRNA_ID)){
    ncRNA_TSS_df<-rbind(ncRNA_TSS_df, dejesus_tss_df[i,])
  }
}

nrow(ncRNA_TSS_df)
#209

View(ncRNA_TSS_df)
