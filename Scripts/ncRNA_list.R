library(dplyr)

#add annotations from Arnvig, 2014 to dataframe
arnvig_table <-
  read.delim(
    "Arnvig_ncRNA_2014.csv",
    sep = ",",
    header = TRUE,
    stringsAsFactors = F
  )
head(arnvig_table)

ncRNA_list <- arnvig_table$New.annotation
length(ncRNA_list)

start <- arnvig_table$X5...RACE.
stop <- arnvig_table$X3.
strand <- arnvig_table$Strand
verified_by <- arnvig_table$Verified.by
TSS <- arnvig_table$X5...TSS.
old_name <- arnvig_table$Other
references <- arnvig_table$Reference

arnvig_df <-
  data.frame(
    ncRNA_list,
    start,
    stop,
    strand,
    verified_by,
    TSS,
    old_name,
    references,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
colnames(arnvig_df) <-
  c("annotation",
    "start",
    "stop",
    "strand",
    "verified",
    "TSS",
    "aka",
    "ref")
head(arnvig_df)

## see below where I did using two statements and gsub (this is better)
arnvig_df$strand <- lapply(arnvig_df$strand, function(x) {
  ifelse(x == "F", "+", "-")
})
View(arnvig_df)

# filter for ncRNAs with start and stop
useful<-arnvig_df[arnvig_df$start != "—" & arnvig_df$stop != "—", ]
useful
unique_useful<-useful[which(!(useful$annotation %in% ncRNA_TSS_df$sRNA_ID)),]
unique_useful

# add these 3 to df
distance = abs(as.numeric(arnvig_df$start[16])-as.numeric(arnvig_df$TSS[16]))
distance
ncRv13661_row <- c("H37Rv", arnvig_df$annotation[16], arnvig_df$start[16], arnvig_df$stop[16],
                    arnvig_df$strand[16], 1,
                    arnvig_df$TSS[16], "", c("Arnvig, 2009, DiChiara,2010"))
ncRNA_tss_df
ncRNA_TSS_df<-rbind(ncRNA_TSS_df, ncRv13661_row)

View(ncRNA_TSS_df)

# annotations from Gerrick, 2018
gerrick_table <-
  read.csv(
    "Gerrick_sRNAs.csv",
    sep = ",",
    header = T,
    stringsAsFactors = F
  )
head(gerrick_table)
gerrick_df <- select(gerrick_table, 1, 3, 4, 2)
ref2 <- rep(c("Gerrick, 2018"))
gerrick_df <- cbind(gerrick_df, ref2)
colnames(gerrick_df) <-
  c("sRNA_ID", "start", "stop", "strand", "reference")
head(gerrick_df)

View(gerrick_df)



#############################################################################

# import shell and cortes TSS coordinates

cortes_table <-
  read.delim(
    "cortes_TSS_exp.csv",
    header = T,
    sep = ',',
    comment.char = "#",
    stringsAsFactors = F
  )

cortes_df1 <- select(cortes_table, 1, 2)
View(cortes_df1)

# change to "+" or "-" for strand
cortes_df1 <- data.frame(lapply(cortes_df1, function(x) {
  gsub("F", "+", x)
}))
cortes_df1 <- data.frame(lapply(cortes_df1, function(x) {
  gsub("R", "-", x)
}))
cortes_df1$Genome.position<-as.numeric(cortes_df1$Genome.position)
nrow(cortes_df1)

# cortes starved tss list
cortes_table2<-
  read.delim(
    "cortes_TSS_starv.csv",
    header=T,
    sep=",",
    comment.char="#",
    stringsAsFactors = F
  )

cortes_df2<-select(cortes_table2, 1, 2)
cortes_df2$Strand<- lapply(cortes_df2$Strand, function(x) {
  ifelse(x == "F", "+", "-")
})
cortes_df2$Genome.position<-cortes_df2$Genome.position
View(cortes_df2)
nrow(cortes_df2)

# # see if we can keep expression data
# condition <- rep(c("starv"))
# cortes2_test <- cbind(cortes_df2, condition)
# head(cortes2_test)
# condition <- rep(c("exp"))
# cortes_test <- cbind(cortes_df, condition)
# head(cortes_test)


# bind two df and select distinct(unique) coordinates (missing condition data)
cortes_df3<-rbind(cortes_df1, cortes_df2)
nrow(cortes_df3)
cortes_df3<-distinct(cortes_df3)
nrow(cortes_df3)
head(cortes_df3)
cortes_df3$Genome.position<-as.numeric((cortes_df3$Genome.position))

# write table of combined cortes TSSs
write.table(cortes_df3, "comb_cortesTSS_srna.txt", row.names = F, quote = F)


shell_table<-
  read.csv(
    "All_TSSs_shell.csv", 
    header=T, 
    stringsAsFactors = F)

shell_df<-select(shell_table, 3,2)
colnames(shell_df)<-c("Genome.position", "Strand")
head(shell_df)
nrow(shell_df)

shell_cortes_tss <- rbind(cortes_df3, shell_df)
nrow(shell_cortes_tss)
#[1] 13272
shell_cortes_tss <- distinct(shell_cortes_tss)
nrow(shell_cortes_tss)
#[1] 11494
# write table of combined TSSs from both cortes and shell (unique)
write.table(shell_cortes_tss, "shell_cortes_srna_tss.txt", quote = F, row.names = F)

##########################################################################

library(GenomicRanges)
# create a genomic range out of gerrick data and see if tss present (vector?)
# make another column for tss, y or n

Gerrick_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(gerrick_df$start, end = gerrick_df$stop),
  strand   = Rle(strand(gerrick_df$strand)),
  sRNA_ID  = gerrick_df$sRNA_ID
  )

#change df to have end coordinates for tss sites?
cortes_df3$end<-cortes_df3$Genome.position
head(cortes_df3)
shell_df$end<-shell_df$Genome.position

cortes_gr <-GRanges(
  seqnames = "H37Rv",
  ranges = IRanges(cortes_df3$Genome.position, end = cortes_df3$end),
  strand = Rle(strand(cortes_df3$Strand))
)

shell_gr <- GRanges(
  seqnames = "H37Rv",
  ranges = IRanges(shell_df$Genome.position, end = shell_df$end),
  strand = Rle(strand(shell_df$Strand))
)

################################################################################

## find overlap of each to tss within 20 bp of start
tss_overlaps <-
  findOverlaps(
    Gerrick_gr,
    cortes_gr,
    type = "start",
    maxgap = 20,
    ignore.strand = FALSE
  )
length(tss_overlaps)
queryHits(tss_overlaps)
subjectHits(tss_overlaps)
gerrick_hits <- as.data.frame(Gerrick_gr[queryHits(tss_overlaps)])
cortes_hits  <- as.data.frame(cortes_gr[subjectHits(tss_overlaps)])

# strict overlap--don't use
#gerrick_overlap<-Gerrick_gr[Gerrick_gr %over% cortes_gr]
#gerrick_with_TSS<-c(gerrick_overlap$sRNA_ID)
#print(gerrick_with_TSS)

# create extra column for TSS
gerrick_df$Cortes_TSS<-matrix("N/A", nrow(gerrick_df))
head(gerrick_df)
# make df with sRNA id and corresponding tss site
tss_pairs<-cbind(gerrick_hits$sRNA_ID, cortes_hits$start)
head(tss_pairs)

# insert TSS into df
for (i in 1:nrow(tss_pairs)){
  gerrick_df[which(gerrick_df$sRNA_ID == tss_pairs[i,1]), 6] <- tss_pairs[i,2]
}
View(gerrick_df)

# do same for shell tss sites
tss_overlaps_shell <-
  findOverlaps(
    Gerrick_gr,
    shell_gr,
    type = "start",
    maxgap = 20,
    ignore.strand = FALSE
  )
length(tss_overlaps_shell)
queryHits(tss_overlaps_shell)
gerrick_hits_shell <-
  as.data.frame(Gerrick_gr[queryHits(tss_overlaps_shell)])
shell_hits <- as.data.frame(shell_gr[subjectHits(tss_overlaps_shell)])
tss_pairs_shell <- cbind(gerrick_hits_shell$sRNA_ID, shell_hits$start)
head(tss_pairs_shell)

gerrick_df$Shell_TSS<-matrix("N/A", nrow(gerrick_df))
for (i in 1:nrow(tss_pairs_shell)){
  gerrick_df[which(gerrick_df$sRNA_ID == tss_pairs_shell[i,1]), 7] <-
    tss_pairs_shell[i,2]
}
View(gerrick_df)
  
# make table of all ncRNAs
write.table(gerrick_df, "gerrick_ncRNAs.txt", quote = FALSE, row.names = FALSE)
# make table of only ncRNAs with TSS
only_tss<- subset(gerrick_df, Cortes_TSS!="N/A" | Shell_TSS!="N/A")  
View(only_tss)
nrow(only_tss)
write.table(only_tss, "gerrick_with_TSS.txt", quote = FALSE, row.names = FALSE)


###########################################################################

## find number of nucleotides to closest TSS site
## distanceToNearest: Returns the distance for each range in x to its nearest 
## neighbor in the subject.
## For distanceToNearest, a Hits object with a column for the query index (queryHits), 
## subject index (subjectHits) and the distance between the pair.

## need to find distance from START of each gerrick ncRNA--either direction?
# make new gr object with single nucleotide start/end for gerrick

# need to make different df for one for each strand since stop==start on - strand
gerrick_fwd_df<-gerrick_df[gerrick_df$strand == "+", ]
nrow(gerrick_fwd_df)
gerrick_rev_df<-gerrick_df[gerrick_df$strand == "-", ]
nrow(gerrick_rev_df)

## make dataframe for each strand using only start coordinate
G_fwd_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(gerrick_fwd_df$start, end = gerrick_fwd_df$start),
  strand   = Rle(strand(gerrick_fwd_df$strand)),
  sRNA_ID  = gerrick_fwd_df$sRNA_ID
)

## use stop coordinate for reverse strand
G_rev_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(gerrick_rev_df$stop, end = gerrick_rev_df$stop),
  strand   = Rle(strand(gerrick_rev_df$strand)),
  sRNA_ID  = gerrick_rev_df$sRNA_ID
)

# tss_nearest_cortes<-
#   distanceToNearest(G_short_gr, cortes_gr, select=c("arbitrary"))
# head(tss_nearest_cortes)
# head(subjectHits(tss_nearest_cortes))
# length(subjectHits(tss_nearest_cortes))
# dist_to_tss<-as.data.frame(mcols(tss_nearest_cortes))
# head(dist_to_tss)
# cortes_tss_hits<-as.data.frame(cortes_gr[subjectHits(tss_nearest_cortes)])
# nrow(cortes_tss_hits)
# TSS_Cortes<-as.data.frame(cortes_tss_hits$start)
# colnames(TSS_Cortes)<-c("TSS")
# dist_to_cortes<-NULL
# dist_to_cortes <-
#   as.data.frame(G_short_gr[queryHits(tss_nearest_cortes)])
# dist_to_cortes<- cbind(dist_to_cortes, dist_to_tss, TSS_Cortes)
# dist_to_cortes$end<-gerrick_df$stop
# View(dist_to_cortes)

## add all TSS together and use that for more complete list
# make combination TSS gr object

total_tss_df<-rbind(cortes_df3, shell_df)
total_tss_df<-distinct(total_tss_df)
nrow(total_tss_df)
head(total_tss_df)

total_tss_gr <- GRanges(
  seqnames = "H37Rv",
  ranges = IRanges(total_tss_df$Genome.position, end = total_tss_df$Genome.position),
  strand = Rle(strand(total_tss_df$Strand))
)

# find distance to nearest tss for fwd strand
tss_nearest_fwd<-
  distanceToNearest(G_fwd_gr, total_tss_gr, select=c("arbitrary"))
head(tss_nearest_fwd)
head(subjectHits(tss_nearest_fwd))
length(subjectHits(tss_nearest_fwd))

#df column for distances
dist_to_tss<-NULL
dist_to_tss<-as.data.frame(mcols(tss_nearest_fwd))
head(dist_to_tss)
#extract df column for TSS coordinate from hits info
fwd_tss_hits<- NULL
fwd_tss_hits<-as.data.frame(total_tss_gr[subjectHits(tss_nearest_fwd)])
nrow(fwd_tss_hits)
TSS_fwd<-NULL
TSS_fwd<-as.data.frame(fwd_tss_hits$start)
colnames(TSS_fwd)<-c("TSS")

# make dataframe
f_dist_to_total<-NULL
f_dist_to_total <-
  as.data.frame(G_fwd_gr[queryHits(tss_nearest_fwd)])
f_dist_to_total<- cbind(f_dist_to_total, dist_to_tss, TSS_fwd)
f_dist_to_total$end<-gerrick_fwd_df$stop
#remove column 'width'
f_dist_to_total$width<-NULL
# reorder columns
f_dist_to_total<-f_dist_to_total[,c(1,5,2,3,4,6,7)]
View(f_dist_to_total)  

# # add column for source of TSS
# tss_source<-NULL
# tss_source<-as.data.frame(matrix(0, nrow=nrow(f_dist_to_total)), ncol=1)
# colnames(tss_source)[colnames(tss_source) == "V1"] <- "tss_source"
# head(tss_source)
# f_dist_to_total<-cbind(dist_to_total, tss_source)
# head(f_dist_to_total)

#determine source of TSS

## add source in
# for (i in 1:nrow(dist_to_total)){
#   if (dist_to_total$TSS[i] %in% cortes_df3$Genome.position & 
#       dist_to_total$TSS[i] %in% shell_df$Genome.position){
#       dist_to_total$tss_source[i] <- c("Cortes, Shell")
#   }else if (dist_to_total$TSS[i] %in% shell_df$Genome.position){
#     dist_to_total$tss_source[i] <- c("Shell")
#   }else{
#     dist_to_total$tss_source[i] <- c("Cortes")
#   }
# }
# 
# View(dist_to_total)
# #write table
# write.table(dist_to_total, "Gerrick_ncRNA_TSS_positions.txt", quote = F, row.names = F)

## way better way to avoid loops and use ifelse on all data (same exact results)
## also, don't need to create column first and cbind--creates as it goes along


#f_dist_to_total$tss_source <- NULL

f_dist_to_total$tss_source<- with(f_dist_to_total,
              ifelse (TSS %in% cortes_df3$Genome.position &
                TSS %in% shell_df$Genome.position, c("Cortes, Shell"),
              ifelse (TSS %in% shell_df$Genome.position, 
                c("Shell"), c("Cortes"))))

head(f_dist_to_total)

## do all again with reverse strand

# find distance to nearest tss for rev strand
tss_nearest_rev<-
  distanceToNearest(G_rev_gr, total_tss_gr, select=c("arbitrary"))
head(tss_nearest_rev)
head(subjectHits(tss_nearest_rev))
length(subjectHits(tss_nearest_rev))

#df column for distances
dist_to_tss<-NULL
dist_to_tss<-as.data.frame(mcols(tss_nearest_rev))
head(dist_to_tss)
#extract df column for TSS coordinate from hits info
tss_hits<- NULL
tss_hits<-as.data.frame(total_tss_gr[subjectHits(tss_nearest_rev)])
nrow(tss_hits)
TSS_rev<-NULL
TSS_rev<-as.data.frame(tss_hits$start)
colnames(TSS_rev)<-c("TSS")

# make dataframe
r_dist_to_total<-NULL
r_dist_to_total <-
  as.data.frame(G_rev_gr[queryHits(tss_nearest_rev)])
r_dist_to_total<- cbind(r_dist_to_total, dist_to_tss, TSS_rev)
head(r_dist_to_total)
head(gerrick_rev_df)
# change start of reverse strands to start from original 
# (used stop for gr)
r_dist_to_total$start<-gerrick_rev_df$start
#remove column 'width'
r_dist_to_total$width<-NULL
# reorder columns
r_dist_to_total<-r_dist_to_total[,c(1,5,2,3,4,6,7)]
View(r_dist_to_total)  
## all distances seem to be short by 1. should i add one to each 
## distance on rev strand or is this the convention? 

# add source col
r_dist_to_total$tss_source<- with(r_dist_to_total,
                                  ifelse (TSS %in% cortes_df3$Genome.position &
                                            TSS %in% shell_df$Genome.position, c("Cortes, Shell"),
                                          ifelse (TSS %in% shell_df$Genome.position, 
                                                  c("Shell"), c("Cortes"))))

head(r_dist_to_total)


## combine both fwd and rev into one dataframe, include column for source of ncRNA?

ncRNA_TSS_df<-rbind(f_dist_to_total, r_dist_to_total)
head(ncRNA_TSS_df)
length(ncRNA_TSS_df[ncRNA_TSS_df$distance<100,2])
#102

## add citation

ncRNA_TSS_df$citation<-rep(c("Gerrick, 2018"))
head(ncRNA_TSS_df)
nrow(ncRNA_TSS_df)
# 189

## add dejesus predictions to larger list of ncRNAs

for (i in 1:nrow(dejesus_tss_df)){
  if (!(dejesus_tss_df$sRNA_ID[i] %in% ncRNA_TSS_df$sRNA_ID)){
    ncRNA_TSS_df<-rbind(ncRNA_TSS_df, dejesus_tss_df[i,])
  }
}

nrow(ncRNA_TSS_df)
#209

# reorder by start
#ncRNA_TSS_df<-ncRNA_TSS_df[order(ncRNA_TSS_df$start),] 
## not sure this is good idea because messes with row numbers and 
## makes it hard to index? also doesn't really put in numeric order?
## put back to original
test<-ncRNA_TSS_df
ncRNA_TSS_df$index <- as.numeric(row.names(ncRNA_TSS_df))
ncRNA_TSS_df<-ncRNA_TSS_df[order(ncRNA_TSS_df$index), ]
head(ncRNA_TSS_df)

# find any repeat sRNAs in gerrick and other sources

base_list<-ncRNA_TSS_df[grep("ncRv[0-9]+$", ncRNA_TSS_df$sRNA_ID),2]
head(base_list)

#which(base_list == "ncRv1072")
#base_list[8]

repeats<-NULL
for (i in 1:length(base_list)){
  if (length(grep(base_list[i], ncRNA_TSS_df$sRNA_ID)) > 1){
    repeats<-c(repeats,i)
  }
}
repeats
# [1]   3   4   6  14  39  95  96  97  99 103 106
base_list[repeats]
#[1] "ncRv0897"  "ncRv1072"  "ncRv1181"  "ncRv11315" "ncRv10128"
#[6] "ncRv0186"  "ncRv12220" "ncRv10243" "ncRv13661" "ncRv10609"
#[11] "ncRv10071"

# manually curated the duplicates to make sure longest possible transcript
# most of these were just same transcript shifted one or two bp, 
# some were actually on opposite strand, so had same base name. 


# write table
write.table(ncRNA_TSS_df, "ncRNAs_and_TSS.txt", quote = FALSE, row.names = FALSE)




####################################################################


# can i make a 'histogram' to show distribution of distances from TSS?

close<-ncRNA_TSS_df$distance[which(ncRNA_TSS_df$distance<200)]
close<-sort(close)

length(close)
#138
max(close)
par(lwd=1)
dist_to_total$distance

plot(close,
     main=c("ncRNAs"), 
     type="p",
     xlim=c(0,length(close)),
     #xlab=c("ncRNA"),
     ylab = c("distance to nearest TSS"),
     #pch = 16 
     cex = .5
     )

old_par<-par()
par()
#divide data into 50bp regions for only 500bp
bins<-seq(0,200,10)
distances<-cut(close, bins)
transform(table(distances))

# histogram of frequency of distance to TSS
hist(close, 
     breaks=bins, 
     xlab="distance to TSS", 
     main="frequency of distance to TSS",
     ylim = c(0,100),
     xlim = c(0,200),
     col = "grey"
)

# violin plot
#install.packages("vioplot")
library(vioplot)
vioplot(close, ylab = "distance to TSS", names="predicted", xlab="ncRNAs with TSS w/in 200bp", col="lightblue")
# not helpful when i include all the data
vioplot(ncRNA_TSS_df$distance)



par(mfrow=c(1, 1))

plot(density(close), main="distance to TSS", col="green")


###############################################################################

## using pairs to find intersection of overlapping range (TSS site in hit)
## this restricts to just when TSS inside ncRNA 
p<-findOverlapPairs(Gerrick_gr, cortes_gr)
pintersect(p)
p
length(p)
cortes_tss<-as.data.frame(pintersect(p))
cortes_tss

gerrick_df$TSS<-matrix("N/A", nrow(gerrick_df))
head(gerrick_df)




common_arnvig<- unlist(lapply(ncRNA_df$annotation, function(x) {
  x %in% dist_to_total$sRNA_ID
}))
# show which ncRNAs are in common between gerrick and arnvig
ncRNA_df$annotation[common_arnvig]
# 5 in common:
# "ncRv11147c" "ncRv11264c" "ncRv11733"  "ncRv12659"  "ncRv13059"

ncRNA_df$annotation[!common_arnvig] ##shows what is not in common

head(dist_to_total)
head(ncRNA_df)

## could add to a df but is it worth it? just definitely test the arnvig ones, as well

combo<-dist_to_total

for (i in 1:nrow(ncRNA_df)){
  if (!(ncRNA_df$annotation[i] %in% dist_to_total$sRNA_ID[i])){
    r<-with(ncRNA_df, c("H37Rv", annotation[i], start[i], stop[i], strand[i], "", TSS[i], ""))
    combo<-rbind(combo, r)
  }
}
View(combo)
