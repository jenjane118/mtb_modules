## make list of ncRNAs with TSS within 10bp of start

# read in table from file "~/git/mtb_modules/ncRNA_verified.txt"
library(dplyr)

ncrna_table <-
  read.csv(
    "ncRNA_verified.txt",
    sep = " ",
    header = T,
    stringsAsFactors = F
  )
head(ncrna_table)
nrow(ncrna_table)
#204 predicted from rna-seq

# new ncrna list with tss within 10bp
srna_tss<-ncrna_table[which(ncrna_table$distance < 11),]
nrow(srna_tss)
srna_tss
# 67 sRNAs predicted with TSS within 10 bp
t<-c(1,2,3,4,5,7)
write.table(srna_tss[,t], "Data/gerrick_with_tss.txt", quote = FALSE, row.names = FALSE)


# lets see how many have no blots?
arnvig_table <-
  read.delim(
    "Data/Arnvig_ncRNA_2014.csv",
    sep = ",",
    header = TRUE,
    stringsAsFactors = F
  )

arnvig_srnas<-arnvig_table$New.annotation
srna_tss[which(srna_tss$sRNA_ID %in% arnvig_srnas),]
# only 3 have been verified with northern blot (?) (ncRv11264c, ncRv11733-DrrS, ncRv12659)

# 82 diff expressed in at 1+ conditions
de_cond_1<-read.delim(
    "Data/diff_ex_srnas_gerrick.csv",
    sep = ",",
    header = TRUE,
    stringsAsFactors = F
  )
# how many have TSS?
de_tss<-srna_tss[which(srna_tss$sRNA_ID %in% de_cond_1$name),]
nrow(de_tss)
# 63 are differentially expressed in at least one condition and have TSS within 10nt
srna_tss[!which(srna_tss$sRNA_ID %in% de_cond_1$name),2]
# all of 67 TSS including sRNAs are also dif expressed

# these 3 were diff expressed in 3+ conditions (ncRv11803, ncRv11846, and ncRv12659) in Gerrick
check<-c("ncRv11803", "ncRv11846", "ncRv12659")
ncrna_table[which(ncrna_table$sRNA_ID %in% check),]
# these have TSS within 1 nt

library(rtracklayer)
RNA_filter<-list(type="ncRNA")
myco_rnas<-readGFF('ref_seqs/mycobrowser_H37Rv_gff_v3.gff', filter=RNA_filter)
colnames(myco_rnas)
head(myco_rnas)
nrow(myco_rnas)

p<-c(3,4,5,7,14,15)
write.table(myco_rnas[,p], "Data/mycobrowser_ncrnas.txt", quote = FALSE, row.names = FALSE)


# read in entire gff file and filter with rtracklayer
#my_cols<-c("seqid","source","type","start","end")
RNA_filter<-list(type="ncRNA")
# column 'Locus'=MTB000115' and 'Name' = ncRv10666, 'Comments' =supported by sRNA-Seq in H37Rv
# compare with 'Name' to list of diff expressed sRNAs from gerrick (de_cond_1)
f<-de_cond_1[which(de_cond_1$name %in% myco_rnas),]
f
nrow(f)
#42 are in list on mycobrowser (must use DeJesus list?)
# compare with list of tss verified (de_tss)
t<-de_tss[which(de_tss$sRNA_ID %in% myco_rnas$Name), 1:4]
t
nrow(t)
colnames(myco_rnas)

# find which ones have TSS in first 10 nt
# import table
tss_all<-read.delim("Data/shell_cortes_srna_tss.txt", delim='\t', header = TRUE, stringsAsFactors = F )
head(tss_all)
# find which ones have TSS in first 10 nt
# import table
tss_all<-read.delim("Data/shell_cortes_srna_tss.txt", sep=' ', header = TRUE, stringsAsFactors = F )
head(tss_all)
View(myco_rnas)
# how many of srna-seq are verified with TSS at start??
# make myco_rnas into genomic ranges obj for each strand
myco_fwd_df<-myco_rnas[myco_rnas$strand == "+", ]
nrow(myco_fwd_df)
myco_rev_df<-myco_rnas[myco_rnas$strand == "-", ]
nrow(myco_rev_df)

f_mycobrowser_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(myco_fwd_df$start, end = myco_fwd_df$start),
  strand   = Rle(strand(myco_fwd_df$strand)),
  sRNA_ID  = myco_fwd_df$Name,
  locus    = myco_fwd_df$Locus
)

# for rev strand
r_mycobrowser_gr <- GRanges(
  seqnames = "H37Rv",
  ranges   = IRanges(myco_rev_df$end, end = myco_rev_df$end),
  strand   = Rle(strand(myco_rev_df$strand)),
  sRNA_ID  = myco_rev_df$Name,
  locus    = myco_rev_df$Locus
)
tss_nearest_fwd<-
  distanceToNearest(f_mycobrowser_gr, total_tss_gr, select=c("arbitrary"))
head(tss_nearest_fwd)
head(subjectHits(tss_nearest_fwd))
length(subjectHits(tss_nearest_fwd))
dist_to_tss<-NULL
dist_to_tss<-as.data.frame(mcols(tss_nearest_fwd))
head(dist_to_tss)
fwd_tss_hits<- NULL
fwd_tss_hits<-as.data.frame(total_tss_gr[subjectHits(tss_nearest_fwd)])
nrow(fwd_tss_hits)
TSS_fwd<-NULL
TSS_fwd<-as.data.frame(fwd_tss_hits$start)
colnames(TSS_fwd)<-c("TSS")
head(dist_to_tss)
#extract df column for TSS coordinate from hits info
fwd_tss_hits<- NULL
fwd_tss_hits<-as.data.frame(total_tss_gr[subjectHits(tss_nearest_fwd)])
fwd_tss_hits
TSS_fwd<-NULL
TSS_fwd<-as.data.frame(fwd_tss_hits$start)
colnames(TSS_fwd)<-c("TSS")
# make dataframe
f_dist_to_total<-NULL
f_dist_to_total <- as.data.frame(f_mycobrowser_gr[queryHits(tss_nearest_fwd)])
f_dist_to_total<- cbind(f_dist_to_total, dist_to_tss, TSS_fwd)
f_dist_to_total$end<-myco_fwd_df$stop
#remove column 'width'
f_dist_to_total$width<-NULL
View(f_dist_to_total)

# find distance to nearest tss for rev strand
tss_nearest_rev<-
  distanceToNearest(r_mycobrowser_gr, total_tss_gr, select=c("arbitrary"))
head(tss_nearest_rev)
head(subjectHits(tss_nearest_rev))
length(subjectHits(tss_nearest_rev))
# put in df
#df column for distances
dist_to_tss<-NULL
dist_to_tss<-as.data.frame(mcols(tss_nearest_rev))
head(dist_to_tss)
#extract df column for TSS coordinate from hits info
rev_tss_hits<- NULL
rev_tss_hits<-as.data.frame(total_tss_gr[subjectHits(tss_nearest_rev)])
nrow(rev_tss_hits)
TSS_rev<-NULL
TSS_rev<-as.data.frame(rev_tss_hits$start)
colnames(TSS_rev)<-c("TSS")
# make dataframe
r_dist_to_total<-NULL
r_dist_to_total <- as.data.frame(r_mycobrowser_gr[queryHits(tss_nearest_rev)])
r_dist_to_total<- cbind(r_dist_to_total, dist_to_tss, TSS_rev)
r_dist_to_total$end<-myco_rev_df$stop
#remove column 'width'
r_dist_to_total$width<-NULL
View(r_dist_to_total)
#put together
myco_TSS_df<-rbind(f_dist_to_total, r_dist_to_total)
head(myco_TSS_df)
#order by start
ordered_myco_tss<-myco_TSS_df[order(myco_TSS_df$start),]
nrow(ordered_myco_tss)
#how many within 10nt
tss_verified_myco<-ordered_myco_tss[ordered_myco_tss$distance<11,]
myco_names<-tss_verified_myco$sRNA_ID
nrow(tss_verified_myco)
#38 within 10 nt of TSS

myco_names
#write to file
write.table(tss_verified_myco, "Data/mycobrowser_with_TSS.txt", quote = FALSE, row.names = FALSE)

#assemble table of all 'supported' ncrnas, add column about TSS
head(myco_rnas)
myco<-as.data.frame(myco_rnas)
head(tss_verified_myco)
library(dplyr)

master_df<- select(myco, 'Name', 'Locus', 'start', 'end', 'strand')
head(master_df)
# if start in myco_rnas
for (i in 1: nrow(master_df)){
  if (master_df$Name[i] %in% tss_verified_myco$sRNA_ID){
    master_df$TSS[i]<-'TSS'
  }
  else{
    master_df$TSS[i]<-''
  }
}
nrow(master_df[which(master_df$TSS != ''),])
# 38 TSS found in list
# add in whether northern/rlm-race

northern_race<-readGFF("Data/verified_myco.txt")
for (i in 1: nrow(master_df)){
  if (master_df$Name[i] %in% northern_race$Name){
    master_df$verified[i]<-"RLM-RACE"
  }
  else{
    master_df$verified[i]<-''
  }
}

View(master_df)

# file of all ncRNAs verified either by TSS mapping or experimental
a<-master_df[which(master_df$TSS !='' | master_df$verified !=''),]
nrow(a)
write.table(master_df[which(master_df$TSS !='' | master_df$verified !=''),],
            "Data/verified_list_review.txt", quote = FALSE, row.names = FALSE)



##############################################################

