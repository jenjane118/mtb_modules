# long srnas

# run this on each dataset to see how many long srnas

library(rtracklayer)

annot <- import.gff3("~/git/mtb_modules/output_BH_07_12/PRJNA390669_12.gff3")
pred.sRNA <- subset(annot, type == "putative_sRNA")
pred.UTR <- subset(annot, type == "putative_UTR")
w<-width(pred.sRNA)
mean(w)
#522
min(w)
#47
max(w)
#6098
long_srnas<-as.data.frame(pred.sRNA[which(width(pred.sRNA)>1000),2:3])
write.table(long_srnas, "long_srnas.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
nrow(long_srnas)
#45

annot <- import.gff3("~/git/mtb_modules/output_BH_07_12/PRJNA327080_15.gff3")
pred.sRNA <- subset(annot, type == "putative_sRNA")
pred.UTR <- subset(annot, type == "putative_UTR")
w<-width(pred.sRNA)
mean(w)
#351
min(w)
#84
max(w)
#4130
long_srnas<-as.data.frame(pred.sRNA[which(width(pred.sRNA)>1000),2:3])
#write.table(long_srnas, "long_srnas.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
nrow(long_srnas)
#10
long_srnas

annot <- import.gff3("~/git/mtb_modules/output_BH_08_12/PRJNA278760_22.gff3")
pred.sRNA <- subset(annot, type == "putative_sRNA")
pred.UTR <- subset(annot, type == "putative_UTR")
w<-width(pred.sRNA)
mean(w)
#198
min(w)
#33
max(w)
#1287
long_srnas<-as.data.frame(pred.sRNA[which(width(pred.sRNA)>1000),2:3])
#write.table(long_srnas, "long_srnas.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
nrow(long_srnas)
#2

long_srnas
