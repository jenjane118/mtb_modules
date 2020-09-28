library(dplyr)

arnvig_table<-read.delim("Arnvig_ncRNA_2014.csv", sep=",", header=TRUE, stringsAsFactors=F)
head(arnvig_table)



ncRNA_list<-arnvig_table$New.annotation
length(ncRNA_list)
ncRNA_list

gerrick_table<-read.delim("Gerrick2018_pred_sRNAs.csv", sep=",", header=T, stringsAsFactors=FALSE)
head(gerrick_table)


for (i in 1:length(gerrick_table$name)){
  if (!(gerrick_table$name[i] %in% ncRNA_list)){
    ncRNA_list<-c(ncRNA_list, gerrick_table$name[i])
  }
}

length(ncRNA_list)
ncRNA_list

dejesus_table<-read.delim("dejesus_predicted_sRNA_list.csv", header=T, sep=",", comment.char="#")
head(dejesus_table)


for (i in 1:length(dejesus_table$sRNA.ID..)){
  if (!(dejesus_table$sRNA.ID..[i] %in% ncRNA_list)){
    ncRNA_list<-c(ncRNA_list, dejesus_table$sRNA.ID..[i])
  }
}

length(ncRNA_list)

