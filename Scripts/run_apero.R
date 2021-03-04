# run apero
# https://github.com/Simon-Leonard/APERO

install.packages("reshape2")
install.packages("snowfall")
devtools::install_github("Simon-Leonard/APERO")

library(Rsamtools)
library(reshape2)
library(snowfall)

library(APERO)

# Load annotation file (needs to be in .ptt form?)
#ptt=read.csv("E:/APERO/NC_016810.ptt",sep="\t",skip=2,header=T,stringsAsFactors = F)

# write ptt (protein table) for use 
#Location    Strand    Length    PID    Gene    Synonym    Code    COG    Product

# load csv from ncbi (tabular protein file)
h37<-read.csv("ref_seqs/proteins_166_159857.csv", header = TRUE, stringsAsFactors = F)
head(h37)
h37rv<-NULL

#h37_ptt<-data.frame(paste(h37$Start, h37$Stop, sep=".."), h37$Strand, h37$Length, 
#                h37$Protein.product, h37$Locus.tag, h37$Locus, " ", " ", 
#                h37$Protein.Name)
#colnames(h37_ptt)<-c("Location", "Strand", "Length", "PID", "Gene", "Synonym",
#                     "Code", "COG", "Product")
#head(h37_ptt)
#write.table(h37_ptt, "H37Rv.ptt", sep="\t", quote = F, row.names = F)

#load annotation file

h37_ptt<-read.csv("ref_seqs/H37Rv.ptt", sep="\t",header=T,stringsAsFactors = F)
head(h37_ptt)

#wmax= min accepted width of start peak
# min_dist = minimum dist between separated start sites

# 5'end detection
res<-APERO_start_detection(work_dir = "~/bam_files/PRJNA327080_15/", bam_name = "SRR3725585_sorted.bam",
                          ptt_file = h37_ptt, wmax = 10, min_dist = 10, enrichment = 0.1, min_read_number = 0, genome_size = 4411532
)

write.table(res, file="apero_5prime_SRR3725585.txt", sep="\t", quote=F)

head(res)
nrow(res)


#3'end detection
res2<-APERO_end_detection(work_dir = "~/bam_files/PRJNA327080_15/", start_table = res, mTEX_bam = "SRR3725585_sorted.bam",
                         readthrough_proportion = 0.01, Fmin=NA, thread_number = 8, genome_size = 4411532, ptt_file=h37_ptt
)


bam_file<-"SRR3725585_sorted.bam"
basename(bam_file)
tools::file_path_sans_ext(bam_file)
bam_name<-sub('_sorted.bam', '', bam_file) 
bam_name
paste("/d/in16/u/sj003/APERO/apero_3prime_", bam_name, ".txt", sep="")
print("finished")
