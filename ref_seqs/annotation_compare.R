## compare ncRNAs in genbank, mycobrowser, ena/ncbi annotation files

# get genbank ncRNA entries

BiocManager::install("genbankr")
browseVignettes("genbankr")
library(genbankr)

h37rv<-readGenBank("GCF_000195955.2_ASM19595v2_genomic.gbff")
head(h37rv)

of<-otherFeatures(h37rv)
head(of)
ncrna<-of[which(of$type=="ncRNA"),]
ncrna_df<-data.frame(ncrna$locus_tag, ncrna$gene, ncrna$product, ncrna$note)
head(ncrna_df)
nrow(ncrna_df)
View(ncrna_df)

# get ncbi ncRNA entries

library(rtracklayer)
RNA_filter<-list(type="ncRNA")
h37gff<-readGFF("MtbH37RvNC_000962.3.gff")
ncbi_df<-as.data.frame(h37gff)
colnames(h37gff)
ncbi_df<-select(ncbi_df, 3,16,18,19,20, 26)
head(ncbi_df)
ncbi_ncrna<-ncbi_df[which(ncbi_df$gene_biotype=="ncRNA"),]
View(ncbi_ncrna)
#20
# all the same as in genbank
# others as 'sequence features' including 'fragment of sRNA'
sf<-ncbi_df[which(ncbi_df$type=="sequence_feature"),]
View(sf)
u<-h37gff[which(!h37gff$locus_tag %in% ncrna_df$ncrna.locus_tag), ]

# ena
ena_filter<-list(type="ncRNA_gene")
h37gff2<-readGFF("Mtb_h37rv.ASM19595v2_AL123456.3.gff3", filter=ena_filter)
nrow(h37gff2)
#30
ena_names<-h37gff2$Name
ena_names[! ena_names %in% ncrna_df$ncrna.gene]
#[1] "mcr19"     "mpr5"      "ncrMT1234" "mpr6"      "mcr5"      "mpr11"     "mpr12"    
#[8] "mpr17"     "mpr18"     "ncrMT3949"

# mycobrowser overlap
myco<-readGFF("mycobrowser_H37Rv_gff_v3.gff", filter=RNA_filter)
colnames(myco)
type_myco<-myco[which(myco$type=="ncRNA"),]
nrow(type_myco)
#92
v_ncrna<-myco[which( myco$functional_Category=="stable RNAs"),]
nrow(v_ncrna)
v_ncrna$Name
#30
dejesus_myco<-myco[which(myco)]

#There seem to be 20 genbank/ncbi/ena ncRNAs
#10 more include northern blot found 'fragments' of srnas
#mycobrowser includes 62 gerrick predictions from RNAseq
