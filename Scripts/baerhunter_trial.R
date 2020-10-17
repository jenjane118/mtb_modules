## using baerhunter
## using: https://github.com/irilenia/baerhunter/blob/master/vignettes/baerhunter.Rmd


install.packages("devtools")
## answer "Yes"
library(devtools)

## need dependencies 'Rsubread', 'DESeq2'
## Skipping 5 packages not available: 
## DESeq2, Rsubread, Rsamtools, GenomicAlignments, IRanges
library(GenomicAlignments)
library(IRanges)
library(Rsamtools)

BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("Rsubread")
library(Rsubread)

devtools::install_github("irilenia/baerhunter")
library(baerhunter)

## with vignettes (not working)
#devtools::install_github("irilenia/baerhunter", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), force=TRUE)
#vignette(all=FALSE)
# not present

## Updating the genome annotation using the RNA-seq signal
## We start by calling baerhunter's 
## *feature_file_editor()* function to predict intergenic elements 
## (sRNAs and UTRs) based on the RNA-seq signal and existing annotation 
## (available as a gff3 file).
## The bam files used in this analysis are reduced versions of the originals 
## covering only the first 10,000 positions in the genome.  

# create a directory to hold the output files..
if (!(dir.exists("./output_BH_5_10/"))) { dir.create("./output_BH_5_10/")
}
feature_file_editor(bam_directory=system.file("extdata/", package="baerhunter"),  
                    original_annotation_file="Mycobacterium_tuberculosis_h37rv.ASM19595v2.40.chromosome.Chromosome.gff3",
                    annot_file_dir = system.file("extdata/", package="baerhunter"),
                    output_file="output_BH_5_10/mtb_5_10.gff3", 
                    original_sRNA_annotation="ncRNA", 
                    low_coverage_cutoff=5, 
                    min_sRNA_length=40, 
                    high_coverage_cutoff=10, 
                    min_UTR_length=50, 
                    paired_end_data=FALSE, 
                    strandedness="stranded")

## Counting reads against the new annotation produced by baerhunter
## Once new putative ncRNA features are added to the genome annotation, we use the *count_features()* function to count reads against both the original and newly annotated features.

count_features(bam_dir = system.file("extdata", package="baerhunter"),
               annotation_dir= "output_BH_5_10/", 
               annotation_file = "mtb_5_10.gff3", 
               output_dir = "output_BH_5_10/",
               chromosome_alias_file = system.file("extdata","chromosome.txt", package="baerhunter") , 
               target_features = c("gene", "putative_sRNA", "putative_UTR"), 
               strandedness = "stranded", 
               is_paired_end= FALSE)

## Filtering predictions by level of expression
## Occasionally, it is preferred to filter out low-expressed transcripts, 
## both because of noise in the RNA-seq data and because features with low 
## expression are unlikely to be useful for further downstream analysis. 
## The *tpm_flag_filtering()* function is used here to keep only putative 
## ncRNAs with higher expression. 

## To filter transcripts by expression, we calculate first TPM 
## (transcripts per million) values for each feature of interest using 
## the *tpm_normalisation()* function, then add expression level flags 
## to the annotation file and finally filter out transcripts with flags 
## corresponding to lower expression values.

# Filter out low expression putative sRNAs
# Calculate TPM values
tpm<- tpm_normalisation(count_table="output_BH_5_10/putative_sRNA_Counts.csv", 
                        complete_gff="output_BH_5_10/mtb_5_10.gff3", 
                        target_feature="putative_sRNA", 
                        output_file="output_BH_5_10/putative_sRNA_TPM.csv")
#repeat for putative UTRs
tpm <- tpm_normalisation(count_table="output_BH_5_10/putative_UTR_Counts.csv", 
                         complete_gff="output_BH_5_10/mtb_5_10.gff3", 
                         target_feature="putative_UTR", 
                         output_file="output_BH_5_10/putative_UTR_TPM.csv")
# produce a single file of TPM values for all non-coding RNAs
system("( cat output_BH_5_10/putative_sRNA_TPM.csv ; tail -n +2 output_BH_5_10/putative_UTR_TPM.csv ; ) | cat > output_BH_5_10/putative_ncRNA_TPM.csv")

# flag features according to their TPM values
tpm_flagging( tpm_data="output_BH_5_10/putative_ncRNA_TPM.csv", 
              complete_annotation="output_BH_5_10/mtb_5_10.gff3", 
              output_file="output_BH_5_10/flagged_mtb_5_10.gff3")
# filter features according to their flags
tpm_flag_filtering( complete_annotation_file="output_BH_5_10/flagged_mtb_5_10.gff3", 
                    target_flag="high_expression_hit", 
                    target_features=c("putative_sRNA", "putative_UTR"), 
                    output_file="output_BH_5_10/filtered_high_expression_sRNA_mtb_5_10.gff3")
