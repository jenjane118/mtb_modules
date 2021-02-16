# test bh with only parameter change (changed bug back)
# 3/12/2020
#first load libraries
library(IRanges)
library(GenomicAlignments)
library(Rsamtools)
library(devtools)

remove.packages("baerhunter")
devtools::load_all("~/R_projects/baerhunter")
load_all("~/R_projects/baerhunter")
library(baerhunter)

# then redefine functions on feature_file_editor.R

#call function for test bams (1)

feature_file_editor(bam_directory="~/bam_files/PRJNA390669_12/",
                    original_annotation_file = "ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3",
                    annot_file_dir = ".",
                    output_file = "test2t.3.12.gff3",
                    original_sRNA_annotation = "ncRNA", 
                    min_sRNA_length=30, 
                    min_UTR_length=50,
                    paired_end_data=TRUE, 
                    strandedness ="reversely_stranded")


