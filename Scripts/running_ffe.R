## run feature file editor from baerhunter 
## re-define functions from inside script
## 07/12/20


library(GenomicAlignments)
library(IRanges)
library(Rsamtools)
library(devtools)
library(Rsubread)
# to get local edited versions:
devtools::load_all("~/R_projects/baerhunter")
# to get versions from my cloned bh repo:
#devtools::install_github("jenjane118/baerhunter", dependencies=FALSE)
library(baerhunter)

## functions to get functions to re-define
source("~/R_projects/baerhunter/R/feature_file_editor.R")

#source_url("https://github.com/jenjane118/baerhunter/blob/master/R/feature_file_editor.R")
# this doesn't work: error in source(temp_file)

feature_file_editor(bam_directory="~/bam_files/PRJNA390669_12/",
                    original_annotation_file = "ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3",
                    annot_file_dir = ".",
                    output_file = "test.3.12.gff3",
                    original_sRNA_annotation = "ncRNA", 
                    min_sRNA_length=30, 
                    min_UTR_length=50,
                    paired_end_data=TRUE, 
                    strandedness ="reversely_stranded")
