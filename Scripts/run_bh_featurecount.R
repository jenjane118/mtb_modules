# usage: bh_count_features.R

# script to use count_features function to count reads for newly annotated ncRNA and UTRs in 
# multiple datasets

library(devtools)
#devtools::install_github("jenjane118/baerhunter", dependencies=FALSE)
# local edits
devtools::load_all("/d/in16/u/sj003/baerhunter")
library(baerhunter)
library(Rsamtools)
library(Rsubread)

## functions to get functions to re-define (from local version)
#source("/d/in16/u/sj003/baerhunter/R/feature_file_editor.R")

# create new directory for results

if (!(dir.exists("./output_BH_15_03/"))) { dir.create("./output_BH_15_03/")
}

# Read in list of datasets
directory_list <- scan("dataset_list.txt", what="", sep="\n")

for (i in 1: length(directory_list)){
  dataset_name<-sub("/d/in19/u/zchayyt/sample_accessions/", "",    directory_list[i])
  dataset_name<-sub("/BWA_mem/", "", dataset_name)
  output_directory<-paste("/d/in16/u/sj003/output_BH_15_03/", dataset_name)
  
  print(output_directory)
  
  # run feature_file_editor to create .gff3 file with new annotations
  #count_features(bam_dir=directory_list[i],
  #                  annotation_dir= "/d/in16/u/sj003/",
  #                  annotation_file="combined_gffs_01_03.gff3",
  #                  output_dir=output_directory,
  #                  chromosome_alias_file="/d/in16/u/sj003/chromosome.txt"
  #                  target_features=c("gene","putative_sRNA", "putative_UTR")
  #                  strandedness="reversely_stranded",
  #                  is_paired_end=TRUE
  #                  )
}

