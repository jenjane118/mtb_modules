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
library(DESeq2)
BiocManager::install("Rsubread")
library(Rsubread)

devtools::install_github("irilenia/baerhunter", dependencies = FALSE)
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
feature_file_editor(bam_directory=,  
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



### functions inside feature_file_editor
# runs for each strand
peak_union_calc <- function(bam_location = "/PRJNA327080_15", target_strand="+", low_coverage_cutoff=300, high_coverage_cutoff=450,  peak_width=40, paired_end_data = TRUE, strandedness  = "reversely-stranded") {
  ## Find all BAM files in the directory.
  bam_files <- list.files(path = bam_location, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
  
  
  # use GenomicAlignments (dependent on Rsamtools) to identify mapped sequence start/end
  # from .bam file
  
  peak_union <- IRanges()
  ## Go over each BAM file to extract coverage peaks for a target strand and 
  # gradually build a union of all peak sets.
  for (f in bam_files) {
    ## Read a BAM file in accordance with its type and select only the reads aligning to a target strand.
    strand_alignment <- c()
    if (paired_end_data == FALSE & strandedness  == "stranded") {
      file_alignment <- readGAlignments(f)
      strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]
    } else if (paired_end_data == TRUE & strandedness  == "stranded") {
      # positively stranded, strandMode = 1 
      file_alignment <- readGAlignmentPairs(f, strandMode = 1)
      strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]
    } else if (paired_end_data == FALSE & strandedness  == "reversely_stranded") {
      file_alignment <- readGAlignments(f)
      relevant_strand <- c()
      if (target_strand=="+") {
        relevant_strand <- "-"
      } else {
        relevant_strand <- "+"
      }
      strand_alignment <- file_alignment[strand(file_alignment)==relevant_strand,]
    } else if (paired_end_data == TRUE & strandedness  == "reversely_stranded") {
      # reversely stranded, strandMode = 2
      file_alignment <- readGAlignmentPairs(f, strandMode = 2)
      strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]
    }
    
    ## Create a strand coverage vector and extract it
    # IRanges function 'coverage' counts the number of ranges over each position.
    strand_cvg <- coverage(strand_alignment)
    list_components <- names(strand_cvg)
    target <- c()
    if (length(list_components)==1) {
      target <- list_components
    } else {
      return(paste("Invalid BAM file:",f, sep = " "))
    }
    
    ## at this point we could determine coverage quantiles and use 5% (?) threshold 
    ## to establish parameters for selecting boundaries of expression peaks
    
    #viewApply is IRanges method
    
    ## Cut the coverage vector to obtain the expression peaks with the coverage above the low cut-off values.
    peaks <- slice(strand_cvg[[target]], lower = low_coverage_cutoff, includeLower=TRUE)
    ## Examine the peaks for the stretches of coverage above the high cut-off. The stretches have to be a defined width.
    ## this uses intrinsic function (below) peak_analysis
    test <- viewApply(peaks, function(x) peak_analysis(x,high_coverage_cutoff,peak_width))
    ## Select only the peaks that satisfy the high cut-off condition.
    selected_peaks <- peaks[lapply(test, function(x) !is.null(x))==TRUE]
    ## Convert peak coordinates into IRanges.
    peaks_IRange <- IRanges(start = start(selected_peaks), end = end(selected_peaks))
    ## Calculate the peak union in with the previous peak sets.
    peak_union <- union(peak_union,peaks_IRange)
  }
  return(peak_union)
}


#' Peak checking for the second coverage threshold and width.
#' 
#' This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
#' 
#' @param View_line A line from a RleViews object.
#' @param high_cutoff An integer indicating the high coverage threshold value.
#' @param min_sRNA_length An integer indicating the minimum sRNA length (peak width).
#' 
#' @return Returns a RleViews line if it satisfies conditions.
#' 
#' @export
peak_analysis <- function(View_line, high_cutoff, min_sRNA_length) {
  ## This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
  cvg_string <- as.vector(View_line)
  target_peak <- which(cvg_string>high_cutoff)
  if (length(target_peak)>min_sRNA_length) {
    return(View_line)
  }
}

