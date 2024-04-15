#' Install and load required packages
#' 
#' This code chunk installs and loads the necessary packages for the analysis.
#' 
#' @import devtools
#' @import MissensePathoR
#' @import EnhancedVolcano
#' @import dplyr
#' @import edgeR
#'
devtools::install_github("Lola-W/MissensePathoR", build_vignettes = TRUE)
library(MissensePathoR)
library(dplyr)
library(edgeR)

#' Prepare variant data
#'
#' This function prepares the variant data by renaming columns and creating a new 'CHROM' column.
#' 
#' @param processed_gene_whole_genome_disease A data frame containing variant information.
#'
#' @return A data frame with renamed columns and a new 'CHROM' column.
#'
#' @examples
#' variantsample <- prepare_variant_data(processed_gene_whole_genome_disease)
prepare_variant_data <- function(processed_gene_whole_genome_disease) {
  processed_gene_whole_genome_disease %>%  
    rename(POS = pos, REF = ref, ALT = alt) %>%  
    mutate(CHROM = paste0("chr", chr))
}

#' Summarize pathogenicity scores
#'
#' This function summarizes pathogenicity scores by disease and cell subclass.
#'
#' @param prediction A data frame with predicted pathogenicity scores.
#'
#' @return None. The function prints summary statistics.
#'
#' @examples
#' summarize_patho_scores(prediction)
summarize_patho_scores <- function(prediction) {
  scoreSummary(prediction, category = "disease")
  scoreSummary(prediction, category = "Subclass")
  diseasePrediction <- prediction %>% rename(group = disease)
  classSummary(diseasePrediction) 
}

