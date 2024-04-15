#' @title Install and load required packages
#'
#' @description This function checks if the required packages are installed and installs them if necessary.
#' It then loads the required libraries for the analysis.
#'
#' @return None
#'
#' @examples
#' install_and_load_packages()
#'
#' @import BiocManager
#' @import cowplot
#' @import tidyverse
#' @import MutationalPatterns
#' @import VariantAnnotation
#' @import GenomicRanges
#' @import BSgenome.Hsapiens.UCSC.hg38
#'
#' @export
install_and_load_packages <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  required_packages <- c(
    "MutationalPatterns",
    "BSgenome.Hsapiens.UCSC.hg38",
    "VariantAnnotation",
    "GenomicRanges",
    "tidyverse",
    "cowplot"
  )
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }
  
  library(MutationalPatterns)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(tidyverse)
  library(cowplot)
  
  theme_set(theme_cowplot())
}

#' @title Define cell classes
#'
#' @description This function defines a named vector of cell classes.
#'
#' @return A named vector of cell classes
#'
#' @examples
#' cell_classes <- define_cell_classes()
#'
#' @export
define_cell_classes <- function() {
  c(
    "L2/3 IT" = "Glutamatergic", "L4 IT" = "Glutamatergic", "L5/6 NP" = "Glutamatergic",
    "L5 ET" = "Glutamatergic", "L5 IT" = "Glutamatergic", "L6b" = "Glutamatergic",
    "L6 CT" = "Glutamatergic", "L6 IT" = "Glutamatergic", "L6 IT Car3" = "Glutamatergic",
    "Chandelier" = "GABAergic", "Lamp5" = "GABAergic", "Lamp5 Lhx6" = "GABAergic",
    "Pax6" = "GABAergic", "Sncg" = "GABAergic", "Sst" = "GABAergic",
    "Sst Chodl" = "GABAergic", "Vip" = "GABAergic", "Pvalb" = "GABAergic",
    "Astrocyte" = "Non-neuronal", "Oligodendrocyte" = "Non-neuronal",
    "OPC" = "Non-neuronal", "Endothelial" = "Non-neuronal",
    "Microglia-PVM" = "Non-neuronal", "VLMC" = "Non-neuronal"
  )
}

#' @title Perform mutational spectrum analysis
#'
#' @description This function performs mutational spectrum analysis on the provided GRangesList.
#'
#' @param grl A GRangesList object containing the genomic ranges and subclass information.
#' @param ref_genome The reference genome to use for the analysis (default: "BSgenome.Hsapiens.UCSC.hg38").
#'
#' @return A list containing the type occurrences, subclasses, class vector, and unnamed class vector.
#'
#' @examples
#' mutational_spectrum_analysis(grl)
#'
#' @export
mutational_spectrum_analysis <- function(grl, ref_genome = "BSgenome.Hsapiens.UCSC.hg38") {
  type_occurrences <- mut_type_occurrences(grl, ref_genome)
  subclasses <- names(grl)
  class_vector <- cell_classes[subclasses]
  class_vector_unnamed <- unname(class_vector)
  
  list(
    type_occurrences = type_occurrences,
    subclasses = subclasses,
    class_vector = class_vector,
    class_vector_unnamed = class_vector_unnamed
  )
}

#' @title Perform subset analysis
#'
#' @description This function performs mutational spectrum analysis on a subset of the provided GRangesList.
#'
#' @param grl A GRangesList object containing the genomic ranges and subclass information.
#' @param ref_genome The reference genome to use for the analysis (default: "BSgenome.Hsapiens.UCSC.hg38").
#'
#' @return A matrix containing the subset type occurrences for the selected subclasses.
#'
#' @examples
#' subset_analysis(grl)
#'
#' @export
subset_analysis <- function(grl, ref_genome = "BSgenome.Hsapiens.UCSC.hg38") {
  grSubset <- c(grl$Astrocyte, grl$`L2/3 IT`, grl$Pvalb)
  grlSubset <- split(grSubset, grSubset$Subclass)
  subset_type_occurrence <- mut_type_occurrences(grlSubset, ref_genome)
  subset_type_occurrence[c("Astrocyte", "L2/3 IT", "Pvalb"), ]
}

# Main analysis
install_and_load_packages()
cell_classes <- define_cell_classes()
processed_whole <- readRDS("~/projects/metadata/processed_whole_filtered_0.4.rds")

gr <- GRanges(
  seqnames = Rle(processed_whole$chr),
  ranges = IRanges(start = processed_whole$pos, end = processed_whole$pos),
  strand = Rle(rep("*", nrow(processed_whole))),
  REF = processed_whole$ref,
  ALT = processed_whole$alt,
  Subclass = processed_whole$Subclass
)
GenomeInfoDb::genome(gr) <- "hg38"
grl <- split(gr, gr$Subclass)

mutational_spectrum_results <- mutational_spectrum_analysis(grl)
subset_type_occurrence <- subset_analysis(grl)