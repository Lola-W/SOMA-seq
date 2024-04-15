# Dependencies
install_packages <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("EnhancedVolcano")
  
  needed_packages <- c("edgeR", "EnhancedVolcano", "ggplot2", "gprofiler2", "tidyverse", 
                       "gridExtra", "nebula", "dplyr")
  
  new_packages <- needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
  if (length(new_packages) > 0) {
    install.packages(new_packages, dependencies = TRUE)
  }
  
  sapply(needed_packages, require, character.only = TRUE)
}

#' Perform mutation analysis using NEBULA-HL
#'
#' This function performs mutation analysis using the NEBULA-HL tool. It calculates the sum of mutational
#' frequencies for each gene within pseudobulks, aggregated for each cell type group in each donor.
#' The function uses a single-cell differential expression analysis tool, NEBULA-HL, for its accuracy
#' in estimating the cell-level overdispersion.
#'
#' @param counts_matrix A matrix containing the mutation counts.
#' @param meta_cell_info A data frame containing the metadata for each cell.
#' @param sample_info A data frame containing the sample information.
#'
#' @return The result of the NEBULA-HL analysis.
#'
#' @import edgeR EnhancedVolcano ggplot2 gprofiler2 tidyverse gridExtra nebula dplyr
#'
#' @examples
#' result <- perform_mutation_analysis(counts_matrix, meta_cell_info, sample_info)
#'
perform_mutation_analysis <- function(counts_matrix, meta_cell_info, sample_info) {
  
  # Filtering and processing meta_cell_info
  filtered_meta_cell_info <- meta_cell_info %>%
    filter(barcode %in% processed_whole$barcode)
  
  processed_meta_cell_info <- filtered_meta_cell_info %>%
    group_by(donor_id, Subclass) %>%
    summarise(total_read_count = sum(read_count),
              distinct_barcodes = n_distinct(barcode),
              .groups = 'drop') %>%
    mutate(total_cell_count = distinct_barcodes)
  
  combined_data <- processed_meta_cell_info %>%
    left_join(sample_info[, c("external_donor_name", "ROI", "Sex", "Age", "Disease")], 
              by = c("donor_id" = "external_donor_name"))
  
  meta_cell_info <- subset(combined_data, select = -c(distinct_barcodes))
  
  # Calculating barcode count sum
  barcode_count_sum <- meta_cell_info %>%
    dplyr::group_by(Subclass, donor_id) %>%
    summarise(total_read_count = sum(total_read_count), .groups = "drop") %>%
    dplyr::mutate(sample_name = paste0(donor_id, "-", Subclass))
  
  # Data cleaning and preparation
  data <- data.table::setDT(enrichData_subclass[!is.na(enrichData_subclass$Subclass),])
  counts_matrix <- data.table::dcast(data[!is.na(data$ensembl_gene),], ensembl_gene ~ sample_name, value.var = "count", fun.aggregate = sum)
  counts_matrix <- data.frame(counts_matrix, row.names = 1)
  keep <- rowSums(counts_matrix) > 0
  counts_matrix_filtered <- counts_matrix[keep, ]
  counts_matrix <- data.matrix(counts_matrix_filtered)
  
  donor_ids <- sapply(colnames(counts_matrix_filtered), function(x) {
    sub("^([^\\.]+\\.[^\\.]+\\.[^\\.]+).*", "\\1", x)
  })
  group_info <- match(donor_ids, samples$sample_name)
  disease_group <- samples$group[group_info]
  cell_types <- gsub("^([^\\.]+\\.[^\\.]+\\.[^\\.]+).", "", colnames(counts_matrix_filtered))
  
  pred_data <- data.frame(disease_group, cell_type = cell_types)
  pred_data$disease_group <- factor(pred_data$disease_group)
  pred_data$cell_type <- factor(pred_data$cell_type)
  pred_data$cell_type <- relevel(pred_data$cell_type, "Pvalb") # MODIFY IT: could change the cell type to desired refer
  
  df <- model.matrix(~ disease_group + cell_type, data = pred_data)
  
  # Performing NEBULA-HL analysis
  counts_matrix <- counts_matrix[, colnames(counts_matrix) %in% barcode_off$donor_id]
  donor_ids <- colnames(counts_matrix)
  re_subclass <- nebula(counts_matrix, donor_ids, pred = df, cpc = 0, mincp = 3, ncore = 5)
  
  return(re_subclass)
}


#' Automate rank calculation and saving for GSEA
#'
#' This function calculates the rank for each gene based on the log fold change and p-value columns,
#' orders the genes by rank, and saves the ranked list as a .rnk file for GSEA analysis.
#'
#' @param data A data frame containing the log fold change and p-value columns for each comparison.
#' @param output_dir The directory where the output .rnk files will be saved. Default is the current working directory.
#'
#' @return None. The function saves the ranked list files in the specified output directory.
#'
#' @examples
#' automate_rank_and_save(res, "~/projects/metadata/enrichment/rnk_files/cpc")
#'
#' @import dplyr
#' @export
automate_rank_and_save <- function(data, output_dir = getwd()) {
  data <- data[complete.cases(data),]
  base_names <- gsub("logFC_", "", names(data)[grepl("logFC_", names(data))])[-1]
  
  for (base in base_names) {
    logFC_col <- paste0("logFC_", base)
    p_col <- paste0("p_", base)
    
    if (!logFC_col %in% names(data) | !p_col %in% names(data)) {
      next 
    }
    
    data$rank <- (-log10(data[[p_col]])) * sign(data[[logFC_col]])
    ranked_list <- data[order(data$rank, decreasing = T),c("gene", "rank")]
    
    output_filename <- file.path(output_dir, paste0("/", base, ".rnk"))
    write.table(ranked_list, file = output_filename, sep = "\t", 
                quote = FALSE, row.names = FALSE)
  }
}

