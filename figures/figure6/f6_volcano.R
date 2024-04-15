#' Plot Volcano Plots for NEBULA Results
#'
#' This function generates a volcano plot for differential expression data, highlighting significant genes
#' based on log fold change and p-value thresholds. It is tailored for NEBULA analysis results, using
#' gene symbols for labeling and adjusting plot dimensions automatically to accommodate data ranges.
#'
#' @param dataset Data frame containing gene expression results with columns for genes, fold change, and p-values.
#' @param FCname Name of the column in `dataset` containing log fold changes.
#' @param pName Name of the column in `dataset` containing p-values.
#' @param typeName Descriptive name for the condition or type, used in the plot title.
#'
#' @return A ggplot object representing the volcano plot.
#'
#' @import ggplot2
#' @import EnhancedVolcano
#'
#' @examples
#' volcano_plot <- plot_volcano_nebula(data, "logFC_Epilepsy", "p_Value", "Epilepsy vs Control")
#' print(volcano_plot)
plot_volcano_nebula <- function(dataset, FCname, pName, typeName = "Disease") {
  # Prepare inputs
  sample_volcano <- data.frame(
    ensembl_Symbol = dataset$gene,  # Use the 'gene' column for gene identifiers
    LogFC = dataset[, FCname],  # Log fold change
    P.Val = dataset[, pName]  # Unadjusted p-values
  )
  
  # Filter out infinite and NA values
  sample_volcano <- na.omit(sample_volcano)
  sample_volcano <- sample_volcano[!is.infinite(sample_volcano$LogFC) & !is.infinite(sample_volcano$P.Val),]
  
  # Ensure there's enough data to plot
  if (nrow(sample_volcano) < 2) {
    stop("Not enough data to plot.")
  }
  
  # Adjust plot limits
  xlim_vals <- range(sample_volcano$LogFC, na.rm = TRUE) + c(-1, 1)
  ylim_vals <- range(-log10(sample_volcano$P.Val), na.rm = TRUE) + c(-1, 1)
  
  # Generate volcano plot
  EnhancedVolcano(sample_volcano,
                  lab = sample_volcano$ensembl_Symbol,
                  x = 'LogFC',
                  y = 'P.Val',
                  xlim = xlim_vals,
                  ylim = ylim_vals,
                  pCutoff = 0.05,
                  FCcutoff = 1.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75,
                  legendPosition = 'right',
                  labSize = 3.0,
                  legendLabSize = 12,
                  title = paste0(typeName, " vs Control"),
                  subtitle = "P-value cutoff = 0.05, Log Fold Change cutoff = 1")
}

#' Plot Mutational Rate for Specified Gene, use for EDA
#'
#' This function calculates and visualizes the mutational rate of a specified gene across different subclasses.
#' The rate is normalized by the logarithm of the total read count, providing insights into the mutational burden.
#'
#' @param enrichData Data frame of enriched gene data including counts.
#' @param barcode_count_sum Data frame summarizing read counts by barcode.
#' @param gene_name Name of the gene to analyze.
#'
#' @return A ggplot object representing the mutational rate as a boxplot.
#'
#' @import ggplot2
#' @importFrom dplyr filter left_join summarise mutate
#' @importFrom scales log10
#'
#' @examples
#' mutational_plot <- plot_mutational_rate(enrichData, barcode_count_sum, "Gene1")
#' print(mutational_plot)
plot_mutational_rate <- function(enrichData, barcode_count_sum, gene_name) {
  # Filter data for the specified gene
  enrichData_filtered <- enrichData %>%
    dplyr::filter(ensembl_gene == gene_name)
  
  # Merge with barcode counts to calculate mutational rates
  merged_data <- enrichData_filtered %>%
    dplyr::left_join(barcode_count_sum, by = c("donor_id", "Subclass", "sample_name")) %>%
    dplyr::mutate(mutational_rate = count / log10(total_read_count + 1))  # Adjustment to avoid log(0)
  
  # Generate and return plot
  ggplot(merged_data, aes(x = Subclass, y = mutational_rate, fill = Subclass)) +
    geom_boxplot() +
    labs(title = paste("Mutational Rate for Gene", gene_name), 
         x = "Subclass", 
         y = "Mutational Rate (count/log10(read count))") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")
}
