#' Perform mixed effect model and create a heatmap of pairwise comparisons
#'
#' This function performs a mixed effect model using the provided data, extracts the pairwise comparison results,
#' and creates a heatmap visualization of the p-values from the comparisons. The heatmap is annotated with cell class information.
#'
#' @param meta_cell_info Meta cell information data frame.
#' @param processed_whole Processed whole data.
#' @param sample_info Sample information data frame.
#' @param my_color Color palette for the heatmap.
#'
#' @return A ComplexHeatmap object representing the pairwise comparison heatmap.
#'
#' @import ComplexHeatmap
#' @importFrom viridis viridis
#'
#' @examples
#' result <- perform_mixed_effect_model(meta_cell_info, processed_whole, sample_info, my_color)
#' model <- result$model
#' heatmap <- create_pairwise_comparison_heatmap(model)
#' draw(heatmap)
create_pairwise_comparison_heatmap <- function(meta_cell_info, processed_whole, sample_info, my_color) {
  result <- perform_mixed_effect_model(meta_cell_info, processed_whole, sample_info, my_color)
  model <- result$model
  # Testing pairwise comparisons using Tukey's method
  tukey_test <- glht(model, linfct = mcp(Subclass = "Tukey"))
  tukey_test_sum <- summary(tukey_test) # Summary of the test results
  confint(tukey_test)  # Adjusted p-values for the pairwise comparisons

  # Extracting comparison and p-value
  comp <- tukey_test_sum$linfct 
  p_values <- tukey_test_sum$test$pvalues

  # Converting p-values to a matrix format suitable for heatmap
  pairs <- strsplit(rownames(comp), " - ")  
  subclasses <- unique(unlist(pairs))
  p_matrix <- matrix(NA, nrow = length(subclasses), ncol = length(subclasses), dimnames = list(subclasses, subclasses))

  for (i in 1:length(p_values)) {
    pair <- pairs[[i]]
    p_matrix[pair[1], pair[2]] <- p_values[i]
    p_matrix[pair[2], pair[1]] <- p_values[i]  
  }

  p_matrix[is.na(p_matrix)] <- 1
  p_matrix_log <- -log10(p_matrix)

  # Create an ordered list of subclasses based on cell class categories
  cell_class_order <- c("GABAergic","Glutamatergic", "Non-neuronal") 
  ordered_subclasses <- names(cell_classes)[order(match(cell_classes, cell_class_order))]

  # Ensure the matrix rows and columns are in this order
  p_matrix_ordered <- p_matrix[ordered_subclasses, ordered_subclasses]

  cell_class_annotations <- factor(cell_classes[ordered_subclasses], levels = cell_class_order)

  # Create annotations for rows and columns
  col_annotation <- HeatmapAnnotation(cell_class = cell_class_annotations, which = "col", 
                                      col = list(cell_class = c(Glutamatergic = "#bca2cd", GABAergic = "#f5928a", "Non-neuronal" = "#c9bd2e")))
  row_annotation <- HeatmapAnnotation(cell_class = cell_class_annotations, which = "row",
                                      col = list(cell_class = c(Glutamatergic = "#bca2cd", GABAergic = "#f5928a", "Non-neuronal" = "#c9bd2e")))
  
  p_matrix_log_ordered <- -log10(p_matrix_ordered)
  
  # Create the heatmap with annotations
  ht <- Heatmap(p_matrix_log_ordered, name = "-log10(p-value)", col = viridis::viridis(256),
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = TRUE, show_column_names = TRUE, 
                row_names_side = "left", column_names_side = "bottom", top_annotation = col_annotation)
  
  return(ht)
}