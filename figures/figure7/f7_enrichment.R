#' Visualize top enriched gene sets from GSEA results
#' 
#' This function takes the output from GSEA analysis and creates a scatter plot of the top enriched gene sets,
#' with the size of the points representing the gene set size, color representing the FDR q-value, 
#' and the x-axis representing the enrichment score.
#'
#' @param gsea_result The output data frame from GSEA analysis.
#' @param title The title of the plot.
#' @param n_pathways The number of top pathways to plot. Default is 10.
#'
#' @return A ggplot object representing the top enriched gene sets.
#'
#' @examples
#' plot_gsea_results(up_reg, "Top Mutation Enriched Gene sets in the Epilepsy Group")
#' plot_gsea_results(neg_reg, "Top Mutation Enriched Gene sets in the Tumor Group", n_pathways = 9)
#' 
#' @import ggplot2
#' @import cowplot
#' @export
plot_gsea_results <- function(gsea_result, title, n_pathways = 10) {
  
  gsea_result <- gsea_result[order(gsea_result$FDR.q.val), ]
  
  res <- data.frame(
    pathway = gsub("\\%.*", "", gsea_result$NAME),
    size = gsea_result$SIZE,
    FDR = gsea_result$FDR.q.val,
    es = gsea_result$ES
  )
  
  ggplot(res[1:n_pathways,]) + 
    geom_point(aes(
      x = es, 
      color = FDR,
      y = pathway,
      size = size)) +
    theme(axis.title.x = element_text(),
          axis.title.y = element_text()) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Enrichment Score",
         color = "FDR", 
         size = "Gene set size", 
         y = NULL,
         title = title) +
    theme_cowplot()
}