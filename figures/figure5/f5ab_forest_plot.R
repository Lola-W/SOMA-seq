#' Create a forest plot of mixed-effect model estimates
#'
#' This function creates a forest plot of the mixed-effect model estimates for log10 mutation count.
#'
#' @param model Mixed-effect model object.
#' @param my_color Named vector of colors for each subclass.
#' @param cell_classes Named vector of cell classes for each subclass.
#' @param class_color Named vector of colors for each cell class.
#'
#' @return A ggplot object representing the forest plot.
#'
#' @examples
#' forest_plot <- create_forest_plot(model, my_color, cell_classes, class_color)
#'
#' @import jtools
#' @import cowplot
#'
#' @export
create_forest_plot <- function(model, my_color, cell_classes, class_color) {
  cell_classes_mapped <- cell_classes[names(my_color)[-1]]
  class_colors_mapped <- class_color[cell_classes_mapped]
  
  jtools::plot_summs(model,
                     model.names = c("Mixed-effect Model"),
                     legend.title = "Mixed-effect Model",
                     point.shape = TRUE) +
    labs(title = "Mixed-effect Model Estimates for log10 Mutation Count",
         y = "Cell Type") +
    theme_cowplot()
}