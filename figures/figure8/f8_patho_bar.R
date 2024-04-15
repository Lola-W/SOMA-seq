#' Visualize the distribution of pathogenicity classes across groups
#'
#' This function plots the distribution of pathogenicity classes (benign, possibly damaging, probably damaging) across disease groups using the MissensePathoR package and the ggplot2 and cowplot packages for additional formatting [[1]](https://poe.com/citation?message_id=157442654852&citation=1)[[2]](https://poe.com/citation?message_id=157442654852&citation=2).
#'
#' @param prediction A data frame containing the predicted pathogenicity scores and disease information.
#'
#' @return A ggplot object visualizing the distribution of pathogenicity classes across disease groups.
#'
#' @examples
#' diseasePrediction <- prediction %>% rename(group = disease)
#' classVis(diseasePrediction) +
#'   theme_cowplot() +
#'   labs(
#'     title = "Distribution of Pathogenicity Classes Across Disease Groups",
#'     y = "Pathogenicity Class",
#'     x = "Disease Group"
#'   )
#'
#' @import MissensePathoR
#' @import ggplot2
#' @import cowplot
#' @import dplyr
#'
#' @export
visualize_pathogenicity_classes <- function(prediction) {
  diseasePrediction <- prediction %>% rename(group = disease)
  
  classVis(diseasePrediction) +
    theme_cowplot() +
    labs(
      title = "Distribution of Pathogenicity Classes Across Disease Groups",
      y = "Pathogenicity Class",
      x = "Disease Group"
    )
}