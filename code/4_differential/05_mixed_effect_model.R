#' Perform mixed-effect model analysis
#'
#' This function performs mixed-effect model analysis on the given data using the lme4 and lmerTest packages.
#'
#' @param meta_cell_info Data frame containing metadata information for each cell.
#' @param processed_whole Data frame containing processed variant data.
#' @param sample_info Data frame containing sample information.
#' @param my_color Named vector of colors for each subclass.
#'
#' @return A list containing the mixed-effect model object and the child data frame used for modeling.
#'
#' @examples
#' result <- perform_mixed_effect_model(meta_cell_info, processed_whole, sample_info, my_color)
#' model <- result$model
#'
#' @import lme4
#' @import lmerTest
#' @import dplyr
#'
#' @export
perform_mixed_effect_model <- function(meta_cell_info, processed_whole, sample_info, my_color = NA) {
  # Filter and process data
  filtered_meta_cell_info <- meta_cell_info %>%
    filter(barcode %in% processed_whole$barcode)
  
  processed_meta_cell_info <- filtered_meta_cell_info %>%
    group_by(donor_id, Subclass) %>%
    summarise(total_read_count = sum(read_count),
              distinct_barcodes = n_distinct(barcode),
              .groups = 'drop') %>%
    mutate(log_read_count = total_read_count,
           log_cell_count = distinct_barcodes)
  
  combined_data <- processed_meta_cell_info %>%
    left_join(sample_info[, c("external_donor_name", "ROI", "Sex", "Age", "Disease")], 
              by = c("donor_id" = "external_donor_name"))
  
  meta_cell_info <- subset(combined_data, select = -c(total_read_count, distinct_barcodes))
  
  processed_whole$freq <- 1
  child <- processed_whole %>%
    # filter(ref == "C" & alt == "T" & !is.na(Subclass)) %>%
    group_by(Subclass, donor_id) %>% 
    summarise(total_freq = sum(freq)) %>%
    mutate(total_freq = log10(total_freq + 0.001))
  
  child <- inner_join(meta_cell_info[!is.na(meta_cell_info$Subclass), ], child, by = c("Subclass", "donor_id")) %>% 
    mutate(log_read_count = log10(log_read_count),
           log_cell_count = log10(log_cell_count),
           Age = log10(Age))
  
  # Refactor base on provided color order
  if (!is.na(my_color)){
  child$Subclass <- factor(child$Subclass, levels = names(my_color))
  }
  
  # Fit the mixed-effect model
  ## MODIFY IT: with your desired model
  model <- lmer(total_freq ~ Subclass + log_read_count + log_cell_count + Age + Sex + Disease + (1|donor_id), data = child)
  
  return(list(model = model, child = child))
}