#' Perform EDA plots on processed mutation data
#'
#' This script generates various exploratory data analysis (EDA) plots to visualize mutation frequencies and loads across cell types and donors.
#'
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#' @import ggpubr
#' @import gridExtra
#'
#' @param processed_whole Processed mutation data frame
#' @param filtered_summary Filtered summary data frame
#' @param my_color Named vector of colors for subclasses
#' @param class_color Named vector of colors for classes
#'
#' @return Multiple plots visualizing mutation frequencies and loads
# Data preparation
processed_whole$freq <- 1
child <- processed_whole %>%
  group_by(Subclass, donor_id, disease) %>%
  filter(ref == "C" & alt == "T") %>%
  summarise(total_freq = sum(freq))

child <- inner_join(filtered_summary[!is.na(filtered_summary$Subclass), ], child, by = c("Subclass", "donor_id"))

# Plotting setup
child$Subclass <- factor(child$Subclass, levels = names(my_color))

# Boxplot of mutation frequency by subclass
freq_sub <- child %>%
  ggplot(aes(x = Subclass, y = total_freq/total_barcodes, fill = Subclass)) +
  geom_boxplot() +
  scale_fill_manual(values = my_color) +
  labs(title = "Mutation Frequency Distribution by Cell Type (Subclass)",
       x = "Cell Type (Subclass)",
       y = "Mutation Count / Cell Count",
       fill = "Subclass") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot of mutation load by subclass 
load_sub <- child %>%
  ggplot(aes(x = Subclass, y = log10(total_freq)/log10(total_read_count), fill = Subclass)) +
  geom_boxplot() +
  theme_cowplot() +
  scale_fill_manual(values = my_color) +
  labs(title = "Mutation Load by Cell Type (Subclass)",
       x = "Cell Type (Subclass)",
       y = "log10(Mutation Count) / log10(Read Count)",
       fill = "Subclass") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Scatter plot of read count vs mutation count  
f3a1 <- child %>%
  ggplot(aes(x = total_read_count, y = total_freq)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE) +
  ggpubr::stat_cor(method = "pearson") +
  labs(title = "Read Count vs. Mutation Count",
       x = "Read Count",
       y = "Mutation Count") +
  theme_cowplot()

# Scatter plot of cell count vs mutation count
f3a2 <- child %>%
  ggplot(aes(x = total_barcodes, y = total_freq)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE) +
  ggpubr::stat_cor(method = "pearson") +
  labs(title = "Cell Count vs. Mutation Count",
       x = "Cell Count", 
       y = "Mutation Count") +
  theme_cowplot()

gridExtra::grid.arrange(f3a2, f3a1, ncol = 2)  

# Bar plot of mutation frequency by donor in L2/3 IT cells
#' @rdname mutation_frequency_donor_l23it
#' @export
child %>%
  filter(Subclass == "L2/3 IT") %>%
  ggplot(aes(x = reorder(donor_id, -total_freq/total_barcodes), y = total_freq/total_barcodes, fill = disease)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Epilepsy" = "#68bd45", "Tumor and epilepsy" = "#7e7e7e", "Tumor" = "#8150a0"),
                    name = "Disease Status") +
  labs(title = "Mutation Frequency by Donor in L2/3 IT Cells",
       x = "Donor",
       y = "Mutation Count / Cell Count") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# Boxplot of percentage of mutated L2/3 IT cells  
filtered_summary %>%
  filter(Subclass == "L2/3 IT") %>%
  ggplot(aes(x = Subclass, y = sm_barcodes/total_barcodes, fill = Subclass)) +
  geom_boxplot() +
  labs(title = "Percentage of Mutated L2/3 IT Cells",
       y = "Mutated Cell / Total Cell",
       fill = "Subclass") +
  theme_cowplot() +
  scale_fill_manual(values = my_color)

# Class level analysis  
cell_classes <- c("L2/3 IT" = "Glutamatergic", "L4 IT" = "Glutamatergic", ...)

child_class <- child %>%
  mutate(Class = cell_classes[Subclass]) %>%
  group_by(Class, donor_id, disease) %>%
  summarise(total_freq = sum(total_freq),
            total_read_count = sum(total_read_count),
            total_barcodes = sum(total_barcodes), 
            sm_read_count = sum(sm_read_count),
            sm_barcodes = sum(sm_barcodes),
            .groups = "drop")

# Boxplot of mutation load by class
load_class <- child_class %>%
  ggplot(aes(x = Class, y = log10(total_freq)/log10(total_read_count), fill = Class)) +
  geom_boxplot() +
  scale_fill_manual(values = class_color) +
  labs(title = "Mutation Load by Cell Type (Class)",
       x = "Cell Type (Class)",
       y = "log10(Mutation Count) / log10(Read Count)",
       fill = "Class") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplot of mutation frequency by class  
freq_class <- child_class %>%
  ggplot(aes(x = Class, y = total_freq/total_barcodes, fill = Class)) +
  geom_boxplot() +
  scale_fill_manual(values = class_color) +
  labs(title = "Mutation Frequency Distribution by Cell Type (Class)",
       x = "Cell Type (Class)",
       y = "Mutation Count / Cell Count",
       fill = "Class") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))