### ALPHA DIVERSITY ANALYSIS ###
library(ggplot2)
library(vegan)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)
library(packcircles)
library(cowplot)

# load in the rarefied phyloseq object
ms_rare <- final_filtered_phyloseq_rare



# Make data frame containing diversity measurements and treatment status
alpha_div <- ms_rare %>%
  estimate_richness(measures = c("Shannon", "Simpson", "Observed", "Chao1"))

alpha_div$treatment_status <- factor(sample_data(ms_rare)$treatment_status)

# Comparisons for plots
my_comparisons <- list( c("Control", "Treated"), c("Treated", "Untreated"), c("Control", "Untreated"))




# BOXPLOT
box_plot <- ggboxplot(alpha_div, x = "treatment_status", y = "Shannon",
               color = "treatment_status", palette =c("#7570b3", "#1b9e77", "#d95f02"),
               add = "jitter", shape = "treatment_status") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = 6) +
  
  labs(y = "Shannon Index", x = " ", color = "Treatment Status", shape = "Treatment Status") +
  theme_classic()




# VIOLIN PLOT
violin_plot <- ggviolin(alpha_div, x = "treatment_status", y = "Shannon", fill = "treatment_status",
                        palette =c("#7570b3", "#1b9e77", "#d95f02"),
                        add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = c(4.9, 5.15, 5.4)) +
  stat_compare_means(label.y = 6) +
  
  labs(y = "Shannon Index", x = " ", fill = "Treatment Status") +
  theme_classic()




# BEESWARM PLOT
bee_plot <- ggplot(alpha_div, aes(x = treatment_status, y = Shannon, color = treatment_status)) +
  geom_beeswarm(size = 1.9, alpha = 1, cex = 1.5) + 
  scale_color_manual(values = c("#7570b3", "#1b9e77", "#d95f02")) + 
  scale_shape_identity(guide = "none") + 
  labs(
    y = "Shannon Index", 
    x = " ", 
    color = "Treatment Status"
  ) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  stat_compare_means(label.y = 6) +
  theme_classic()



alpha_combined_plot <- plot_grid(box_plot, violin_plot, bee_plot, box_plot2, labels = c("A", "B", "C"), ncol = 2)
ggsave("alpha_plot_types.png", plot = alpha_combined_plot, width = 9, height = 7)








