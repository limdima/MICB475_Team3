### BETA DIVERSITY ANALYSIS ###
ms_rare <- final_filtered_phyloseq_rare

custom_colors <- c("Control" = "#7570b3", "Treated" = "#1b9e77", "Untreated" = "#d95f02")


# compute distance matrices
unifrac_dist <- distance(ms_rare, method = "unifrac", weighted = TRUE) #weighted unifrac!
bray_dist <- distance(ms_rare, method = "bray")


# perform PCoA
unifrac_pcoa <- ordinate(ms_rare, method = "PCoA", distance = unifrac_dist)
bray_pcoa <- ordinate(ms_rare, method = "PCoA", distance = bray_dist)

# statistical testing
meta_data <- data.frame(sample_data(ms_rare))
unifrac_permanova <- adonis2(unifrac_dist ~ meta_data$treatment_status)
bray_permanova <- adonis2(bray_dist ~ meta_data$treatment_status)


# plot PCoA
unifrac_plot <- plot_ordination(ms_rare, unifrac_pcoa, color = "treatment_status") + 
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  labs(
    x = paste("PC1: 48% Variance"),
    y = paste("PC2: 13.2% Variance"), 
    col = "Treatment Status"
  ) +
  theme(
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10)
  ) +
  stat_ellipse(type = "t", level = 0.95) +
  annotate("text", x = -0.05, y = 0.15, 
           label = paste("PERMANOVA p-value = 0.002"), size = 4, hjust = 0)


bray_plot <- plot_ordination(ms_rare, bray_pcoa, color = "treatment_status") + 
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  labs(
    x = paste("PC1: 8.2% Variance"),
    y = paste("PC2: 6.9% Variance"), 
    col = "Treatment Status"
  ) +
  theme(
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10)
  ) +
  stat_ellipse(type = "t", level = 0.95) +
  annotate("text", x = -0.05, y = 0.4, 
           label = paste("PERMANOVA p-value = 0.001"), size = 4, hjust = 0)


# PERMANOVA RESULTS #
# weighted unifrac: R2 is 0.01868, p-value is 0.002
# bray: R2 is 0.01198, p-value is 0.001

beta_combined_plot <- plot_grid(unifrac_plot, bray_plot, labels = c("E", "F"), ncol = 2)
ggsave("beta_plots.png", plot = beta_combined_plot, width = 13, height = 6)
