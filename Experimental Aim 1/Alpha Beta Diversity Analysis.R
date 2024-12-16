#### ALPHA DIVERSITY ANALYSIS ####
library(ggplot2)
library(vegan)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)
library(packcircles)
library(cowplot)
library(phyloseq)

# load in the rarefied phyloseq object
load("final_filtered_phyloseq_rare.RData")
ms_rare <- final_filtered_phyloseq_rare



# Make data frame containing diversity measurements and treatment status
alpha_div <- ms_rare %>%
  estimate_richness(measures = c("Shannon", "Simpson", "Observed", "Chao1"))

alpha_div$disease_course <- factor(sample_data(ms_rare)$disease_course)

alpha_div$disease_course <- factor(alpha_div$disease_course, levels = c("RRMS", "SPMS", "PPMS"))



# Comparisons for plots
my_comparisons <- list(c("RRMS", "SPMS"),
                       c("SPMS", "PPMS"),
                       c("PPMS", "RRMS"))


# VIOLIN PLOT
# Simpson metric plot
simpson <- ggviolin(alpha_div, x = "disease_course", y = "Simpson", fill = "disease_course",
                        palette =c("#d95f02", "#1b9e77","#7570b3"),
                        add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = c(1.043, 1.07, 1.10)) + 
  
  labs(y = "Simpson Index", x = " ", fill = "Disease Course") +
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12)
  )


# shannon metric plot
shannon <- ggviolin(alpha_div, x = "disease_course", y = "Shannon", fill = "disease_course",
                    palette =c("#d95f02", "#1b9e77","#7570b3"),
                    add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = c(5.31, 5.49, 5.70)) + 
  
  labs(y = "Shannon Index", x = " ", fill = "Disease Course") +
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12)
  )

# Chao1 metric plot
chao1 <- ggviolin(alpha_div, x = "disease_course", y = "Chao1", fill = "disease_course",
                                                 palette =c("#d95f02", "#1b9e77","#7570b3"),
                                                 add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = c(234, 243, 254)) + 
     
  labs(y = "Chao1", x = " ", fill = "Disease Course") +
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12)
  )

# observed metric plot
observed <- ggviolin(alpha_div, x = "disease_course", y = "Observed", fill = "disease_course",
                         palette =c("#d95f02", "#1b9e77","#7570b3"),
                         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = c(234, 243, 254)) + 
  
  labs(y = "Observed", x = " ", fill = "Disease Course") +
  theme_classic()  +
  
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12)
  )


# combine alpha plots
alpha_combined_plot <- plot_grid(shannon, simpson, chao1, observed, labels = c("A", "B", "C", "D"), ncol = 4)


################################################################################
#### BETA DIVERSITY ANALYSIS ####

custom_colors <- c("SPMS" = "#1b9e77", "PPMS" = "#7570b3", "RRMS" = "#d95f02")


# compute distance matrices
unifrac_dist <- distance(ms_rare, method = "unifrac", weighted = TRUE) #weighted unifrac
bray_dist <- distance(ms_rare, method = "bray")


# perform PCoA
unifrac_pcoa <- ordinate(ms_rare, method = "PCoA", distance = unifrac_dist)
bray_pcoa <- ordinate(ms_rare, method = "PCoA", distance = bray_dist)

# statistical testing (with PERMANOVA)
meta_data <- data.frame(sample_data(ms_rare))

unifrac_permanova <- adonis2(unifrac_dist ~ meta_data$disease_course)
bray_permanova <- adonis2(bray_dist ~ meta_data$disease_course)


# plot PCoA
# unifrac metric
unifrac_plot <- plot_ordination(ms_rare, unifrac_pcoa, color = "disease_course") + 
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  
  labs(
    x = paste("PC1: 45% Variance"),
    y = paste("PC2: 15% Variance"),
    col = "Treatment Status"
  )  +
  
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)
  ) +
  
  stat_ellipse(type = "t", level = 0.95) +  #add ellipse
  annotate("text", x = -0.05, y = 0.15, 
           label = paste(""), size = 4, hjust = 0)

# bray Curtis metric
bray_plot <- plot_ordination(ms_rare, bray_pcoa, color = "disease_course") + 
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  
  labs(
    x = paste("PC1: 8.4% Variance"),
    y = paste("PC2: 6.3% Variance"),
    col = "Treatment Status"
  )   +
  
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)
  ) +
  
  stat_ellipse(type = "t", level = 0.95) +
  annotate("text", x = -0.05, y = 0.4, 
           label = paste(""), size = 4, hjust = 0)

# PERMANOVA RESULTS #
# weighted unifrac: R2 is 0.01868, p-value is 0.002
# bray: R2 is 0.01198, p-value is 0.001


#legend <- get_legend(simpson) //retrieve legend from simpson plot (Alpha Diversity Visual.R)


beta_combined_plot <- plot_grid(bray_plot, unifrac_plot, legend, labels = c("E", "F", " "), ncol = 3)

# combine alpha and beta diversity plots
combined_plot <- plot_grid(alpha_combined_plot, 
                           beta_combined_plot,
                           ncol = 1)

ggsave("disease_course_alpha_beta_take3.png", plot = combined_plot, width = 20, height = 12)
