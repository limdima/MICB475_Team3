# Alpha and Beta Analysis
# Note, final_filtered_ms_phyloseq and final_filtered_phyloseq_rare are available on Github, 
# You can load these objects in and start this script at the Alpha Diversity step

# libraries
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)

# load in phyloseq object
load("final_filtered_ms_phyloseq.RData")

# Rarefy phyloseq object
ms_phyloseq_rare <- rarefy_even_depth(final_filtered_ms_phyloseq, rngseed = 4222, sample.size = 9488)



#### FILTERING: RAREFIED PHYLOSEQ ####
# make rarefied phyloseq object into a data frame:
phyloseq_rare_df <- as(sample_data(ms_phyloseq_rare), "data.frame")

# filter for pairs
filtered_phyloseq_rare_df <- phyloseq_rare_df %>%
  mutate(sample_id = rownames(phyloseq_rare_df)) %>% #keep a column for sample_ids
  
  group_by(household) %>%
  filter(n() == 2) %>%
  ungroup()

# prune the rarefied phyloseq to only have sample ids in filtered phyloseq dataframe
final_filtered_phyloseq_rare <- prune_samples(filtered_phyloseq_rare_df$sample_id, ms_phyloseq_rare)

# save filtered phyloseq rarefied
save(final_filtered_phyloseq_rare, file = "final_filtered_phyloseq_rare.RData")



### ALPHA DIVERSITY ANALYSIS ###
# alpha (Shannon) diversity for all samples
ms_rare <- final_filtered_phyloseq_rare
shannon_div <- estimate_richness(ms_rare, measures = "Shannon")

# add a column for treatment status

shannon_div$treatment_status <- sample_data(ms_rare)$treatment_status

# create a box-plot to visualize diversity metrics

shan_vs_treatment <- ggplot(shannon_div, aes(x=treatment_status, y=Shannon)) +
  geom_boxplot() +
  labs(title = "Shannon Diversity vs Treatment Status", 
       x = "Treatment Status",
       y = "Shannon Diversity Index") +
  theme_minimal()

plot(shan_vs_treatment)
ggsave("shannon_div.png", plot = shan_vs_treatment, width = 6, height = 4)

### BETA DIVERSITY ANALYSIS ###
# beta (weighted unifrac) diversity for all samples

w_unifrac_dist <- UniFrac(ms_rare, weighted=TRUE)

# PCoA on weighted unifrac distance

pcoa_wu <- ordinate(ms_rare, method = "PCoA", distance = w_unifrac_dist)

# plot PCoA

w_unifrac_plot <- plot_ordination(ms_rare, pcoa_wu, color = "treatment_status") +
  labs(col = "Treatment Status",
       title = "PCoA for Treatment Status")

plot(w_unifrac_plot)
ggsave("w_unifrac_plot.png", plot = w_unifrac_plot, width = 6, height = 4)

# Statistical tests:

# Alpha diversity (Shannon) using Kruskal-Wallis test:

meta_data <- data.frame(sample_data(ms_rare))
shannon_data <- cbind(meta_data, shannon_div$Shannon)

kruskal.test(shannon_div$Shannon ~ treatment_status, data = shannon_data)

# p-value = 0.7057 > 0.05, the difference in alpha diversity is not significant in individual treatment statuses.

# Beta diversity (weighted unifrac) using PERMANOVA:

permanova_result <- adonis2(w_unifrac_dist ~ meta_data$treatment_status)
permanova_result

# p-value = 0.002 < 0.05, the difference in beta diversity is significant between treatment statuses.
# However, the R2 value is 0.01868, indicating that only a 1.87% of the variation can be attributed to treatment status.

# testing git