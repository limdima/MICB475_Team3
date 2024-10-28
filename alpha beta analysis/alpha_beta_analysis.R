# Alpha and Beta Analysis

# libraries
library(phyloseq)
library(ggplot2)

# load in phyloseq object
load("ms_phyloseq_rare.RData")

# ALPHA DIVERSITY ANALYSIS
# alpha (Shannon) diversity for all samples
shannon_div <- estimate_richness(ms_phyloseq_rare, measures = "Shannon")

# add a column for treatment status
shannon_div$treatment_status <- sample_data(ms_phyloseq_rare)$treatment_status

# create a box-plot to visualize diversity metrics
shan_vs_treatment <- ggplot(shannon_div, aes(x=treatment_status, y=Shannon)) +
  geom_boxplot() +
  labs(title = "Shannon Diversity vs Treatment Status", 
       x = "Treatment Status",
       y = "Shannon Diversity Index") +
  theme_minimal()

plot(shan_vs_treatment)
ggsave("shannon_div.png", plot = shan_vs_treatment, width = 6, height = 4)

# BETA DIVERSITY ANALYSIS
# beta (weighted unifrac) diversity for all samples
w_unifrac_dist <- distance(ms_phyloseq_rare, method = "wunifrac")

# PCoA on weighted unifrac distance
pcoa_wu <- ordinate(ms_phyloseq_rare, method = "PCoA", distance = w_unifrac_dist)

# plot PCoA
w_unifrace_plot <- plot_ordination(ms_phyloseq_rare, pcoa_wu, color = "treatment_status") +
  labs(col = "Treatment Status",
       title = "PCoA for Treatment Status")

plot(w_unifrace_plot)
ggsave("w_unifrac_plot.png", plot = w_unifrac_plot, width = 6, height = 4)

