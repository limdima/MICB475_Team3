# Alpha and Beta Analysis

# libraries
library(phyloseq)
library(ggplot2)

# load in phyloseq object
ms_ps <- load("ms_phyloseq.RData")


# alpha (Shannon) diversity for all samples
shannon_div <- estimate_richness(ms_phyloseq, measures = "Shannon")


shannon_div$treatment_status <- sample_data(ms_phyloseq)$treatment_status


shan_vs_treatment <- ggplot(shann0n, aes(x=treatment_status, y=Shannon)) +
  geom_boxplot() +
  labs(title = "Shannon Diversity vs Treatment Status", 
       x = "Disease Type",
       y = "Shannon Diversity Index") +
  theme_minimal()
