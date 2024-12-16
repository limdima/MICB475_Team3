# Alpha and Beta Analysis
# Note, final_filtered_ms_phyloseq and final_filtered_phyloseq_rare are available on Github, 

# libraries
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)

# load in phyloseq object
load("phyloseq_objects/final_filtered_ms_phyloseq.RData")

# Rarefy phyloseq object
ms_phyloseq_rare <- rarefy_even_depth(final_filtered_ms_phyloseq, rngseed = 4222, sample.size = 9488)



#### FILTERING: RAREFIED PHYLOSEQ ####
# make rarefied phyloseq object into a data frame:
phyloseq_rare_df <- as(sample_data(ms_phyloseq_rare), "data.frame")

# filter for control (alpha/beta analysis compares MS sub-types only)
filtered_phyloseq_rare_df <- phyloseq_rare_df %>%
  mutate(sample_id = rownames(phyloseq_rare_df)) %>% # Keep a column for sample_ids
  filter(disease_course %in% c("SPMS", "RRMS", "PPMS")) %>% # Keep only specific disease courses
  ungroup()

# prune the rarefied phyloseq to only have sample ids in filtered phyloseq data frame
final_filtered_phyloseq_rare <- prune_samples(filtered_phyloseq_rare_df$sample_id, ms_phyloseq_rare)

# save filtered phyloseq rarefied
save(final_filtered_phyloseq_rare, file = "../phyloseq/final_filtered_phyloseq_rare.RData")
