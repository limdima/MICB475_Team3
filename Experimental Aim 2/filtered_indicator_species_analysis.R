# Conducting Indicator Species Analysis

# Load packages for Indicator Species
library(ape)
library(tidyverse)
library(indicspecies)

# Load non-rarefied phyloseq object
phyloseq <- load("final_filtered_ms_phyloseq.RData")

# Run indicator species analysis using the disease_course variable
isa_output <- multipatt(t(otu_table(final_filtered_ms_phyloseq)), cluster = sample_data(final_filtered_ms_phyloseq)$'disease_course')

# Look at results
summary(isa_output)

# Extract taxonomy table
taxtable <- tax_table(ms_phyloseq) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
res <- isa_output$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)

# Save results
save(res, file = "res.RData")

# View results
View(res)