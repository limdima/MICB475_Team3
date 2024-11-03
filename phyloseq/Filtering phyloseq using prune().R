# Load in packages:
library(phyloseq)
library(ape)
library(tidyverse)

# Load in the ms metadata,  feature table, taxonomy file, and rooted tree
metafp <- "corrected_ms_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "ms-mit-chlor-freq-filtered-table.txt"
otu <- read_delim(file =otufp, delim="\t", skip=1)

taxfp <- "ms-taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "ms-tree.nwk"
phylotree <- read.tree(phylotreefp)

# Adjust files to be read into a phyloseq object. Make the phyloseq object.
# feature-table;
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$"#OTU ID"
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# metadata:
meta_df <- as.data.frame(meta[,-1])
rownames(meta_df) <- meta$"sample-id"
META <- sample_data(meta_df)

# taxonomy file:
tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(col=Taxon, sep="; ", into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$"Feature ID"
TAX <- tax_table(tax_mat)

# create and save phyloseq object:
unfiltered_ms_phyloseq <- phyloseq(OTU, META, TAX, phylotree)

# filter phyloseq object (probiotics and special diets), filter for two samples per household
filtered_metadata <- meta %>%
  filter(probiotics == 0) %>%        
  filter(diet_no_special_needs == 1) %>%
  filter(eating_disorder == 0) %>%
  
  group_by(`household`) %>% 
  filter(n() == 2) %>% 
  ungroup()

# Update the phyloseq object, keeping only samples in filtered metadata
# Filters to 515 total samples
final_filtered_ms_phyloseq <- 
  prune_samples(rownames(sample_data(unfiltered_ms_phyloseq)) %in% filtered_metadata$`sample-id`, unfiltered_ms_phyloseq)

# save phyloseq object:
save(final_filtered_ms_phyloseq, file = "final_filtered_ms_phyloseq.RData")

