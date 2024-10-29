# load in packages:

library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)

# load in the sequence table and the metadata table:

otufp <- "ms_export/ms-mit-chlor-freq-filtered-table.txt"
otu <- read_delim(file =otufp, delim="\t", skip=1)

metafp <- "ms_export/corrected_ms_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

# transpose the otu object so the rows become columns and vice versa - this is done so it can be semi_joined against meta later on

otu_transposed <- as.data.frame(t(otu[, -1]))    # transposes the data but excludes the first column as that will be col names
colnames(otu_transposed) <- otu[[1]]     # first column in the original otu table is for the OTU IDs, now as the new col names for the transposed table
otu_transposed <- cbind("sample-id" = rownames(otu_transposed), otu_transposed)   # makes a new column "sample-id" with the original header names from otu, for downstream joining
rownames(otu_transposed) <- NULL   # removes the indexing row names


# filter the metadata to get rid of variables that we deemed as potential confounders
# We think mainly probiotics, diet, and eating disorders will contribute a larger confounding effect on the gut microbiome

filter_meta <- meta %>%
  filter(probiotics == 0) %>%        
  filter(diet_no_special_needs == 1) %>%
  filter(eating_disorder == 0) %>%
  
  group_by(`household`) %>% # samples are paired by household, one MS and one control
  filter(n() == 2) %>% # filters out any rows that do not contain two rows per household
  ungroup()

# after filtering there is less samples in the metadata (522) than OTUs (915)
# semi_join the metadata table and the otu table to filter out samples that we won't be including:
# semi_join will keep only the rows in meta with a match in the otu table, without actually merging the two dataframes

metadata_trimmed <- semi_join(filter_meta, otu_transposed, by = 'sample-id')
# returns 515 samples, meaning some metadata without corresponding OTUs were also removed

otu_filter <- semi_join(otu_transposed, filter_meta, by = 'sample-id')
# returns the same corresponding 515 samples as in the metadata


# re-transpose the otu table so that it is ready to be put into a phyloseq object, similar to before
otu_trimmed <- as.data.frame(t(otu_filter))
colnames(otu_trimmed) <- otu_filter[[1]]
otu_trimmed <- cbind("#OTU ID" = rownames(otu_trimmed), otu_trimmed)   # makes a new column "#OTU ID" as in the original OTU table
rownames(otu_trimmed) <- NULL

# Load in the rest of the files needed for the phyloseq object:

taxfp <- "ms_export/ms-taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "ms_export/ms-tree.nwk"
phylotree <- read.tree(phylotreefp)

# Adjust files for a phyloseq readable object, then make the phyloseq object:

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

# phyloseq object:

ms_phyloseq <- phyloseq(OTU, META, TAX, phylotree)

# Rarefy the phyloseq object. Sample depth is chosen by analysis done in qiime2.

ms_rare <- rarefy_even_depth(ms_phyloseq, rngseed = 1, sample.size = 9488)

# ALPHA DIVERSITY ANALYSIS
# alpha (Shannon) diversity for all samples

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

# BETA DIVERSITY ANALYSIS
# beta (weighted unifrac) diversity for all samples

w_unifrac_dist <- UniFrac(ms_rare, weighted=TRUE)

# PCoA on weighted unifrac distance

pcoa_wu <- ordinate(ms_rare, method = "PCoA", distance = w_unifrac_dist)

# plot PCoA

w_unifrace_plot <- plot_ordination(ms_rare, pcoa_wu, color = "treatment_status") +
  labs(col = "Treatment Status",
       title = "PCoA for Treatment Status")

plot(w_unifrace_plot)
ggsave("w_unifrac_plot.png", plot = w_unifrac_plot, width = 6, height = 4)

# Statistical tests:

