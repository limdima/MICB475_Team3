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

# transpose the otu object so the rows become columns and vice versa:

otu_transposed <- as.data.frame(t(otu))
otu_transposed <- cbind("sample-id" = rownames(otu_transposed), otu_transposed)
colnames(otu_transposed) <- otu_transposed[1,]
colnames(otu_transposed)[1] <- "sample-id"
otu_transposed <- otu_transposed[-1,]

# right_join the metadata table and the otu table to filter out unwanted metadata columns:

combined <- right_join(meta, otu_transposed)

# separate the dataframes:

metadata_trimmed <- combined[, 1:60]

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

# Statistical tests:

