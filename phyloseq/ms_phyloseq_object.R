# Create phyloseq object

# Load in packages:
library(phyloseq)
library(ape)
library(tidyverse)

# Load in the ms metadata,  feature table, taxonomy file, and rooted tree
metafp <- "phyloseq/ms_export/corrected_ms_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "phyloseq/ms_export/ms-mit-chlor-freq-filtered-table.txt"
otu <- read_delim(file =otufp, delim="\t", skip=1)

taxfp <- "phyloseq/ms_export/ms-taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "phyloseq/ms_export/ms-tree.nwk"
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
ms_phyloseq <- phyloseq(OTU, META, TAX, phylotree)
save(ms_phyloseq, file = "phyloseq/ms_phyloseq.RData")


# rarefy samples (for diversity analyses)
ms_phyloseq_rare <- rarefy_even_depth(ms_phyloseq, rngseed = 4222, sample.size = 9488)
save(ms_phyloseq_rare, file = "phyloseq/ms_phyloseq_rare.RData")
