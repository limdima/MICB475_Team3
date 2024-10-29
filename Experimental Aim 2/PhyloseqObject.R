# Load packages for making the Phyloseq object

library(phyloseq)
library(tidyverse)

metafp <- "meta_ms.csv"
meta <- read_delim(metafp, delim=",")

otufp <- "ms-table.txt"
otu <- read_delim(file =otufp, delim="\t", skip=1)

taxfp <- "ms-taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "ms-tree.nwk"
phylotree <- read.tree(phylotreefp)

otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$"#OTU ID"
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

meta_df <- as.data.frame(meta[,-1])
rownames(meta_df) <- meta$"sample-id"
META <- sample_data(meta_df)

tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(col=Taxon, sep="; ", into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$"Feature ID"
TAX <- tax_table(tax_mat)

# Make phyloseq object
ms_phyloseq <- phyloseq(OTU, META, TAX, phylotree)

save(ms_phyloseq, file="ms_phyloseq.RData")