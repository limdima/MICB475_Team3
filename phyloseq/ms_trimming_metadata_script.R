# load in packages:

library(tidyverse)
library(phyloseq)
library(ape)

# load in the sequence table and the metadata table:

otufp <- "ms_export/ms-mit-chlor-freq-filtered-table.txt"
otu <- read_delim(file =otufp, delim="\t", skip=1)

metafp <- "ms_export/corrected_ms_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

# transpose the otu object so the rows become columns and vice versa

otu_transposed <- as.data.frame(t(otu))
colnames(otu_transposed) <- rownames(otu)
otu_transposed <- cbind("sample-id" = rownames(otu_transposed), otu_transposed)
colnames(otu_transposed) <- otu_transposed[1,]
colnames(otu_transposed)[colnames(otu_transposed) == "#OTU ID"] <- "sample-id"
otu_transposed_test <- otu_transposed[-1,]

# will simplify this code later


# right_join the metadata table and the otu table to filter out unwanted metadata columns:

combined <- right_join(meta, otu_transposed)

# separate the dataframes:

metadata_trimmed <- combined[, 1:60]
