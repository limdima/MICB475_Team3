# load in packages:

library(tidyverse)
library(phyloseq)
library(ape)

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

# metadata and otu table both now have 515 samples each after filtering for the variables we will include for analysis

save(metadata_trimmed, file = "metadata_trimmed")
save(otu_trimmed, file = "otu_trimmed")
