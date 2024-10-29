## DESEQ CODE STARTS AT LINE 71 ##
# Lines 4-70: metadata trimming and creating phyloseq object

# load in packages:
library(tidyverse)
library(phyloseq)
library(ape)
library(DESeq2)


# Load in the ms metadata,  feature table, taxonomy file, and rooted tree
otufp <- "ms-mit-chlor-freq-filtered-table.txt"
otu <- read_delim(file =otufp, delim="\t", skip=1)

metafp <- "corrected_ms_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

taxfp <- "ms-taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "ms-tree.nwk"
phylotree <- read.tree(phylotreefp)


### using Cyrus' trimming code
# transpose the otu object so the rows become columns and vice versa
otu_transposed <- as.data.frame(t(otu))
otu_transposed <- cbind("sample-id" = rownames(otu_transposed), otu_transposed)
colnames(otu_transposed) <- otu_transposed[1,]
colnames(otu_transposed)[1] <- "sample-id"
otu_transposed <- otu_transposed[-1,]


# right_join the metadata table and the otu table to filter out unwanted metadata columns:
combined <- right_join(meta, otu_transposed)


# separate the dataframes:
metadata_trimmed <- combined[, 1:60]


### using Natalia's phyloseq code with trimmed metatdata
# Adjust files to be read into a phyloseq object. Make the phyloseq object.
# feature-table;
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$"#OTU ID"
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)


# metadata:
meta_df <- as.data.frame(metadata_trimmed[,-1])
rownames(meta_df) <- metadata_trimmed$"sample-id"
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
save(ms_phyloseq, file = "ms_phyloseq.RData")


### Dima's DESeq2 code
# loading in non-rarefied phyloseq object
load("ms_export/ms_phyloseq.RData")


#converting zeros to ones to avoid zeros error and converting phyloseq object to DESeq2 object
ms_plus1 <- transform_sample_counts(ms_phyloseq, function(x) x+1)
ms_deseq <- phyloseq_to_deseq2(ms_plus1, ~`disease_course`)
DESEQ_ms <- DESeq(ms_deseq)

#results RRMS vs PPMS (RRMS is set as reference)
RRMS_vs_PPMS_res <- results(DESEQ_ms, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("disease_course","PPMS","RRMS"))

#results RRMS vs SPMS (RRMS is set as reference)
RRMS_vs_SPMS_res <- results(DESEQ_ms, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("disease_course","SPMS","RRMS"))

#results PPMS vs SPMS (SPMS is set as reference)
PPMS_vs_SPMS_res <- results(DESEQ_ms, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("disease_course","PPMS","SPMS"))


# Volcano plots: effect size vs significance (see what is set as reference above)
RRMS_vs_PPMS_vol_plot <- RRMS_vs_PPMS_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

RRMS_vs_SPMS_vol_plot <- RRMS_vs_SPMS_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

PPMS_vs_SPMS_vol_plot <- PPMS_vs_SPMS_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))


# Preparing data tables (*sigASVs) for the bar plots
RRMS_vs_PPMS_sigASVs <- RRMS_vs_PPMS_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(RRMS_vs_PPMS_sigASVs)

RRMS_vs_SPMS_sigASVs <- RRMS_vs_SPMS_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(RRMS_vs_SPMS_sigASVs)

PPMS_vs_SPMS_sigASVs <- PPMS_vs_SPMS_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(PPMS_vs_SPMS_sigASVs)


# Retrieving ASV names
RRMS_vs_PPMS_sigASVs_vec <- RRMS_vs_PPMS_sigASVs %>%
  pull(ASV)

RRMS_vs_SPMS_sigASVs_vec <- RRMS_vs_SPMS_sigASVs %>%
  pull(ASV)

PPMS_vs_SPMS_sigASVs_vec <- PPMS_vs_SPMS_sigASVs %>%
  pull(ASV)


# Pruning phyloseq file for significant ASVs (*_vs_*_sigASVs variables are rewritten)
RRMS_vs_PPMS_DESeq <- prune_taxa(RRMS_vs_PPMS_sigASVs_vec,ms_phyloseq)
RRMS_vs_PPMS_sigASVs <- tax_table(RRMS_vs_PPMS_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(RRMS_vs_PPMS_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

RRMS_vs_SPMS_DESeq <- prune_taxa(RRMS_vs_SPMS_sigASVs_vec,ms_phyloseq)
RRMS_vs_SPMS_sigASVs <- tax_table(RRMS_vs_SPMS_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(RRMS_vs_SPMS_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

PPMS_vs_SPMS_DESeq <- prune_taxa(PPMS_vs_SPMS_sigASVs_vec,ms_phyloseq)
PPMS_vs_SPMS_sigASVs <- tax_table(PPMS_vs_SPMS_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(PPMS_vs_SPMS_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


# Bar plots
RRMS_vs_PPMS_bar_plot <- ggplot(RRMS_vs_PPMS_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

RRMS_vs_SPMS_bar_plot <- ggplot(RRMS_vs_SPMS_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

PPMS_vs_SPMS_bar_plot <- ggplot(PPMS_vs_SPMS_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


# saving volcano plots to the local computer
ggsave(filename="RRMS_vs_PPMS_vol_plot.png",RRMS_vs_PPMS_vol_plot)
ggsave(filename="RRMS_vs_SPMS_vol_plot.png",RRMS_vs_SPMS_vol_plot)
ggsave(filename="PPMS_vs_SPMS_vol_plot.png",PPMS_vs_SPMS_vol_plot)


# saving bar plots to the local computer
ggsave(filename="RRMS_vs_PPMS_bar_plot.png",RRMS_vs_PPMS_bar_plot)
ggsave(filename="RRMS_vs_SPMS_bar_plot.png",RRMS_vs_SPMS_bar_plot)
ggsave(filename="PPMS_vs_SPMS_bar_plot.png",PPMS_vs_SPMS_bar_plot)
