# load in packages:
library(tidyverse)
library(phyloseq)
library(DESeq2)


# loading in non-rarefied phyloseq object containing trimmed metadata (created by using ms_trimming_metada_script.R and ms_phyloseq_object.R)
load("ms_export/ms_phyloseq.RData")


#converting zeros to ones to avoid zeros error and converting phyloseq object to DESeq2 object
ms_plus1 <- transform_sample_counts(ms_phyloseq, function(x) x+1)
ms_deseq <- phyloseq_to_deseq2(ms_plus1, ~`disease_course`)
DESEQ_ms <- DESeq(ms_deseq)

#results RRMS vs PPMS 
RRMS_vs_PPMS_res <- results(DESEQ_ms, tidy=TRUE, 
               #(RRMS is set as reference)
               contrast = c("disease_course","PPMS","RRMS"))

#results RRMS vs SPMS 
RRMS_vs_SPMS_res <- results(DESEQ_ms, tidy=TRUE, 
               #(RRMS is set as reference)
               contrast = c("disease_course","SPMS","RRMS"))

#results PPMS vs SPMS 
PPMS_vs_SPMS_res <- results(DESEQ_ms, tidy=TRUE, 
               #(SPMS is set as reference)
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
