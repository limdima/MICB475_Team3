# load in packages:
library(tidyverse)
library(phyloseq)
library(DESeq2)


# loading in non-rarefied phyloseq object
load("ms_phyloseq.RData")


#converting zeros to ones to avoid zeros error and converting phyloseq object to DESeq2 object
ms_plus1 <- transform_sample_counts(ms_phyloseq, function(x) x+1)
ms_deseq <- phyloseq_to_deseq2(ms_plus1, ~`disease_course`)
DESEQ_ms <- DESeq(ms_deseq)

                                    
#results MS subgroup vs control
RRMS_res <- results(DESEQ_ms, tidy=TRUE, 
               #control is set as reference
               contrast = c("disease_course","PPMS","Control"))

SPMS_res <- results(DESEQ_ms, tidy=TRUE, 
                            #control is set as reference
               contrast = c("disease_course","SPMS","Control"))

PPMS_res <- results(DESEQ_ms, tidy=TRUE, 
                            #control is set as reference
               contrast = c("disease_course","PPMS","Control"))


# Volcano plots: effect size vs significance (control is set as reference)
RRMS_vol_plot <- RRMS_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

SPMS_vol_plot <- SPMS_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

PPMS_vol_plot <- PPMS_res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))


# Preparing data tables (*sigASVs) for the bar plots
RRMS_sigASVs <- RRMS_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(RRMS_sigASVs)

SPMS_sigASVs <- SPMS_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(SPMS_sigASVs)

PPMS_sigASVs <- PPMS_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(PPMS_sigASVs)


# Retrieving ASV names
RRMS_sigASVs_vec <- RRMS_sigASVs %>%
  pull(ASV)

SPMS_sigASVs_vec <- SPMS_sigASVs %>%
  pull(ASV)

PPMS_sigASVs_vec <- PPMS_sigASVs %>%
  pull(ASV)


# Pruning phyloseq file for significant ASVs (*_vs_*_sigASVs variables are rewritten)
RRMS_DESeq <- prune_taxa(RRMS_sigASVs_vec,ms_phyloseq)
RRMS_sigASVs <- tax_table(RRMS_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(RRMS_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

SPMS_DESeq <- prune_taxa(SPMS_sigASVs_vec,ms_phyloseq)
SPMS_sigASVs <- tax_table(SPMS_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(SPMS_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

PPMS_DESeq <- prune_taxa(PPMS_sigASVs_vec,ms_phyloseq)
PPMS_sigASVs <- tax_table(PPMS_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(PPMS_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


# Bar plots
RRMS_bar_plot <- ggplot(RRMS_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

SPMS_bar_plot <- ggplot(SPMS_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

PPMS_bar_plot <- ggplot(PPMS_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


# saving volcano plots to the local computer
ggsave(filename="RRMS_vol_plot.png",RRMS_vol_plot)
ggsave(filename="SPMS_vol_plot.png",SPMS_vol_plot)
ggsave(filename="PPMS_vol_plot.png",PPMS_vol_plot)


# saving bar plots to the local computer
ggsave(filename="RRMS_bar_plot.png",RRMS_bar_plot)
ggsave(filename="SPMS_bar_plot.png",SPMS_bar_plot)
ggsave(filename="PPMS_bar_plot.png",PPMS_bar_plot)
