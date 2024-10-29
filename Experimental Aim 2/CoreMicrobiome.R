# load packages for Core Microbiome 
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# load the phyloseq object
load("phyloseq/ms_phyloseq.RData")

# Convert object to relative abundance
phyloseq_rel <- transform_sample_counts(ms_phyloseq, fun = function(x) x/sum(x))


### Subset dataset into treatment and control groups
RRMS_stat <- subset_samples(phyloseq_rel, disease_course == "RRMS")
SPMS_stat <- subset_samples(phyloseq_rel, disease_course == "SPMS")
PPMS_stat<- subset_samples(phyloseq_rel, disease_course == "SPMS")

### Set a prevalence threshold and abundance threshold
RRMS_ASVs <- core_members(RRMS_stat, detection=0, prevalence = 0)
SPMS_ASVs <- core_members(SPMS_stat, detection=0, prevalence = 0)
PPMS_ASVs <- core_members(PPMS_stat, detection=0, prevalence = 0)
# Keep all ASVs including rare ASVs (detection = 0% relative abundance)
# Find ASVs present in any samples (prevalence = 0), trying to tease out what is unique to each group

### Make a Venn-diagram
venn_pd <- ggVennDiagram(x=list(RRMS = RRMS_ASVs, SPMS = SPMS_ASVs, PPMS = PPMS_ASVs)) +
  scale_x_continuous(expand = expansion(mult = .2)) # this line helps fit the labels into the plot border

venn_pd

### save venn diagram
ggsave("Experimental Aim 2/venn_ms_groups.png", plot = venn_pd)

# check what taxa we have
# watch out the following can take a long time and the graph needs to be modified to fit everything
#prune_taxa(PPMS_ASVs,ms_phyloseq) %>%
#  plot_bar(fill = "Genus") +
#  facet_wrap(.~ disease_course, scales = "free")

