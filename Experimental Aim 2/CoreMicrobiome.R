# load packages for Core Microbiome 
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(RColorBrewer)

# load the phyloseq object
load("phyloseq/final_filtered_ms_phyloseq.RData")

# Convert object to relative abundance
phyloseq_rel <- transform_sample_counts(final_filtered_ms_phyloseq, fun = function(x) x/sum(x))


### Subset dataset into treatment and control groups
RRMS_stat <- subset_samples(phyloseq_rel, disease_course == "RRMS")
SPMS_stat <- subset_samples(phyloseq_rel, disease_course == "SPMS")
PPMS_stat<- subset_samples(phyloseq_rel, disease_course == "SPMS")
control_stat <- subset_samples(phyloseq_rel, disease_course == "Control")

### Set a prevalence threshold and abundance threshold
RRMS_ASVs <- core_members(RRMS_stat, detection=0.0001, prevalence = 0.25)
SPMS_ASVs <- core_members(SPMS_stat, detection=0.0001, prevalence = 0.25)
PPMS_ASVs <- core_members(PPMS_stat, detection=0.0001, prevalence = 0.25)
control_ASVs <- core_members(control_stat, detection=0.0001, prevalence = 0.25)
# Exclude very rare ASVs (detection = 0.001% relative abundance)
# Find ASVs present in 33% of samples (prevalence = 0.25), to account for outliers/ASVs unique to only 1 person

### Make a Venn-diagram
venn_pd <- ggVennDiagram(x=list(RRMS = RRMS_ASVs, SPMS = SPMS_ASVs, PPMS = PPMS_ASVs, Control = control_ASVs)) +
  scale_x_continuous(expand = expansion(mult = .2)) # this line helps fit the labels into the plot border

venn_pd

### save venn diagram
ggsave("Experimental Aim 2/venn_ms_groups.png", plot = venn_pd)


# As a preliminary check, check what taxa we have by a heat map
# (not 100% sure how the code works but copied from core microbiome tutorial websites)
# https://microbiome.github.io/tutorials/CoremicrobiotaAmplicon.html

prevalences <- seq(0.25, 1, .1) # set prevalence threshold from 0.33 to 1, in increments of 0.1  
detections <- seq(0.0001, 0.05, 0.002) 
# detection threshold from 0.0001 to get rid of low-abundance ASVs, 
# up to 0.05 (to just view what is rare/common), in increments of 0.002

phyloseq_rel.gen <- aggregate_taxa(phyloseq_rel, "Genus") # changes the ASVs to Genus names

p.core <- plot_core(phyloseq_rel.gen, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .33) + 
  xlab("Detection Threshold (Relative Abundance (%))")

print(p.core) 

