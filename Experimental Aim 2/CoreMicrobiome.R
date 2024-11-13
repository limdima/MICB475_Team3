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
ggsave("Experimental Aim 2/Visualizations/venn_ms_groups.png", plot = venn_pd)






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





### Check intersection between groups of interest (i.e. not shared with control)

only_R <- setdiff(RRMS_ASVs, union(SPMS_ASVs, union(PPMS_ASVs, control_ASVs))) # the 7 ASVs only in RRMS

intersection_RSP <- intersect(RRMS_ASVs, union(PPMS_ASVs, SPMS_ASVs)) %>%
  setdiff(control_ASVs) # the 5 ASVs between RRMS, SPMS, PPMS

intersection_SP <- intersect(PPMS_ASVs, SPMS_ASVs) %>% 
  setdiff(union(control_ASVs, RRMS_ASVs))   # The 21 ASVs found between SPMS and PPMS

only_C <- setdiff(control_ASVs, union(SPMS_ASVs, union(PPMS_ASVs, RRMS_ASVs))) # the 12 ASVs only in Control

### Use the above intersection to prune the phyloseq object
# For RRMS only
subset_RRMS <- prune_taxa(only_R, phyloseq_rel) # Keeps the 7 ASVs that are unique to RRMS samples
taxa_RRMS <- as.data.frame(tax_table(subset_RRMS)) # take taxa names and Abundance values (otu table) into dataframes
abund_RRMS <- as.data.frame(otu_table(subset_RRMS))

taxa_RRMS$ASV <- rownames(taxa_RRMS)   # set the row names to an actual column 
abund_RRMS$ASV <- rownames(abund_RRMS)

merge_RRMS_df <- merge(taxa_RRMS, abund_RRMS, by = "ASV")  # merge the dataframes together
merge_RRMS_tidy <- merge_RRMS_df %>% gather(key = "SampleID", value = "RelAbundance", 
                                            - ASV, -Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  mutate(GenusSpecies = paste(Genus, Species, sep = "_"))


summary_RRMS <- merge_RRMS_tidy %>%
  group_by(GenusSpecies) %>%
  summarize(MeanAbundance = mean(RelAbundance, na.rm = TRUE))  # Summarizes the 7 ASVs (Genus') and the mean relative abundance values


# For in RRMS, SPMS, PPMS only
subset_MS <- prune_taxa(intersection_RSP, phyloseq_rel) # Keeps the 5 ASVs that are unique to RRMS samples
taxa_MS <- as.data.frame(tax_table(subset_MS)) # take taxa names and Abundance values (otu table) into dataframes
abund_MS <- as.data.frame(otu_table(subset_MS))

taxa_MS$ASV <- rownames(taxa_MS)   # set the row names to an actual column 
abund_MS$ASV <- rownames(abund_MS)

merge_MS_df <- merge(taxa_MS, abund_MS, by = "ASV")  # merge the dataframes together
merge_MS_tidy <- merge_MS_df %>% gather(key = "SampleID", value = "RelAbundance", 
                                            - ASV, -Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  mutate(GenusSpecies = paste(Genus, Species, sep = "_"))

summary_MS <- merge_MS_tidy %>%
  group_by(GenusSpecies) %>%
  summarize(MeanAbundance = mean(RelAbundance, na.rm = TRUE))  
# Summarizes the 5 ASVs (Genus') and the mean relative abundance values in the MS groups only\


# For in SPMS, PPMS only (progressive forms of MS)
subset_PMS <- prune_taxa(intersection_SP, phyloseq_rel) # Keeps the 21 ASVs that are unique to RRMS samples
taxa_PMS <- as.data.frame(tax_table(subset_PMS)) # take taxa names and Abundance values (otu table) into dataframes
abund_PMS <- as.data.frame(otu_table(subset_PMS))

taxa_PMS$ASV <- rownames(taxa_PMS)   # set the row names to an actual column 
abund_PMS$ASV <- rownames(abund_PMS)

merge_PMS_df <- merge(taxa_PMS, abund_PMS, by = "ASV")  # merge the dataframes together
merge_PMS_tidy <- merge_PMS_df %>% gather(key = "SampleID", value = "RelAbundance", 
                                        - ASV, -Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  mutate(GenusSpecies = paste(Genus, Species, sep = "_"))

summary_PMS <- merge_PMS_tidy %>%
  group_by(GenusSpecies) %>%
  summarize(MeanAbundance = mean(RelAbundance, na.rm = TRUE))  
# Only 19 rows appear
# Some ASVs correspond to the same genus_species? e.g. g__Peptoniphilus_NA

### Plots for the summary tables

# ASVs in RRMS only
plot_unique_RRMS <- ggplot(summary_RRMS, aes(x = reorder(GenusSpecies, MeanAbundance), y= MeanAbundance)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to RRMS",
       x = "Species",
       y = "Average relative abundance")

ggsave("Experimental Aim 2/Visualizations/RRMS_Unique_Species.png", plot = plot_unique_RRMS)


# ASVs in all MS groups 
plot_unique_MS <- ggplot(summary_MS, aes(x = reorder(GenusSpecies, MeanAbundance), y= MeanAbundance)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to MS groups",
       x = "Species",
       y = "Average relative abundance")

ggsave("Experimental Aim 2/Visualizations/MS_Unique_Species.png", plot = plot_unique_MS)


plot_unique_PMS <- ggplot(summary_PMS, aes(x = reorder(GenusSpecies, MeanAbundance), y= MeanAbundance)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to Progressive MS",
       x = "Species",
       y = "Average relative abundance")

ggsave("Experimental Aim 2/Visualizations/PMS_Unique_Species.png", plot = plot_unique_PMS)

