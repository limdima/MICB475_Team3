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
PPMS_stat<- subset_samples(phyloseq_rel, disease_course == "PPMS")
control_stat <- subset_samples(phyloseq_rel, disease_course == "Control")

### Set a prevalence threshold and abundance threshold
RRMS_ASVs <- core_members(RRMS_stat, detection=0.0001, prevalence = 0.25)
SPMS_ASVs <- core_members(SPMS_stat, detection=0.0001, prevalence = 0.25)
PPMS_ASVs <- core_members(PPMS_stat, detection=0.0001, prevalence = 0.25)
control_ASVs <- core_members(control_stat, detection=0.0001, prevalence = 0.25)
# Exclude very rare ASVs (detection = 0.01% relative abundance)
# Find ASVs present in 33% of samples (prevalence = 0.25), to account for outliers/ASVs unique to only 1 person

### Make a Venn-diagram
venn_pd <- ggVennDiagram(x=list(RRMS = RRMS_ASVs, SPMS = SPMS_ASVs, PPMS = PPMS_ASVs, Control = control_ASVs),
                         label_size = 3.5,
                         label_alpha = 0) +
  scale_x_continuous(expand = expansion(mult = .2)) # this line helps fit the labels into the plot border

# Manually change colors here:
venn_pd$layers[[1]]$mapping <- aes(fill = name)

zero <- "white"
onefive <- "#7ce8ff"
sixten <- "#55d0ff"
eleventwenty <- "#0091bc"
thirtyplus <- "#00609a"
sharedcontrol <- "#b8b8b8"

venn_pd_test <- venn_pd + scale_fill_manual(values = c('PPMS/Control' = sharedcontrol, 
                                       'PPMS'= eleventwenty, 
                                       'SPMS' = eleventwenty, 
                                       'SPMS/Control' = sharedcontrol,
                                       'RRMS/PPMS' = onefive,
                                       'RRMS/SPMS' = onefive,
                                       'RRMS/SPMS/Control' = sharedcontrol,
                                       'RRMS/PPMS/Control' = "#858d90",
                                       'RRMS' = oneten,
                                       'RRMS/SPMS/PPMS' = oneten,
                                       'SPMS/PPMS/Control' = sharedcontrol,
                                       'Control' = sharedcontrol,
                                       'SPMS/PPMS' = sixten,
                                       'RRMS/Control' = "#858d90",
                                       'RRMS/SPMS/PPMS/Control' = "#4c4c4c")) +
  theme(legend.position = "none")
venn_pd_test


### save venn diagram
ggsave("Experimental Aim 2/Visualizations/venn_ms_groups.png", plot = venn_pd_test)




### Check intersection between groups of interest (i.e. not shared with control)

only_R <- setdiff(RRMS_ASVs, union(SPMS_ASVs, union(PPMS_ASVs, control_ASVs))) # the 3 ASVs only in RRMS

only_S <- setdiff(SPMS_ASVs, union(RRMS_ASVs, union(PPMS_ASVs, control_ASVs))) # the 13 ASVs only in SPMS

only_P <- setdiff(PPMS_ASVs, union(RRMS_ASVs, union(SPMS_ASVs, control_ASVs))) # the 19 ASVs only in PPMS

intersection_RSP <- intersect(RRMS_ASVs, PPMS_ASVs) %>%
  intersect(SPMS_ASVs) %>%
  setdiff(control_ASVs) # the 3 ASVs between RRMS, SPMS, PPMS

intersection_SP <- intersect(PPMS_ASVs, SPMS_ASVs) %>% 
  setdiff(union(control_ASVs, RRMS_ASVs))   # The 8 ASVs found between SPMS and PPMS (PMS)



only_C <- setdiff(control_ASVs, union(SPMS_ASVs, union(PPMS_ASVs, RRMS_ASVs))) # the 12 ASVs only in Control

### Use the above intersection to prune the phyloseq object
# For RRMS only
subset_RRMS <- prune_taxa(only_R, phyloseq_rel) # Keeps the 3 ASVs that are unique to RRMS samples
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
  summarize(MeanAbundancePct = mean(RelAbundance, na.rm = TRUE)*100)  # Convert the relative abundance values to percentage
# Summarizes the 3 ASVs (Genus') and the mean relative abundance values


# For SPMS only
subset_SPMS <- prune_taxa(only_S, phyloseq_rel) # Keeps the 13 ASVs that are unique to SPMS samples
taxa_SPMS <- as.data.frame(tax_table(subset_SPMS)) # take taxa names and Abundance values (otu table) into dataframes
abund_SPMS <- as.data.frame(otu_table(subset_SPMS))

taxa_SPMS$ASV <- rownames(taxa_SPMS)   # set the row names to an actual column 
abund_SPMS$ASV <- rownames(abund_SPMS)

merge_SPMS_df <- merge(taxa_SPMS, abund_SPMS, by = "ASV")  # merge the dataframes together
merge_SPMS_tidy <- merge_SPMS_df %>% gather(key = "SampleID", value = "RelAbundance", 
                                          - ASV, -Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  mutate(GenusSpecies = paste(Genus, Species, sep = "_"))

summary_SPMS <- merge_SPMS_tidy %>%
  group_by(GenusSpecies) %>%
  summarize(MeanAbundancePct = mean(RelAbundance, na.rm = TRUE)*100)  # Convert the relative abundance values to percentage
# Summarizes the 13 ASVs (Genus') and the mean relative abundance values

# For PPMS only
subset_PPMS <- prune_taxa(only_P, phyloseq_rel) # Keeps the 19 ASVs that are unique to SPMS samples
taxa_PPMS <- as.data.frame(tax_table(subset_PPMS)) # take taxa names and Abundance values (otu table) into dataframes
abund_PPMS <- as.data.frame(otu_table(subset_PPMS))

taxa_PPMS$ASV <- rownames(taxa_PPMS)   # set the row names to an actual column 
abund_PPMS$ASV <- rownames(abund_PPMS)

merge_PPMS_df <- merge(taxa_PPMS, abund_PPMS, by = "ASV")  # merge the dataframes together
merge_PPMS_tidy <- merge_PPMS_df %>% gather(key = "SampleID", value = "RelAbundance", 
                                            - ASV, -Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  mutate(GenusSpecies = paste(Family, Genus, Species, sep = "_"))

summary_PPMS<- merge_PPMS_tidy %>%
  group_by(GenusSpecies) %>%
  summarize(MeanAbundancePct = mean(RelAbundance, na.rm = TRUE)*100)
# Summarizes the 17 ASVs (Genus') and the mean relative abundance values
# Note that there are 19 ASVs, but 3 of them correspond to g__Christensenellaceae_R-7_group


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
  summarize(MeanAbundancePct = mean(RelAbundance, na.rm = TRUE)*100)   # Convert the relative abundance values to percentage
# Summarizes the 5 ASVs (Genus') and the mean relative abundance values in the MS groups only\

# For in SPMS, PPMS progressive MS
subset_PMS <- prune_taxa(intersection_SP, phyloseq_rel) # Keeps the 5 ASVs that are unique to RRMS samples
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
  summarize(MeanAbundancePct = mean(RelAbundance, na.rm = TRUE)*100)   # Convert the relative abundance values to percentage
# Summarizes the 5 ASVs (Genus') and the mean relative abundance values in the MS groups only\

# For control group only
subset_C <- prune_taxa(only_C, phyloseq_rel) # Keeps the 21 ASVs that are unique to RRMS samples
taxa_C <- as.data.frame(tax_table(subset_C)) # take taxa names and Abundance values (otu table) into dataframes
abund_C <- as.data.frame(otu_table(subset_C))

taxa_C$ASV <- rownames(taxa_C)   # set the row names to an actual column 
abund_C$ASV <- rownames(abund_C)

merge_C_df <- merge(taxa_C, abund_C, by = "ASV")  # merge the dataframes together
merge_C_tidy <- merge_C_df %>% gather(key = "SampleID", value = "RelAbundance", 
                                      - ASV, -Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
  mutate(GenusSpecies = paste(Genus, Species, sep = "_"))

summary_C <- merge_C_tidy %>%
  group_by(GenusSpecies) %>%
  summarize(MeanAbundancePct = mean(RelAbundance, na.rm = TRUE)*100)  # Convert the relative abundance values to percentage
# 12 ASVs unique to control


### Plots for the summary tables

# ASVs in RRMS only
plot_unique_RRMS <- ggplot(summary_RRMS, aes(x = reorder(GenusSpecies, MeanAbundancePct), y= MeanAbundancePct)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to RRMS",
       x = "Species",
       y = "Average relative abundance (%)") +
  theme_classic()
ggsave("Experimental Aim 2/Visualizations/RRMS_Unique_Species.png", plot = plot_unique_RRMS)

# ASVs in SPMS only
plot_unique_SPMS <- ggplot(summary_SPMS, aes(x = reorder(GenusSpecies, MeanAbundancePct), y= MeanAbundancePct)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to SPMS",
       x = "Species",
       y = "Average relative abundance (%)") +
  theme_classic()

ggsave("Experimental Aim 2/Visualizations/SPMS_Unique_Species.png", plot = plot_unique_SPMS)

plot_unique_PPMS <- ggplot(summary_PPMS, aes(x = reorder(GenusSpecies, MeanAbundancePct), y= MeanAbundancePct)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to PPMS",
       x = "Species",
       y = "Average relative abundance (%)") +
  theme_classic()

ggsave("Experimental Aim 2/Visualizations/PPMS_Unique_Species.png", plot = plot_unique_PPMS)

# ASVs in all MS groups 
plot_unique_MS <- ggplot(summary_MS, aes(x = reorder(GenusSpecies, MeanAbundancePct), y= MeanAbundancePct)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to MS groups",
       x = "Species",
       y = "Average relative abundance (%)")+
  theme_classic()

ggsave("Experimental Aim 2/Visualizations/MS_Unique_Species.png", plot = plot_unique_MS)

# ASVs in progressive MS only
plot_unique_PMS <- ggplot(summary_PMS, aes(x = reorder(GenusSpecies, MeanAbundancePct), y= MeanAbundancePct)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to Progressive MS",
       x = "Species",
       y = "Average relative abundance (%)")+
  theme_classic()

ggsave("Experimental Aim 2/Visualizations/PMS_Unique_Species.png", plot = plot_unique_PMS)


# ASVs in control group only
plot_unique_control <- ggplot(summary_C, aes(x = reorder(GenusSpecies, MeanAbundancePct), y= MeanAbundancePct)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Mean Relative abundance per ASV unique to Control",
       x = "Species",
       y = "Average relative abundance (%)")+
  theme_classic()


ggsave("Experimental Aim 2/Visualizations/Control_Unique_Species.png", plot = plot_unique_control)




## Re-labelling plots and making into lollipop: ###

### RRMS LOLLIPOP ###
                                        
plot_unique_RRMS <- ggplot(summary_RRMS, aes(x = reorder(GenusSpecies, MeanAbundancePct), y = MeanAbundancePct)) +
  geom_segment(aes(xend = GenusSpecies, y = 0, yend = MeanAbundancePct), 
               color = "#d95f02", linewidth = 1.5) +  # Increase line thickness
  geom_point(color = "#d95f02", size = 5) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() + 
  labs(title = " ",
       x = "Species unique to RRMS",
       y = "Average relative abundance (%)") +
  scale_x_discrete(labels = c("g__Incertae_Sedis_s__[Clostridium]_leptum" = expression(italic("                     Clostridium leptum")),
                              "g__Family_XIII_AD3011_group_s__uncultured_organism" = "                     Unidentified",
                              "NA_NA" = "                     Unidentified"
  )) +
  theme_classic() +
  theme(
    # Change background color to lighter purple
    axis.text.y = element_text(size = 12), 
    axis.text.x = element_text(size = 10),  
    axis.title.x = element_text(size = 14),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title font size
    panel.background = element_rect(fill = "#ffebcc", color = "transparent"),
  )

### SPMS LOLLIPOP ###
                                        
plot_unique_SPMS <- ggplot(summary_SPMS, aes(x = reorder(GenusSpecies, MeanAbundancePct), y= MeanAbundancePct)) +
  geom_segment(aes(xend = GenusSpecies, y = 0, yend = MeanAbundancePct), 
               color = "#1b9e77", linewidth = 1.5) +  # Increase line thickness
  geom_point(color = "#1b9e77", size = 5) + 
  coord_flip() + 
  labs(title = " ",
       x = "Species unique to SPMS",
       y = "Average relative abundance (%)")+
  
  scale_x_discrete(labels = c("g__Tyzzerella_NA" = expression(italic("Coprococcus phoceensis")),
                              "g__Holdemanella_s__uncultured_bacterium" = expression(italic("Holdemanella porci")),
                              "g__Peptoniphilus_NA" = expression(italic("Peptoniphilus senegalensis")),
                              
                              "g__Blautia_NA" = expression(italic("Blautia stercoris")),
                              "g__Ezakiella_NA" = expression(italic("Ezakiella coagulans")),
                              "g__Porphyromonas_NA" = "Unidentified",
                              "g__S5-A14a_s__uncultured_Anaerovorax" = expression(italic("Casaltella massiliensis")),
                              "NA_NA" = expression(italic("Clostridium hylemonae")),
                              "g__Peptoniphilus_s__Peptoniphilus_lacrimalis" = expression(italic("Peptoniphilus lacrimalis")),
                              "g__Porphyromonas_s__Porphyromonas_asaccharolytica" = expression(italic("Porphyromonas asaccharolytica")),
                              
                              "g__Tyzzerella_s__unidentified" = "Unidentified",
                              "g__Colidextribacter_NA" = expression(italic("Flavonifractor plautii")),
                              
                              "g__Incertae_Sedis_s__uncultured_organism" = expression(italic("Acutalibacter muris"))
  )) +
  theme_classic() +
  theme(
    # Change background color to lighter purple
    panel.background = element_rect(fill = "#ebfaeb", color = "transparent"),
    axis.text.y = element_text(size = 12), 
    axis.text.x = element_text(size = 10),  
    axis.title.x = element_text(size = 14),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 14),
  ) +
  scale_y_continuous(expand = c(0, 0))
                                       
### PPMS LOLLIPOP ###
                                        
plot_unique_PPMS <- ggplot(summary_PPMS, aes(x = reorder(GenusSpecies, MeanAbundancePct), y = MeanAbundancePct)) +
  geom_segment(aes(xend = GenusSpecies, y = 0, yend = MeanAbundancePct), 
               color = "#7570b3", linewidth = 1.5) +  # Increase line thickness
  geom_point(color = "#7570b3", size = 5) + 
  coord_flip() + 
  labs(title = " ",
       x = "Species unique to PPMS",
       y = "Average relative abundance (%)") +
  scale_x_discrete(labels = c("f__Prevotellaceae_g__Prevotella_NA" = expression(italic("Segatella copri")),
                              "f__Ruminococcaceae_g__CAG-352_s__uncultured_bacterium" = expression(italic("Ruminococcoides bili")),
                              "f__Bacteroidaceae_g__Bacteroides_s__Bacteroides_plebeius" = "Unidentified ",
                              "f__Enterobacteriaceae_NA_NA" = expression(italic("Alistipes indistinctus")),
                              "f__Peptostreptococcaceae_g__Terrisporobacter_s__uncultured_bacterium" = expression(italic("Terrisporobacter mayombei")),
                              "f__Lachnospiraceae_g__[Ruminococcus]_torques_group_NA" = expression(italic("Ruminococcus faecis")),
                              "f__Lachnospiraceae_g__Blautia_s__Ruminococcus_sp." = "Unidentified",
                              "f__Izemoplasmatales_g__Izemoplasmatales_s__uncultured_organism" = expression(italic("Blautia faecis")),
                              "f__Christensenellaceae_g__Christensenellaceae_R-7_group_NA" = "Unidentified",
                              "f__Oscillospiraceae_NA_NA" = expression(italic("Intestinimonas butyriciproducens")),
                              "f__Lachnospiraceae_g__Frisingicoccus_NA" = expression(italic("Frisingicoccus caecimuris")),
                              "f__Oscillospiraceae_g__Colidextribacter_s__uncultured_Clostridia" = expression(italic("Colidextribacter massiliensis")),
                              "f__UCG-010_g__UCG-010_s__metagenome" = "Unidentified",
                              "f__Eggerthellaceae_g__Adlercreutzia_NA" = "Unidentified",
                              "f__Ruminococcaceae_g__Negativibacillus_s__uncultured_bacterium" = "Unidentified",
                              "f__UCG-010_g__UCG-010_s__gut_metagenome" = expression(italic("Adlercreutzia equolifaciens")),
                              "f__UCG-010_g__UCG-010_s__uncultured_bacterium" = expression(italic("Negativibacillus massiliensis"))
  )) +
  theme_classic() +
  theme(
    # Change background color to lighter purple
    panel.background = element_rect(fill = "#E6E6FA", color = "transparent"),
    axis.text.y = element_text(size = 12), 
    axis.text.x = element_text(size = 10),  
    axis.title.x = element_text(size = 14),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 14),
  ) +
  scale_y_continuous(expand = c(0, 0))


### Create combined plot ###
                                        
abundance_combined_plot <- ggdraw(
  plot_grid( 
    plot_unique_SPMS, 
    plot_unique_PPMS, 
    plot_unique_RRMS,
    labels = c("A", "B", "C"), 
    ncol = 2
  )
)

ggsave("Experimental Aim 2/Visualizations/abund_unique_spms_ppms_rrms", plot = abundance_combined_plot, width = 15, height = 8)
