##### Install packages #####
# Start by installing all necessary packages when asked if you want to install
# from source, please just type Yes in the terminal below

# If you don't have BiocManager, here is the code to install it
# A lot of you probably already have this so you can skip
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# Use the above list to install all the packages using a for loop
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
# when asked if you want to update all, some or none please type "n" for none

# After installing all of its above dependencies, install ggpicrust2
install.packages("ggpicrust2")

#### Load packages ####
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(dplyr)


#### Import files and preparing tables ####
#Importing the pathway PICrust2
abundance_file <- "Experimental Aim 3/pathway_abundance.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1)
abundance_data  =as.data.frame(abundance_data)
colnames(abundance_data)[1] <- "pathways"

#Import your metadata file, no need to filter yet
metadata <- read_delim("Experimental Aim 3/corrected_ms_metadata.tsv")

#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$disease_course),]

#Example Looking at subject number
#If you have multiple variants, filter your metadata to include only 2 at a time
metadata_SPMS_vs_control <- filter(metadata, disease_course == "SPMS" | disease_course == "Control")
metadata_PPMS_vs_control <- filter(metadata, disease_course == "PPMS" | disease_course == "Control")
metadata_RRMS_vs_control <- filter(metadata, disease_course == "RRMS" | disease_course == "Control")
metadata_SPMS_vs_PPMS <- filter(metadata, disease_course == "SPMS" | disease_course == "PPMS")
metadata_SPMS_vs_RRMS <- filter(metadata, disease_course == "SPMS" | disease_course == "RRMS")
metadata_PPMS_vs_RRMS <- filter(metadata, disease_course == "PPMS" | disease_course == "RRMS")


#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = metadata_SPMS_vs_control$'sample-id'
sample_names = append(sample_names, "pathways")
abundance_data_filtered = abundance_data[,colnames(abundance_data) %in% sample_names] #This step is the actual filtering
# 512 samples in abundance data filtered for SPMS and control


#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata = metadata[metadata$`sample-id` %in% abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq ####
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathways"), 
                                        metadata = metadata, group = "disease_course", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.005 <- abundance_daa_results_df %>% filter(p_values < 0.005)  # Filtered for p vales < 0.005 (more stringent, filters out more pathways)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.005,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.005)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathways %in% feature_with_p_0.005$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description

#this line will change for each dataset. 512 represents the number of samples in the filtered abundance table for SPMS/Control
abundance_desc = abundance_desc[,-c(513:ncol(abundance_desc))] 

# Generate a heatmap
heatmap_SPMS_control <- pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), 
                metadata = metadata, 
                group = "disease_course")

ggsave("Experimental Aim 3/Visualizations/heatmap SPMS vs Control.png", plot = heatmap_SPMS_control)


# NEED TO TROUBLESHOOT: 
# Seemingly because in the abundance_data_filtered, there are some cells that have 0.0000e+00 which is messing with the PCA function
# To circumvent, log transform then removed Inf values (though this may remove a lot of rows)

# Generate pathway PCA plot
meta_rename = metadata %>% dplyr::rename(`sample_name` = `sample-id`)
colnames_abundance_data_filtered <- abundance_data_filtered %>% column_to_rownames("pathways")
colnames_abundance_data_filtered <-  log(colnames_abundance_data_filtered)
colnames_abundance_data_filtered <- colnames_abundance_data_filtered[apply(colnames_abundance_data_filtered, 1, function(x) !any(x == -Inf)), ]

pathway_pca(abundance = colnames_abundance_data_filtered, 
            metadata = meta_rename, 
            group = "disease_course")


### TO DO IF PREFERED: 
# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")

# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata, "subject")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change

sig_res <- sig_res[order(sig_res$log2FoldChange),]
ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")
