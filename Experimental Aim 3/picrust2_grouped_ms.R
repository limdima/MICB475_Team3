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
# 916 samples (+1 column for pathway names)

#Import your metadata file, no need to filter yet
metadata <- read_delim("Experimental Aim 3/corrected_ms_metadata.tsv")

#Remove NAs for your column of interest in this case subject
metadata <- metadata[!is.na(metadata$disease_course),]


#If you have multiple variants, filter your metadata to include only 2 at a time

metadata_SPMS_vs_control <- filter(metadata, disease_course == "SPMS" | disease_course == "Control")
metadata_PPMS_vs_control <- filter(metadata, disease_course == "PPMS" | disease_course == "Control")
metadata_RRMS_vs_control <- filter(metadata, disease_course == "RRMS" | disease_course == "Control")


#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names_SPMS = metadata_SPMS_vs_control$'sample-id'
sample_names_SPMS = append(sample_names_SPMS, "pathways")
abundance_data_filtered_SPMS = abundance_data[,colnames(abundance_data) %in% sample_names_SPMS] #This step is the actual filtering
# 512 samples in abundance data filtered for SPMS and control


sample_names_PPMS = metadata_PPMS_vs_control$'sample-id'
sample_names_PPMS = append(sample_names_PPMS, "pathways")
abundance_data_filtered_PPMS = abundance_data[,colnames(abundance_data) %in% sample_names_PPMS] #This step is the actual filtering
# 514 samples in abundance data filtered for PPMS and control


sample_names_RRMS = metadata_RRMS_vs_control$'sample-id'
sample_names_RRMS = append(sample_names_RRMS, "pathways")
abundance_data_filtered_RRMS = abundance_data[,colnames(abundance_data) %in% sample_names_RRMS] #This step is the actual filtering
# 804 samples in abundance data filtered for RRMS and control


#Removing individuals with no data that caused a problem for pathways_daa()  (No change in the abundance data filtered)
abundance_data_filtered_SPMS =  abundance_data_filtered_SPMS[, colSums(abundance_data_filtered_SPMS != 0) > 0]
abundance_data_filtered_PPMS =  abundance_data_filtered_PPMS[, colSums(abundance_data_filtered_PPMS != 0) > 0]
abundance_data_filtered_RRMS =  abundance_data_filtered_RRMS[, colSums(abundance_data_filtered_RRMS != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered_SPMS) = NULL
rownames(abundance_data_filtered_PPMS) = NULL
rownames(abundance_data_filtered_RRMS) = NULL


#verify samples in metadata match samples in abundance_data
abun_samples_SPMS = rownames(t(abundance_data_filtered_SPMS[,-1])) #Getting a list of the SPMS sample names in the newly filtered abundance data
metadata_SPMS = metadata[metadata$`sample-id` %in% abun_samples_SPMS,] #making sure the filtered metadata only includes these samples
# 512 sample IDs


abun_samples_PPMS = rownames(t(abundance_data_filtered_PPMS[,-1])) # Same but for PPMS
metadata_PPMS = metadata[metadata$`sample-id` %in% abun_samples_PPMS,] 
# 514 sample IDs


abun_samples_RRMS = rownames(t(abundance_data_filtered_RRMS[,-1])) # Same but for RRMS
metadata_RRMS = metadata[metadata$`sample-id` %in% abun_samples_RRMS,] 
# 804 sample IDs


#### DESEq ####
#Perform pathway DAA using DESEQ2 method
#SPMS vs Control
abundance_daa_results_df_SPMS <- pathway_daa(abundance = abundance_data_filtered_SPMS %>% column_to_rownames("pathways"), 
                                        metadata = metadata_SPMS, 
                                        group = "disease_course", 
                                        reference = "Control",
                                        daa_method = "DESeq2")

#PPMS vs Control
abundance_daa_results_df_PPMS <- pathway_daa(abundance = abundance_data_filtered_PPMS %>% column_to_rownames("pathways"), 
                                        metadata = metadata_PPMS, 
                                        group = "disease_course", 
                                        reference = "Control",
                                        daa_method = "DESeq2")
#RRMS vs Control
abundance_daa_results_df_RRMS <- pathway_daa(abundance = abundance_data_filtered_RRMS %>% column_to_rownames("pathways"), 
                                        metadata = metadata_RRMS, 
                                        group = "disease_course", 
                                        reference = "Control",
                                        daa_method = "DESeq2")


# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df_SPMS <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df_SPMS, 
                                                       ko_to_kegg = FALSE)

metacyc_daa_annotated_results_df_PPMS <- pathway_annotation(pathway = "MetaCyc", 
                                                            daa_results_df = abundance_daa_results_df_PPMS, 
                                                            ko_to_kegg = FALSE)

metacyc_daa_annotated_results_df_RRMS <- pathway_annotation(pathway = "MetaCyc", 
                                                            daa_results_df = abundance_daa_results_df_RRMS, 
                                                            ko_to_kegg = FALSE)



# Generating a bar plot representing log2FC from the custom deseq2 function
# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("Experimental Aim 3/DESeq2_function.R")

# Run the function on your own data
res_SPMS =  DEseq2_function(abundance_data_filtered_SPMS, metadata_SPMS, "disease_course")
res_SPMS$feature =rownames(res_SPMS)
res_SPMS_desc = inner_join(res_SPMS,metacyc_daa_annotated_results_df_SPMS, by = "feature")
#res_SPMS_desc = res_SPMS_desc[, -c(8:13)] # Don't really need to remove the columns, but could for clarity
# View(res_SPMS_desc)

res_PPMS =  DEseq2_function(abundance_data_filtered_PPMS, metadata_PPMS, "disease_course")
res_PPMS$feature =rownames(res_PPMS)
res_PPMS_desc = inner_join(res_PPMS,metacyc_daa_annotated_results_df_PPMS, by = "feature")
# View(res_PPMS_desc)


res_RRMS =  DEseq2_function(abundance_data_filtered_RRMS, metadata_RRMS, "disease_course")
res_RRMS$feature =rownames(res_RRMS)
res_RRMS_desc = inner_join(res_RRMS,metacyc_daa_annotated_results_df_RRMS, by = "feature")
# View(res_RRMS_desc)


# Filter to only include significant pathways, can be changed but RRMS has a lot between 0.5-1 fold change
sig_res_SPMS = res_SPMS_desc %>%
  filter(pvalue < 0.05) %>%         # p < 0.05
  filter(abs(log2FoldChange) > 1) # Also filtered by Log2fold change

sig_res_PPMS = res_PPMS_desc %>%
  filter(pvalue < 0.05) %>%         # p < 0.05
  filter(abs(log2FoldChange) > 1) # Also filtered by Log2fold change

sig_res_RRMS = res_RRMS_desc %>%
  filter(pvalue < 0.05) %>%         # p < 0.05
  filter(abs(log2FoldChange) > 1) #largest log fold change is about 1


combined_sig_res <- rbind(
    sig_res_SPMS,
    sig_res_PPMS,
    sig_res_RRMS)

combined_sig_res <- combined_sig_res[order(combined_sig_res$log2FoldChange),]

foldchange_pathways <- ggplot(data = combined_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = group1))+
  geom_bar(stat = "identity", position =  position_dodge(width = 0.8, preserve = "single")) + 
  theme_bw() +
  scale_fill_manual(values = c("#9994D6", "#d95f02", "#1b9e77")) +
  labs(x = "Log2 Fold Change Relative to Control", y="Pathways", fill = "MS Group")

ggsave("Experimental Aim 3/Visualizations/DESEQ Pathways relative to control.png", plot = foldchange_pathways)



#######################################################################################################################
# Below is some junk that may be useful in the future idk - AH



# Can check if this aligns with the pathway_daa() function
# Not super useful since can't really filter by log2 fold change the same way
# Combine all the abundance results into one dataframe 

combined_pathway_daa_results_df <- rbind(
  abundance_daa_results_df_SPMS,
  abundance_daa_results_df_PPMS,
  abundance_daa_results_df_RRMS)

metacyc_daa_annotated_results_df_combined <- pathway_annotation(pathway = "MetaCyc", 
                                                            daa_results_df = combined_pathway_daa_results_df, 
                                                            ko_to_kegg = FALSE)
# Filter p-values to only significant ones
significant_features <- combined_pathway_daa_results_df %>% 
  filter(p_values < 0.005)  # Filtered for p values < 0.05 - Why doesn't Evelyn use p adjusted in her tutorial?


significant_features_with_desc <- inner_join(significant_features, metacyc_daa_annotated_results_df_combined, by = "feature")
significant_features_with_desc <- significant_features_with_desc %>%
    select(-ends_with(".y")) %>%
    rename_with(~str_remove(., "\\.x$"), ends_with(".x")) %>%
    group_by(feature, group1) %>%
    slice_head(n = 1) %>%   # this removes the duplicate rows by selecting the first row in each "grouped" feature+group
    ungroup()



### Code for the heatmap or PCA won't really work as they only take 2 variables at a time
# Even when combining individual pathway_daa() runs, heat map and PCA doesn't work with 
# all three disease groups together (duplicate row names are not allowed)

#Changing the pathway column to description for the results 
#significant_features_with_desc <- inner_join(significant_features, metacyc_daa_annotated_results_df, by = "feature")

# inner join creates .x and .y columns (this occured as there are many-to-many comparisons)
# Need to remove .y columns (which are from metacyc_daa_annotated_results_df, not from feature_with_p_0.005)
# Then remove duplicate rows in .x columns
# Make sure # of rows match significant_features (201 entries)

#significant_features_with_desc <- significant_features_with_desc %>%
#  select(-ends_with(".y")) %>%
#  rename_with(~str_remove(., "\\.x$"), ends_with(".x")) %>%
#  group_by(feature, group1) %>%
#  slice_head(n = 1) %>%   # this removes the duplicate rows by selecting the first row in each "grouped" feature+group
#  ungroup()
# 201 entries, all features have the same name as significant_features and descriptions are added
# Check number of unique features
#n_distinct(feature_desc$feature) # 140 unique features

#significant_features_with_desc$feature <- significant_features_with_desc$description # Replaces the feature row with the description
#significant_features_with_desc <- significant_features_with_desc[,c(1:7)] # get rid of the description column as it's replaced the feature column already
#colnames(significant_features_with_desc) <- colnames(significant_features) # makes sure the colnames are the same as before


# Filter abundance_data for the significant features
#abundance_data_significant <- abundance_data
#colnames(abundance_data_significant)[1] <- "feature"  # Changes the first column name to feature for join step below
#abundance_data_significant <- abundance_data_significant[abundance_data_significant$feature %in% significant_features_with_desc$feature, ]


#Changing the pathway column to description for the abundance table
#abundance <- abundance_data %>% 
#  filter(pathways %in% significant_features$feature) # 140 rows for features (same as above, 140 unique features)
#colnames(abundance)[1] <- "feature"   # Changes the first column name to feature for join step below
#abundance_desc <- inner_join(abundance, metacyc_daa_annotated_results_df, by = "feature") # 3 times more rows, for
#abundance_desc$feature <- abundance_desc$description
#abundance_desc <- abundance_desc %>% filter(p_adjust < 0.005)

#this line will change for each dataset. 512 represents the number of samples in the filtered abundance table for SPMS/Control
# abundance_desc <- abundance_desc[,-c(918:ncol(abundance_desc))]  # same number of columns (916 samples + 1 features column) as abundance








