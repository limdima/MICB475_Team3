library(dplyr)
library(tidyverse)
library(readxl)
# install.packages("readxl")
set.seed(1234) # set random seed

meta_untrimmed <- read_excel("ms_metadata.xlsx")
manifest_untrimmed <- read.delim("ms_manifest.tsv", sep="\t")

# remove duplicate rows

meta <- distinct(meta_untrimmed) 
  
# filter out probiotics, (etc.), and also any unpaired samples after the filtering

filter_meta <- meta %>%
  filter(probiotics == 0) %>%
  filter(diet_no_special_needs == 1) %>%
  filter(eating_disorder == 0) %>%
  
  group_by(`household`) %>% # samples are paired by household, one MS and one control
  filter(n() == 2) %>% # filters out any rows that do not contain two rows per household
  ungroup()

# get spms, rrms, and ppms samples

SPMS_filtered <- filter(filter_meta, disease_course=="SPMS") 
RRMS_filtered <- filter(filter_meta, disease_course=="RRMS")
PPMS_filtered <- filter(filter_meta, disease_course=="PPMS")

filter_meta %>% group_by(disease_course) %>% summarise(n = n())

# Changing the column names in manifest_untrimmed to match meta

colnames(manifest_untrimmed)  <- c("sample-id","absolute-filepath")

# Combining both data frames together with left join to filter out deleted rows in meta, then separating the data frames again

combined <- left_join(meta, manifest_untrimmed)
manifest <- combined[c('sample-id', 'absolute-filepath')]

write.csv(filter_meta, "meta_ms.csv", row.names = FALSE)
