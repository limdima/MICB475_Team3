library(dplyr)
library(tidyverse)
library(readxl)
# install.packages("readxl")

meta_untrimmed <- read_excel("ms_metadata.xlsx")

# remove duplicate rows

meta <- distinct(meta_untrimmed) 
  
# filter out probiotics, (etc.), and also any unpaired samples after the filtering

filter_meta <- filter(meta, probiotics == 0) %>%
# add more filters here for variables we don't need  # filter() %>%
  
  group_by(`household`) %>% # samples are paired by household, one MS and one control
  filter(n() == 2) %>% # filters out any rows that do not contain two rows per household
  ungroup()
