# Title: Metagenomic Sample Distance Analysis
# Author: Feargal Ryan
# Date: 2025-02-14
# GitHub: https://github.com/feargalr/
#
# Description:
# This script processes metagenomic samples  and calculates pairwise sample 
# distances using Aitchison, Bray-Curtis, and Manhattan distance metrics. 
#
# Libraries
library(phyloseq)
library(vegan)
library(microbiome)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# First select input data. As an example we will select public data from
# the curatedMetagenomicData library and pick a dataset with samples from
# 3 countries.

library(curatedMetagenomicData)
allmeta_data <- sampleMetadata
table(allmeta_data$age_category,allmeta_data$body_site)

subset.df = allmeta_data[allmeta_data$body_site %in% c("milk","vagina"),]







VatanenT_2016_meta <- allmeta_data[allmeta_data$study_name == "AsnicarF_2017", ]
rownames(VatanenT_2016_meta) <- VatanenT_2016_meta$sample_id

# Load Count Data
VatanenT_2016cmd_object <- curatedMetagenomicData(
  pattern = "AsnicarF_2017.relative_abundance",
  counts = TRUE,
  dryrun = FALSE
)
VatanenT_2016counts.df <- assay(VatanenT_2016cmd_object[[1]])

# Filter Data
VatanenT_2016counts.df <- VatanenT_2016counts.df[apply(VatanenT_2016counts.df > 1000, 1, sum) > 4, ]
VatanenT_2016md.df <- allmeta_data[allmeta_data$sample_id %in% colnames(VatanenT_2016counts.df), ]
rownames(VatanenT_2016md.df) <- VatanenT_2016md.df$sample_id

# Assign Species Names
spec_vec <- gsub("s__", "", sapply(rownames(VatanenT_2016counts.df), function(y) strsplit(y, "\\|")[[1]][7]))
rownames(VatanenT_2016counts.df) <- spec_vec

# Create Phyloseq Object
OTU <- otu_table(VatanenT_2016counts.df, taxa_are_rows = TRUE)
map_data <- sample_data(VatanenT_2016_meta)
physeq <- phyloseq(OTU, map_data)

######################################################################################################
# Function to Compute and Plot Distance Metrics
compute_and_plot_distance <- function(physeq, method, title, y_label) {
  dist_mat <- vegdist(t(otu_table(physeq)), method = method)
  
  # Convert distance matrix to data frame
  dist_matrix <- as.matrix(dist_mat)
  samples <- rownames(dist_matrix)
  
  dist_df <- expand.grid(sample1 = samples, sample2 = samples, stringsAsFactors = FALSE) %>%
    filter(sample1 < sample2) %>%
    mutate(distance = dist_matrix[cbind(sample1, sample2)]) %>%
    left_join(VatanenT_2016_meta %>% select(sample_id, country) %>% rename(country1 = country), 
              by = c("sample1" = "sample_id")) %>%
    left_join(VatanenT_2016_meta %>% select(sample_id, country) %>% rename(country2 = country), 
              by = c("sample2" = "sample_id")) %>%
    mutate(pair_type = ifelse(country1 == country2, "Within-country", "Between-country"))
  
  # Plot results
  ggplot(dist_df, aes(x = pair_type, y = distance)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    labs(title = title, x = "Comparison Type", y = y_label) +
    stat_compare_means()
}
######################################################################################################

# Compute and Plot Aitchison Distance
ps_clr <- microbiome::transform(physeq, "clr")
compute_and_plot_distance(ps_clr, "euclidean", "Pairwise Aitchison Distances", "Aitchison Distance")

# Compute and Plot Bray-Curtis Distance
compute_and_plot_distance(physeq, "bray", "Pairwise Bray-Curtis Distances", "Bray-Curtis Distance")

# Compute and Plot Manhattan Distance
compute_and_plot_distance(physeq, "manhattan", "Pairwise Manhattan Distances", "Manhattan Distance")
