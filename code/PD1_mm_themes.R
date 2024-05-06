#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/23/2024
# Description: PD1 Kluc analysis colors and themes
#==============================================================================#

library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(plyr)
library(circlize)
library(googlesheets4)

#==============================================================================
# Colors and themes
#==============================================================================

# ggplot theme
manuscript_theme <- theme(text = element_text(size = 6),
                          axis.text.x = element_text(size = 6),
                          axis.text.y = element_text(size = 6),  
                          axis.title.x = element_text(size = 6),
                          axis.title.y = element_text(size = 6))

# Colors for plotting
# Define colors for each level of categorical variables

# Clusters
clusters <- as.factor(c(0, seq(19)))
cluster_col <- colorRampPalette(brewer.pal(11, "Paired"))(length(clusters))
names(cluster_col) <- levels(clusters)

# Cell types
gs4_deauth()

PD1_tables  <- gs4_get("https://docs.google.com/spreadsheets/d/1iWYBouwQlQboI-rwiujC0QKJ6lq9XeTffbKm2Nz8es0/edit?usp=sharing")
sheet_names(PD1_tables)
celltype_annot <- read_sheet(PD1_tables, sheet = "Cluster annotations")
head(celltype_annot)
length(unique(celltype_annot$annotation))

PD1_Kluc_celltypes <- sort(celltype_annot$annotation)
PD1_Kluc_celltype_col <- colorRampPalette(brewer.pal(11, "Spectral"))(length(PD1_Kluc_celltypes))
names(PD1_Kluc_celltype_col) <- PD1_Kluc_celltypes

# Colors for each treatment group and time point
treatment_day_col <- c("NEO_12" = "#85428a",
                       "NEO_16" = "#bb7abf",
                       "CTRL_12" = "#41697a",
                       "CTRL_16" = "#86b9cf",
                       "ADJ_16" = "#678c23")

treatment_col <- c("NEO" = "#85428a",
                   "CTRL" = "#41697a",
                   "ADJ" = "#678c23")

