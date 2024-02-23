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
# https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing
#gs4_deauth()
#tumor_tables  <- gs4_get("https://docs.google.com/spreadsheets/d/1ApwXjEVtpPB87al6q3ab8TKvZYJTh3iNH1cuO-A_OoU/edit?usp=sharing")
#sheet_names(tumor_tables)
#celltype_annot <- read_sheet(tumor_tables, sheet = "Cluster annotations")
#head(celltype_annot)
#length(unique(celltype_annot$annotation))
#
#tumor_celltype_col <- celltype_annot$color_fig1
#names(tumor_celltype_col) <- celltype_annot$annotation


