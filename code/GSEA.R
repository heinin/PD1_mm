#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/15/2024
# Description: GSEA pathway analysis on the Kluc tumor CD45+ scRNAseq data
#==============================================================================#

# Packages and environment variables
suppressPackageStartupMessages({
  #library(cli)
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(tidyverse)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  #library(workflowr)
  library(googlesheets4)
  library(escape)})

setwd("/home/hnatri/PD1_mm/")
set.seed(9999)
options(ggrepel.max.overlaps = Inf)

# Colors, themes, cell type markers, and plot functions
source("/home/hnatri/PD1_mm/code/utilities.R")
source("/home/hnatri/PD1_mm/code/PD1_mm_themes.R")
source("/home/hnatri/PD1_mm/code/CART_plot_functions.R")

# Importing data
seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/PD1_mm_Seurat/PD1_mm_Seurat_merged.Rds")
seurat_data$Treatment_Day <- paste0(seurat_data$Treatment, "_", seurat_data$Day)

### Running GSEA using R/escape
# Importing gene sets and focusing on canonical KEGG, RACTOME, and BIOCARTA
# pathways
GS <- getGeneSets(species = "Homo sapiens", library ="C2")
GS_CANONICAL <- GS[grep("KEGG|REACTOME|BIOCARTA", names(GS),
                        ignore.case = TRUE)]

res <- enrichIt(obj = seurat_data,
                gene.sets = GS_CANONICAL,
                groups = 50, cores = 4)
seurat_data <- AddMetaData(seurat_data, res)

saveRDS(seurat_data, "/scratch/hnatri/CART/PD1_mm_Seurat_GSEA.Rds")
saveRDS(res, "/scratch/hnatri/CART/KlucCD45pos_GSEA_res.rds")

