#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/03/2024
# Description: CellChat on Kluc CD45+ scRNA-seq
#==============================================================================#

#==============================================================================
# Import libraries
#==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(CellChat)
  library(googlesheets4)
  library(plyr)})

setwd("/home/hnatri/PD1_mm/")
set.seed(9999)

#==============================================================================
# Import data
#==============================================================================

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/PD1_mm_Seurat/PD1_mm_Seurat_GSEA.Rds")

# Updating annotations
gs4_deauth()
markers_annotations  <- gs4_get("https://docs.google.com/spreadsheets/d/1iWYBouwQlQboI-rwiujC0QKJ6lq9XeTffbKm2Nz8es0/edit?usp=sharin#g")
sheet_names(markers_annotations)
annotations <- read_sheet(markers_annotations, sheet = "Cluster annotations")

seurat_data$celltype <- mapvalues(seurat_data$snn_res.0.5,
                                  from = annotations$snn_res.0.5,
                                  to = annotations$annotation)

seurat_data$celltype <- factor(seurat_data$celltype,
                               levels = sort(as.character(unique(seurat_data$celltype))))

Idents(seurat_data) <- seurat_data$celltype

#==============================================================================
# Creating inputs
#==============================================================================

Idents(seurat_data) <- seurat_data$celltype

# Inputs for CellChat:
# Log-normalized counts and cell labels (cluster/celltype)
# Creating CC objects for CD3_high/low, merging for comparative analysis, as
# well as an object for all samples together
all_counts <- LayerData(seurat_data, assay = "RNA", layer = "data")
all_CC <- createCellChat(object = all_counts,
                         meta = seurat_data@meta.data,
                         group.by = "celltype")

## CD3_high
NEO_seurat <- subset(seurat_data, subset = Treatment_Day == "NEO_16")
NEO_counts <- LayerData(NEO_seurat, assay = "RNA", layer = "data")
NEO_CC <- createCellChat(object = NEO_counts, meta = NEO_seurat@meta.data, group.by = "celltype")

ADJ_seurat <- subset(seurat_data, subset = Treatment_Day == "ADJ_16")
ADJ_counts <- GetAssayData(ADJ_seurat, assay = "RNA", layer = "data")
ADJ_CC <- createCellChat(object = ADJ_counts, meta = ADJ_seurat@meta.data, group.by = "celltype")

# Adding metadata
NEO_CC <- addMeta(NEO_CC, meta = NEO_seurat@meta.data)
NEO_CC <- setIdent(NEO_CC, ident.use = "celltype")
ADJ_CC <- addMeta(ADJ_CC, meta = ADJ_seurat@meta.data)
ADJ_CC <- setIdent(ADJ_CC, ident.use = "celltype")

NEO_CC@idents <- factor(NEO_CC@idents, levels = sort(unique(NEO_CC@idents)))
ADJ_CC@idents <- factor(ADJ_CC@idents, levels = sort(unique(ADJ_CC@idents)))

#==============================================================================
# Running CellChat
#==============================================================================

# Setting the ligand-receptor database
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# Use a subset of CellChatDB for cell-cell communication analysis
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

CC_objects <- list("all_CC" = all_CC,
                   "NEO_CC" = NEO_CC,
                   "ADJ_CC" = ADJ_CC)

CC_objects <- lapply(CC_objects, function(xx){
  # Set the used database in the object
  xx@DB <- CellChatDB.use
  
  # Preprocessing the expression data for cell-cell communication analysis
  # Subset the expression data of signaling genes for saving computation cost
  xx <- subsetData(xx) # This step is necessary even if using the whole database
  
  xx <- identifyOverExpressedGenes(xx)
  xx <- identifyOverExpressedInteractions(xx)
  
  # Project gene expression data onto PPI (Optional: when running it, USER should
  # set `raw.use = FALSE` in the function `computeCommunProb()` in order to use
  # the projected data)
  xx <- projectData(xx, PPI.human)
  
  ## Inference of cell-cell communication network
  # Compute the communication probability and infer cellular communication network
  xx <- computeCommunProb(xx) # raw.use = FALSE
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  xx <- filterCommunication(xx, min.cells = 10)
  
  # Infer the cell-cell communication at a signaling pathway level
  unique(xx@idents)
  
  xx <- computeCommunProbPathway(xx)
  
  # Calculate the aggregated cell-cell communication network
  xx <- aggregateNet(xx)
  
  xx
})

# Merging objects for comparative analysis
#CC_objects_compare <- CC_objects[c("CR_SD_CC", "PD_CC")]
CC_objects_compare <- CC_objects[c("NEO_CC", "ADJ_CC")]
CC_merged_object <- mergeCellChat(CC_objects_compare, add.names = names(CC_objects_compare))

saveRDS(CC_objects_compare, "/tgen_labs/banovich/BCTCSF/PD1_mm_Seurat/NEO_ADJ_16_CellChat_compare.rds")
saveRDS(CC_merged_object, "/tgen_labs/banovich/BCTCSF/PD1_mm_Seurat/NEO_ADJ_16_CellChat_merged_object.rds")
