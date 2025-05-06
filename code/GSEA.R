#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/15/2024
# Description: GSEA pathway analysis on the Kluc tumor CD45+ scRNAseq data
#==============================================================================#

# Packages and environment variables
suppressPackageStartupMessages({
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
  library(escape)
  library(GSVA)
  library(MatrixGenerics)})

setwd("/home/hnatri/PD1_mm/")
set.seed(9999)
options(ggrepel.max.overlaps = Inf)

# Colors, themes, cell type markers, and plot functions
source("/home/hnatri/PD1_mm/code/utilities.R")
source("/home/hnatri/PD1_mm/code/PD1_mm_themes.R")
source("/home/hnatri/PD1_mm/code/CART_plot_functions.R")

#==============================================================================#
# Modified GSVA functions
#==============================================================================#

.GS.check <- function(gene.sets) {
  if(is.null(gene.sets)) {
    stop("Please provide the gene.sets you would like to use for 
            the enrichment analysis")
  }
  egc <- gene.sets
  if(inherits(egc, what = "GeneSetCollection")){
    egc <- GSEABase::geneIds(egc) # will return a simple list, 
    #which will work if a matrix is supplied to GSVA
  }
  return(egc)
}

.cntEval <- function(obj, 
                     assay = "RNA", 
                     type = "counts") {
  if (inherits(x = obj, what = "Seurat")) {
    cnts <- obj@assays[[assay]][type]
  } else if (inherits(x = obj, what = "SingleCellExperiment")) {
    pos <- ifelse(assay == "RNA", "counts", assay) 
    if(assay == "RNA") {
      cnts <- assay(obj,pos)
    } else {
      cnts <- assay(altExp(obj), pos)
    }
  } else {
    cnts <- obj
  }
  cnts <- cnts[rowSums2(cnts) != 0,]
  return(cnts)
}

.split_data.matrix <- function(matrix, chunk.size=1000) {
  ncols <- dim(matrix)[2]
  nchunks <- (ncols-1) %/% chunk.size + 1
  
  split.data <- list()
  min <- 1
  for (i in seq_len(nchunks)) {
    if (i == nchunks-1) {  #make last two chunks of equal size
      left <- ncols-(i-1)*chunk.size
      max <- min+round(left/2)-1
    } else {
      max <- min(i*chunk.size, ncols)
    }
    split.data[[i]] <- matrix[,min:max]
    min <- max+1    #for next chunk
  }
  return(split.data)
}

.gsva.setup <- function(data, egc) {
  params.to.use <- gsvaParam(exprData = data,
                             geneSets = egc,
                             kcdf = "Poisson")
  return(params.to.use)
}

is_seurat_object <- function(obj) inherits(obj, "Seurat")
is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
is_seurat_or_se_object <- function(obj) {
  is_seurat_object(obj) || is_se_object(obj)
}

.mapGeneSetsToFeatures <- function(gsets, features) {
  
  ## Aaron Lun's suggestion at
  ## https://github.com/rcastelo/GSVA/issues/39#issuecomment-765549620
  gsets2 <- CharacterList(gsets)
  mt <- match(gsets2, features)
  mapdgenesets <- as.list(mt[!is.na(mt)])
  
  if (length(unlist(mapdgenesets, use.names=FALSE)) == 0)
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")
  
  mapdgenesets
}

split_rows <- function (matrix, chunk.size = 1000) 
{
  nrows <- dim(matrix)[1]
  if(is.vector(matrix)){
    nrows <- length(matrix)
  }
  nchunks <- (nrows - 1)%/%chunk.size + 1
  split.data <- list()
  min <- 1
  for (i in seq_len(nchunks)) {
    if (i == nchunks - 1) {
      left <- nrows - (i - 1) * chunk.size
      max <- min + round(left/2) - 1
    }
    else {
      max <- min(i * chunk.size, nrows)
    }
    split.data[[i]] <- matrix[min:max,]
    min <- max + 1
  }
  return(split.data)
}

split_vector <- function (vector, chunk.size = 1000) 
{
  
  nrows <- length(vector)
  nchunks <- (nrows - 1)%/%chunk.size + 1
  split.data <- list()
  min <- 1
  for (i in seq_len(nchunks)) {
    if (i == nchunks - 1) {
      left <- nrows - (i - 1) * chunk.size
      max <- min + round(left/2) - 1
    }
    else {
      max <- min(i * chunk.size, nrows)
    }
    split.data[[i]] <- vector[min:max]
    min <- max + 1
  }
  return(split.data)
}

performNormalization <- function(sc.data,
                                 enrichment.data = NULL,
                                 assay = "escape",
                                 gene.sets = NULL,
                                 make.positive = FALSE,
                                 scale.factor = NULL,
                                 groups = NULL) {
  if(!is.null(assay)) {
    if(is_seurat_object(sc.data)) {
      assay.present <- assay %in% names(sc.data@assays)
    } else if (is_se_object(sc.data)) {
      assay.present <- assay %in% assays(sc.data)
    }
  } else {
    assay.present <- FALSE
  }
  
  .pull.Enrich <- function(sc, enrichment.name) {
    if (inherits(sc, "Seurat")) {
      values <- Matrix::t(sc[[enrichment.name]]["data"])
    } else if (inherits(sc, "SingleCellExperiment")) {
      if(length(assays(altExp(sc))) == 1) {
        values <- t(assay(altExps(sc)[[enrichment.name]]))
      }
    }
  }
  
  if(is_seurat_or_se_object(sc.data) & !is.null(assay) & assay.present) {
    enriched <- .pull.Enrich(sc.data, assay)
  } else {
    enriched <- enrichment.data
  }
  
  if(!is.null(scale.factor) & length(scale.factor) != dim(sc.data)[2]) {
    stop("If using a vector as a scale factor, please ensure the length matches the number of cells.")
  }
  
  #Getting the gene sets that passed filters
  egc <- .GS.check(gene.sets)
  names(egc) <- str_replace_all(names(egc), "_", "-")
  egc <- egc[names(egc) %in% colnames(enriched)]
  
  #Isolating the number of genes per cell expressed
  if(is.null(groups)){
    chunks <- dim(enriched)[[1]]
  }else{
    chunks <- min(groups, dim(enriched)[[1]])
  }
  
  if (is.null(scale.factor)) {
    cnts <- .cntEval(sc.data, assay = "RNA", type = "counts")
    print("Calculating features per cell...")
    egc.sizes <- lapply(egc, function(x){
      scales<-unname(Matrix::colSums(cnts[which(rownames(cnts) %in% x),]!=0))
      scales[scales==0] <- 1
      scales
    })
    egc.sizes <- split_rows(do.call(cbind,egc.sizes), chunk.size=chunks)
    rm(cnts)
  }else{
    egc.sizes <- split_vector(scale.factor, chunk.size=chunks)
  }
  enriched <- split_rows(enriched, chunk.size=chunks)
  
  print("Normalizing enrichment scores per cell...")
  #Dividing the enrichment score by number of genes expressed
  
  enriched <- mapply(function(scores, scales){
    scores/scales
  }, enriched, egc.sizes, SIMPLIFY = FALSE)
  
  enriched <- do.call(rbind, enriched)
  if(make.positive){
    enriched <- apply(enriched, 2, function(x){
      x+max(0, -min(x))
    })
  }
  if(is_seurat_or_se_object(sc.data)) {
    if(is.null(assay)) {
      assay <- "escape"
    }
    sc.data <- .adding.Enrich(sc.data, enriched, paste0(assay, "_normalized"))
    return(sc.data)
  } else {
    return(enriched)
  }
}

escape.matrix <- function(input.data, 
                          gene.sets = NULL, 
                          method = "ssGSEA", 
                          groups = 1000, 
                          min.size = 5,
                          normalize = FALSE,
                          make.positive = FALSE,
                          BPPARAM = SerialParam(),
                          ...) {
  egc <- .GS.check(gene.sets)
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts")
  egc.size <- lapply(egc, function(x) length(which(rownames(cnts) %in% x)))
  if (!is.null(min.size)){
    remove <- unname(which(egc.size < min.size))
    if(length(remove) > 0) {
      egc <- egc[-remove]
      egc.size <- egc.size[-remove]
      if(length(egc) == 0) {
        stop("No gene sets passed the minimum length - please reconsider the 'min.size' parameter")
      }
    }
  }
  
  scores <- list()
  splits <- seq(1, ncol(cnts), by=groups)
  print(paste('Using sets of', groups, 'cells. Running', 
              length(splits), 'times.'))
  split.data <- .split_data.matrix(matrix=cnts, chunk.size=groups)
  
  for (i in seq_along(splits)) {
    last <- min(ncol(cnts), i+groups-1)
    if(method == "GSVA") {
      parameters <- .gsva.setup(split.data[[i]], egc)
    } else if (method == "ssGSEA") {
      parameters <- .ssGSEA.setup(split.data[[i]], egc)
    }
    if(method %in% c("ssGSEA", "GSVA")) {
      a <- suppressWarnings(gsva(param = parameters, 
                                 verbose = FALSE,
                                 BPPARAM = BPPARAM))
    }
    scores[[i]] <- a
  }
  
  # Unable to calculate scores for some pathways for some sets
  complete_pathways <- lapply(scores, function(xx){
    rownames(xx)
  })
  
  complete_pathways <- Reduce(intersect, complete_pathways)
  length(complete_pathways)
  
  complete_GS_CANONICAL <- GS_CANONICAL[complete_pathways]
  
  scores_complete <- lapply(scores, function(xx){
    xx[complete_pathways,]
  })
  
  scores_complete <- do.call(cbind, scores_complete)
  output <- t(as.matrix(scores_complete))
  
  gene.sets_complete <- gene.sets[scores_complete]
  
  #Normalize based on dropout
  #if(normalize) {
  #  output <- performNormalization(sc.data = input.data,
  #                                 enrichment.data = output,
  #                                 #assay = NULL,
  #                                 gene.sets = scores_complete,
  #                                 make.positive = make.positive,
  #                                 groups = groups)
  #}
  return(output)
}

runEscape <- function(input.data, 
                      gene.sets = NULL, 
                      method = "ssGSEA", 
                      groups = 1000, 
                      min.size = 5,
                      normalize = FALSE,
                      make.positive = FALSE,
                      new.assay.name = "escape",
                      BPPARAM = BiocParallel::SerialParam(),
                      ...) {
  #.checkSingleObject(input.data)
  enrichment <- escape.matrix(input.data = input.data,
                              gene.sets = gene.sets,
                              method = method,
                              groups = groups,
                              min.size = min.size,
                              BPPARAM = BPPARAM)
  
  input.data <- .adding.Enrich(input.data, enrichment, new.assay.name)
  return(input.data)
}

.adding.Enrich <- function(sc, enrichment, enrichment.name) {
  if (inherits(sc, "Seurat")) {
    if (as.numeric(substr(sc@version,1,1)) == 5) {
      new.assay <- suppressWarnings(CreateAssay5Object(
        data = as.matrix(t(enrichment))))
    } else {
      new.assay <- suppressWarnings(CreateAssayObject(
        data = as.matrix(t(enrichment))))
    }
    
    suppressWarnings(sc[[enrichment.name]] <- new.assay)
  } else if (inherits(sc, "SingleCellExperiment")) {
    altExp(sc, enrichment.name) <- SummarizedExperiment(assays = t(enrichment))
    names(assays(altExp(sc, enrichment.name))) <- enrichment.name
  }
  return(sc)
}


#==============================================================================#
# Import data and run escape/GSVA
#==============================================================================#

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/PD1_mm_Seurat/PD1_mm_Seurat_merged.Rds")
seurat_data$Treatment_Day <- paste0(seurat_data$Treatment, "_", seurat_data$Day)

# Updating annotations
gs4_deauth()
markers_annotations  <- gs4_get("https://docs.google.com/spreadsheets/d/1iWYBouwQlQboI-rwiujC0QKJ6lq9XeTffbKm2Nz8es0/edit?usp=sharin#g")
sheet_names(markers_annotations)
annotations <- read_sheet(markers_annotations, sheet = "Cluster annotations")

seurat_data$celltype <- mapvalues(seurat_data$snn_res.0.5,
                                  from = annotations$snn_res.0.5,
                                  to = annotations$annotation)

seurat_data$celltype <- factor(seurat_data$celltype,
                               levels = sort(unique(as.character(seurat_data$celltype))))

DefaultAssay(seurat_data) <- "RNA"

# C2 = curated gene sets,
# H = Hallmark
GS <- getGeneSets(species = "Mus musculus", library = c("C2", "H"))
GS_CANONICAL <- GS[grep("KEGG|REACTOME|BIOCARTA|HALLMARK", names(GS),
                        ignore.case = TRUE)]

escape_res <- runEscape(input.data = seurat_data,
                        method = "GSVA", 
                        new.assay.name = "escapeGSVA",
                        gene.sets = GS_CANONICAL, 
                        groups = 500,
                        min.size = 5)

escape_res <- performNormalization(sc.data = escape_res,
                                   assay = "escapeGSVA", 
                                   gene.sets = GS_CANONICAL,
                                   scale.factor = escape_res$nFeature_RNA,
                                   #groups = 500,
                                   make.positive = T)

saveRDS(escape_res, "/tgen_labs/banovich/BCTCSF/PD1_mm_Seurat/PD1_mm_Seurat_GSEA.Rds")



