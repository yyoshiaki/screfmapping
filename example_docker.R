## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(Azimuth)
library(patchwork)
library(tidyverse)
library(sctransform)
library(Matrix)

setwd('/home/rstudio/autoimmune_10x')
source('/screfmapping/ref_mapping_seuratobj.R')
source('/screfmapping/utils_seurat.R')

## ----Load　reference---------------------------------------------------------------------------
# Reference for CD4T cells are included in Docker repo. No need to modify this section for CD4T analysis.
# Azimuth
reference <- LoadReference(path = "/screfmapping/data/Azimuth/human_pbmc_v1.0.0")
# Symphony
load("/screfmapping/data/ref_Reference_Mapping_20220525.RData")
file.copy(from = '/screfmapping/data/cache_symphony_sct.uwot', 
  to = '/home/rstudio/autoimmune_10x/cache_symphony_sct.uwot')

## ----parameter setting (change here)---------------------------------------------------------------
project.name <- "example"
prefix <- paste0("./output/", project.name, "/", project.name)
dir.create(paste0("./output/", project.name), recursive = T)

## ----reading data (change here)--------------------------------------------------------------------
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/filtered_gene_bc_matrices/hg19/")
q <- CreateSeuratObject(counts = pbmc.data,
                        project = project.name,
                        assay = "RNA",
                        min.cells = 3,
                        min.features = 200)
# It may need the following lines if nFeature_RNA and nCount_RNA are not calculated.
q$nFeature_RNA <- colSums(GetAssayData(q, layer = "counts") > 0)
q$nCount_RNA <- colSums(GetAssayData(q, layer = "counts"))
q@meta.data['log_umi'] <- log10(q$nCount_RNA)

## ----extraction of CD4T-----------------------------------------------------------------------------
# if your seurat object contains only CD4T, skip this step
q <- extract_cells_seuratobj(q, reference, prefix)
write.csv(q@meta.data[, c('predicted.celltype.l1.score', 'predicted.celltype.l1')],
      file = paste0(prefix, '_Azimuth.csv'))

# load extracted CD4T
query_obj <- readRDS(paste0(prefix, "_CD4T_AssayData.rds"))
query_obj <- CreateSeuratObject(counts = query_obj,
                        project = project.name,
                        assay = "RNA",
                        min.cells = 3,
                        min.features = 200)
# It may need the following lines if nFeature_RNA and nCount_RNA are not calculated.
query_obj$nFeature_RNA <- colSums(GetAssayData(query_obj, layer = "counts") > 0)
query_obj$nCount_RNA <- colSums(GetAssayData(query_obj, layer = "counts"))
query_obj@meta.data['log_umi'] <- log10(query_obj$nCount_RNA)

## ----run Symphony-----------------------------------------------------------------------------
reference_mapping_seuratobj(ref, query_obj, prefix)

## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
