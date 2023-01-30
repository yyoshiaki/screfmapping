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

## ----Load_the_reference---------------------------------------------------------------------------
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

## ----extraction of CD4T-----------------------------------------------------------------------------
# if your seurat object contains only CD4T, skip this step
extract_cells_seuratobj(q, reference, prefix)

# load extracted CD4T
query_obj <- readRDS(paste0(prefix, "_CD4T_AssayData.rds"))
query_obj <- CreateSeuratObject(counts = query_obj,
                        project = project.name,
                        assay = "RNA",
                        min.cells = 3,
                        min.features = 200)

## ----run Symphony-----------------------------------------------------------------------------
reference_mapping_seuratobj(ref, query_obj, prefix)

## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()