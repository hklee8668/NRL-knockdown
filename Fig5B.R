# ==============================================================================
# Title: scRNA-seq Analysis Pipeline for Fig 5B
# Paper: Lee, H., Jang, H., Chae, JB. et al. In vivo efficacy of NRL knockdown with cell-penetrating siRNA in retinal degeneration. Sci Rep (2025)
# Description:
#   This script performs standard scRNA-seq analysis including:
#   Preprocessing, Normalization, Cell Cycle Regression, Clustering, and Visualization.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(RColorBrewer) # For color
library(ggplot2)


# ------------------------------------------------------------------------------
# 2. Download data & create Seurat objects
# ------------------------------------------------------------------------------
base.dir <- "/path/to/data"
sample.file = file.path(base.dir, "final_Dimension_15_Resolution_0.6_celltype.RDS")
sample.rds <- readRDS(sample.file)


# ------------------------------------------------------------------------------
# 3. Visualization
# ------------------------------------------------------------------------------
# order
celltype <- c("Microglia", "Muller_glia", "Astrocyte", "Endothelial_cell", "Pericyte", "Bipolar_cell", "Cone", "Rod", "Melanocyte", "Amacrine_cell", "Horizontal_cell", "Retinal_ganglion_cell", "Unknown")
sample.rds@meta.data$celltype <- factor(sample.rds@meta.data$celltype, levels = celltype)

# color
y.getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
y.celltype.color <- y.getPalette(length(unique(sample.rds@meta.data$celltype)))

DimPlot(object = sample.rds, reduction = "umap", group.by = "celltype", pt.size = 0.0005, cols = y.celltype.color)
