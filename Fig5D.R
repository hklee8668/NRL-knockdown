# ==============================================================================
# Title: scRNA-seq Analysis Pipeline for Fig 5D
# Paper: Lee, H., Jang, H., Chae, JB. et al. In vivo efficacy of NRL knockdown with cell-penetrating siRNA in retinal degeneration. Sci Rep (2025)
# Description:
#   This script performs standard scRNA-seq analysis including:
#   Preprocessing, Normalization, Cell Cycle Regression, Clustering, and Visualization.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(stringr)      # For str_to_title()
library(RColorBrewer) # For color
library(ggplot2)


# ------------------------------------------------------------------------------
# 2. Download data
# ------------------------------------------------------------------------------
base.dir <- "/path/to/data"
sample.file = file.path(base.dir, "Photoreceptor_subclustering_final_Dimension_10_Resolution_0.4.RDS")
sample.rds <- readRDS(sample.file)

# ------------------------------------------------------------------------------
# 3. Visualization
# ------------------------------------------------------------------------------

# color for seurat_clusters
x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
x.cluster.color <- x.getPalette(length(unique(sample.rds@meta.data$seurat_clusters)))

DimPlot(object = sample.rds, reduction = "umap", label = T, pt.size = 1, label.size = 15, cols = x.cluster.color, repel = TRUE)

# color for celltype
Cone.cluster <- c(2)
Cone_like <- c(4)
none.cluster <- c(5)

sample.rds$trajectory_celltype <- ifelse(sample.rds$seurat_clusters %in% Cone.cluster, "Cone", ifelse(sample.rds$seurat_clusters %in% Cone_like, "Cone_like", ifelse(sample.rds$seurat_clusters %in% none.cluster, "Unknown", "Rod")))
sample.rds$trajectory_celltype <- factor(sample.rds$trajectory_celltype, levels = c("Rod", "Cone", "Cone_like", "Unknown"))

DimPlot(object = sample.rds, reduction = "umap", group.by = "sub_celltype", pt.size = 1, cols = c("#F999B7", "#C4DFAA", "gray")) + 
theme(text = element_text(size = 30), axis.text = element_text(size = 30), legend.key.size = unit(1.1, 'cm'), plot.title = element_blank(), legend.text = element_text(size=30)) +
guides(color = guide_legend(override.aes = list(size=7)))

