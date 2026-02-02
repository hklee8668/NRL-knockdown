# ==============================================================================
# Title: scRNA-seq Analysis Pipeline for Fig 5F
# Paper: Lee, H., Jang, H., Chae, JB. et al. In vivo efficacy of NRL knockdown with cell-penetrating siRNA in retinal degeneration. Sci Rep (2025)
# Description:
#   Performs single-cell trajectory analysis using Monocle 2.
#   Infers pseudotime and visualizes cell differentiation grouping by clusters.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(monocle)
library(RColorBrewer)
library(ggplot2)


# ------------------------------------------------------------------------------
# 2. Data Preparation (Seurat -> Monocle)
# ------------------------------------------------------------------------------
base.dir <- "/path/to/data"
sample.file = file.path(base.dir, "Photoreceptor_subclustering_final_Dimension_10_Resolution_0.4.RDS")
sample.rds <- readRDS(sample.file)

x.data <- as(as.matrix(sample.rds@assays$RNA@data), 'sparseMatrix')
x.data <- x.data[rowSums(x.data) > 0,]

x.pd <- new('AnnotatedDataFrame', data = sample.rds@meta.data)
x.fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(x.data), row.names = row.names(x.data)))

x.cds <- newCellDataSet(x.data, phenoData = x.pd, featureData = x.fd, expressionFamily = negbinomial.size())


# ------------------------------------------------------------------------------
# 3. Preprocessing & Gene Selection
# ------------------------------------------------------------------------------
x.cds <- estimateSizeFactors(x.cds)
x.cds <- estimateDispersions(x.cds)

x.cds <- detectGenes(x.cds, min_expr = 0.1)

x.markers <- FindAllMarkers(object = sample.rds, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", assay="RNA", slot="data")
x.markers <- subset(x.markers, p_val_adj < 0.05)

x.order.genes <- unique(as.character(x.markers$gene))

x.cds <- setOrderingFilter(x.cds, x.order.genes)


# ------------------------------------------------------------------------------
# 4. Dimensionality Reduction & Ordering (DDRTree)
# ------------------------------------------------------------------------------
temp.cds <- reduceDimension(cds = x.cds, max_components = 4, method = 'DDRTree', norm_method = 'log')
temp.cds <- orderCells(temp.cds, reverse = TRUE)


# ------------------------------------------------------------------------------
# 5. Visualization
# ------------------------------------------------------------------------------
# pseudotime (Fig 5F, left)
plot_cell_trajectory(temp.cds, color_by = "Pseudotime", cell_size = 4) + scale_x_reverse()

# cluster (Fig 5F, right)
x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
x.cluster.color <- x.getPalette(length(unique(sample.rds@meta.data$seurat_clusters)))

plot_cell_trajectory(temp.cds, color_by = "seurat_clusters", cell_size = 2) + scale_color_manual(values = x.cluster.color) + scale_x_reverse()

