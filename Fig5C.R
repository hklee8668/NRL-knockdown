# ==============================================================================
# Title: scRNA-seq Analysis Pipeline for Fig 5C
# Paper: Lee, H., Jang, H., Chae, JB. et al. In vivo efficacy of NRL knockdown with cell-penetrating siRNA in retinal degeneration. Sci Rep (2025)
# Description:
#   CPM normalization, Hierarchical clustering, and Visualization.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(readxl)
library(ggplot2)
library(reshape2)     # For melt()
library(stringr)      # For str_to_title()
library(RColorBrewer)
library(dplyr)        # For summarise()


# ------------------------------------------------------------------------------
# 2. Data preparation
# ------------------------------------------------------------------------------
base.dir <- "/path/to/data"
sample.file = file.path(base.dir, "final_Dimension_15_Resolution_0.6_celltype.RDS")
sample.rds <- readRDS(sample.file)

marker.file = file.path(base.dir, "marker_gene.xlsx")
marker.df = read_excel(marker.file)

celltype.df <- read.table(file.path(base.dir, "celltype.txt"), header = T)

# Calculate average expression (CPM) per 'seurat_clusters'
meta.df = sample.rds@meta.data

rna.mat = sample.rds@assays$RNA@counts
cpm.mat = RelativeCounts(rna.mat, scale.factor = 1e6)
cpm.mat = cpm.mat[rowSums(cpm.mat) > 0, ]

cluster.array = as.character(unique(meta.df$seurat_clusters))
for(i in 1:length(cluster.array)) {
  tmp.cluster = cluster.array[i]
  tmp.cells = rownames(subset(meta.df, seurat_clusters == tmp.cluster))
  tmp.cpm.df = as.data.frame(cpm.mat[,tmp.cells])
  tmp.bulk.cpm.df = data.frame(rowSums(tmp.cpm.df)/ncol(tmp.cpm.df))
  colnames(tmp.bulk.cpm.df) = tmp.cluster
  if(i == 1) {
    bulk.cpm.df = tmp.bulk.cpm.df
  } else {
    bulk.cpm.df = cbind(bulk.cpm.df, tmp.bulk.cpm.df)
  }
}


# ------------------------------------------------------------------------------
# 3. Prepare plot data (Z-score calculation)
# ------------------------------------------------------------------------------
# Expected excel columns: "Genes", "Types"
Genes = c(unique(marker.df$Genes))
marker.groups = c(unique(marker.df$Types))

z.plot.df = data.frame()
for(i in 1:length(marker.groups)) {
  tmp.group = marker.groups[i]
  markers = subset(marker.df, Types == tmp.group)$Genes

  markers <- str_to_title(markers) # mouse
  
  marker.cpm.df = bulk.cpm.df[markers,]
  marker.cpm.df$Marker = markers
  
  plot.df = melt(marker.cpm.df)
  colnames(plot.df) = c("Marker", "Cluster", "CPM")
  plot.df$Group = tmp.group
  
  max.cluster = max(as.numeric(as.character(plot.df$Cluster)))
  plot.df$Cluster = as.factor(plot.df$Cluster)
  plot.df$Cluster = factor(plot.df$Cluster, level = c(0:max.cluster))
  plot.df$log2CPM = log2(plot.df$CPM + 1)
  
  # Z normalize
  tmp.z.plot.df = data.frame()
  genes.array = unique(plot.df[,1])
  for(k in 1:length(genes.array)) {
    tmp.marker = genes.array[k]
    tmp.plot.df = subset(plot.df, Marker == tmp.marker)
    tmp.plot.df$Z_CPM = (tmp.plot.df$CPM - mean(tmp.plot.df$CPM))/sd(tmp.plot.df$CPM)
    tmp.plot.df$Z_log2CPM = (tmp.plot.df$log2CPM - mean(tmp.plot.df$log2CPM))/sd(tmp.plot.df$log2CPM)
    tmp.z.plot.df = rbind(tmp.z.plot.df, tmp.plot.df)
  }
  z.plot.df = rbind(z.plot.df, tmp.z.plot.df)
}

b.array = which(is.na(z.plot.df[,3]))
if(length(b.array) > 0) {
  z.plot.df = z.plot.df[-b.array,]
}


# ------------------------------------------------------------------------------
# 4. Hierarchical clustering & factor reordering
# ------------------------------------------------------------------------------
clust.df = data.frame()
z.plot.df$GM = paste0(z.plot.df$Group, ":", z.plot.df$Marker)
marker.array = unique(z.plot.df$Marker)
for(i in 1:length(marker.array)) {
  tmp.marker = marker.array[i]
  cpm.array = subset(z.plot.df, Marker == tmp.marker)$Z_log2CPM
  cluster.array = subset(z.plot.df, Marker == tmp.marker)$Cluster
  tmp.clust.df = as.data.frame(t(data.frame(cpm.array)))
  colnames(tmp.clust.df) = cluster.array
  rownames(tmp.clust.df) = tmp.marker
  tmp.clust.df = tmp.clust.df[,order(colnames(tmp.clust.df))]
  b.array = grep("\\.", colnames(tmp.clust.df))
  if(length(b.array) > 0) {
    tmp.clust.df = tmp.clust.df[,-b.array]
  }
  clust.df = rbind(clust.df, tmp.clust.df)
}
row.ord = hclust(dist(clust.df, method = "euclidean"), method = "ward.D2" )$order
col.ord = hclust(dist(as.data.frame(t(clust.df)), method = "euclidean"), method = "ward.D2")$order

# Set factor levels
z.plot.df$Marker = as.factor(z.plot.df$Marker)

z.plot.df$Cluster = as.factor(z.plot.df$Cluster)
z.plot.df$Cluster = factor(z.plot.df$Cluster, levels = colnames(clust.df)[col.ord])

celltype <- c("Microglia", "Muller_glia", "Astrocyte", "Endothelial_cell", "Pericyte", "Bipolar_cell", "Cone", "Rod", "Melanocyte", "Amacrine_cell", "Horizontal_cell", "Retinal_Ganglion_cell")
z.plot.df$Group = factor(z.plot.df$Group, levels = celltype)

z.plot.df$celltype <- "celltype"
for(i in unique(z.plot.df$Cluster)){
        z.plot.df[z.plot.df$Cluster == i,]$celltype <- celltype.df[celltype.df$seurat_clusters == i,]$celltype
}
z.plot.df$celltype <- factor(z.plot.df$celltype, levels = c(celltype, "Unknown"))


# ------------------------------------------------------------------------------
# 7. Visualization
# ------------------------------------------------------------------------------
ggplot() +
geom_tile(data = z.plot.df, aes(x = Marker, y = Cluster, fill = Z_log2CPM)) +
scale_fill_gradient2(low = "white", high = "red", midpoint = 0) +
facet_nested(vars(celltype), vars(Group), scales = "free") +
theme_classic(base_size = 35) +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text = element_text(), legend.title = element_blank(), strip.background = element_blank(), plot.title = element_blank(), axis.line = element_line(size = 2), strip.text.x = element_blank(),strip.text.y = element_text(angle = 0)) +
guides(fill = guide_colourbar(barwidth = 3, barheight = 15)) +
guides(color = guide_legend(override.aes = list(size = 20)))
