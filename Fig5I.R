# ==============================================================================
# Title: scRNA-seq Analysis Pipeline for Fig 5I
# Paper: Lee, H., Jang, H., Chae, JB. et al. In vivo efficacy of NRL knockdown with cell-penetrating siRNA in retinal degeneration. Sci Rep (2025)
# Description:
#   Visualizes expression of Rod and Cone marker genes across specific cell types
#   (Rod, Cone_like, Cone) using Violin plots.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(reshape2)     # For melt()
library(ggpubr)       # For stat_compare_means()
library(ggh4x)        # For facet_nested()


# ------------------------------------------------------------------------------
# 2. Data Preparation
# ------------------------------------------------------------------------------
base.dir <- "/path/to/data"
sample.file = file.path(base.dir, "Photoreceptor_subclustering_final_Dimension_10_Resolution_0.4.RDS")
sample.rds <- readRDS(sample.file)

# set cluster 4 = "Cone_like"
sample.rds$trajectory_celltype <- NA
sample.rds$trajectory_celltype[sample.rds$seurat_clusters %in% c(0, 1, 2)] <- "Rod"
sample.rds$trajectory_celltype[sample.rds$seurat_clusters %in% c(4)]       <- "Cone_like"
sample.rds$trajectory_celltype[sample.rds$seurat_clusters %in% c(3, 5)]    <- "Cone"

sample.rds$trajectory_celltype <- factor(sample.rds$trajectory_celltype, levels = c("Rod", "Cone_like", "Cone"))


# ------------------------------------------------------------------------------
# 3. Define Marker Genes
# ------------------------------------------------------------------------------
Rod_marker_genes = c("Cnga1", "Cngb1", "Gnat1", "Gnb1", "Nr2e3", "Nrl", "Pde6a", "Reep6", "Rho")
Cone_marker_genes = c("Arr3", "Ccdc136", "Cnga3", "Cngb3", "Gnat2", "Gngt2", "Opn1mw", "Opn1sw", "Pde6h", "Sall3")


# ------------------------------------------------------------------------------
# 4. make plot dataframe
# ------------------------------------------------------------------------------
data <- FetchData(sample.rds, vars = c(c(Rod_marker_genes, Cone_marker_genes), "trajectory_celltype"), slot = "data")
melt_data <- melt(data)

melt_data$log2 <- log2(melt_data$value)

sub.melt_data <- melt_data[melt_data$value > 0,]
sub.melt_data <- subset(sub.melt_data, variable %in% Rod_marker_genes) # Rod_marker_genes & Cone_marker_genes


# ------------------------------------------------------------------------------
# 5. Visualization
# ------------------------------------------------------------------------------
ggplot(sub.melt_data, aes(variable, value, fill = trajectory_celltype)) +
geom_violin(trim = FALSE, size = 0.1) +
ylim(-1, 7.5) +
scale_fill_manual(values = c("#F999B7", "#F3FDE8", "#C4DFAA")) + # "Rod" = "#F999B7", "Cone_like" = "#F3FDE8", "Cone" = "#C4DFAA"
geom_point(stat = 'summary', fun = 'mean', color = 'red', size = 1, position = position_dodge(width = 1)) +
ylab("log2(Expression)") +
facet_nested(~variable, scales = "free_x") +
theme_bw() +
theme(axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.background = element_rect(size = 0.1), strip.text.x = element_text(face = "bold.italic", size = 16, colour = "black"), strip.placement = "outside", strip.switch.pad.grid = unit('0.25', "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = "white", size = 0.1)) +
stat_compare_means(label = "p.signif", size = 10, label.x.npc = "left", label.x = 0.7, label.y = 6.5, color = 'black')
