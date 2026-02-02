# ==============================================================================
# Title: scRNA-seq Analysis Pipeline for Fig 5E
# Paper: Lee, H., Jang, H., Chae, JB. et al. In vivo efficacy of NRL knockdown with cell-penetrating siRNA in retinal degeneration. Sci Rep (2025)
# Description:
#   Calculates and visualizes Log2 Fold Change of Rod and Cone marker genes
#   between Retina_NRL (Knockdown) and Retina_PBS (Control).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(readxl)
library(stringr) # For str_to_title


# ------------------------------------------------------------------------------
# 2. Download data
# ------------------------------------------------------------------------------
base.dir <- "/path/to/data"
sample.file = file.path(base.dir, "Photoreceptor_subclustering_final_Dimension_10_Resolution_0.4.RDS")
sample.rds <- readRDS(sample.file)

siNRL.rds <- sample.rds[,sample.rds$orig.ident == 'si_NRL']
PBS.rds <- sample.rds[,sample.rds$orig.ident == 'PBS']

# ------------------------------------------------------------------------------
# 3. Load markers
# ------------------------------------------------------------------------------
marker.file = file.path(base.dir, "Photoreceptor_marker_gene.xlsx")
marker.df = read_excel(marker.file, col_names = FALSE)

Rod_marker_genes = str_to_title(as.character(unlist(marker.df[1, ])))
Cone_marker_genes = str_to_title(as.character(unlist(marker.df[2, ])))


# ------------------------------------------------------------------------------
# 4. Calculate log2(fold change)
# ------------------------------------------------------------------------------
# For Cone marker plot, change 'Rod_marker_genes' to 'Cone_marker_genes'
foldchange.df <- log2(rowMeans(exp(siNRL.rds@assays$RNA@data[intersect(rownames(siNRL.rds), str_to_title(Rod_marker_genes)),,drop = FALSE]))/rowMeans(exp(PBS.rds@assays$RNA@data[intersect(rownames(PBS.rds), str_to_title(Rod_marker_genes)),,drop = FALSE]))) 
foldchange.df <- melt(foldchange.df)

foldchange.df$gene <- rownames(foldchange.df)
foldchange.df <- foldchange.df[order(foldchange.df$value),]
foldchange.df$gene <- factor(foldchange.df$gene, levels = c(foldchange.df$gene))


# ------------------------------------------------------------------------------
# 4. Visualization
# ------------------------------------------------------------------------------
ggplot(data = foldchange.df, aes(x = gene, y = value, fill = value < 0)) +
geom_bar(stat="identity", width = 0.7) + 
scale_fill_manual(values = c("#F999B7", "gray")) + # #F999B7"(Rod_marker_genes), "#C4DFAA"(Cone_marker_genes)
geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.5) +
geom_hline(yintercept = -0.5, linetype = "dashed", linewidth = 0.5) +
geom_hline(yintercept = 0, linewidth = 1) +
geom_text(aes(label = foldchange.df$gene), vjust = -0.2, size = 4) +
theme_minimal(base_size = 20) +
ylab("log2 (Fold Change)") +
guides(fill = FALSE) +
theme(axis.text.y = element_text(size = 20, color = "black"), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.y = element_line(size = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
