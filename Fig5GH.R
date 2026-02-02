# ==============================================================================
# Title: scRNA-seq Analysis Pipeline for Fig 5G, H
# Paper: Lee, H., Jang, H., Chae, JB. et al. In vivo efficacy of NRL knockdown with cell-penetrating siRNA in retinal degeneration. Sci Rep (2025)
# Description:
#   Fig. 5G: Box plot of Nrl expression (Log2 normalized)
#   Fig. 5H: Bar plot of cell type proportion ratio (Nrl Knockdown vs Control)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(reshape2)     # For melt()
library(RColorBrewer)
library(ggplot2)
library(ggpubr)       # For stat_compare_means
library(dplyr)

# ------------------------------------------------------------------------------
# 2. Data Preparation
# ------------------------------------------------------------------------------
base.dir <- "/path/to/data"
sample.file = file.path(base.dir, "Photoreceptor_subclustering_final_Dimension_10_Resolution_0.4.RDS")
sample.rds <- readRDS(sample.file)


# ==============================================================================
# Figure 5G: Nrl expression boxplot
# ==============================================================================
data <- FetchData(sample.rds, vars = c("Nrl", "seurat_clusters"), slot = "data")
melt_data <- melt(data)
melt_data$log2 <- log2(melt_data$value)

# color
x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
x.cluster.color <- x.getPalette(length(unique(sample.rds@meta.data$seurat_clusters)))

# draw plot
ggplot(melt_data, aes(x = seurat_clusters, y = log2, fill = seurat_clusters)) +
geom_boxplot(outlier.stroke = 0.1) +
scale_fill_manual(values = x.cluster.color) +
ggtitle("Nrl") +
ylab("log2(Expression)") +
xlab("Cluster") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(face = "italic", hjust = 0.5), axis.text.x = element_text(size = 17, angle = 0, hjust = 0.5, color = 'black'), axis.text.y = element_text(size = 20, color = 'black'), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
geom_point(stat = 'summary', fun = 'mean', color = 'red', size = 2, position = position_dodge(width = 1)) +
stat_summary(fun = mean, geom = "text", col = "red", fontface = "bold", size = 5, aes(y = stage(log2, after_stat = 2.5), label = round(..y.., digits = 2)), position = position_dodge(width = 1)) +
stat_compare_means(size = 4, label.x.npc = "left", label.x = 1.2, label.y = 3, color = 'black')


# ==============================================================================
# Figure 5H: Ratio by ratio
# ==============================================================================
x.count.df <- data.frame()
for (temp.id in unique(sample.rds$orig.ident)){
  print(temp.id)
  
  temp.df <- subset(x = sample.rds, orig.ident == temp.id)
  
  temp.count.df <- data.frame(Samples = unique(temp.df$orig.ident), table(temp.df$seurat_clusters))
  colnames(temp.count.df) <- c("Samples", "Clusters", "Counts")
  
  x.count.df <- rbind(x.count.df, temp.count.df)
}

sum_siNRL <- sum(x.count.df[x.count.df$Samples == "si_NRL",]$Counts)
sum_PBS <- sum(x.count.df[x.count.df$Samples == "PBS",]$Counts)

# Ratio
ratio.df <- data.frame(Clusters = c(0:5), Percentage = subset(x.count.df, Samples == "si_NRL")$Counts/subset(x.count.df, Samples == "PBS")$Counts)
ratio.df <- ratio.df[order(ratio.df$Percentage),]
order.ratio.df <- ratio.df[order(ratio.df$Percentage, decreasing = F),]

x.count.df$Clusters <- factor(x.count.df$Clusters, levels = order.ratio.df$Clusters)

# Ratio by ratio
ratio.2.df <- x.count.df
ratio.2.df$Ratio <- c(subset(ratio.2.df, Samples == "si_NRL")$Counts/sum_siNRL, subset(ratio.2.df, Samples == "PBS")$Counts/sum_PBS)

ratiobyratio.df <- data.frame(Clusters = c(0:5), Ratio = subset(ratio.2.df, Samples == "si_NRL")$Ratio/subset(ratio.2.df, Samples == "PBS")$Ratio)

ratiobyratio.df <- ratiobyratio.df[order(ratiobyratio.df$Ratio),]
ratiobyratio.df <- ratiobyratio.df[order(ratiobyratio.df$Ratio, decreasing = F),]

ratiobyratio.df$Clusters <- factor(ratiobyratio.df$Clusters, levels = ratiobyratio.df$Clusters)

order.cluster.color <- c('0' = c(x.cluster.color))
names(order.cluster.color) <- c(0:5)

ggplot(data = ratiobyratio.df, aes(x = Clusters, y = Ratio, fill = Clusters)) +
geom_bar(stat = "identity") +
scale_fill_manual(values = order.cluster.color) +
xlab("Clusters") +
ylab("si-NRL ratio / PBS ratio") +
geom_hline(yintercept = 2, linetype = "dashed", linewidth = 2, color = "gray") +
theme_classic(base_size = 50) + theme(legend.position = "none")
