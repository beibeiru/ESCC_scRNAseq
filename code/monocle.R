library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

source("path.R")

sc.integrated <- readRDS(file = paste0(dataPath,"sc.integrated.rds"))


cds <- as.cell_data_set(sc.integrated)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
g <- wrap_plots(p1, p2)

ggsave("t2.jpg", g, width = 40, height = 20, dpi=300, units = "cm")

cds <- learn_graph(cds)
g <-  plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

ggsave("t3.jpg", g, width = 40, height = 40, dpi=300, units = "cm")

