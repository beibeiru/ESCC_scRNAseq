library(Seurat)
library(SingleCellExperiment)
library(slingshot)

batches <- c("H0","H6","H48","D77","D169m")

outputPath <- "../results/timePoints5/"
dir.create(outputPath)
integrationPath <- paste0(outputPath,"integration/")
dir.create(integrationPath)
singlePath <- paste0(outputPath,"single/")
dir.create(singlePath)

sc.integrated <- readRDS(file = paste0(outputPath,"sc.integrated.rds"))

rdt <- sc.integrated @ reductions $ tsne @ cell.embeddings
rdu <- sc.integrated @ reductions $ umap @ cell.embeddings
cl <- sc.integrated @ active.ident

sim <- SingleCellExperiment(assays = List(counts = as.matrix(sc.integrated@assays$RNA@counts)))
reducedDims(sim) <- SimpleList(umap = rdu , tsne = rdt)
colData(sim)$our <- cl


sim <- slingshot(sim, clusterLabels = 'our', reducedDim = 'tsne')

library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
plotcol[is.na(plotcol)] <- "grey"

jpeg(paste0(integrationPath,"corrected_cluster_pseudotime.jpg"))
plot(reducedDims(sim)$tsne, col = plotcol, pch=16, asp = 1, main="red -> blue")
lines(SlingshotDataSet(sim), lwd=2, col='black')
dev.off()
