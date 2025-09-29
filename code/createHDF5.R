source("path.R")

sc.integrated <- readRDS(file = paste0(dataPath,"sc.integrated.rds"))

barcodes <- sc.integrated @ assays $ RNA @ counts @ Dimnames [[2]]
barcodes_new <- c()
for(i in 1:length(batches))
{
	barcodes_sub <- barcodes[grepl(batches[i],barcodes)]
	barcodes_sub <- sapply(strsplit(barcodes_sub,"_",fixed=T),function(x) return(paste0(x[2])))
	barcodes_sub <- gsub(1,i,barcodes_sub)
	barcodes_new <- c(barcodes_new,barcodes_sub)
}



coord <- sc.integrated@reductions$tsne@cell.embeddings
coord <- cbind(barcodes_new,coord)
colnames(coord) <- c("Barcode","X Coordinate","Y Coordinate")

write.csv(coord, file = paste0("../results/5TimePoints/integration/cloupe/tSNE_Coordinates.csv"),row.names=F,quote=F)

coord <- sc.integrated@reductions$umap@cell.embeddings
coord <- cbind(barcodes_new,coord)
colnames(coord) <- c("Barcode","X Coordinate","Y Coordinate")

write.csv(coord, file = paste0("../results/5TimePoints/integration/cloupe/UMAP_Coordinates.csv"),row.names=F,quote=F)


clust <- cbind(barcodes_new,paste0("cluster ",sc.integrated@ meta.data$ seurat_clusters))
colnames(clust) <- c("Barcode","Cluster")

write.csv(clust, file = paste0("../results/5TimePoints/integration/cloupe/clusters.csv"),row.names=F,quote=F)


timep <- cbind(barcodes_new,sc.integrated@ meta.data$ orig.ident)
colnames(timep) <- c("Barcode","Time")

write.csv(timep, file = paste0("../results/5TimePoints/integration/cloupe/timePoints.csv"),row.names=F,quote=F)
