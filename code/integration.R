source("path.R")

rm(sc.combined)
rm(sc.list)

statMat <- data.frame()

for(i in 1:length(batches))
{
	sc.data <- Read10X(data.dir = paste0("../data/",batches[i],"/filtered_feature_bc_matrix/"))
	colnames(sc.data) <- paste0(batches[i],"_",colnames(sc.data))

	sc <- CreateSeuratObject(counts = sc.data, project = batches[i])
	
	sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
	
	
	g <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggplot2::ggsave(paste0(QCPath,batches[i],".jpg"), g, width = 30, height = 10, dpi=200, units = "cm")

	statMat[batches[i],1] <- dim(sc)[1]
	statMat[batches[i],2] <- dim(sc)[2]
	
	sc <- subset(sc, subset = nFeature_RNA >= 200 & nFeature_RNA <= 8000 & nCount_RNA >= 1000 & percent.mt <= 15)
	
	g <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggplot2::ggsave(paste0(QCPath,batches[i],"_filtered.jpg"), g, width = 30, height = 10, dpi=200, units = "cm")

	statMat[batches[i],3] <- dim(sc)[2]
	statMat[batches[i],4] <- paste0(round((statMat[batches[i],3]/statMat[batches[i],2])*100),"%")

	if(!exists("sc.combined"))
	{
		sc.combined <- sc
	}else{
		sc.combined <- merge(sc.combined, y = sc, project = batches[i])
	}
	
	if(!exists("sc.list"))
	{
		sc.list <- list()
	}
	sc.list[[batches[i]]] <- sc
}

names(statMat) <- c("geneNumber","cellNumber","cellNumberAfterFiltering","cellPercentageAfterFiltering")
write.table(statMat,paste0(outputPath,"sc_stat.txt"),quote=F,sep="\t")

# sc.combined

sc.combined <- NormalizeData(sc.combined, normalization.method = "LogNormalize", scale.factor = 10000)
sc.combined <- FindVariableFeatures(sc.combined, selection.method = "vst", nfeatures = 2000)

sc.combined <- ScaleData(sc.combined)

sc.combined <- RunPCA(sc.combined, npcs = pcsTop, verbose = FALSE)
sc.combined <- RunTSNE(sc.combined, dims = 1:pcsTop, check_duplicates = FALSE)
sc.combined <- RunUMAP(sc.combined, dims = 1:pcsTop)

sc.combined @ meta.data $ orig.ident <- factor(sc.combined @ meta.data $ orig.ident,levels=batches)

for(reduc in c("umap","tsne"))
{
	g <- DimPlot(sc.combined, reduction = reduc, group.by = "orig.ident")
	ggplot2::ggsave(paste0(uncorrectedPath,reduc,".jpg"), g, width = 18, height = 16, dpi=300, units = "cm")
}

saveRDS(sc.combined, file = paste0(dataPath,"sc.combined.rds"))


# sc.integrated

for (i in 1:length(sc.list)) {
    sc.list[[i]] <- NormalizeData(sc.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
    sc.list[[i]] <- FindVariableFeatures(sc.list[[i]], selection.method = "vst", nfeatures = 000)
}

sc.anchors <- FindIntegrationAnchors(object.list = sc.list, dims = 1:pcsTop)
sc.integrated <- IntegrateData(anchorset = sc.anchors, dims = 1:pcsTop)

sc.integrated <- ScaleData(sc.integrated, verbose = FALSE)

sc.integrated <- RunPCA(sc.integrated, npcs = pcsTop, verbose = FALSE)
sc.integrated <- RunTSNE(sc.integrated, dims = 1:pcsTop, check_duplicates = FALSE)
sc.integrated <- RunUMAP(sc.integrated, dims = 1:pcsTop)

sc.integrated @ meta.data $ orig.ident <- factor(sc.integrated @ meta.data $ orig.ident,levels=batches)

for(reduc in c("umap","tsne"))
{
	g <- DimPlot(sc.integrated, reduction = reduc, group.by = "orig.ident")
	ggplot2::ggsave(paste0(dimensionReductionPath,reduc,"_group.by_batch.jpg"), g, width = 18, height = 16, dpi=300, units = "cm")

	g <- DimPlot(sc.integrated, reduction = reduc, group.by = "orig.ident", split.by="orig.ident")
	ggplot2::ggsave(paste0(dimensionReductionPath,reduc,"_group.by_batch_split.by_batch.jpg"), g, width = 50, height = 10, dpi=300, units = "cm")

}

saveRDS(sc.integrated, file = paste0(dataPath,"sc.integrated.rds"))
