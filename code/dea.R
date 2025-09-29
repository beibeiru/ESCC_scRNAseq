library(limma)
library(Seurat)

counts2CPM <- function(matrix.data)
{
	matrix.data <- t(t(matrix.data)*1e6/colSums(matrix.data))
	
	matrix.data
}

timePoint <- "H0"
sc.data <- Read10X(data.dir = paste0("../data/",timePoint,"/filtered_feature_bc_matrix/"))
colnames(sc.data) <- paste0(timePoint,"_",colnames(sc.data))

sc.matrix.data <- as.matrix(sc.data)
sc.matrix.data <- sc.matrix.data[rowSums(sc.matrix.data)>0,]

sc.matrix.data <- counts2CPM(sc.matrix.data)

sc.matrix.data <- round(log2(sc.matrix.data/100+1),1)

sc.matrix.anno <- read.csv("../results/integration_corrected_label.csv",row.names=1,as.is=T)
sc.matrix.anno <- sc.matrix.anno[sc.matrix.anno[,1]%in%timePoint,]



				
				TT <- as.numeric(sc.matrix.anno[,2]%in%c(0,4,5,7,8))
				WT <- as.numeric(sc.matrix.anno[,2]%in%c(1,3,6))
				
				design <- cbind(TT,WT)
				fit <- lmFit(sc.matrix.data,design)
				cont.matrix <- makeContrasts(TTvsWT=TT-WT,levels=design)
				fit2 <- contrasts.fit(fit, cont.matrix)
				fit2 <- eBayes(fit2)
				res <- topTable(fit2,coef=1,number=nrow(sc.matrix.data))
				
				write.csv(res,"../results/H0_differentialExpression_c04578_vs_c136.csv",quote=F)
				
				
				
				
				
DefaultAssay(sc.integrated) <- "RNA"
sc.integrated <- NormalizeData(sc.integrated, verbose = FALSE)

geneList <- readLines("geneList.txt")
geneList <- geneList[geneList%in%rownames(sc.integrated)]

#for(gene in geneList)
for(gene in "GSTM3")
{
	p1 <- FeaturePlot(sc.integrated, c(gene))
	ggplot2::ggsave(paste0("../results/marker/",gene,".jpg"), p1, width = 15, height = 15, dpi=200, units = "cm")
}

#p2 <- DimPlot(sc.integrated, reduction = "tsne", label = TRUE)
#ggplot2::ggsave(paste0("../results/integration_corrected_cluster.jpg"), p2, width = 23, height = 16, dpi=200, units = "cm")


		