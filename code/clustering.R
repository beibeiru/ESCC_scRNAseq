source("path.R")

sc.integrated <- readRDS(file = paste0(dataPath,"sc.integrated.rds"))

sc.integrated <- FindNeighbors(sc.integrated, dims = 1:pcsTop)
sc.integrated <- FindClusters(sc.integrated, resolution = 0.6)

for(reduc in c("umap","tsne"))
{
	g <- DimPlot(sc.integrated, reduction = reduc, group.by = "seurat_clusters", label=TRUE)
	ggplot2::ggsave(paste0(dimensionReductionPath,reduc,"_group.by_cluster.jpg"), g, width = 18, height = 16, dpi=300, units = "cm")

	g <- DimPlot(sc.integrated, reduction = reduc, group.by = "seurat_clusters", split.by="seurat_clusters", ncol=4)
	ggplot2::ggsave(paste0(dimensionReductionPath,reduc,"_group.by_cluster_split.by_cluster.jpg"), g, width = 30, height = 27, dpi=300, units = "cm")
	
	g <- DimPlot(sc.integrated, reduction = reduc, group.by = "seurat_clusters", split.by="orig.ident", label=TRUE)
	ggplot2::ggsave(paste0(dimensionReductionPath,reduc,"_group.by_cluster_split.by_batch.jpg"), g, width = 50, height = 10, dpi=300, units = "cm")
}

cluster2group <- list(
	decreased=c("0","3","4","6"),
	increased=c("1","2","5","8"),
	peaked=c("7","10","",""),
	stable=c("9","11","",""))

write.table(as.data.frame(cluster2group),file=paste0(clusterPath,"cluster2group.csv"), quote=F,sep=",",row.names=F)

temp <- as.character(sc.integrated @ meta.data $ seurat_clusters)
temp[temp%in%cluster2group[[1]]] <- names(cluster2group)[1]
temp[temp%in%cluster2group[[2]]] <- names(cluster2group)[2]
temp[temp%in%cluster2group[[3]]] <- names(cluster2group)[3]
temp[temp%in%cluster2group[[4]]] <- names(cluster2group)[4]

sc.integrated @ meta.data $ manual_groups <- factor(temp)

saveRDS(sc.integrated, file = paste0(dataPath,"sc.integrated.rds"))


# cluster batch stat
batch_clusters <- data.frame(
	TP=sc.integrated @ meta.data $ orig.ident,
	CL=sc.integrated @ meta.data $ seurat_clusters,
	stringsAsFactors=F)

rownames(batch_clusters) <- colnames(sc.integrated)

write.csv(batch_clusters,paste0(clusterPath,"cluster.csv"),quote=F)


CLtable <- c(table(batch_clusters$CL))

library(ggplot2)
library(dplyr)

# Create Data
fracMat <- data.frame(
  group=names(CLtable),
  value=CLtable,
  stringsAsFactors=FALSE
)

fracMat[[1]] <- factor(fracMat[[1]],levels=fracMat[[1]])

# Compute the position of labels
fracMat <- fracMat %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(fracMat$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  mutate(lname = paste0(group," (",value,", ",round(prop,1),"%)") )

fracMat[[5]] <- factor(fracMat[[5]],levels=rev(fracMat[[5]]))

# Basic piechart
g <- ggplot(fracMat, aes(x="", y=prop, fill=lname)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  labs(fill = "Cluster")+
  #theme(legend.position="none") +
  geom_text(aes(y = ypos, label = group), color = "white", size=6)

ggsave(paste0(clusterPath,"cluster_fraction_in_all_timePoints.jpg"), g, width = 16, height = 12, dpi=300, units = "cm")

##############

fracMat <- tapply(batch_clusters$CL,batch_clusters$TP,function(x) table(x)/length(x))
fracMat <- fracMat[batches]

fracMat <- matrix(unlist(fracMat),nrow=length(batches),byrow=T)
rownames(fracMat) <- batches
colnames(fracMat) <- as.character(sort(unique(batch_clusters$CL)))

write.csv(t(fracMat),paste0(clusterPath,"cluster_fraction_in_each_timePoint.csv"),quote=F)

fracMat.m <- reshape2::melt(fracMat)
fracMat.m[[1]] <- factor(fracMat.m[[1]], levels = batches)
fracMat.m[[2]] <- as.character(fracMat.m[[2]])
fracMat.m[[2]] <- factor(fracMat.m[[2]], levels = colnames(fracMat))
fracMat.m[[3]] <- round(fracMat.m[[3]],3)

library(ggplot2)

#band plot

g <- ggplot(fracMat.m, aes(x = Var1, y = value, group=Var2, shape=Var2))+
  #geom_line(size=1.5)+
  geom_area(aes(fill=Var2),color="white")+
  xlab("Time Points")+
  ylab("% Cells")+
  labs(fill = "Cluster")+
  theme_bw() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_text(size = 14,colour = "black")
  )
ggsave(paste0(clusterPath,"cluster_fraction_in_each_timePoint.jpg"), g, width = 15, height = 12, dpi=300, units = "cm")



for(i in 1:length(batches))
{
	sc.data <- Read10X(data.dir = paste0("../data/",batches[i],"/filtered_feature_bc_matrix/"))
	colnames(sc.data) <- paste0(batches[i],"_",colnames(sc.data))

	sc <- CreateSeuratObject(counts = sc.data, project = batches[i])
	
	sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
	sc <- subset(sc, subset = nFeature_RNA >= 200 & nFeature_RNA <= 8000 & nCount_RNA >= 1000 & percent.mt <= 15)
	
	sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
	
	sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
	
	sc <- ScaleData(sc)
	
	sc <- RunPCA(sc, features = VariableFeatures(object = sc))
	
	sc <- RunTSNE(sc, dims = 1:pcsTop, check_duplicates = FALSE)


	sc.matrix.anno <- read.csv(paste0(clusterPath,"cluster.csv"),row.names=1,as.is=T)
	sc.matrix.anno <- sc.matrix.anno[sc.matrix.anno[,1]%in%batches[i],]



	clustering10x <- read.csv(paste0("../data/",batches[i],"/Graph-based.csv"),as.is=T)
	clustering10x[[2]] <- gsub("Cluster ", "", clustering10x[[2]])
	
	clustering10x[,1] <- paste0(batches[i],"_",clustering10x[,1])
	clustering10x <- clustering10x[clustering10x[,1]%in%rownames(sc.matrix.anno),]
	
	head(Idents(object = sc))
	sc <- SetIdent(object = sc, value=clustering10x[[2]])
	sc @ active.ident <- factor(sc @ active.ident, levels=as.character(sort(as.numeric(unique(clustering10x[[2]])))))
	
	g <- DimPlot(sc, reduction = "tsne", label = TRUE, pt.size = 0.6) 


	ggplot2::ggsave(paste0(singlePath,batches[i],"_cluster_10x.jpg"), g, width = 16, height = 16, dpi=200, units = "cm")


	
	sc <- SetIdent(object = sc, value=sc.matrix.anno[[2]])
	
	sc @ active.ident <- factor(sc @ active.ident, levels=as.character(sort(unique(sc.matrix.anno[[2]]))))
	
	g <- DimPlot(sc, reduction = "tsne", label = TRUE, pt.size = 0.6) 


	ggplot2::ggsave(paste0(singlePath,batches[i],"_cluster_our.jpg"), g, width = 16, height = 16, dpi=200, units = "cm")

}

