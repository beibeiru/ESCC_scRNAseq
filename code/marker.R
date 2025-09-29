source("path.R")
library(clusterProfiler)

sc.integrated <- readRDS(file = paste0(dataPath,"sc.integrated.rds"))

DefaultAssay(sc.integrated) <- "RNA"

iMarkers <- c("CD44","CST6","C19orf33","TACSTD2","S100A14","RHOD","TM4SF1")
iMarkers <- c("UCHL1","CYBA","CSRP2","BEX4","BEX2","IFITM2","TMEM98","AKR1B1")
iMarkers <- c("GPRC5A","MAL2","KRT7","TINAGL1","ARHGDIB","UCA1")

g <- FeaturePlot(sc.integrated, iMarkers, reduction = "tsne")
ggplot2::ggsave(paste0(markerPath,"iMarkers_6.jpg"), g, width = 20, height = 30, dpi=200, units = "cm")

g <- FeaturePlot(sc.integrated, iMarkers, reduction = "tsne", split.by = "orig.ident")
ggplot2::ggsave(paste0(markerPath,"iMarkers_6_tsne.jpg"), g, width = 50, height = 60, dpi=200, units = "cm")

g <- FeaturePlot(sc.integrated, iMarkers, reduction = "umap", split.by = "orig.ident")
ggplot2::ggsave(paste0(markerPath,"iMarkers_6_umap.jpg"), g, width = 50, height = 60, dpi=200, units = "cm")


as.matrix(sc.integrated@assays$RNA@data) -> expr
expr_0 <- expr[,x=="H0"]

geneList <- c(
	"CD44","CST6","C19orf33","TACSTD2","S100A14","RHOD","TM4SF1",
	"UCHL1","CYBA","CSRP2","BEX4","BEX2","IFITM2","TMEM98","AKR1B1",
	"GPRC5A","MAL2","KRT7","TINAGL1","ARHGDIB","UCA1"
	)

expr_0_sub <- expr_0[rownames(expr_0)%in%geneList,]

write.csv(expr_0_sub,paste0(markerPath,"iMarkers.csv"),quote=F)









iMarkers1 <- c("CD44","CST6","C19orf33","TACSTD2","S100A14","RHOD","TM4SF1")
iMarkers2 <- c("UCHL1","CYBA","CSRP2","BEX4","BEX2","IFITM2","TMEM98","AKR1B1")

iMarkers <- c(iMarkers1,iMarkers2)

smmy<-data.frame()
as.matrix(sc.integrated@assays$RNA@data) -> normMat
as.character(sc.integrated@meta.data$orig.ident) -> anno

for(Marker in iMarkers)
{
	for(batch in batches)
	{
		temp <- normMat[Marker,anno==batch]
		
		smmy[Marker,batch] <- sum(temp>0)/length(temp)
	
	}


}

write.table(smmy,paste0(markerPath,"iMarkers_percentage.txt"),quote=F)


clusters <- as.numeric(levels(sc.integrated @ meta.data $ seurat_clusters))

for(i in clusters)
{
	markerPath.s <- paste0(markerPath,"cluster_",i,"/")
	dir.create(markerPath.s)
	
	markers <- FindMarkers(sc.integrated, ident.1 = i)
	markers_pos <- markers[markers[,2]>0,]
	markers_neg <- markers[markers[,2]<0,]
	
	write.csv(markers,paste0(markerPath.s,"cluster_",i,"_marker.csv"),quote=F)
	
	topMarkers <- rownames(markers_pos)[1:15]
	
	g <- FeaturePlot(sc.integrated, topMarkers, reduction = "tsne", ncol=5)
	ggplot2::ggsave(paste0(markerPath.s,"cluster_",i,"_marker_pos_top.jpg"), g, width = 50, height = 30, dpi=200, units = "cm")

	markerList <- list('pos'=rownames(markers_pos),'neg'=rownames(markers_neg))

	for(pn in c("pos","neg"))
	{
		symbol2id <- bitr(markerList[[pn]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		
		kk <- enrichKEGG(
			gene = symbol2id[,2], 
			organism = 'hsa', 
			qvalueCutoff = 0.05)
		
		for(xx in 1:dim(kk)[1])
		{
			geneIDs <- unlist(strsplit(kk[xx,"geneID"],"/"))
			geneSymbols <- symbol2id[match(geneIDs,symbol2id[,2]),1]
			geneIDs_geneSymbols <- paste0(geneIDs,"_",geneSymbols)
			
			kk@result$geneID[xx] <- paste0(geneIDs_geneSymbols,collapse="/")
		}
		
		write.csv(kk,paste0(markerPath.s,"cluster_",i,"_marker_",pn,"_enrichKEGG.csv"),quote=F)
		
		for(gosub in c("BP","MF","CC"))
		{
			ego <- enrichGO(
				gene = symbol2id[,2],
				OrgDb = org.Hs.eg.db,
				ont = gosub,
				qvalueCutoff = 0.05,
			   readable = TRUE)
		
			write.csv(ego,paste0(markerPath.s,"cluster_",i,"_marker_",pn,"_enrichGO_",gosub,".csv"),quote=F)
		}
		
		for(MSigDB in c("h.all","c2.cgp","c3.all","c6.all","c7.all"))
		{
			gmt <- read.gmt(paste0("../data/",MSigDB,".v7.1.entrez.gmt"))

			egmt <- enricher(symbol2id[,2], TERM2GENE=gmt)
			
			for(xx in 1:dim(egmt)[1])
			{
				geneIDs <- unlist(strsplit(egmt[xx,"geneID"],"/"))
				geneSymbols <- symbol2id[match(geneIDs,symbol2id[,2]),1]
				geneIDs_geneSymbols <- paste0(geneIDs,"_",geneSymbols)
				
				egmt@result$geneID[xx] <- paste0(geneIDs_geneSymbols,collapse="/")
			}
		
			write.csv(egmt,paste0(markerPath.s,"cluster_",i,"_marker_",pn,"_enrichMSigDB_",MSigDB,".csv"),quote=F)
		}
		
	}
}

summaryTable <- data.frame()
for(i in clusters)
{
	markerPath.s <- paste0(markerPath,"cluster_",i,"/")
	dir.create(markerPath.s)
	
	markers <- FindMarkers(sc.integrated, ident.1 = i)
	markers_pos <- markers[markers[,2]>0,]
	markers_neg <- markers[markers[,2]<0,]
	
	write.csv(markers,paste0(markerPath.s,"cluster_",i,"_marker.csv"),quote=F)
	
	topMarkers <- rownames(markers_pos)[1:15]
	
	g <- FeaturePlot(sc.integrated, topMarkers, reduction = "tsne", ncol=5)
	ggplot2::ggsave(paste0(markerPath.s,"cluster_",i,"_marker_pos_top.jpg"), g, width = 50, height = 30, dpi=200, units = "cm")

	markerList <- list('pos'=rownames(markers_pos),'neg'=rownames(markers_neg))

	for(pn in c("pos"))
	{
		symbol2id <- bitr(markerList[[pn]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

		for(MSigDB in c("h.all"))
		{
			gmt <- read.gmt(paste0("../data/",MSigDB,".v7.1.entrez.gmt"))

			egmt <- enricher(symbol2id[,2], TERM2GENE=gmt, pvalueCutoff = 2, qvalueCutoff = 2, minGSSize = 1,maxGSSize = 10000)
			
			for(xx in 1:dim(egmt)[1])
			{
				geneIDs <- unlist(strsplit(egmt[xx,"geneID"],"/"))
				geneSymbols <- symbol2id[match(geneIDs,symbol2id[,2]),1]
				geneIDs_geneSymbols <- paste0(geneIDs,"_",geneSymbols)
				
				egmt@result$geneID[xx] <- paste0(geneIDs_geneSymbols,collapse="/")
			}
			
			summaryTable[egmt[,"ID"],as.character(i)] <- -log10(egmt[,"qvalue"])
		}
		
	}
}

summaryTable <- as.matrix(summaryTable)
rownames(summaryTable) <- gsub("HALLMARK_","",rownames(summaryTable))
summaryTable[is.na(summaryTable)] <- 0


library(ComplexHeatmap)
tiff(paste0(markerPath,"cluster_marker_",pn,"_enrichMSigDB_",MSigDB,".tif"),width=8,height=12,units="in",res=300)
Heatmap(
    summaryTable, 
    
    col=c("blue","white","yellow","red"),
    
    name = "-Log10(qvalue)", 
    cluster_rows = TRUE,
    cluster_columns = TRUE,
        
    show_row_names = TRUE,
    show_column_names = TRUE,
    
    column_names_rot = 0
)
dev.off()



g <- FeaturePlot(sc.integrated, c("CD44"), reduction = "tsne")
ggplot2::ggsave(paste0(markerPath,"iMarkers_stem.jpg"), g, width = 10, height = 10, dpi=200, units = "cm")

g <- VlnPlot(sc.integrated, c("CD44"))
ggplot2::ggsave(paste0(markerPath,"iMarkers_stem_vln.jpg"), g, width = 25, height = 10, dpi=200, units = "cm")











groups <- as.character(levels(sc.integrated @ meta.data $ manual_groups))

for(i in groups)
{
	markerPath.s <- paste0(markerPath,"group_",i,"/")
	dir.create(markerPath.s)
	
	markers <- FindMarkers(sc.integrated, ident.1 = i, group.by = 'manual_groups')
	markers_pos <- markers[markers[,2]>0,]
	markers_neg <- markers[markers[,2]<0,]
	
	write.csv(markers,paste0(markerPath.s,"group_",i,"_marker.csv"),quote=F)
	
	topMarkers <- rownames(markers_pos)[1:15]
	
	g <- FeaturePlot(sc.integrated, topMarkers, reduction = "tsne", ncol=5)
	ggplot2::ggsave(paste0(markerPath.s,"group_",i,"_markers_pos_top.jpg"), g, width = 50, height = 30, dpi=200, units = "cm")
	
	markerList <- list('pos'=rownames(markers_pos),'neg'=rownames(markers_neg))

	for(pn in c("pos","neg"))
	{
		symbol2id <- bitr(markerList[[pn]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		
		kk <- enrichKEGG(
			gene = symbol2id[,2], 
			organism = 'hsa', 
			qvalueCutoff = 0.05)
		
		for(xx in 1:dim(kk)[1])
		{
			geneIDs <- unlist(strsplit(kk[xx,"geneID"],"/"))
			geneSymbols <- symbol2id[match(geneIDs,symbol2id[,2]),1]
			geneIDs_geneSymbols <- paste0(geneIDs,"_",geneSymbols)
			
			kk@result$geneID[xx] <- paste0(geneIDs_geneSymbols,collapse="/")
		}
		
		write.csv(kk,paste0(markerPath.s,"group_",i,"_marker_",pn,"_enrichKEGG.csv"),quote=F)
		
		for(gosub in c("BP","MF","CC"))
		{
			ego <- enrichGO(
				gene = symbol2id[,2],
				OrgDb = org.Hs.eg.db,
				ont = gosub,
				qvalueCutoff = 0.05,
			   readable = TRUE)
		
			write.csv(ego,paste0(markerPath.s,"group_",i,"_marker_",pn,"_enrichGO_",gosub,".csv"),quote=F)
		}
		
		for(MSigDB in c("h.all","c2.cgp","c3.all","c6.all","c7.all"))
		{
			gmt <- read.gmt(paste0("../data/",MSigDB,".v7.1.entrez.gmt"))

			egmt <- enricher(symbol2id[,2], TERM2GENE=gmt)
			
			for(xx in 1:dim(egmt)[1])
			{
				geneIDs <- unlist(strsplit(egmt[xx,"geneID"],"/"))
				geneSymbols <- symbol2id[match(geneIDs,symbol2id[,2]),1]
				geneIDs_geneSymbols <- paste0(geneIDs,"_",geneSymbols)
				
				egmt@result$geneID[xx] <- paste0(geneIDs_geneSymbols,collapse="/")
			}
		
			write.csv(egmt,paste0(markerPath.s,"group_",i,"_marker_",pn,"_enrichMSigDB_",MSigDB,".csv"),quote=F)
		}
	}
}


