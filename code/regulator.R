source("path.R")
library(SCENIC)

sc.integrated <- readRDS(file = paste0(dataPath,"sc.integrated.rds"))

DefaultAssay(sc.integrated) <- "RNA"

cellInfo <- data.frame(seuratCluster=Idents(sc.integrated))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)


  setwd('/home/data/results/workspace/nCoV/SCENIC')
  regulonAUC <- importAUCfromText(file.path("auc_mtx.csv"))
  regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
  macro.cellinfo = read.delim2("/home/data/results/workspace/nCoV/SCENIC/3-macrophage-SCENIC-meta.csv",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
  rownames(macro.cellinfo) = macro.cellinfo$ID
  
  regulonActivity_byCellType <- sapply(split(rownames(macro.cellinfo), macro.cellinfo$celltype),
                                       function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
  
  write.table(regulonActivity_byCellType_Scaled,file='3-SCENIC-heatmap.txt',row.names = TRUE,quote = FALSE,sep='\t')
  
  pdf(file="3-SCENIC-heatmap.pdf", width = 6, height = 13)
  pheatmap::pheatmap(regulonActivity_byCellType_Scaled, fontsize_row = 5, fontsize_col = 12, 
                     color=colorRampPalette(c("#00a9ff","white","#F8766D"))(100), breaks=seq(-3, 3, length.out = 100),
                     treeheight_row=10, treeheight_col=10, border_color=NA)
  dev.off()