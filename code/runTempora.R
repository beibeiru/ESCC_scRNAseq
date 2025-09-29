source("path.R")
library(Seurat)
library(Tempora)

sc.integrated <- readRDS(file = paste0(dataPath,"sc.integrated.rds"))

sc_tempora <- ImportSeuratObject(
	sc.integrated, 
	assayType = "RNA",
	clusters = "seurat_clusters",
	timepoints = "orig.ident",
	timepoint_order = c("H0", "H6", "H48", "D77", "D145m")
	)

sc_tempora@ cluster.metadata$ label <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11")

sc_tempora <- CalculatePWProfiles(
	sc_tempora, 
	#gmt_path = "Human_GOBP_AllPathways_no_GO_iea_September_01_2021_symbol.gmt",
	gmt_path = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Pathways/Human_Reactome_September_01_2021_symbol.gmt",
	method="gsva", min.sz = 5, max.sz = 1000, parallel.sz = 1
	)
                
sc_tempora <- BuildTrajectory(
	sc_tempora, 
	n_pcs = 10, difference_threshold = 0.1, loadings = 0.5
	)	
	
	
	
	library(igraph)
	
	sc_tempora -> object
	edge_graph <- igraph::graph_from_data_frame(d=object@trajectory, vertices = object@cluster.metadata, directed = T)
    l <- igraph::layout_with_sugiyama(edge_graph, layers = object@cluster.metadata$Cluster_time_score, maxiter = 1000)

	pdf(paste0("xxxx.pdf"))
	
	colours <- c("green","skyblue","grey","orange","red")
 	plot.igraph(edge_graph, ylim=c(-1,1), ylab = "Inferred time", layout = l$layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(object@meta.data$Timepoints)))), vertex.frame.color="white",
                   vertex.label.color="black", edge.lty = E(edge_graph)$type)
    legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border = "black")
    axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
      
    dev.off()