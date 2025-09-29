
dataPath <- "/data/rub2/project/SpaCE/results/ST_data/TCGA/"

cancers <- list.files(paste0(dataPath,"mRNAseq_gene_exp_symbol"))
	
read.gene.exp <- function(cancer)
{
	gene_exp_path <- paste0(dataPath,"mRNAseq_gene_exp_symbol/")

	gene_exp_link <- paste0(gene_exp_path,cancer)
	gene_exp <- read.csv(gene_exp_link,as.is=T,row.names=1,sep="\t")
	
	gene_exp
}

extract.samples <- function(data_matrix,type)
{
	sample_type <- as.numeric(substr(colnames(data_matrix),14,15))
	if(type=="T")
	{
		TS <- sample_type>=1&sample_type<=9
		new_matrix <- data_matrix[,TS]
	}
	if(type=="N")
	{
		NS <- sample_type>=10&sample_type<=19
		new_matrix <- data_matrix[,NS,drop=F]
	}
	return(new_matrix)
}

extract.patients <- function(data_matrix)
{
	data_matrix <- data_matrix[,sort(colnames(data_matrix)),drop=F]
	data_matrix <- data_matrix[,!duplicated(substr(colnames(data_matrix),1,12)),drop=F]
	colnames(data_matrix) <- substr(colnames(data_matrix),1,12)
	return(data_matrix)
}

data_version <- "20160128"

read.clinical <- function(cancer)
{
	clinical_start  <- paste0(dataPath,"/clinical/gdac.broadinstitute.org_")
	clinical_middle <- paste0(".Merge_Clinical.Level_1.",data_version,"00.0.0/")
	clinical_end    <- ".merged_only_clinical_clin_format.txt"
	
	clinical_link <- paste0(clinical_start,cancer,clinical_middle,cancer,clinical_end)
	clinical <- read.csv(clinical_link,header=T,row.names=1,sep="\t",as.is=T)
	names(clinical) <- gsub("-",".",toupper(clinical["patient.bcr_patient_barcode",]),fixed=T)
	clinical
}

generate.clinical.surv.matrix <- function(clinical)
{
	ind_keep <- grep("vital_status",rownames(clinical))
	vital <- as.matrix(clinical[ind_keep,])
	vital[which(vital=="")] <- NA
	vital_status <- c()
	for (j in 1:dim(vital)[2])
	{
		if(sum(is.na(vital[,j])) < dim(vital)[1])
		{
			vital_status <- c(vital_status,ifelse("dead" %in% vital[,j],1,0))
		} else {
			vital_status <- c(vital_status,NA)
		}
	}
	
	ind_keep <- grep("days_to_death",rownames(clinical))
	death <- data.matrix(clinical[ind_keep,])
	days_to_death <- c()
	for (j in 1:dim(death)[2])
	{
		if(sum(is.na(death[,j])) < dim(death)[1])
		{
			days_to_death <- c(days_to_death,max(as.numeric(death[,j]),na.rm=T))
		} else {
			days_to_death <- c(days_to_death,NA)
		}
	}
	
	ind_keep <- grep("days_to_last_followup",rownames(clinical))
	followup <- data.matrix(clinical[ind_keep,])
	days_to_last_followup <- c()
	for (j in 1:dim(followup)[2])
	{
		if( clinical["patient.follow_ups.follow_up.vital_status",j] %in% c("dead", "deceased") |
			clinical["patient.vital_status",j] %in% c("dead", "deceased") )
		{
			days_to_last_followup <- c(days_to_last_followup,NA)
		} else {
			if(sum(is.na(followup[,j])) < dim(followup)[1])
			{
				days_to_last_followup <- c(days_to_last_followup,max(as.numeric(followup[,j]),na.rm=T))
			} else {
				days_to_last_followup <- c(days_to_last_followup,NA)
			}
		}
	}

	clinical_surv <- rbind(vital_status,days_to_death,days_to_last_followup)
	colnames(clinical_surv) <- names(clinical)
	clinical_surv <- clinical_surv[,!is.na(clinical_surv[1,])]
	clinical_surv <- clinical_surv[,!(clinical_surv[1,]==0&!is.na(clinical_surv[2,]))]
	
	days_to_death_or_fup <- c()
	for (j in 1:dim(clinical_surv)[2])
	{
		days_to_death_or_fup <- c(days_to_death_or_fup,ifelse(clinical_surv[1,j]==1,clinical_surv[2,j],clinical_surv[3,j]))
	}
	
	clinical_surv <- rbind(clinical_surv,days_to_death_or_fup)
	clinical_surv <- clinical_surv[,!is.na(clinical_surv[4,])]
	clinical_surv <- clinical_surv[,clinical_surv[4,]>0]

	clinical_surv
}

source("path.R")


gmt <- list()
for(i in c(0:11))
{
	markerPath.s <- paste0(markerPath,"cluster_",i,"/")
	
	markers <- read.csv(paste0(markerPath.s,"cluster_",i,"_marker.csv"),row.names=1,as.is=T)
	markers_pos <- markers[markers[,2]>0,]
	
	gmt[[as.character(i)]] <- rownames(markers_pos)
}



cancers      <- c("ACC","BLCA","BRCA","BRCA.LumA","BRCA.Basal","BRCA.Normal","BRCA.LumB","BRCA.Her2","CESC","CHOL",
                  "COAD","ESCA",
                  "GBM","HNSC","HNSC.HPVpos","HNSC.HPVneg","KICH",
                  "KIRC","KIRP","LGG","LIHC",
                  "LUAD","LUSC","MESO","OV","PAAD",
                  "PCPG","PRAD","READ","SARC","SKCM","SKCM.primary","SKCM.metastatic",
                  "STAD","THCA",
                  "UCEC","UCS","UVM")

cancers      <- c("ACC","BLCA","BRCA","CESC",
                  "COAD","ESCA",
                  "GBM","HNSC","KICH",
                  "KIRC","KIRP","LGG","LIHC",
                  "LUAD","LUSC","MESO","OV","PAAD",
                  "PRAD","READ","SARC","SKCM",
                  "STAD","THCA",
                  "UCEC","UCS","UVM")
                  
cancers      <- c("BLCA","BRCA","CESC",
                  "COAD","ESCA",
                  "GBM","HNSC","KICH",
                  "KIRC","KIRP","LGG","LIHC",
                  "LUAD","LUSC","OV","PAAD",
                  "PRAD","READ","SARC","SKCM",
                  "STAD","THCA",
                  "UCEC")

survival_table <- matrix(NA,ncol=12,nrow=length(cancers) )
colnames(survival_table) <- as.character( c(0:11) )
rownames(survival_table) <- cancers


for(cancer in cancers)
{		
		if(cancer%in%c("HNSC.HPVpos","HNSC.HPVneg"))
		{
			cdata <- read.gene.exp("HNSC")
		}else if(cancer%in%c("SKCM.primary","SKCM.metastatic")){
			cdata <- read.gene.exp("SKCM")
			if(cancer=="SKCM.primary"){
				cdata <- cdata[,!grepl(".06",colnames(cdata),fixed=TRUE)]
			}else{
				cdata <- cdata[, grepl(".06",colnames(cdata),fixed=TRUE)]
			}
		}else if(cancer%in%c("BRCA.LumA","BRCA.Basal","BRCA.Normal","BRCA.LumB","BRCA.Her2")){
			cdata <- read.gene.exp("BRCA")
		}else{
			cdata <- read.gene.exp(cancer)
		}
		
		cdata_T <- extract.samples(cdata,"T")
		cdata_P <- extract.patients(cdata_T)
		
		if(cancer%in%c("HNSC.HPVpos","HNSC.HPVneg"))
		{
			clinical <- read.clinical("HNSC")
			clinical <- clinical[,!(clinical["patient.hpv_status_by_ish_testing",]%in%"positive"&
				clinical["patient.hpv_status_by_p16_testing",]%in%"negative")]
			clinical <- clinical[,!(clinical["patient.hpv_status_by_ish_testing",]%in%"negative"&
				clinical["patient.hpv_status_by_p16_testing",]%in%"positive")]
				
			if(cancer=="HNSC.HPVpos"){
				clinical <- clinical[,(clinical["patient.hpv_status_by_ish_testing",]%in%"positive"|
					clinical["patient.hpv_status_by_p16_testing",]%in%"positive")]
			}else{
				clinical <- clinical[,(clinical["patient.hpv_status_by_ish_testing",]%in%"negative"|
					clinical["patient.hpv_status_by_p16_testing",]%in%"negative")]
			}
		}else if(cancer%in%c("SKCM.primary","SKCM.metastatic")){
			clinical <- read.clinical("SKCM")
		}else if(cancer%in%c("BRCA.LumA","BRCA.Basal","BRCA.Normal","BRCA.LumB","BRCA.Her2")){
			clinical <- read.clinical("BRCA")
			patientsSub <- rownames(subtypes.BRCA)[subtypes.BRCA[,"TCGA.Subtype"]%in%cancer]
			clinical <- clinical[,colnames(clinical)%in%gsub("-",".",patientsSub)]
		}else{
			clinical <- read.clinical(cancer)
		}
		
		clinical_surv <- generate.clinical.surv.matrix(clinical)
		clinical_surv <- rbind(clinical_surv,Age=unlist(clinical["patient.age_at_initial_pathologic_diagnosis",colnames(clinical_surv)]))
		clinical_surv <- rbind(clinical_surv,Gender=unlist(clinical["patient.gender",colnames(clinical_surv)]))

		olp_patients <- intersect(colnames(cdata_P),colnames(clinical_surv))
		cdata_P_olp <- cdata_P[,olp_patients]
		clinical_surv_olp <- clinical_surv[,olp_patients]


		gsvaScore <- GSVA::gsva(as.matrix(cdata_P_olp),gmt)

	
		for(item in rownames(gsvaScore))
		{
			neww <- as.data.frame(cbind(t(clinical_surv_olp),mes=gsvaScore[item,]))
			neww[,1] <- as.numeric(as.character(neww[,1]))
			neww[,2] <- as.numeric(as.character(neww[,2]))
			neww[,3] <- as.numeric(as.character(neww[,3]))
			neww[,4] <- as.numeric(as.character(neww[,4]))
			neww[,5] <- as.numeric(as.character(neww[,5]))
			neww[,7] <- as.numeric(as.character(neww[,7]))
			
			library(survival)
			
			if(sum(is.na(neww[,5]))==dim(neww)[1]&sum(is.na(neww[,6]))==dim(neww)[1])
			{
				coxmodel_fit <- coxph(Surv(days_to_death_or_fup, vital_status) ~ mes, data = neww)
			}
			
			if(sum(is.na(neww[,5]))!=dim(neww)[1]&sum(is.na(neww[,6]))==dim(neww)[1])
			{
				coxmodel_fit <- coxph(Surv(days_to_death_or_fup, vital_status) ~ mes + Age, data = neww)
			}
			
			if(sum(is.na(neww[,5]))==dim(neww)[1]&sum(is.na(neww[,6]))!=dim(neww)[1])
			{
				if(length(unique(neww[,6]))==2)
				{	
					coxmodel_fit <- coxph(Surv(days_to_death_or_fup, vital_status) ~ mes + Gender, data = neww)
				}else{
					coxmodel_fit <- coxph(Surv(days_to_death_or_fup, vital_status) ~ mes, data = neww)
				}
			}
			
			if(sum(is.na(neww[,5]))!=dim(neww)[1]&sum(is.na(neww[,6]))!=dim(neww)[1])
			{
				if(length(unique(neww[,6]))==2)
				{	
					coxmodel_fit <- coxph(Surv(days_to_death_or_fup, vital_status) ~ mes + Age + Gender, data = neww)
				}else{
					coxmodel_fit <- coxph(Surv(days_to_death_or_fup, vital_status) ~ mes + Age , data = neww)
				}
			}
			
			coxmodel_obj <- summary(coxmodel_fit);
			
			survival_table[cancer,item] <- coxmodel_obj$coefficients["mes",4]
			

		}


}

survival_table <- as.data.frame(survival_table)

xx1 <- rownames(survival_table)

xx1[xx1=="BLCA"] <- "Bladder"
xx1[xx1=="BRCA"] <- "Breast"
xx1[xx1=="CESC"] <- "Cervical and endocervical"
xx1[xx1=="COAD"] <- "Colon"
xx1[xx1=="ESCA"] <- "Esophageal"
xx1[xx1=="GBM"] <- "Glioblastoma"
xx1[xx1=="HNSC"] <- "Head and Neck"
xx1[xx1=="KICH"] <- "Kidney chromophobe"
xx1[xx1=="KIRC"] <- "Kidney renal clear"
xx1[xx1=="KIRP"] <- "Kidney renal papillary"
xx1[xx1=="LGG"] <- "Glioma"
xx1[xx1=="LIHC"] <- "Liver"
xx1[xx1=="LUAD"] <- "Lung adeno."
xx1[xx1=="LUSC"] <- "Lung squamous"
xx1[xx1=="OV"] <- "Ovarian"
xx1[xx1=="PAAD"] <- "Pancreatic"
xx1[xx1=="PRAD"] <- "Prostate"
xx1[xx1=="READ"] <- "Rectum"
xx1[xx1=="SARC"] <- "Sarcoma"
xx1[xx1=="SKCM"] <- "Melanoma"
xx1[xx1=="STAD"] <- "Stomach"
xx1[xx1=="THCA"] <- "Thyroid"
xx1[xx1=="UCEC"] <- "Endometrial"

xx1 -> rownames(survival_table)



library(ComplexHeatmap)
library(circlize)
tiff(paste0(markerPath,"cluster_marker_pos_survival.tif"),width=8,height=12,units="in",res=300)
Heatmap(
    survival_table, 
    
	col=colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
	 
    name = "Z-score", 
    cluster_rows = TRUE,
    cluster_columns = TRUE,
        
    show_row_names = TRUE,
    show_column_names = TRUE,
    
    column_names_rot = 0
)
dev.off()



