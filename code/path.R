library(Seurat)

batches <- c("H0","H6","H48","D77","D145m")

dataPath <- "../data/5"

outputPath <- "../results/5TimePoints/"
dir.create(outputPath)

QCPath <- paste0(outputPath,"QC/")
dir.create(QCPath)

uncorrectedPath <- paste0(outputPath,"uncorrected/")
dir.create(uncorrectedPath)

integrationPath <- paste0(outputPath,"integration/")
dir.create(integrationPath)

singlePath <- paste0(outputPath,"single/")
dir.create(singlePath)

dimensionReductionPath <- paste0(integrationPath,"dimensionReduction/")
dir.create(dimensionReductionPath)

clusterPath <- paste0(integrationPath,"cluster/")
dir.create(clusterPath)

markerPath <- paste0(integrationPath,"marker/")
dir.create(markerPath)


pcsTop <- 50
