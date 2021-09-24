gse <- getGEO("GSE151229", GSEMatrix = TRUE)
show(gse)
dim(pData(gse[1]))
head(pData(gse[[1]])[, 1:4])
