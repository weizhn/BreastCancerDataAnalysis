library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(tidyr)
library(gghalves)
library(tidyverse)
sample = ""
rds = "tissue.symbol.rds"
print(sample)
print(rds)
out = ""

print(paste("SampleID： ",sample ,sep = " "))
wk = paste0(out,sample,"/")
if(!dir.exists(wk)){
dir.create(wk)
print("0")    
} else{
print("1")    
}
setwd(wk)
st <- readRDS(rds)    
reference <- readRDS('reference.rds')
print(table(reference@cell_types))
coords <- data.frame(x = st@meta.data$x,y = st@meta.data$y)
rownames(coords) = rownames(st@meta.data)
names(coords) <- c('x', 'y')
counts_matrix <- GetAssayData(st, assay = "Spatial", slot = "counts")
puck <- SpatialRNA(coords=coords, counts=counts_matrix)
myRCTD <- create.RCTD(puck, reference, counts_MIN =0,UMI_min=0,max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
save(myRCTD,file= paste0(wk,sample,"_cellbin_myRCTD.Rdata"))
results <- myRCTD@results
norm_weights <- normalize_weights(results$weights)
saveRDS(norm_weights, file = paste0(sample, "_cellbin.norm_weights.rds"))