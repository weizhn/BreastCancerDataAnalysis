library(patchwork)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(vctrs)
library(cluster)
library(clustree)
library(harmony)

setwd("")
path <- 'DATA/'
samples <- c('A3787T3','A3788T2','A3789T2','A3790T3','A3791T3','A3792T3','A3793T4','A3794T4','A4108T2','A4109T4')
second <- '/filter'
for (sam in samples) {
    print(sam)
    data.dir <- paste0(path, sam, second)
    exp <- Read10X(
        data.dir = data.dir,
        gene.column = 2,
        cell.column = 1,
        unique.features = TRUE,
        strip.suffix = FALSE
    )
    rownames(exp) <- gsub("_", "-", rownames(exp))
    object <- CreateSeuratObject(
        counts = exp,
        project = sam,
        min.cells = 10,
        min.features = 200
    )
    Idents(object) <- sam
    object@meta.data$Sample <- sam
    save(object, file = paste0(sam, '_ob.rda')) 
}

setwd("")
load('A3787T3_ob.rda')
A3787T3 <- object
A3787T3$Type <- 'T3'
A3787T3$Sample <- 'A3787T3'

load('A3788T2_ob.rda')
A3788T2 <- object
A3788T2$Type <- 'T2'
A3788T2$Sample <- 'A3788T2'

load('A3789T2_ob.rda')
A3789T2 <- object
A3789T2$Type <- 'T2'
A3789T2$Sample <- 'A3789T2'

load('A3790T3_ob.rda')
A3790T3 <- object
A3790T3$Type <- 'T3'
A3790T3$Sample <- 'A3790T3'

load('A3791T3_ob.rda')
A3791T3 <- object
A3791T3$Type <- 'T3'
A3791T3$Sample <- 'A3791T3'

load('A3792T3_ob.rda')
A3792T3 <- object
A3792T3$Type <- 'T3'
A3792T3$Sample <- 'A3792T3'

load('A3793T4_ob.rda')
A3793T4 <- object
A3793T4$Type <- 'T4'
A3793T4$Sample <- 'A3793T4'

load('A3794T4_ob.rda')
A3794T4 <- object
A3794T4$Type <- 'T4'
A3794T4$Sample <- 'A3794T4'

load('A4108T2_ob.rda')
A4108T2 <- object
A4108T2$Type <- 'T2'
A4108T2$Sample <- 'A4108T2'

load('A4109T4_ob.rda')
A4109T4 <- object
A4109T4$Type <- 'T4'
A4109T4$Sample <- 'A4109T4'

filtob <- function(object,nfeamin,nfeamax,nCountmin,nCountmax,permt,Hb,Rpsl){
    obnew <- PercentageFeatureSet(object ,"^MT-", col.name = "percent.mt")
    obnew <- PercentageFeatureSet(obnew ,"^HB", col.name = "percent.Hb")
    obnew <- PercentageFeatureSet(obnew ,"^RP[SL]", col.name = "percent.Rpsl")
    obnew <- obnew[,obnew$nFeature_RNA > nfeamin & 
                    obnew$nFeature_RNA < nfeamax &
                    obnew$nCount_RNA > nCountmin &
                    obnew$nCount_RNA < nCountmax &
                    obnew$percent.mt < permt & 
                    obnew$percent.Hb < Hb &
                    obnew$percent.Rpsl < Rpsl]
    obnew <- SCTransform(obnew, assay = "RNA", verbose = TRUE,vars.to.regress='percent.mt',vst.flavor = 'v1')
    sampName <- obnew$Sample 
    newname <- paste0(sampName,'_',colnames(obnew))
    obnew <- RenameCells(obnew,new.names=newname)   
    return(obnew)
}

A3787T3 <- filtob(A3787T3,200,7500,200,40000,30,15,15)
A3788T2 <- filtob(A3788T2,200,7500,200,40000,30,15,15)
A3789T2 <- filtob(A3789T2,200,7500,200,40000,30,15,15)
A3790T3 <- filtob(A3790T3,200,7500,200,40000,30,15,15)
A3791T3 <- filtob(A3791T3,200,7500,200,40000,30,15,15)
A3792T3 <- filtob(A3792T3,200,7500,200,40000,30,15,15)
A3793T4 <- filtob(A3793T4,200,7500,200,40000,30,15,15)
A3794T4 <- filtob(A3794T4,200,7500,200,40000,30,15,15)
A4108T2 <- filtob(A4108T2,200,7500,200,40000,30,15,15)
A4109T4 <- filtob(A4109T4,200,7500,200,40000,30,15,15)
save(A3787T3,A3788T2,A3789T2,A3790T3,A3791T3,A3792T3,A3793T4,A3794T4,A4108T2,A4109T4,file='filtob.rda')

object_list <- list(A3787T3,A3788T2,A3789T2,A3790T3,A3791T3,A3792T3,A3793T4,A3794T4,A4108T2,A4109T4)
CANmerge <- merge(object_list[[1]], y = object_list[2:10], 
                  add.cell.ids = c('A3787T3','A3788T2','A3789T2','A3790T3','A3791T3','A3792T3','A3793T4','A3794T4','A4108T2','A4109T4'))

setwd('')
gene_counts_per_cell <- CANmerge$nFeature_RNA
gene_counts_summary <- summary(gene_counts_per_cell)
gene_counts_summary.show <- data.frame(Quantile=names(gene_counts_summary),
                                       summary=as.numeric(gene_counts_summary))
write.table(gene_counts_summary.show,file='CANgene_counts_summary.show',sep = '\t',quote = F,row.names = F)

data <- data.frame(gene_counts_per_cell)
colnames(data) <- 'gene_counts_per_cell'
data$Cell <- rownames(data)
png('CANCellCount.png',width=300,height=200,unit='mm',res=300)
ggplot(data,aes(gene_counts_per_cell))+
       geom_histogram(binwidth = 50,fill='blue',color='blue')+
       ggtitle("Distribution of Gene counts")+
       xlab("Gene counts")+ylab("Cell number")+theme_bw()+
       theme(axis.text=element_text(size=30),
             title=element_text(size=30),
             axis.title=element_text(size=30))+
       theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()

Idents(CANmerge) <- CANmerge$Sample
p <- VlnPlot(CANmerge, features = c("nFeature_RNA", "nCount_RNA",
                                "percent.Hb",'percent.mt','percent.Rpsl'), 
                                ncol = 3,pt.size = 0.1,raster = FALSE)
p1 <- lapply(list(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]]),function(x){
      return(x+theme(axis.title.x=element_blank(),
                    axis.text=element_text(size=20),
                    title=element_text(size=20)))
})        
p2 <- plot_grid(plotlist=p1,ncol=3)
ggsave(p2,filename='merge_VlnPlotForHbMt_CAN.pdf',width=15,height=12)

Idents(CANmerge) <- CANmerge$Sample
pcor <- FeatureScatter(CANmerge,feature1='nFeature_RNA',
                    feature2 = 'nCount_RNA',raster = FALSE)+
                    theme(axis.title.x=element_blank(),
                    axis.text=element_text(size=20),
                    axis.text.x=element_text(angle=45,hjust=1,vjust=1),
                    title=element_text(size=20),
                    legend.title=element_text(size=20),
                    legend.text=element_text(size=20))
ggsave(pcor,filename='pcorCAN.pdf',width=6.5,height=6)

DefaultAssay(CANmerge) <- 'RNA'
CANmerge <- NormalizeData(CANmerge)
CANmerge <- FindVariableFeatures(CANmerge,nfeatures=3000)
VariableFeatures <- VariableFeatures(CANmerge)[1:3000]
CANmerge <- ScaleData(CANmerge, features = VariableFeatures,vars.to.regress = c("percent.mt", "percent.Rpsl", "percent.Hb"))
CANmerge <- RunPCA(CANmerge)
save(CANmerge,file='CANmerge.rda')

setwd('/public/home/weizhn/Guochunming/scRNA_Progess_0826_DLMTRRNA')
load('/public/home/weizhn/Guochunming/scRNA_Progess_0826_DLMTRRNA/CANmerge.rda')
options(future.globals.maxSize = 1000000 * 1024^2)
CANInte <- IntegrateLayers(
  object = CANmerge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE,kmeans_init_nstart=20, kmeans_init_iter_max=2000
)
print('Joint')
CANInte[["RNA"]] <- JoinLayers(CANInte[["RNA"]])
save(CANInte,file='CANInteJoint.rda')

data <- data.frame(Matrix::colSums(CANInte))
colnames(data) <- 'library_size'
data$Cell <- rownames(data)
png('distri_lib_size_CAN.png',width=300,height=200,unit='mm',res=300)
ggplot(data,aes(library_size))+
       geom_histogram(binwidth = 50,fill='blue',color='blue')+
       ggtitle("Distribution of library size")+
       xlab("library size")+ylab("Cell number")+theme_bw()+
       theme(axis.text=element_text(size=30),
             title=element_text(size=30),
             axis.title=element_text(size=30))+
       theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()

CANfinal <- FindNeighbors(CANInte, reduction = 'harmony', dims = 1:30)
CANfinal <- FindClusters(CANfinal, resolution = seq(0.1,1,0.1))
CANfinal <- RunUMAP(CANfinal, reduction = 'harmony', dims = 1:30, reduction.name = 'harmony_umap')
CANfinal <- RunTSNE(CANfinal, reduction = 'harmony', dims = 1:30, reduction.name = 'harmony_tsne')
save(CANfinal,file='CANfinal_harmony.rda')
CANmeta <- CANfinal@meta.data
CANmeta$Cell <- rownames(CANmeta)
CANexp <- CANfinal@assays$RNA$data

setwd('')
myreso <- seq(0.1,1,0.1)
for (reso in myreso){
    cluster.name <- paste0('RNA_snn_res.',reso)
    p1 <- DimPlot(
    CANfinal,label=TRUE,repel=FALSE,raster=FALSE,
    reduction = 'harmony_tsne',
    group.by = c("Sample", cluster.name),
    combine = TRUE, label.size = 2
    )
    ggsave(p1,filename=paste0('reso',reso,'_inteTSNE.png'),width=15,height=6)

    p2 <- DimPlot(
    CANfinal,label=TRUE,repel=FALSE,raster=FALSE,
    reduction = 'harmony_umap',
    group.by = c("Sample", cluster.name),
    combine = TRUE, label.size = 2
    )
    ggsave(p2,filename=paste0('reso',reso,'_inteUMAP.png'),width=15,height=6)
}


library(Seurat)
load('/public/home/weizhn/Guochunming/scRNA_Progess_0826_DLMTRRNA/CANfinal_harmony.rda')
load('/public/home/weizhn/Guochunming/scRNA_Progess_0826_DLMTRRNA/new_UPDATE/CANfinal_UpdateType.rda')
data <- CANfinal
genes_to_plot <- c('SPARC', 'COL4A1', 'TAGLN', 'CALD1', 'FBLN1', 'RARRES2', 'CCDC80', 'SFRP2','VCAN')
output_dir <- ""
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

for (gene in genes_to_plot) {
  if (gene %in% rownames(GetAssayData(data))) {
    p <- FeaturePlot(data, features = gene, reduction = "harmony_umap",raster=FALSE) +
      scale_color_gradient(low = "lightgrey", high = "blue") +
      theme_minimal() +
      ggtitle(paste0(gene, " Expression"))

    ggsave(
      filename = file.path(output_dir, paste0("umap_", gene, ".png")),
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )

    umap_coords <- Embeddings(data, "harmony_umap")
    gene_expression <- FetchData(data, vars = gene)
    combined_data <- cbind(umap_coords, gene_expression)
    write.csv(
      combined_data,
      file.path(output_dir, paste0("umap_", gene, "_expression.csv")),
      row.names = TRUE
    )
    
    message(paste0("succeed ", gene, " result"))
  } else {
    message(paste0("no ", gene, "exist"))
  }
}

library(Seurat)
library(ggplot2)

load('CANfinal_harmony.rda')
data <- CANfinal

genes_to_plot <- c('CD247','ENSG00000290067','ITK','BCL11B','TXK','LINC00402','LINC00861','LINC01550','CD3E','PRKCQ',
                    'FCRL1','LINC00926','BANK1','FCRL2','CLEC17A','CD19','ENSG00000302378','ENSG00000258760','COL19A1','BLK',
                    'ENSG00000290801','ANKRD30A','ENSG00000286208','IGF1R','ENSG00000310520','ENSG00000278022','ENSG00000304952','LINC01087','ENSG00000307340','CPB1',
                    'POSTN','FN1','COL1A1','FAP','COL1A2','SULF1','BGN','ISLR','VCAN','COL3A1',
                    'ENSG00000253108','SLC22A3','ENSG00000289842','ADRA1A','SLC26A7','CCL21','KEL','ZNF804B','ANKRD29','ENSG00000301611',
                    'MS4A4E','SIGLEC1','FPR3','ITGAX','CSF1R','MS4A6A','CPVL','SLC1A3','MS4A14','TLR2',
                    'CEMIP','PROX1','MARCO','MYCT1','EGFL7','FLRT2','MUC3A','SEMA3D','BMPER','STOX2',
                    'TREML1','LAMP3','ENSG00000231873','CCL22','NCCRP1','WNT5B','CD1E','WFDC21P','ST3GAL6','ENSG00000303924',
                    'AQP1','BTNL9','CD34','PLVAP','FLT1','TPO','MECOM','PCDH17','VWF','ADGRL4',
                    'RGS13','ASPM','RRM2','TOP2A','MKI67','BUB1B','NCAPG','ANLN','ENSG00000225885','EPS15-AS1',
                    'ENSG00000240040','IGKC','IGLC2','IGHA1','JCHAIN','IGLL5','ENSG00000294816','IGHG4','IGHG3','IGHG1',
                    'ENSG00000290592','CLEC4C','ENSG00000291038','LINC01478','ENSG00000288703','CUX2','ENSG00000231272','LINC01226','PACSIN1','ENSG00000293411',
                    'ENSG00000290007','AQP7','TRARG1','PLIN4','ADIPOQ','GPD1','PLIN1','CIDEA','GPAM','SLC19A3')

resolution <- "RNA_snn_res.0.1"

output_dir <- "/public/home/weizhn/Guochunming/scRNA_Progess_0826_DLMTRRNA/touying/T"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!resolution %in% colnames(data@meta.data)) {
  stop(paste("reso", resolution, "error", 
             paste(colnames(data@meta.data), collapse = ", ")))
}

Idents(data) <- data@meta.data[[resolution]]

all_clusters <- levels(Idents(data))
numeric_clusters <- as.numeric(as.character(all_clusters))
ordered_clusters <- as.character(sort(numeric_clusters))

valid_genes <- genes_to_plot[genes_to_plot %in% rownames(GetAssayData(data))]
if (length(valid_genes) == 0) {
  stop("no gene")
}

dot_data <- DotPlot(
  object = data,
  features = valid_genes,
  scale = TRUE,
  scale.by = "radius"
)$data

plot_data <- dot_data %>%
  mutate(
    features.plot = factor(features.plot, levels = valid_genes),
    id = factor(id, levels = ordered_clusters)
  ) %>%
  arrange(id, features.plot)

custom_dot_plot <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(
    low = "lightgrey", 
    high = "blue",
    name = "average\n(Z-score)"
  ) +
  scale_size(
    name = "cell\n%",
    range = c(0, 8), 
    limits = c(0, 100)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    x = "gene",
    y = "Cluster",
    title = paste("express - reso:", resolution)
  ) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  )

pdf_width <- min(30, max(8, length(valid_genes) * 0.85))
pdf_height <- min(15, max(6, length(ordered_clusters) * 0.7))

pdf_file <- file.path(output_dir, paste0("dotplot_resolution_", gsub("\\.", "_", resolution), ".pdf"))
pdf(pdf_file, width = pdf_width, height = pdf_height)
print(custom_dot_plot)
dev.off()

write.csv(
  plot_data,
  file.path(output_dir, paste0("dotplot_data_resolution_", gsub("\\.", "_", resolution), ".csv")),
  row.names = FALSE
)


#allmarer##############
load('CANfinal_harmony.rda')
setwd('')
reso <- 0.1
Idents(CANfinal) <- CANfinal@meta.data[,'RNA_snn_res.0.1']
allmarkers <- FindAllMarkers(CANfinal, only.pos = FALSE, min.pct = 0.2, 
                        logfc.threshold = 0.1,assay='RNA',slot='data')
save(allmarkers,file=paste0('allmarkerharmony.reso',reso,'.rda'))
write.csv(allmarkers,file=paste0('allmarkerharmony.reso',reso,'.csv'),row.names=F,quote=T)
write.table(allmarkers,'allmarkerharmony.txt',sep='\t',row.names=F,quote=F)
allmarkerfilt <- allmarkers[allmarkers$p_val_adj<0.05,]
save(allmarkerfilt,file=paste0('allmarkerfiltharmony.reso',reso,'.rda'))
write.csv(allmarkerfilt,file=paste0('allmarkerfiltharmony.reso',reso,'.csv'),row.names=F,quote=T)
write.table(allmarkerfilt,'allmarkerfiltharmony.txt',sep='\t',row.names=F,quote=F)
table(CANfinal$RNA_snn_res.0.1)

#####################
library(openxlsx)
load('CANfinal_UpdateType.rda')
setwd('')
result_table <- table(CANfinal$Sample, CANfinal$Celltype)
result_df <- as.data.frame(result_table)
wb <- createWorkbook()
addWorksheet(wb, "Cross Table")
writeData(wb, "Cross Table", result_df)
saveWorkbook(wb, "Sample_celltype.xlsx", overwrite = TRUE)

result_table <- table(CANfinal$Type, CANfinal$Celltype)
result_df <- as.data.frame(result_table)
wb <- createWorkbook()
addWorksheet(wb, "Cross Table")
writeData(wb, "Cross Table", result_df)
saveWorkbook(wb, "Type_celltype.xlsx", overwrite = TRUE)


load('CANfinal_harmony.rda')
setwd('')
CANmeta <- CANfinal@meta.data
CANmeta$Celltype <- CANmeta$RNA_snn_res.0.1
CANmeta$Celltype <- gsub('^0$','T/NK cell',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^1$','B cell',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^2$','Epithelial',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^3$','Fib2',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^4$','Macrophage',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^5$','Fib1',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^6$','Plasma cell',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^7$','Epithelial',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^8$','EC',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^9$','LEC',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^10$','cDCs',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^11$','pDCs',CANmeta$Celltype)
CANmeta$Celltype <- gsub('^12$','Adipocyte',CANmeta$Celltype)
CANfinal$Celltype <- CANmeta[colnames(CANfinal),'Celltype']
saveRDS(CANfinal, file = "CANfinal_UpdateType.rds")
save(CANfinal,file='CANfinal_UpdateType.rda')
setwd('')
load('CANfinal_UpdateType.rda')

genes_to_plot <- c('ERBB4','IGF1','IGF1R','NRG1')
p <- VlnPlot(object = CANfinal,features = genes_to_plot,
  group.by = "Celltype", 
  split.by='Type',
  stack = TRUE,
  flip = TRUE,
  cols = c('#39C5BB','#EE0000','#9999FF')) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5),
        axis.title.x=element_blank(),
        strip.text.y = element_text(size = 16,face='italic'))
ggsave(p,filename='cortexGeneheatdot.pdf',width=8,height=10)
##################################################################################
setwd('')
load('CANfinal_UpdateType.rda')
skipped_types <- c()
minpct <- 0.1
logfc <- 0.1
Idents(CANfinal) <- CANfinal$Celltype 
DEGs <- data.frame()

for (type in unique(CANfinal$Celltype)) {
  subset_data <- subset(CANfinal, subset = Celltype == type)

  cell_counts <- table(subset_data$Type)
  wt_count <- cell_counts["T4"]
  cko_count <- cell_counts["T3"]

  if (is.na(wt_count) || is.na(cko_count) || wt_count < 3 || cko_count < 3) {
    message("skip ", type, ": T4=", wt_count, " cells, T3=", cko_count, " cells")
    skipped_types <- c(skipped_types, type)
    next 
  }
  
  DEGtmp <- FindMarkers(
    object = subset_data,
    ident.1 = "T4", 
    ident.2 = "T3",
    group.by = 'Type',
    only.pos = FALSE,
    min.pct = minpct,
    test.use = "wilcox",
    logfc.threshold = logfc,
    assay = 'RNA',
    slot = 'data'
  )

  DEGtmp$Celltype <- type
  DEGtmp$gene <- rownames(DEGtmp)

  DEGs <- rbind(DEGs, DEGtmp)
}

if (length(skipped_types) > 0) {
  message("\n skip:")
  message(paste(skipped_types, collapse = ", "))
} else {
  message("DEG over")
}

save(DEGs, file='DEG_CANfinal_T4T3.rda')
write.table(DEGs, file='DEG_CANfinal_T4T3.txt', row.names=F, quote=F, sep='\t')
DEGfilt <- DEGs[DEGs$p_val_adj < 0.05,]
save(DEGfilt, file='DEGfilt_CANfinal_T4T3.rda')
write.table(DEGfilt, file='DEGfilt_CANfinal_T4T3.txt', row.names=F, quote=F, sep='\t')

DEGsta <- table(DEGfilt$Celltype)
DEGsta <- data.frame(Celltype=names(DEGsta), number=as.vector(DEGsta))
write.table(DEGsta, file='DEGsta_CANfinal_T4T3.txt', row.names=F, quote=F, sep='\t')
##################################################################
setwd('')
load('CANfinal_UpdateType.rda')

celltype1 <- "Fib1"
celltype2 <- "Fib2"

skipped_types <- c()
minpct <- 0.1
logfc <- 0.1

Idents(CANfinal) <- CANfinal$Celltype

cell_counts <- table(CANfinal$Celltype)
count1 <- cell_counts[celltype1]
count2 <- cell_counts[celltype2]

if (is.na(count1) || is.na(count2) || count1 < 3 || count2 < 3) {
  message("skip: ", celltype1, "=", count1, " cells, ", celltype2, "=", count2, " cells")
  skipped_types <- c(celltype1, celltype2)
} else {

  DEGs <- FindMarkers(
    object = CANfinal,
    ident.1 = celltype1, 
    ident.2 = celltype2,
    only.pos = FALSE,
    min.pct = minpct,
    test.use = "wilcox",
    logfc.threshold = logfc,
    assay = 'RNA',
    slot = 'data'
  )

  DEGs$comparison <- paste0(celltype1, "_vs_", celltype2)
  DEGs$gene <- rownames(DEGs)
}

if (length(skipped_types) > 0) {
  message("\n skip:")
  message(paste(skipped_types, collapse = ", "))
} else {
  message("DEG over: ", celltype1, " vs ", celltype2)
}

if (exists("DEGs")) {
  save(DEGs, file='DEG_CANfinal_celltype_comparison.rda')
  DEGfilt <- DEGs[DEGs$p_val_adj < 0.05,]
  save(DEGfilt, file='DEGfilt_CANfinal_celltype_comparison.rda')
  write.table(DEGfilt, file='DEGfilt_CANfinal_celltype_comparison.txt', row.names=F, quote=F, sep='\t')

  DEGsta <- data.frame(
    comparison = paste0(celltype1, "_vs_", celltype2),
    number = nrow(DEGfilt)
  )
  write.table(DEGsta, file='DEGsta_CANfinal_celltype_comparison.txt', row.names=F, quote=F, sep='\t')
}
##############################################################################
setwd('')
load('DEG_CANfinal_T4T3.rda')

library(ggplot2)
library(ggrepel)
library(dplyr)
library(extrafont)
library(openxlsx)

tryCatch({
  if(! "Arial" %in% fonts()) {
    font_import(prompt = FALSE)
    loadfonts()
  }
}, error = function(e) {
  message("no Arial，use: ", e$message)
})

colors <- c("up" = "red", "down" = "blue", "ns" = "grey")

if (!dir.exists("output")) {
  dir.create("output")
}

for(cell_type in unique(DEGs$Celltype)) {

  cell_data <- DEGs %>% filter(Celltype == cell_type)

  cell_data <- cell_data %>%
    mutate(
      expression = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0.1 ~ "up",
        p_val_adj < 0.05 & avg_log2FC < -0.1 ~ "down",
        TRUE ~ "ns"
      )
    )

  up_genes <- cell_data %>% 
    filter(expression == "up") %>%
    select(gene, avg_log2FC, p_val, p_val_adj, Celltype)
  
  down_genes <- cell_data %>% 
    filter(expression == "down") %>%
    select(gene, avg_log2FC, p_val, p_val_adj, Celltype)

  wb <- createWorkbook()

  addWorksheet(wb, "up")
  writeData(wb, "up", up_genes, rowNames = FALSE)

  addWorksheet(wb, "down")
  writeData(wb, "down", down_genes, rowNames = FALSE)

  safe_cell_type <- gsub("[^[:alnum:]]", "_", cell_type)
  excel_filename <- file.path("output", paste0(safe_cell_type, "_DEGs_filtered_p0.05.xlsx"))
  saveWorkbook(wb, excel_filename, overwrite = TRUE)
  
  message("save ", cell_type, " to ", excel_filename)

  top_up <- cell_data %>%
    filter(expression == "up") %>%
    arrange(p_val_adj) %>%
    head(10)
  
  top_down <- cell_data %>%
    filter(expression == "down") %>%
    arrange(p_val_adj) %>%
    head(10)

  p <- ggplot(cell_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = expression), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", size = 16),
      axis.title = element_text(family = "Arial", face = "bold", size = 12),
      axis.text = element_text(family = "Arial", size = 10),
      legend.title = element_text(family = "Arial", face = "bold"),
      legend.text = element_text(family = "Arial")
    ) +
    labs(
      title = paste("Volcano Plot:", cell_type),
      x = "log2(Fold Change)",
      y = "-log10(Adjusted P-value)"
    ) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")

  if(nrow(top_up) > 0) {
    p <- p + geom_text_repel(
      data = top_up,
      aes(label = gene),
      size = 3,
      color = "red",
      family = "Arial",
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3
    )
  }
  
  if(nrow(top_down) > 0) {
    p <- p + geom_text_repel(
      data = top_down,
      aes(label = gene),
      size = 3,
      color = "blue",
      family = "Arial",
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3
    )
  }

  png_filename <- file.path("output", paste0("volcano_", safe_cell_type, ".png"))
  ggsave(
    filename = png_filename,
    plot = p,
    width = 10,
    height = 8,
    dpi = 300,
    limitsize = FALSE
  )

  message("type: ", cell_type, " (volcano: ", png_filename, ")")
}

all_data <- DEGs %>%
  mutate(
    expression = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.1 ~ "up",
      p_val_adj < 0.05 & avg_log2FC < -0.1 ~ "down",
      TRUE ~ "ns"
    )
  )

all_top_up <- all_data %>%
  filter(expression == "up") %>%
  arrange(p_val_adj) %>%
  head(10)

all_top_down <- all_data %>%
  filter(expression == "down") %>%
  arrange(p_val_adj) %>%
  head(10)

p_all <- ggplot(all_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = expression), alpha = 0.5, size = 1.2) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", size = 18),
    axis.title = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", size = 11),
    legend.title = element_text(family = "Arial", face = "bold"),
    legend.text = element_text(family = "Arial")
  ) +
  labs(
    title = "Volcano Plot: All Cell Types",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-value)"
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

if(nrow(all_top_up) > 0) {
  p_all <- p_all + geom_text_repel(
    data = all_top_up,
    aes(label = gene),
    size = 3.5,
    color = "red",
    family = "Arial",
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3
  )
}

if(nrow(all_top_down) > 0) {
  p_all <- p_all + geom_text_repel(
    data = all_top_down,
    aes(label = gene),
    size = 3.5,
    color = "blue",
    family = "Arial",
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3
  )
}
png_filename_all <- file.path("output", "volcano_all_cell_types.png")
ggsave(
  filename = png_filename_all,
  plot = p_all,
  width = 12,
  height = 9,
  dpi = 300,
  limitsize = FALSE
)


##############################################################################
library(ggplot2)
library(ggrepel)
library(dplyr)
library(extrafont)
library(openxlsx)

tryCatch({
  if(! "Arial" %in% fonts()) {
    font_import(prompt = FALSE)
    loadfonts()
  }
}, error = function(e) {
  message("no Arial，use: ", e$message)
})

colors <- c("up" = "red", "down" = "blue", "ns" = "grey")

if (!dir.exists("output")) {
  dir.create("output")
}

if (!"comparison" %in% colnames(DEGs)) {
  stop("no comparison")
}

for(comp in unique(DEGs$comparison)) {

  comp_data <- DEGs %>% filter(comparison == comp)

  comp_data <- comp_data %>%
    mutate(
      expression = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0.1 ~ "up",
        p_val_adj < 0.05 & avg_log2FC < -0.1 ~ "down",
        TRUE ~ "ns"
      )
    )

  up_genes <- comp_data %>% 
    filter(expression == "up") %>%
    select(gene, avg_log2FC, p_val, p_val_adj, comparison)
  
  down_genes <- comp_data %>% 
    filter(expression == "down") %>%
    select(gene, avg_log2FC, p_val, p_val_adj, comparison)

  wb <- createWorkbook()

  if(nrow(up_genes) > 0) {
    addWorksheet(wb, "up")
    writeData(wb, "up", up_genes, rowNames = FALSE)
  }

  if(nrow(down_genes) > 0) {
    addWorksheet(wb, "down")
    writeData(wb, "down", down_genes, rowNames = FALSE)
  }

  safe_comp <- gsub("[^[:alnum:]]", "_", comp)
  excel_filename <- file.path("output", paste0(safe_comp, "_DEGs_filtered_p0.05.xlsx"))
  saveWorkbook(wb, excel_filename, overwrite = TRUE)
  
  message("saved ", comp, " updown gene ", excel_filename)

  top_up <- comp_data %>%
    filter(expression == "up") %>%
    arrange(p_val_adj) %>%
    head(10)
  
  top_down <- comp_data %>%
    filter(expression == "down") %>%
    arrange(p_val_adj) %>%
    head(10)

  p <- ggplot(comp_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = expression), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", family = "Arial", size = 16),
      axis.title = element_text(family = "Arial", face = "bold", size = 12),
      axis.text = element_text(family = "Arial", size = 10),
      legend.title = element_text(family = "Arial", face = "bold"),
      legend.text = element_text(family = "Arial")
    ) +
    labs(
      title = paste("Volcano Plot:", comp),
      x = "log2(Fold Change)",
      y = "-log10(Adjusted P-value)"
    ) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")

  if(nrow(top_up) > 0) {
    p <- p + geom_text_repel(
      data = top_up,
      aes(label = gene),
      size = 3,
      color = "red",
      family = "Arial",
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3
    )
  }
  
  if(nrow(top_down) > 0) {
    p <- p + geom_text_repel(
      data = top_down,
      aes(label = gene),
      size = 3,
      color = "blue",
      family = "Arial",
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3
    )
  }

  png_filename <- file.path("output", paste0("volcano_", safe_comp, ".png"))
  ggsave(
    filename = png_filename,
    plot = p,
    width = 10,
    height = 8,
    dpi = 300,
    limitsize = FALSE
  )

  message("complete: ", comp, " (volcano: ", png_filename, ")")
}
