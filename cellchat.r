library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("simplify", "igraph")
conflict_prefer("as_data_frame", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("crossing", "tidyr")
conflict_prefer("%--%", "lubridate")

library(CellChat)
library(tidyverse)
library(Seurat)
library(future)
library(NMF)
library(ggplot2)
library(reshape2)
library(patchwork)
library(ggalluvial)
library(ggplot2)
library(ggalluvial)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 100*1024^3)
setwd('')
load('CANfinal_UpdateType.rda')

color_mapping <- c(
  "Fib1" = "#F5D2A8",
  "Fib2" = "#FCED82",
  "T/NK cell" = "#EE934E",
  "B cell" = "#3C77AF",
  "Epithelial" = "#D1352B",
  "Plasma cell" = "#F5CFE4",
  "Macrophage" = "#9B5B33",
  "cDCs" = "#B383B9",
  "pDCs" = "#8FA4AE",
  "EC" = "#AECDE1",
  "LEC" = "#D2EBC8",
  "Adipocyte" = "#7F7F7F"
)

analysis_groups <- c("T2", "T3", "T4", "all")
analysis_group <- "T2"
for(analysis_group in analysis_groups) {
  
  cat("Processing sample:", analysis_group, "\n")

  sample_dir <- paste0('', analysis_group)
  if(!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE)
  }
  setwd(sample_dir)

  Idents(CANfinal) <- "Type"

  if(analysis_group == "all") {
    CANfinal_subset <- subset(CANfinal, idents = c("T2", "T3", "T4"))
    cat("Created 'all' group by merging T2, T3, T4 samples\n")
  } else {
    CANfinal_subset <- subset(CANfinal, idents = analysis_group)
  }

  cellchat <- createCellChat(object = CANfinal_subset,
                           meta = CANfinal_subset@meta.data,
                           group.by = "Celltype")
  
  cat("CellChat object created for", analysis_group, "\n")
  cat("Number of cells:", ncol(CANfinal_subset), "\n")

  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) 
  plan(multisession)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, paste0("net_lr_", analysis_group, ".csv"))
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, paste0("net_pathway_", analysis_group, ".csv"))
  all_pathways <- cellchat@netP$pathways
  cat("Total pathways detected for", analysis_group, ":", length(all_pathways), "\n")
  if(length(all_pathways) > 0) {
    if(length(all_pathways) >= 20) {
      pathways.show <- all_pathways[1:20]
    } else {
      pathways.show <- all_pathways
    }
    cat("Selected pathways for", analysis_group, ":", paste(pathways.show, collapse = ", "), "\n")
    write.table(pathways.show, paste0("selected_pathways_", analysis_group, ".txt"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  } else {
    cat("No pathways detected for", analysis_group, ". Skipping pathway visualizations.\n")
    pathways.show <- NULL
  }

  saveRDS(cellchat, file = paste0("cellchat_processed_", analysis_group, ".rds"))

  groupSize <- as.numeric(table(cellchat@idents))

  celltypes <- levels(cellchat@idents)
  colors_for_plot <- color_mapping[celltypes]
  cellchat <- aggregateNet(cellchat)
  pdf(paste0("Interaction_Count_and_Strength_", analysis_group, ".pdf"), width = 10, height = 5)
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, title.name = "Number of interactions", color.use = colors_for_plot)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, title.name = "Interaction weights/strength", color.use = colors_for_plot)
  dev.off()

  pdf(paste0("Cell_Specific_Interaction_Count_", analysis_group, ".pdf"), width = 12, height = 12)
  mat <- cellchat@net$count
  par(mfrow = c(3,3), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                     arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i],
                     color.use = colors_for_plot)
  }
  dev.off()
  pdf(paste0("Cell_Specific_Interaction_Strength_", analysis_group, ".pdf"), width = 12, height = 12)
  mat <- cellchat@net$weight
  par(mfrow = c(3,3), xpd=T)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                     arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i],
                     color.use = colors_for_plot)
  }
  dev.off()
  
  cat("Available cell types for", analysis_group, ":", paste(celltypes, collapse = ", "), "\n")
  if(!is.null(pathways.show) && length(pathways.show) > 0) {

    pdf(paste0("circle_", analysis_group, ".pdf"), width = 12, height = 12)
    par(mfrow=c(1,1))
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = colors_for_plot)
    dev.off()

    pdf(paste0("chord_", analysis_group, ".pdf"), width = 12, height = 12)
    par(mfrow=c(1,1))
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", color.use = colors_for_plot)
    dev.off()
    pdf(paste0("circle_plots_", analysis_group, ".pdf"), width = 14, height = 10)
    par(mfrow = c(2, 3))
    
    for (pathway in pathways.show) {
      netVisual_aggregate(
        cellchat, 
        signaling = pathway, 
        layout = "circle",
        vertex.label.cex = 1.2,
        color.use = colors_for_plot
      )
      
      title(main = pathway, cex.main = 1.5)
    }
    dev.off()
    pdf(paste0("chord_plots_", analysis_group, ".pdf"), width = 12, height = 12)
    
    for (pathway in pathways.show) {
      plot.new()
      title(main = pathway, cex.main = 2)
      
      netVisual_aggregate(
        cellchat, 
        signaling = pathway, 
        layout = "chord",
        big.gap = 20,
        small.gap = 1,
        color.use = colors_for_plot
      )
    }
    dev.off()
    pdf(paste0("ALL_LR_contribution_", analysis_group, ".pdf"), width = 12, height = 12)
    netAnalysis_contribution(cellchat, signaling = pathways.show)
    dev.off()

    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    pdf(paste0("Signaling_Role_Network_", analysis_group, ".pdf"), width = 15, height = 6)
    netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, font.size = 10, color.use = colors_for_plot)
    dev.off()
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 8, color.heatmap = "Reds")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 8, color.heatmap = "Blues")
    pdf(paste0("Signaling_Pattern_Heatmaps_", analysis_group, ".pdf"), width = 14, height = 10)
    ht1 + ht2
    dev.off()
  }

  outgoing_dir <- paste0(sample_dir, "/outgoing")
  if(!dir.exists(outgoing_dir)) {
    dir.create(outgoing_dir)
  }
  setwd(outgoing_dir)
  
  tryCatch({
    pdf("selectK.pdf", width = 14, height = 10)
    selectK(cellchat, pattern = "outgoing")
    dev.off()
    
    pdf("CommunicationPatterns.pdf", width = 14, height = 14)
    nPatterns = 4
    cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, height = 12)
    dev.off()
    library(ggplot2)
    library(ggalluvial)
    
    pdf("Communicationriver.pdf", width = 14, height = 14)
    netAnalysis_river(cellchat, pattern = "outgoing")
    dev.off()
    pdf("Communicationdot.pdf", width = 14, height = 14)
    netAnalysis_dot(cellchat, pattern = "outgoing")
    dev.off()
    
    out_data_matrix <- as.data.frame(cellchat@netP[["pattern"]][["outgoing"]][["data"]])
    write.csv(out_data_matrix, paste0("NMF_out_data_matrix_", analysis_group, ".csv"))
  }, error = function(e) {
    cat("Error in outgoing analysis:", e$message, "\n")
  })
  incoming_dir <- paste0(sample_dir, "/incoming")
  if(!dir.exists(incoming_dir)) {
    dir.create(incoming_dir)
  }
  setwd(incoming_dir)
  
  tryCatch({
    pdf("selectK.pdf", width = 14, height = 14)
    selectK(cellchat, pattern = "incoming")
    dev.off()
    
    pdf("CommunicationPatterns.pdf", width = 14, height = 14)
    nPatterns = 6
    cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, height = 12)
    dev.off()

    library(ggplot2)
    library(ggalluvial)
    
    pdf("Communicationriver.pdf", width = 14, height = 14)
    netAnalysis_river(cellchat, pattern = "incoming")
    dev.off()
    pdf("Communicationdot.pdf", width = 14, height = 14)
    netAnalysis_dot(cellchat, pattern = "incoming")
    dev.off()
    
    in_data_matrix <- as.data.frame(cellchat@netP[["pattern"]][["incoming"]][["data"]])
    write.csv(in_data_matrix, paste0("NMF_in_data_matrix_", analysis_group, ".csv"))
  }, error = function(e) {
    cat("Error in incoming analysis:", e$message, "\n")
  })
  
  cat("Analysis completed for", analysis_group, "\n")

  rm(cellchat, CANfinal_subset, df.net, df.netp)
  if(exists("out_data_matrix")) rm(out_data_matrix)
  if(exists("in_data_matrix")) rm(in_data_matrix)
  gc()
}