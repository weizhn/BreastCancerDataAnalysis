library(Seurat)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggrepel)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(ggpubr)

setwd('')
load("CANfinal_UpdateType.rda")
expr_matrix <- GetAssayData(CANfinal, assay = "RNA", slot = "data")
expr_matrix <- as.matrix(expr_matrix)

genes_expressed <- rowSums(expr_matrix > 0)
expr_matrix <- expr_matrix[genes_expressed > ncol(expr_matrix) * 0.01, ]

read_list_file <- function(list_file) {
  list_data <- readLines(list_file)
  gene_sets_list <- list()
  total_lines <- length(list_data)

  for (line in list_data) {
    line_clean <- gsub('"', '', line)
    fields <- strsplit(line_clean, ",")[[1]]
    if (length(fields) >= 3) {
      pathway_name <- trimws(fields[2])
      gene_symbol <- trimws(fields[3])

      if (!pathway_name %in% names(gene_sets_list)) {
        gene_sets_list[[pathway_name]] <- character(0)
      }
      gene_sets_list[[pathway_name]] <- c(gene_sets_list[[pathway_name]], gene_symbol)
    }
  }

  gene_sets_list <- lapply(gene_sets_list, unique)
  gene_sets_list <- gene_sets_list[sapply(gene_sets_list, length) >= 5]
  return(gene_sets_list)
}

list_file_path <- "h.all.v2023.1.Hs.symbols.list"
gene_sets <- read_list_file(list_file_path)

gs_sizes <- sapply(gene_sets, length)
for(i in 1:min(5, length(gene_sets))) {
  pathway_name <- names(gene_sets)[i]
  print(paste(pathway_name, ":", length(gene_sets[[pathway_name]]), "个基因"))
}


start_time <- Sys.time()
gsva_gene_sets <- gene_sets

all_genes_in_sets <- unique(unlist(gsva_gene_sets))
genes_in_expr <- rownames(expr_matrix)
matched_genes <- intersect(all_genes_in_sets, genes_in_expr)

filtered_gene_sets <- list()
for(pathway_name in names(gsva_gene_sets)) {
  pathway_genes <- gsva_gene_sets[[pathway_name]]
  filtered_genes <- intersect(pathway_genes, genes_in_expr)
  if(length(filtered_genes) >= 1) {
    filtered_gene_sets[[pathway_name]] <- filtered_genes
  }
}

gsva_scores <- gsva(
  expr = expr_matrix,
  gset.idx.list = filtered_gene_sets,
  method = "gsva",
  kcdf = "Gaussian",
  min.sz = 5,
  max.sz = 500,
  parallel.sz = 2,
  verbose = TRUE
)
end_time <- Sys.time()

write.csv(gsva_scores, file = "GSVA_scores_offline_all_cells_HALLMARK.csv", row.names = TRUE)

group_info <- CANfinal$Type
unique_groups <- unique(group_info)

gsva_mean_by_group <- matrix(NA, nrow = nrow(gsva_scores), ncol = length(unique_groups))
colnames(gsva_mean_by_group) <- unique_groups
rownames(gsva_mean_by_group) <- rownames(gsva_scores)

for (i in 1:length(unique_groups)) {
  group_cells <- which(group_info == unique_groups[i])
  if (length(group_cells) > 0) {
    gsva_mean_by_group[, i] <- rowMeans(gsva_scores[, group_cells, drop = FALSE])
  }
}

write.csv(gsva_mean_by_group, file = "GSVA_scores_offline_mean_by_group_HALLMARK.csv", row.names = TRUE)


design <- model.matrix(~0 + factor(group_info))
colnames(design) <- levels(factor(group_info))
colnames(gsva_scores) <- colnames(expr_matrix)
rownames(design) <- colnames(gsva_scores)

fit <- lmFit(gsva_scores, design)

group_names <- colnames(design)
if (length(group_names) >= 2) {
  contrast_pairs <- combn(group_names, 2)
  contrast_formulas <- apply(contrast_pairs, 2, function(pair) {
    paste(pair[1], "-", pair[2])
  })
  
  contrast_matrix <- makeContrasts(contrasts = contrast_formulas, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  for (i in 1:length(contrast_formulas)) {
    diff_results <- topTable(fit2, coef = i, number = Inf, adjust.method = "BH")
    diff_results <- cbind(Pathway = rownames(diff_results), diff_results)

    contrast_name <- gsub(" - ", "_vs_", contrast_formulas[i])
    write.csv(diff_results, file = paste0("GSVA_differential_pathways_HALLMARK_", contrast_name, ".csv"), row.names = FALSE)
    
    print(paste("Contrast", contrast_formulas[i], "yielded", 
                sum(diff_results$adj.P.Val < 0.05), "significantly differential pathways"))
  }
} else {
  print("Insufficient number of groups")
}

dir.create("GSVA_offline_HALLMARK_plots", showWarnings = FALSE)

pdf_file <- "GSVA_offline_HALLMARK_visualization_results.pdf"
pdf(pdf_file, width = 11, height = 8.5)
top_n_pathways <- 50
if (nrow(gsva_scores) > top_n_pathways) {
  pathway_var <- apply(gsva_scores, 1, var)
  plot_pathways <- names(sort(pathway_var, decreasing = TRUE))[1:min(top_n_pathways, nrow(gsva_scores))]
} else {
  plot_pathways <- rownames(gsva_scores)
}
max_cells_heatmap <- 1000
if (ncol(gsva_scores) > max_cells_heatmap) {
  set.seed(123)
  sampled_cells <- sample(colnames(gsva_scores), max_cells_heatmap)
  heatmap_data <- gsva_scores[plot_pathways, sampled_cells]
  cell_annotations <- data.frame(Group = group_info[sampled_cells])
} else {
  heatmap_data <- gsva_scores[plot_pathways, ]
  cell_annotations <- data.frame(Group = group_info)
}

rownames(cell_annotations) <- colnames(heatmap_data)

group_colors <- rainbow(length(unique(group_info)))
names(group_colors) <- unique(group_info)
annotation_colors <- list(Group = group_colors)

pheatmap(
  heatmap_data,
  main = "HALLMARK Pathway Scores Heatmap (Offline Analysis)",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  annotation_col = cell_annotations,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "row"
)

if (length(plot_pathways) > 0) {
  display_pathways <- head(plot_pathways, min(12, length(plot_pathways)))
  plot_data <- data.frame()
  for (pathway in display_pathways) {
    if (pathway %in% rownames(gsva_scores)) {
      pathway_scores <- gsva_scores[pathway, ]
      plot_data <- rbind(plot_data, 
                         data.frame(
                           Pathway = pathway,
                           Score = pathway_scores,
                           Group = group_info,
                           Cell = names(pathway_scores)
                         ))
    }
  }
  plot_data$Pathway_short <- sapply(as.character(plot_data$Pathway), function(x) {
    if (nchar(x) > 60) {
      paste0(substr(x, 1, 57), "...")
    } else {
      x
    }
  })
  p_box <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, alpha = 0.7) +
    facet_wrap(~ Pathway_short, scales = "free_y", ncol = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      strip.text = element_text(size = 8)
    ) +
    labs(title = "HALLMARK Pathway Score Distribution by Group",
         subtitle = "Top pathways with highest variance",
         x = "Group", y = "GSVA Score") +
    scale_fill_manual(values = group_colors)
  
  print(p_box)
}

if (length(plot_pathways) >= 10) {
  top_pathways_corr <- cor(t(gsva_scores[head(plot_pathways, 20), ]))
  pheatmap(
    top_pathways_corr,
    main = "HALLMARK Pathway-Pathway Correlation (Top 20 Pathways)",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    clustering_method = "average",
    fontsize_row = 8,
    fontsize_col = 8,
    display_numbers = FALSE
  )
}

if (length(unique_groups) >= 3 && nrow(gsva_mean_by_group) >= 6) {
  library(fmsb)
  top_6_pathways <- rownames(gsva_mean_by_group)[1:min(6, nrow(gsva_mean_by_group))]
  radar_data <- as.data.frame(t(gsva_mean_by_group[top_6_pathways, ]))

  max_values <- apply(radar_data, 2, max)
  min_values <- apply(radar_data, 2, min)
  radar_data <- rbind(rep(max(max_values), ncol(radar_data)), 
                      rep(0, ncol(radar_data)), 
                      radar_data)

  op <- par(mar = c(1, 1, 2, 1))
  radarchart(radar_data,
             axistype = 1,
             pcol = group_colors[rownames(radar_data)[3:nrow(radar_data)]],
             pfcol = alpha(group_colors[rownames(radar_data)[3:nrow(radar_data)]], 0.3),
             plwd = 2,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "grey",
             caxislabels = seq(0, max(max_values), length.out = 5),
             cglwd = 0.8,
             vlcex = 0.8,
             title = "HALLMARK Pathway Activity Radar Chart (Top 6 Pathways)")
  legend("topright", 
         legend = rownames(radar_data)[3:nrow(radar_data)],
         bty = "n", pch = 20, col = group_colors[rownames(radar_data)[3:nrow(radar_data)]],
         text.col = "black", cex = 0.8, pt.cex = 1.5)
  par(op)
}

dev.off()