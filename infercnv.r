library(Seurat)
library(dplyr)
library(infercnv)
library(rtracklayer)
setwd('')
load("CANfinal_UpdateType.rda")
analysis_group <- "all"
if (analysis_group != "all") {
  Idents(CANfinal) <- "Type"
  CANfinal_subset <- subset(CANfinal, idents = analysis_group)
} else {
  CANfinal_subset <- CANfinal
}

set.seed(42)
all_cells <- colnames(CANfinal_subset)
sample_size <- round(length(all_cells) * 0.2)
sampled_cells <- sample(all_cells, size = sample_size)
CANfinal_sampled <- subset(CANfinal_subset, cells = sampled_cells)

exprMatrix <- as.matrix(GetAssayData(CANfinal_sampled, slot = 'counts'))
cellAnnota <- subset(CANfinal_sampled@meta.data, select = 'Celltype')
cellAnnota$Celltype <- as.factor(cellAnnota$Celltype)

options(scipen = 100)

output_dir <- paste0('PLASMA_Tdenoise_', analysis_group, '/')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = exprMatrix,
  annotations_file = cellAnnota,
  delim = "\t",
  gene_order_file = "gencode_homo_gene_pos.txt",
  ref_group_names = c("Plasma cell")
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = output_dir,
  cluster_by_groups = FALSE,
  hclust_method = "ward.D2",
  plot_steps = FALSE,
  HMM = TRUE,
  num_threads = 1,
  denoise = TRUE,
  write_expr_matrix = TRUE
)
infercnv_obj_medianfiltered <- infercnv::apply_median_filtering(infercnv_obj)
plot_dir <- paste0(output_dir, 'plot_cnv/')
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
assignInNamespace(
  x = "get_group_color_palette",
  value = function(group_names) {
    cols <- color_mapping[group_names]
    missing <- is.na(cols)
    if (any(missing)) {
      fallback <- grDevices::rainbow(sum(missing))
      cols[missing] <- fallback
    }
    names(cols) <- group_names
    return(cols)
  },
  ns = "infercnv"
)

infercnv::plot_cnv(
  infercnv_obj_medianfiltered,
  out_dir = plot_dir,
  output_filename = paste0('infercnv.median_filtered.', analysis_group),
  x.range = "auto",
  x.center = 1,
  title = paste0("infercnv - ", analysis_group),
  color_safe_pal = FALSE
)