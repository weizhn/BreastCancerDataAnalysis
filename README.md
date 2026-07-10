# BreastCancerDataAnalysis

scRNA.r – Integrates, clusters, and annotates scRNA‑seq data from multiple samples, performs differential expression, and generates marker‑based visualizations.

GSVA_HALLMARKER.r – Computes HALLMARK pathway scores via GSVA on single‑cell expression data, compares pathway activity across groups, and produces differential pathway analyses.

infercnv.r – Infers copy number variations (CNVs) from scRNA‑seq counts using InferCNV, with median filtering and reference‑based denoising.

RCTD.r – Deconvolves spatial transcriptomics spots using RCTD, assigning cell types based on a single‑cell reference and saving normalized weights.

stereoseq_mixed_definition.py – To mitigate the influence of spatial transcriptomic mixing effects caused by limited resolution and data variation, we classified cells lacking sufficient local neighborhood consensus as ‘mixed’. These cells were excluded from downstream cell-type-specific analyses to ensure the robustness of our spatial annotations.

cellchat.r – Constructs and analyses cell‑cell communication networks via CellChat across different sample groups, identifying significant ligand‑receptor pathways and interaction patterns.

phenoptr.py – The nearest neighbor distance from spatial transcriptome spots to epithelial cells was calculated using phenoptr package, and the distance statistics were summarized according to cell type. The specificity of fib2 cells was compared, and the visualization results such as distance distribution, histogram and scatter trend graph were generated.
