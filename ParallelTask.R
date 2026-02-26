# ============================================================
# PARALLEL PIPELINES FOR 10x scRNA-seq (2 samples)
#
# Pipeline A: Filtered -> Seurat -> scDblFinder -> Recluster -> Plots
# Pipeline B: Raw+Filtered -> SoupX -> Seurat -> scDblFinder -> Recluster -> Plots
#
# Outputs:
# - Side-by-side UMAPs (A vs B)
# - Side-by-side FeaturePlots for key markers (A vs B)
# - Side-by-side DotPlots (A vs B)
# - Doublet diagnostics for both
#
# Teaching notes included throughout.
# ============================================================

# ---------------------------
# 0) User paths + extraction
# ---------------------------
zip_file <- "C:/Users/pednet/Desktop/Leipzig Task/Revised/Data_for_task.zip"
out_dir  <- "C:/Users/pednet/Desktop/Leipzig Task/Revised/Data_for_task_unzipped"
res_dir  <- "C:/Users/pednet/Desktop/Leipzig Task/Revised/results_parallel"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

unzip(zip_file, exdir = out_dir)

# Teaching note:
# Always inspect the directory after unzip to confirm where the dataset root is.
# Some zips contain a top folder, others extract files directly.
list.files(out_dir, recursive = FALSE)

# Detect base folder robustly
base1 <- file.path(out_dir, "Data_for_task")
base  <- if (dir.exists(base1)) base1 else out_dir

# ---------------------------
# 1) Libraries
# ---------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# SoupX (CRAN)
if (!requireNamespace("SoupX", quietly = TRUE)) {
  stop("SoupX not installed. Install with: install.packages('SoupX')")
}

# scDblFinder (Bioconductor)
if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
    !requireNamespace("scDblFinder", quietly = TRUE)) {
  stop("Bioconductor packages missing. Run:\n  if(!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n  BiocManager::install(c('SingleCellExperiment','scDblFinder'))")
}
library(SingleCellExperiment)
library(scDblFinder)

# ---------------------------
# 2) Common settings (keep identical between pipelines)
# ---------------------------
set.seed(123)              # controls randomness for clustering/UMAP where used
umap_seed <- 123
dims_use  <- 1:30
npcs_use  <- 50
res_use   <- 0.4

# Marker panels for plotting (target + context)
marker_adipo <- c("ADIPOQ","PLIN1","FABP4","LPL")
marker_aspc  <- c("PDGFRA","DCN","COL1A1","LUM","PI16","CD55")
marker_immune <- c("LST1","C1QC","MRC1","CD163")
marker_endo   <- c("VWF","PECAM1","EMCN","PLVAP")
marker_mural  <- c("RGS5","ACTA2","TAGLN","MYH11")

marker_panel <- c(marker_adipo, marker_aspc, marker_immune, marker_endo, marker_mural)

# ---------------------------
# 3) Paths + Read10X (raw + filtered)
# ---------------------------
s1_filt <- file.path(base, "21L005039", "filtered_feature_bc_matrix")
s1_raw  <- file.path(base, "21L005039", "raw_feature_bc_matrix")

s2_filt <- file.path(base, "21L005043", "filtered_feature_bc_matrix")
s2_raw  <- file.path(base, "21L005043", "raw_feature_bc_matrix")

m1_filt <- Read10X(data.dir = s1_filt)
m2_filt <- Read10X(data.dir = s2_filt)

m1_raw  <- Read10X(data.dir = s1_raw)
m2_raw  <- Read10X(data.dir = s2_raw)

cat("Dims (genes x barcodes)\n")
cat("S1 filtered:", dim(m1_filt), " | raw:", dim(m1_raw), "\n")
cat("S2 filtered:", dim(m2_filt), " | raw:", dim(m2_raw), "\n")

# ---------------------------
# 4) Helper functions
# ---------------------------

# ---- 4A) Add QC metrics
add_qc_metrics <- function(obj) {
  obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
  obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGQZ]")
  return(obj)
}

# ---- 4B) Initial QC filter (broad first pass)
qc_filter <- function(obj) {
  # Teaching note:
  # These are broad thresholds designed to remove obvious low-quality cells and extreme outliers.
  # In adipose/snRNA data, percent.mt is often low, so we keep it permissive.
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
}

# ---- 4C) Minimal Seurat workflow used BEFORE doublet detection
run_seurat_basic <- function(obj, label) {
  # Teaching note:
  # We need a reasonable embedding + clusters so that doublet detection can leverage structure.
  # This is not the final clustering yet.
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  obj <- RunPCA(obj, npcs = npcs_use, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res_use, random.seed = 123, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, seed.use = umap_seed, verbose = FALSE)
  
  return(obj)
}

# ---- 4D) Doublet detection (scDblFinder) + return singlets only
run_doublets_scDblFinder <- function(obj, label) {
  # Teaching note:
  # scDblFinder simulates doublets and classifies cells as singlet/doublet.
  # We run it after initial clustering so it can use neighborhood structure.
  
  sce <- as.SingleCellExperiment(obj)
  
  set.seed(123)
  sce <- scDblFinder(
    sce,
    samples  = sce$sample_id,
    clusters = as.factor(obj$seurat_clusters)
  )
  
  obj$scDblFinder.score <- colData(sce)$scDblFinder.score
  obj$scDblFinder.class <- colData(sce)$scDblFinder.class
  
  # Save summary
  tab <- table(obj$scDblFinder.class, obj$sample_id)
  write.csv(as.data.frame(tab),
            file = file.path(res_dir, paste0(label, "_doublet_table.csv")),
            row.names = FALSE)
  
  # Diagnostics plots
  p_class <- DimPlot(obj, reduction = "umap", group.by = "scDblFinder.class") +
    ggtitle(paste0(label, ": scDblFinder class"))
  p_score <- FeaturePlot(obj, features = "scDblFinder.score") +
    ggtitle(paste0(label, ": scDblFinder score"))
  
  ggsave(file.path(res_dir, paste0(label, "_doublet_umap_class.png")),
         p_class, width = 7, height = 6, dpi = 300)
  ggsave(file.path(res_dir, paste0(label, "_doublet_umap_score.png")),
         p_score, width = 7, height = 6, dpi = 300)
  
  # QC violin by class
  p_vln <- VlnPlot(obj,
                   features = c("nFeature_RNA","nCount_RNA","percent.mt","scDblFinder.score"),
                   group.by = "scDblFinder.class",
                   ncol = 4, pt.size = 0.05) +
    plot_annotation(title = paste0(label, ": QC by scDblFinder class"))
  ggsave(file.path(res_dir, paste0(label, "_doublet_qc_violin.png")),
         p_vln, width = 14, height = 4, dpi = 300)
  
  # Remove doublets
  obj_sing <- subset(obj, subset = scDblFinder.class == "singlet")
  return(obj_sing)
}

# ---- 4E) Final re-run after doublet removal (clean clustering)
run_seurat_final <- function(obj, label) {
  # Teaching note:
  # After removing doublets, we re-run the full workflow to obtain the final clean clusters.
  # Otherwise, clusters/UMAP still reflect cells that are no longer present.
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  obj <- RunPCA(obj, npcs = npcs_use, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res_use, random.seed = 123, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, seed.use = umap_seed, verbose = FALSE)
  
  # Save core UMAPs
  p_samp <- DimPlot(obj, reduction = "umap", group.by = "sample_id") +
    ggtitle(paste0(label, ": UMAP by sample"))
  p_clu  <- DimPlot(obj, reduction = "umap", label = TRUE) +
    ggtitle(paste0(label, ": UMAP by cluster"))
  
  ggsave(file.path(res_dir, paste0(label, "_final_umap_by_sample.png")),
         p_samp, width = 7, height = 6, dpi = 300)
  ggsave(file.path(res_dir, paste0(label, "_final_umap_by_cluster.png")),
         p_clu, width = 7.5, height = 6, dpi = 300)
  
  saveRDS(obj, file = file.path(res_dir, paste0(label, "_final_obj.rds")))
  return(obj)
}

# ---- 4F) SoupX correction per sample (raw + filtered -> corrected filtered counts)
run_soupx_correct <- function(raw_counts, filt_counts, sample_name) {
  # Teaching note:
  # SoupX uses empty droplets (raw) to model the "soup" (ambient RNA),
  # then subtracts estimated contamination from cell-associated droplets (filtered).
  #
  # To estimate contamination well, SoupX benefits from cluster labels of filtered cells.
  # We generate light clustering for this purpose only.
  
  # Temporary Seurat clustering for SoupX guidance
  tmp <- CreateSeuratObject(counts = filt_counts, project = sample_name, min.cells = 3, min.features = 100)
  tmp <- add_qc_metrics(tmp)
  
  tmp <- NormalizeData(tmp, verbose = FALSE)
  tmp <- FindVariableFeatures(tmp, nfeatures = 2000, verbose = FALSE)
  tmp <- ScaleData(tmp, vars.to.regress = "percent.mt", verbose = FALSE)
  tmp <- RunPCA(tmp, npcs = 30, verbose = FALSE)
  tmp <- FindNeighbors(tmp, dims = 1:20, verbose = FALSE)
  tmp <- FindClusters(tmp, resolution = 0.3, random.seed = 123, verbose = FALSE)
  tmp <- RunUMAP(tmp, dims = 1:20, seed.use = umap_seed, verbose = FALSE)
  
  # Build SoupX channel
  sc <- SoupX::SoupChannel(tod = raw_counts, toc = filt_counts, calcSoupProfile = TRUE)
  
  sc <- SoupX::setClusters(sc, setNames(as.character(Idents(tmp)), colnames(tmp)))
  sc <- SoupX::setDR(sc, Embeddings(tmp, "umap"))
  
  # Estimate contamination automatically (usually OK)
  sc <- SoupX::autoEstCont(sc, doPlot = FALSE)
  
  # Optional: cap rho if you ever see over-correction (not necessary with your rho values)
  # rho <- sc$metaData$rho
  # sc  <- SoupX::setContaminationFraction(sc, min(rho, 0.10))
  
  corrected <- SoupX::adjustCounts(sc, roundToInt = TRUE)
  
  # Save rho for reporting
  rho <- sc$metaData$rho
  writeLines(paste0(sample_name, " SoupX rho (estimated contamination): ", unique(rho)),
             con = file.path(res_dir, paste0(sample_name, "_SoupX_rho.txt")))
  
  return(corrected)
}

# ---- 4G) Plotting helpers: FeaturePlot grid and DotPlot
plot_markers_and_save <- function(obj, label) {
  # Teaching note:
  # FeaturePlots show where marker genes are expressed on the UMAP.
  # DotPlots summarize marker expression across clusters (great for PPT).
  
  p_adipo <- FeaturePlot(obj, features = marker_adipo, ncol = 2) +
    plot_annotation(title = paste0(label, ": Adipocyte markers"))
  p_aspc  <- FeaturePlot(obj, features = marker_aspc, ncol = 3) +
    plot_annotation(title = paste0(label, ": ASPC markers"))
  p_imm   <- FeaturePlot(obj, features = marker_immune, ncol = 2) +
    plot_annotation(title = paste0(label, ": Macrophage markers"))
  
  ggsave(file.path(res_dir, paste0(label, "_feature_adipocyte.png")), p_adipo, width = 9, height = 7, dpi = 300)
  ggsave(file.path(res_dir, paste0(label, "_feature_aspc.png")),      p_aspc,  width = 12, height = 7, dpi = 300)
  ggsave(file.path(res_dir, paste0(label, "_feature_macrophage.png")),p_imm,   width = 9, height = 7, dpi = 300)
  
  dp <- DotPlot(obj, features = marker_panel) + RotatedAxis() +
    ggtitle(paste0(label, ": Marker DotPlot (clusters)"))
  ggsave(file.path(res_dir, paste0(label, "_dotplot_markers.png")), dp, width = 14, height = 5, dpi = 300)
  
  invisible(TRUE)
}

# ---------------------------
# 5) PIPELINE A: Filtered + doublets
# ---------------------------
cat("\n=== Running Pipeline A: FILTERED + scDblFinder ===\n")

A_obj1 <- CreateSeuratObject(counts = m1_filt, project = "21L005039", min.cells = 3, min.features = 100)
A_obj2 <- CreateSeuratObject(counts = m2_filt, project = "21L005043", min.cells = 3, min.features = 100)

A_obj1$sample_id <- "21L005039"; A_obj2$sample_id <- "21L005043"
A_obj1 <- add_qc_metrics(A_obj1)
A_obj2 <- add_qc_metrics(A_obj2)

# Pre-filter QC plots (optional quick export)
ggsave(file.path(res_dir, "A_prefilter_violin_S1.png"),
       VlnPlot(A_obj1, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb"), ncol = 5, pt.size = 0.05),
       width = 14, height = 4, dpi = 300)
ggsave(file.path(res_dir, "A_prefilter_violin_S2.png"),
       VlnPlot(A_obj2, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb"), ncol = 5, pt.size = 0.05),
       width = 14, height = 4, dpi = 300)

A_obj1 <- qc_filter(A_obj1)
A_obj2 <- qc_filter(A_obj2)

A <- merge(A_obj1, y = A_obj2, add.cell.ids = c("S1","S2"), project = "Leipzig_A_filtered")

DefaultAssay(A) <- "RNA"
A <- JoinLayers(A) 
A <- run_seurat_basic(A, label = "A_basic")
A_sing <- run_doublets_scDblFinder(A, label = "A")
A_final <- run_seurat_final(A_sing, label = "A_filtered_dblclean")

plot_markers_and_save(A_final, label = "A_filtered_dblclean")


# ---------------------------
# 6) PIPELINE B: SoupX + doublets
# ---------------------------
cat("\n=== Running Pipeline B: RAW+FILTERED -> SoupX + scDblFinder ===\n")

# SoupX correction per sample (returns corrected filtered counts)
m1_sx <- run_soupx_correct(m1_raw, m1_filt, "21L005039")
m2_sx <- run_soupx_correct(m2_raw, m2_filt, "21L005043")

# Build Seurat objects from corrected counts
B_obj1 <- CreateSeuratObject(counts = m1_sx, project = "21L005039", min.cells = 3, min.features = 100)
B_obj2 <- CreateSeuratObject(counts = m2_sx, project = "21L005043", min.cells = 3, min.features = 100)

B_obj1$sample_id <- "21L005039"; B_obj2$sample_id <- "21L005043"
B_obj1 <- add_qc_metrics(B_obj1)
B_obj2 <- add_qc_metrics(B_obj2)

# Pre-filter QC plots (optional)
ggsave(file.path(res_dir, "B_prefilter_violin_S1.png"),
       VlnPlot(B_obj1, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb"), ncol = 5, pt.size = 0.05),
       width = 14, height = 4, dpi = 300)
ggsave(file.path(res_dir, "B_prefilter_violin_S2.png"),
       VlnPlot(B_obj2, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb"), ncol = 5, pt.size = 0.05),
       width = 14, height = 4, dpi = 300)

B_obj1 <- qc_filter(B_obj1)
B_obj2 <- qc_filter(B_obj2)

B <- merge(B_obj1, y = B_obj2, add.cell.ids = c("S1","S2"), project = "Leipzig_B_SoupX")
DefaultAssay(B) <- "RNA"
B <- JoinLayers(B) 
B <- run_seurat_basic(B, label = "B_basic")
B_sing <- run_doublets_scDblFinder(B, label = "B")
B_final <- run_seurat_final(B_sing, label = "B_soupx_dblclean")

plot_markers_and_save(B_final, label = "B_soupx_dblclean")


# ---------------------------
# 7) Side-by-side comparisons (A vs B)
# ---------------------------
cat("\n=== Creating side-by-side comparison plots ===\n")

# UMAP by sample
pA_samp <- DimPlot(A_final, reduction = "umap", group.by = "sample_id") + ggtitle("A: Filtered + dblclean")
pB_samp <- DimPlot(B_final, reduction = "umap", group.by = "sample_id") + ggtitle("B: SoupX + dblclean")
p_samp_side <- pA_samp | pB_samp
ggsave(file.path(res_dir, "COMPARE_umap_by_sample_A_vs_B.png"), p_samp_side, width = 14, height = 6, dpi = 300)

# UMAP by cluster
pA_clu <- DimPlot(A_final, reduction = "umap", label = TRUE) + ggtitle("A: clusters")
pB_clu <- DimPlot(B_final, reduction = "umap", label = TRUE) + ggtitle("B: clusters")
p_clu_side <- pA_clu | pB_clu
ggsave(file.path(res_dir, "COMPARE_umap_by_cluster_A_vs_B.png"), p_clu_side, width = 15, height = 6, dpi = 300)

# FeaturePlots side-by-side (adipocyte markers)
pA_adipo <- FeaturePlot(A_final, features = marker_adipo, ncol = 2) + plot_annotation(title = "A: Adipocyte markers")
pB_adipo <- FeaturePlot(B_final, features = marker_adipo, ncol = 2) + plot_annotation(title = "B: Adipocyte markers")
ggsave(file.path(res_dir, "COMPARE_feature_adipocyte_A_vs_B.png"), (pA_adipo | pB_adipo), width = 18, height = 7, dpi = 300)

# FeaturePlots side-by-side (ASPC markers)
pA_aspc <- FeaturePlot(A_final, features = marker_aspc, ncol = 3) + plot_annotation(title = "A: ASPC markers")
pB_aspc <- FeaturePlot(B_final, features = marker_aspc, ncol = 3) + plot_annotation(title = "B: ASPC markers")
ggsave(file.path(res_dir, "COMPARE_feature_ASPC_A_vs_B.png"), (pA_aspc | pB_aspc), width = 22, height = 7, dpi = 300)

# FeaturePlots side-by-side (macrophage markers)
pA_imm <- FeaturePlot(A_final, features = marker_immune, ncol = 2) + plot_annotation(title = "A: Macrophage markers")
pB_imm <- FeaturePlot(B_final, features = marker_immune, ncol = 2) + plot_annotation(title = "B: Macrophage markers")
ggsave(file.path(res_dir, "COMPARE_feature_macrophage_A_vs_B.png"), (pA_imm | pB_imm), width = 18, height = 7, dpi = 300)

# DotPlots side-by-side
dpA <- DotPlot(A_final, features = marker_panel) + RotatedAxis() + ggtitle("A: DotPlot")
dpB <- DotPlot(B_final, features = marker_panel) + RotatedAxis() + ggtitle("B: DotPlot")
ggsave(file.path(res_dir, "COMPARE_dotplot_A_vs_B.png"), (dpA / dpB), width = 14, height = 10, dpi = 300)

cat("DONE. All outputs saved to:\n", res_dir, "\n")


AverageExpression(A_final, features = c("ADIPOQ","PLIN1","FABP4","LPL","PDGFRA","DCN","COL1A1","LUM","PI16","CD55"))$RNA
AverageExpression(B_final, features = c("ADIPOQ","PLIN1","FABP4","LPL","PDGFRA","DCN","COL1A1","LUM","PI16","CD55"))$RNA
