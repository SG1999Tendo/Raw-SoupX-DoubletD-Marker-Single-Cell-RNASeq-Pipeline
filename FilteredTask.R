# ============================================================
# # ---- Setup and extraction --===
#FILTERED 10x (MEX) -> QC -> CLUSTERING -> DOUBLETS -> CLEAN CLUSTERING
#
# Dataset: 2 patient-derived 10x Genomics scRNA-seq samples
#   - 21L005039
#   - 21L005043
#
# Goal of this script (Part 1):
#   1) Standard preprocessing and quality control (QC)
#   2) Perform doublet detection (scDblFinder) and remove predicted doublets
#   3) Re-run clustering on the cleaned (singlet-only) dataset
#
# NOTE:
# - This script uses FILTERED matrices only (Cell Ranger cell calls).
# - We keep parameters fixed + set seeds for reproducibility.
# - Later (Part 2), we will annotate adipocytes + ASPCs and do pseudobulk.
# ============================================================


# ---------------------------
# 0) User paths
# ---------------------------
zip_file <- "C:/Users/pednet/Desktop/Leipzig Task/Revised/Data_for_task.zip"
out_dir  <- "C:/Users/pednet/Desktop/Leipzig Task/Revised/Data_for_task_unzipped"
res_dir  <- "C:/Users/pednet/Desktop/Leipzig Task/Revised/results_filtered_dblclean"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# ---------------------------
# 1) Unzip and locate data root
# ---------------------------
# Why:
# The task data are provided as a zip archive. We extract it once and then read 10x matrices
# from the resulting folder structure.
unzip(zip_file, exdir = out_dir)

# Some zips contain a top-level folder; others extract directly.
base1 <- file.path(out_dir, "Data_for_task")
base  <- if (dir.exists(base1)) base1 else out_dir

# Quick sanity check: confirm expected folders exist
stopifnot(dir.exists(file.path(base, "21L005039")))
stopifnot(dir.exists(file.path(base, "21L005043")))

# ---------------------------
# 2) Load packages
# ---------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# Doublet detection packages (Bioconductor)
# Install once if missing:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("SingleCellExperiment", "scDblFinder"))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scDblFinder)
})

# ---------------------------
# 3) Reproducibility settings
# ---------------------------
# Why:
# UMAP and clustering can vary run-to-run unless seeds are fixed.
set.seed(123)
umap_seed <- 123
npcs_use  <- 50
dims_use  <- 1:30
res_use   <- 0.4

# ---------------------------
# 4) Define paths and read FILTERED 10x matrices
# ---------------------------
# Why:
# Filtered matrices contain cell-associated barcodes called by Cell Ranger and are the
# standard input for Seurat workflows.
s1_filt <- file.path(base, "21L005039", "filtered_feature_bc_matrix")
s2_filt <- file.path(base, "21L005043", "filtered_feature_bc_matrix")

m1_filt <- Read10X(data.dir = s1_filt)
m2_filt <- Read10X(data.dir = s2_filt)

cat("Dimensions (genes x cells):\n")
cat("21L005039 filtered:", dim(m1_filt), "\n")
cat("21L005043 filtered:", dim(m2_filt), "\n")

# ---------------------------
# 5) Create Seurat objects (one per sample)
# ---------------------------
# Why:
# Keeping samples separate initially lets us inspect QC per sample and avoid mixing issues early.
obj1 <- CreateSeuratObject(counts = m1_filt, project = "21L005039", min.cells = 3, min.features = 100)
obj2 <- CreateSeuratObject(counts = m2_filt, project = "21L005043", min.cells = 3, min.features = 100)

# Store sample identifiers (crucial later for pseudobulk)
obj1$sample_id <- "21L005039"
obj2$sample_id <- "21L005043"

# ---------------------------
# 6) Compute QC metrics
# ---------------------------
# Why:
# QC metrics help remove poor-quality cells and obvious artifacts.
# - nFeature_RNA: number of detected genes (too low = low quality; too high = potential doublet)
# - nCount_RNA: total UMIs (very high = potential doublet)
# - percent.mt: mitochondrial fraction (high = stressed/damaged cells; often low in snRNA)
# - percent.ribo / percent.hb: optional signals of composition/contamination
add_qc_metrics <- function(obj) {
  obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL")
  obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGQZ]")
  return(obj)
}

obj1 <- add_qc_metrics(obj1)
obj2 <- add_qc_metrics(obj2)

# ---------------------------
# 7) Pre-filter QC plots (save for documentation)
# ---------------------------
# Why:
# Showing QC distributions BEFORE filtering makes your filtering decisions defensible.
p_vln1 <- VlnPlot(obj1, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb"),
                  ncol = 5, pt.size = 0.05) + plot_annotation(title = "QC (pre-filter): 21L005039")
p_vln2 <- VlnPlot(obj2, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb"),
                  ncol = 5, pt.size = 0.05) + plot_annotation(title = "QC (pre-filter): 21L005043")

ggsave(file.path(res_dir, "QC_prefilter_violin_21L005039.png"), p_vln1, width = 14, height = 4, dpi = 300)
ggsave(file.path(res_dir, "QC_prefilter_violin_21L005043.png"), p_vln2, width = 14, height = 4, dpi = 300)

p_sc1 <- FeatureScatter(obj1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("21L005039: nCount vs nFeature")
p_sc2 <- FeatureScatter(obj1, feature1 = "nCount_RNA", feature2 = "percent.mt")   + ggtitle("21L005039: nCount vs %MT")
p_sc3 <- FeatureScatter(obj2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("21L005043: nCount vs nFeature")
p_sc4 <- FeatureScatter(obj2, feature1 = "nCount_RNA", feature2 = "percent.mt")   + ggtitle("21L005043: nCount vs %MT")

ggsave(file.path(res_dir, "QC_prefilter_scatter.png"), (p_sc1 | p_sc2) / (p_sc3 | p_sc4), width = 12, height = 10, dpi = 300)

# ---------------------------
# 8) QC filtering (broad first pass)
# ---------------------------
# Why:
# Remove obvious low-quality cells and extreme outliers. We start broad to avoid
# accidentally removing valid adipose cell types.
#
# NOTE:
# Your data show very low percent.mt, so the mito cutoff mainly removes extreme outliers.
qc_filter <- function(obj) {
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
}

obj1_f <- qc_filter(obj1)
obj2_f <- qc_filter(obj2)

cat("Cells retained after QC filtering:\n")
cat("21L005039:", ncol(obj1_f), "of", ncol(obj1), "\n")
cat("21L005043:", ncol(obj2_f), "of", ncol(obj2), "\n")

# ---------------------------
# 9) Merge samples (post-QC)
# ---------------------------
# Why:
# Merge
obj <- merge(obj1_f, y = obj2_f, add.cell.ids = c("S1","S2"), project = "Leipzig_filtered")

# IMPORTANT (Seurat v5): join layers BEFORE any downstream steps
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")

# Now do the initial workflow (this becomes your true baseline UMAP)
set.seed(123)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, random.seed = 123, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, seed.use = 123, verbose = FALSE)

# Save initial UMAPs
p_umap_samp_pre <- DimPlot(obj, reduction = "umap", group.by = "sample_id") +
  ggtitle("Pre-doublet: UMAP by sample")
p_umap_clu_pre  <- DimPlot(obj, reduction = "umap", label = TRUE) +
  ggtitle("Pre-doublet: UMAP by cluster")

ggsave(file.path(res_dir, "UMAP_pre_doublet_by_sample.png"), p_umap_samp_pre, width = 7, height = 6, dpi = 300)
ggsave(file.path(res_dir, "UMAP_pre_doublet_by_cluster.png"), p_umap_clu_pre, width = 8, height = 6, dpi = 300)

# ---------------------------
# 11) Doublet detection (scDblFinder)
# ---------------------------
# Fix Seurat v5 multi-layer issue BEFORE as.SingleCellExperiment()
  # makes "counts", "data", "scale.data" behave normally
# Optional: remove extremely low-count cells (helps scDblFinder stability)
obj <- subset(obj, subset = nCount_RNA >= 200 & nFeature_RNA >= 200)
# sanity check: should now include "counts" and "data" (not counts.<sample>)
Layers(obj[["RNA"]])
# Why:
# Doublets (2 cells captured in one droplet) can create hybrid expression profiles and
# artificially distort clustering/marker detection. We detect and remove them.
sce <- as.SingleCellExperiment(obj)

set.seed(123)
sce <- scDblFinder(
  sce,
  samples  = sce$sample_id,                      # model doublets per sample
  clusters = as.factor(obj$seurat_clusters)       # leverage cluster structure
)

obj$scDblFinder.score <- colData(sce)$scDblFinder.score
obj$scDblFinder.class <- colData(sce)$scDblFinder.class

dbl_tab <- table(obj$scDblFinder.class, obj$sample_id)
print(dbl_tab)
write.csv(as.data.frame(dbl_tab),
          file = file.path(res_dir, "doublet_table_scDblFinder.csv"),
          row.names = FALSE)

# Doublet diagnostic plots
p_dbl_class <- DimPlot(obj, reduction = "umap", group.by = "scDblFinder.class") +
  ggtitle("scDblFinder: singlet vs doublet")
p_dbl_score <- FeaturePlot(obj, features = "scDblFinder.score") +
  ggtitle("scDblFinder score")

ggsave(file.path(res_dir, "doublet_umap_class.png"), p_dbl_class, width = 7, height = 6, dpi = 300)
ggsave(file.path(res_dir, "doublet_umap_score.png"), p_dbl_score, width = 7, height = 6, dpi = 300)

p_dbl_vln <- VlnPlot(obj,
                     features = c("nFeature_RNA","nCount_RNA","percent.mt","scDblFinder.score"),
                     group.by = "scDblFinder.class",
                     ncol = 4, pt.size = 0.05) +
  plot_annotation(title = "QC metrics by scDblFinder class")
ggsave(file.path(res_dir, "doublet_qc_violin.png"), p_dbl_vln, width = 14, height = 4, dpi = 300)

# ---------------------------
# 12) Remove predicted doublets
# ---------------------------
# Why:
# We keep only singlets to produce a cleaner final embedding and more reliable markers.
obj_sing <- subset(obj, subset = scDblFinder.class == "singlet")

cat("Cells retained after doublet removal:", ncol(obj_sing), "\n")

# ---------------------------
# 13) Re-run final Seurat workflow (clean singlet-only object)
# ---------------------------
# Why:
# After removing doublets, the neighborhood graph changes. Re-running prevents clusters
# from being biased by cells that are no longer present.
obj_final <- NormalizeData(obj_sing, verbose = FALSE)
obj_final <- FindVariableFeatures(obj_final, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj_final <- ScaleData(obj_final, vars.to.regress = "percent.mt", verbose = FALSE)
obj_final <- RunPCA(obj_final, npcs = npcs_use, verbose = FALSE)
obj_final <- FindNeighbors(obj_final, dims = dims_use, verbose = FALSE)
obj_final <- FindClusters(obj_final, resolution = res_use, random.seed = 123, verbose = FALSE)
obj_final <- RunUMAP(obj_final, dims = dims_use, seed.use = umap_seed, verbose = FALSE)

# Save final UMAPs
p_umap_samp <- DimPlot(obj_final, reduction = "umap", group.by = "sample_id") +
  ggtitle("FINAL (singlets): UMAP by sample")
p_umap_clu  <- DimPlot(obj_final, reduction = "umap", label = TRUE) +
  ggtitle("FINAL (singlets): UMAP by cluster")

ggsave(file.path(res_dir, "UMAP_FINAL_by_sample.png"), p_umap_samp, width = 7, height = 6, dpi = 300)
ggsave(file.path(res_dir, "UMAP_FINAL_by_cluster.png"), p_umap_clu, width = 8, height = 6, dpi = 300)

# Save final object for downstream annotation + pseudobulk
saveRDS(obj_final, file = file.path(res_dir, "obj_filtered_dblclean_FINAL.rds"))

cat("\nDONE (Part 1).\n")
cat("Outputs saved in:\n", res_dir, "\n")
cat("Next (Part 2): marker discovery -> annotate adipocytes & ASPCs -> pseudobulk.\n")

