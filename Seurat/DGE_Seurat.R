library(Seurat)
library(tidyverse)

setwd("~/Differential-Gene-Expression-Analysis/Seurat/")

# read data in HDF5 format
counts1 <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')

# Preliminary Filters: 
# Read depth: cells /gene > 2 
# Low quality cells: Cells with #genes below user-determined threshold
# doublets / multiplets
seurat <- CreateSeuratObject(counts = counts1$`Gene Expression`, project = "NSCLC", min.cells = 3, min.features = 200)

# percentage of reads mapped to the mitochondrial genome
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")

# QC metrics: Visualize #genes / read depth / #mitchondrial reads across cells 
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

# Visualize feature-feature relationships
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_smooth(method = 'lm')

# filter cells based on number of genes and % mitochondrial reads
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization: 
# 1.Scaling: No. features in each cell are divided by the total counts for that cell and multiplied by the scale.factor(10000). 
# 2.Log transformation
seurat <- NormalizeData(seurat)

# Feature selection: find top 2000 genes with high cell-to-cell variation 
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)

# Visualize 20 most highly variable genes
top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot2

# Data Scaling: Prevent highly expressed genes from dominating the analysis
seurat <- ScaleData(seurat)

# Linear Dimensionality Reduction
seurat <- RunPCA(seurat, npcs = 50)

# Visualize:  To identify primary sources of heterogeneity in a dataset 
DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)

# decide which componenets to include in downstream "differential" analysis?
ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))

# Cluster cells : 
# Get Shared Neighbor Network
seurat <- FindNeighbors(seurat, dims = 1:20)
# Find Communities in network
seurat <- FindClusters(seurat, resolution = 0.1) # low resolution parameter is used to get a broad clustering

# Non-linear dimension reduction for visualization
seurat <- RunUMAP(seurat, dims = 1:20)

# Visualize clusters
DimPlot(seurat, label = TRUE)

# Find differentially expressed features across all clusters
pbmc.markers = FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers, file = "NSCLC_DEGs.csv")
saveRDS(seurat, file = "seurat_20kNSCLC.rds")


