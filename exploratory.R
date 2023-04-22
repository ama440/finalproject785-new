setwd("~/Documents/1-UNC/1-Classes/BIOS785/Final Project")

# Libraries
library(readxl)
library(tidyverse)
library(Seurat)
library(patchwork)
library(data.table)

# Read in data
ctrl1 <- fread("Data/ctrl1.csv", header = T)
ctrl1 <- setDF(ctrl1)
ctrl2 <- fread("Data/ctrl2.csv", header = T)
ctrl2 <- setDF(ctrl2)
ctrl3 <- fread("Data/ctrl3.csv", header = T)
ctrl3 <- setDF(ctrl3)
kd1 <- fread("Data/kdkd1.csv", header = T)
kd1 <- setDF(kd1)
kd2 <- fread("Data/kdkd2.csv", header = T)
kd2 <- setDF(kd2)
kd3 <- fread("Data/kdkd3.csv", header = T)
kd3 <- setDF(kd3)

# Combine datasets; weird order but this is so that the combined df aligns with metadata
df <- cbind(ctrl1, kd2)
df <- cbind(df, ctrl3)
df <- cbind(df, ctrl2)
df <- cbind(df, kd3)
df <- cbind(df, kd1)
rm(ctrl1, ctrl2, ctrl3, kd1, kd2, kd3)

# Set names of combined dataset
rownames(df) <- df$GENE
df <- df %>% select(-GENE)

metadata <- fread("Data/metadata.tsv", header = T)
metadata <- setDF(metadata)
metadata <- metadata[-1,]
rownames(metadata) <- metadata$NAME
metadata <- metadata %>% select(-NAME)

# OPTIONAL: use only 2 kd mice and 2 control mice
keep_df <- c(grep("CTRL1", names(df)),
          grep("CTRL2", names(df)),
          grep("KDKD1", names(df)),
          grep("KDKD2", names(df)))
df <- df[,keep_df]

keep_meta <- c(grep("CTRL1", rownames(metadata)),
               grep("CTRL2", rownames(metadata)),
               grep("KDKD1", rownames(metadata)),
               grep("KDKD2", rownames(metadata)))
metadata <- metadata[keep_meta,]

# ALTERNATIVELY: use only 1 of each
keep_df <- c(grep("CTRL1", names(df)),
             grep("KDKD1", names(df)))
df <- df[,keep_df]

keep_meta <- c(grep("CTRL1", rownames(metadata)),
               grep("KDKD1", rownames(metadata)))
metadata <- metadata[keep_meta,]

# ALTERNATIVELY: use all 6 mice but sample cells randomly
set.seed(123)
# trainset <- sample(seq(1:nrow(boston)), size=0.8*nrow(boston), replace=FALSE)
# train <- boston[trainset,]
# test <- boston[-trainset,]

# Create Seurat object
temp <- as(as.matrix(df), "sparseMatrix")
seu <- CreateSeuratObject(counts=temp, meta.data = metadata)
seu[["percent.mito"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
rm(df, temp)
seu

# Filtering step
seu <- subset(seu, nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 5)
seu

# Define animal sex
orig.ident <- seu[["donor_id"]]
sex <- array(data = "", dim = dim(orig.ident))
sex[which(orig.ident == "CTRL1" | orig.ident == "CTRL2" | orig.ident == "KDKD1" | orig.ident == "KDKD2")] <- "MALE"
sex[which(orig.ident == "CTRL3" | orig.ident == "KDKD3")] <- "FEMALE"
seu$sex <- sex

# Normalize data
seu <- NormalizeData(seu)

# Scale data
# Note: need to find variable features first so that ScaleData can scale/center only them
# Otherwise, the command uses too much memory and fails
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seu), 10)
print(top10)
seu <- ScaleData(seu)

## Perform Data Integration
# Split the dataset into a list of six seurat objects
seu.list <- SplitObject(seu, split.by = "donor_id")

# Normalize and identify variable features for each dataset independently
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seu.list)
anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features)

# This command creates an 'integrated' data assay
seu.combined <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the 
# original unmodified data still resides in the original assay 
DefaultAssay(seu.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seu.combined <- ScaleData(seu.combined, verbose = FALSE)
seu.combined <- RunPCA(seu.combined, npcs = 15, verbose = FALSE)
seu.combined <- RunUMAP(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- FindNeighbors(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- FindClusters(seu.combined, resolution = 0.2)

# PCA
pca <- DimPlot(seu.combined, reduction = "pca", group.by = "donor_id") + 
  ggtitle("PCA Plot")
ggsave("PCA.png", pca)

cluster_df <- data.frame(cbind(Idents(seu.combined), metadata$donor_id))
names(cluster_df) <- c("assigned_cluster", "donor_id")
cluster_df$assigned_cluster <- as.numeric(cluster_df$assigned_cluster)
table(cluster_df$assigned_cluster, cluster_df$donor_id)

# Create UMAP
umap_celltype <- DimPlot(seu.combined, reduction = "umap", group.by = "cell_type__custom") + 
  ggtitle("UMAP Plot")
ggsave("UMAP_celltype.png", umap_celltype)

umap_donor <- DimPlot(seu.combined, reduction = "umap", group.by = "donor_id") + 
  ggtitle("UMAP Plot")
ggsave("UMAP_donor.png", umap_donor)

# Save Seurat object as RDS
saveRDS(seu.combined, "Data/int_ctrl1_kd1.rds")


ElbowPlot(object=seu.combined, ndims = 15)
FeaturePlot(seu.combined, features = c("Braf","Raf1","Mapk1"))
