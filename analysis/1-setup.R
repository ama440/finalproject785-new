setwd("~/Documents/1-UNC/1-Classes/BIOS785/finalproject785/")

# Libraries
library(readxl)
library(tidyverse)
library(Seurat)
library(patchwork)
library(data.table)

# Read in data
ctrl1 <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/ctrl1.csv", header = T)
ctrl1 <- setDF(ctrl1)
ctrl2 <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/ctrl2.csv", header = T)
ctrl2 <- setDF(ctrl2)
ctrl3 <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/ctrl3.csv", header = T)
ctrl3 <- setDF(ctrl3)
kd1 <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/kdkd1.csv", header = T)
kd1 <- setDF(kd1)
kd2 <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/kdkd2.csv", header = T)
kd2 <- setDF(kd2)
kd3 <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/kdkd3.csv", header = T)
kd3 <- setDF(kd3)

# Combine datasets
df <- cbind(ctrl1, ctrl2)
df <- cbind(df, ctrl3)
df <- cbind(df, kd1)
df <- cbind(df, kd2)
df <- cbind(df, kd3)
rm(ctrl1, ctrl2, ctrl3, kd1, kd2, kd3)

# Set names of combined dataset
rownames(df) <- df$GENE
df <- df %>% select(-GENE)

metadata <- fread("Data/metadata.tsv", header = T)
metadata <- setDF(metadata)
metadata <- metadata[-1,]
rownames(metadata) <- metadata$NAME
metadata <- metadata %>% select(-NAME)
metadata$disease_status <- ifelse(metadata$donor_id %in% c("CTRL1", "CTRL2", "CTRL3"), 
                                  "CTRL", "KDKD")

# Subset based on metadata, keeping all podocyte cells to improve power in de analysis
set.seed(123)
keep_ctrl1 <- sample(rownames(metadata %>% filter(cell_type__custom != "Podocyte" & donor_id == "CTRL1")), 3000)
keep_ctrl2 <- sample(rownames(metadata %>% filter(cell_type__custom != "Podocyte" & donor_id == "CTRL2")), 3000)
keep_ctrl3 <- sample(rownames(metadata %>% filter(cell_type__custom != "Podocyte" & donor_id == "CTRL3")), 3000)
keep_kd1 <- sample(rownames(metadata %>% filter(cell_type__custom != "Podocyte" & donor_id == "KDKD1")), 3000)
keep_kd2 <- sample(rownames(metadata %>% filter(cell_type__custom != "Podocyte" & donor_id == "KDKD2")), 3000)
keep_kd3 <- sample(rownames(metadata %>% filter(cell_type__custom != "Podocyte" & donor_id == "KDKD3")), 3000)
keep_pod <- rownames(metadata %>% filter(cell_type__custom == "Podocyte")) # rownames of all podocyte cells
keep <- c(keep_ctrl1, keep_ctrl2, keep_ctrl3, keep_kd1, keep_kd2, keep_kd3, keep_pod)
metadata <- metadata[keep,]
df <- df[,keep] # Extremely important: rownames of metadata must equal colnames of df
summary(rownames(metadata) == colnames(df)) # check whether we have all matches
write.csv(keep, file="cells_to_keep.csv", row.names=F)

# Create Seurat object
temp <- as(as.matrix(df), "sparseMatrix")
seu <- CreateSeuratObject(counts=temp, meta.data = metadata)
seu[["percent.mito"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
rm(df, temp)
seu

# Simple filtering step
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
seu <- subset(seu, percent.mito < .5)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

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
rm(seu)

# Normalize and identify variable features for each dataset independently
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seu.list)
anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features)
rm(seu.list)

# This command creates an 'integrated' data assay
seu.combined <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the 
# original unmodified data still resides in the original assay 
DefaultAssay(seu.combined) <- "integrated"

# Save Seurat object as RDS
saveRDS(seu.combined, "~/Documents/1-UNC/1-Classes/BIOS785/Data/integrated.rds")
