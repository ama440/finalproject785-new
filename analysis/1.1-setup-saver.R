setwd("~/Documents/1-UNC/1-Classes/BIOS785/finalproject785/")

# Libraries
library(readxl)
library(tidyverse)
library(Seurat)
library(patchwork)
library(data.table)

# Read in data
ctrl1.s <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/SAVER_bios785/s_c1_sav.csv", 
               header = T)
ctrl1.s <- setDF(ctrl1.s)
ctrl2.s <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/SAVER_bios785/s_c2_sav.csv",
               header = T)
ctrl2.s <- setDF(ctrl2.s)
ctrl3.s <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/SAVER_bios785/s_c3_sav.csv",
               header = T)
ctrl3.s <- setDF(ctrl3.s)
kd1.s <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/SAVER_bios785/s_k1_sav.csv", 
             header = T)
kd1.s <- setDF(kd1.s)
kd2.s <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/SAVER_bios785/s_k2_sav.csv",
             header = T)
kd2.s <- setDF(kd2.s)
kd3.s <- fread("~/Documents/1-UNC/1-Classes/BIOS785/Data/SAVER_bios785/s_k3_sav.csv",
             header = T)
kd3.s <- setDF(kd3.s)

# Combine datasets
df.s <- cbind(ctrl1.s, ctrl2.s)
df.s <- cbind(df.s, ctrl3.s)
df.s <- cbind(df.s, kd1.s)
df.s <- cbind(df.s, kd2.s)
df.s <- cbind(df.s, kd3.s)
rm(ctrl1.s, ctrl2.s, ctrl3.s, kd1.s, kd2.s, kd3.s)

# Set names of combined dataset
rownames(df.s) <- df.s$V1
df.s <- df.s %>% select(-V1)

metadata <- fread("../Data/metadata.tsv", header = T)
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
df.s <- df.s[,keep] # Extremely important: rownames of metadata must equal colnames of df
summary(rownames(metadata) == colnames(df.s)) # check whether we have all matches
write.csv(keep, file="cells_to_keep.csv", row.names=F)

# Create Seurat object
temp <- as(as.matrix(df.s), "sparseMatrix")
seu.s <- CreateSeuratObject(counts=temp, meta.data = metadata)
seu.s[["percent.mito"]] <- PercentageFeatureSet(seu.s, pattern = "^mt-")
rm(df.s, temp)
seu.s

# Simple filtering step
VlnPlot(seu.s, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
seu.s <- subset(seu.s, percent.mito < .35)
VlnPlot(seu.s, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# Define animal sex
orig.ident <- seu.s[["donor_id"]]
sex <- array(data = "", dim = dim(orig.ident))
sex[which(orig.ident == "CTRL1" | orig.ident == "CTRL2" | orig.ident == "KDKD1" | orig.ident == "KDKD2")] <- "MALE"
sex[which(orig.ident == "CTRL3" | orig.ident == "KDKD3")] <- "FEMALE"
seu.s$sex <- sex

# Normalize data
seu.s <- NormalizeData(seu.s)

# Scale data
# Note: need to find variable features first so that ScaleData can scale/center only them
# Otherwise, the command uses too much memory and fails
seu.s <- FindVariableFeatures(seu.s, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seu.s), 10)
print(top10)
seu.s <- ScaleData(seu.s)

saveRDS(seu.s, "~/Documents/1-UNC/1-Classes/BIOS785/Data/seu_s.rds")

seu.s <- readRDS("~/Documents/1-UNC/1-Classes/BIOS785/Data/seu_s.rds")
## Perform Data Integration
# Split the dataset into a list of six seurat objects
seu.s.list <- SplitObject(seu.s, split.by = "donor_id")
rm(seu.s)

# Normalize and identify variable features for each dataset independently
seu.s.list <- lapply(X = seu.s.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features.s <- SelectIntegrationFeatures(object.list = seu.s.list)
anchors.s <- FindIntegrationAnchors(object.list = seu.s.list, anchor.features = features.s)
saveRDS(anchors.s, "~/Documents/1-UNC/1-Classes/BIOS785/Data/anchors_s.rds")
rm(seu.s.list)

# This command creates an 'integrated' data assay
seu.combined.s <- IntegrateData(anchorset = anchors.s)

# specify that we will perform downstream analysis on the corrected data note that the 
# original unmodified data still resides in the original assay 
DefaultAssay(seu.combined.s) <- "integrated"

# Save Seurat object as RDS
saveRDS(seu.combined.s, "~/Documents/1-UNC/1-Classes/BIOS785/Data/integrated_saver.rds")
