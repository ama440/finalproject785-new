seu.combined <- readRDS("~/Documents/1-UNC/1-Classes/BIOS785/Data/integrated.rds")

# Run the standard workflow for visualization and clustering
seu.combined <- ScaleData(seu.combined, verbose = FALSE)
seu.combined <- RunPCA(seu.combined, npcs = 30, verbose = FALSE)
seu.combined <- RunUMAP(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- RunTSNE(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- FindNeighbors(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- FindClusters(seu.combined, resolution = 0.2)

# PCA
pca <- DimPlot(seu.combined, reduction = "pca", group.by = "donor_id") + 
  ggtitle("PCA Plot")
ggsave("figures/PCA.png", pca)

cluster_df <- data.frame(cbind(Idents(seu.combined), seu.combined$cell_type__custom))
names(cluster_df) <- c("assigned_cluster", "donor_id")
cluster_df$assigned_cluster <- as.numeric(cluster_df$assigned_cluster)
table(cluster_df$assigned_cluster, cluster_df$donor_id)

# Create UMAP
umap_celltype <- DimPlot(seu.combined, reduction = "umap", group.by = "cell_type__custom") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_celltype.png", umap_celltype)

umap_donor <- DimPlot(seu.combined, reduction = "umap", group.by = "donor_id") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_donor.png", umap_donor)

umap_disease <- DimPlot(seu.combined, reduction = "umap", group.by = "disease_status") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_disease.png", umap_disease)

elbow <- ElbowPlot(object=seu.combined, ndims = 30)
ggsave("figures/elbow_plot.png", elbow)
# elbow appears around 15 PCs hence why I used 15 for the umaps 

# Create t-SNE
tsne_celltype <- DimPlot(seu.combined, reduction = "tsne", group.by = "cell_type__custom") + 
  ggtitle("t-SNE Plot")
ggsave("figures/TSNE_celltype.png", tsne_celltype)

DimPlot(seu.combined, reduction = "tsne", group.by = "disease_status") + 
  ggtitle("t-SNE Plot")


### Do the same for SAVER-imputed Seurat object ###
seu.combined.s <- readRDS("~/Documents/1-UNC/1-Classes/BIOS785/Data/integrated_saver.rds")
# Run the standard workflow for visualization and clustering
seu.combined.s <- ScaleData(seu.combined.s, verbose = FALSE)
seu.combined.s <- RunPCA(seu.combined.s, npcs = 30, verbose = FALSE)
seu.combined.s <- RunUMAP(seu.combined.s, reduction = "pca", dims = 1:15)
# seu.combined.s <- RunTSNE(seu.combined.s, reduction = "pca", dims = 1:15)
seu.combined.s <- FindNeighbors(seu.combined.s, reduction = "pca", dims = 1:15)
seu.combined.s <- FindClusters(seu.combined.s, resolution = 0.2)

# PCA
pca.s <- DimPlot(seu.combined.s, reduction = "pca", group.by = "donor_id") + 
  ggtitle("PCA Plot")
ggsave("figures/PCA_saver.png", pca.s)

cluster_df <- data.frame(cbind(Idents(seu.combined.s), seu.combined.s$cell_type__custom))
names(cluster_df) <- c("assigned_cluster", "donor_id")
cluster_df$assigned_cluster <- as.numeric(cluster_df$assigned_cluster)
table(cluster_df$assigned_cluster, cluster_df$donor_id)

# Create UMAP
umap_celltype.s <- DimPlot(seu.combined.s, reduction = "umap", group.by = "cell_type__custom") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_celltype_saver.png", umap_celltype.s)

umap_donor.s <- DimPlot(seu.combined.s, reduction = "umap", group.by = "donor_id") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_donor_saver.png", umap_donor.s)

umap_disease.s <- DimPlot(seu.combined.s, reduction = "umap", group.by = "disease_status") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_disease_saver.png", umap_disease.s)

elbow.s <- ElbowPlot(object=seu.combined.s, ndims = 30)
ggsave("figures/elbow_plot_saver.png", elbow.s)
# elbow appears around 15 PCs hence why I used 15 for the umaps 

# Create t-SNE
tsne_celltype <- DimPlot(seu.combined.s, reduction = "tsne", group.by = "cell_type__custom") + 
  ggtitle("t-SNE Plot")
ggsave("figures/TSNE_celltype_saver.png", tsne_celltype)

DimPlot(seu.combined.s, reduction = "tsne", group.by = "disease_status") + 
  ggtitle("t-SNE Plot")
