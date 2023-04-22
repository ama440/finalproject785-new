# IMPORTANT: switch back to RNA because it is not advisable to use integrated
# data for differential expression analysis
DefaultAssay(seu.combined) <- "RNA"
# To switch back to integrated data, use DefaultAssay(seu.combined) <- "integrated"

# Set levels of seurat object;
# Goal is to compare gene expression between kd and control mice within individual cells
Idents(seu.combined) <- "disease_status"
levels(seu.combined)

# List of cell types
unique(seu.combined$cell_type__custom)

#### Differential Expression Analysis within cell types ####
# Podocyte
marker_genes <- FindMarkers(subset(seu.combined, subset = cell_type__custom == "Podocyte"), 
                             ident.1 = "CTRL", ident.2 = "KDKD")
marker_genes$p_val_adj = p.adjust(marker_genes$p_val, method='fdr')
head(marker_genes %>% arrange(desc(avg_log2FC)), 15)
head(marker_genes, 10)

celltypes <- unique(seu.combined$cell_type__custom)
ngenes <- 15
for (celltype in celltypes) {
  # Find marker genes
  marker_genes <- FindMarkers(subset(seu.combined, subset = cell_type__custom == celltype), 
                              ident.1 = "CTRL", ident.2 = "KDKD")
  # Use Benjamini & Hochberg procedure to adjust p-values; less conservative than default Bonferroni correction
  marker_genes$p_val_adj = p.adjust(marker_genes$p_val, method='fdr')
  # Write dataframes to csv
  if (celltype == "TAL/DCT") {celltype = "TAL-DCT"}
  write.csv(head(marker_genes, ngenes), 
            file.path("markers_by_cell", paste0(celltype, ".csv")))
}

# Read in csv's of interest
pod_markers <- read_csv("markers_by_cell/Podocyte.csv")
pod_markers$...1

# Ridge plots - Visualize single cell expression distributions in each cluster
ridge_pod <- RidgePlot(subset(seu.combined, subset = cell_type__custom == "Podocyte"), 
                       features = pod_markers$...1, ncol = 5)
ggsave("figures/ridge_podocyte.png", ridge_pod, width = 15, height = 9)


Idents(seu.combined) <- "donor_id"
DotPlot(subset(seu.combined, subset = cell_type__custom == "Podocyte"),
               features = c("Mapk1","Raf1","Braf")) + RotatedAxis()
DotPlot(subset(seu.combined, subset = cell_type__custom %in% c("PT-S1", "PT-S2", "PT-S3")),
        features = c("Mapk1","Raf1","Braf")) + RotatedAxis()
