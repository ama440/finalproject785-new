library(GSA, quietly=TRUE, warn.conflicts=FALSE, verbose=FALSE)
# library(devtools)
# install_github("ctlab/fgsea")
library(fgsea)
library(data.table)
library(ggplot2)

KEGGpathways=fgsea::gmtPathways("c2.cp.v6.2.symbols.gmt")
mouse_cp=fgsea::gmtPathways("m2.cp.v2023.1.Mm.symbols.gmt")
# Find mouse pathways here: https://www.gsea-msigdb.org/gsea/msigdb/mouse/genesets.jsp?collection=CP

# Set up object appropriately
DefaultAssay(seu.combined) <- "RNA"
Idents(seu.combined) <- "disease_status"
levels(seu.combined)

# Obtain a ranked list of genes; rank by -log10(pval)*sign(FC)
marker_genes <- FindMarkers(seu.combined, 
                            ident.1 = "CTRL", ident.2 = "KDKD", 
                            logfc.threshold   = -Inf,
                            min.pct           = -Inf,
                            min.cells.feature = 0,
                            min.cells.group   = 0)
write.csv(marker_genes, "marker_genes.csv")

# Because of the way this function is setup, 
# positive avg_logFC indicates higher expression in ctrl mice, and
# negative avg_logFC indicates higher expression in kdkd mice

marker_genes <- read.csv("marker_genes.csv")
marker_genes$p_val_adj <- p.adjust(marker_genes$p_val, method='fdr')
marker_genes$rank_metric <- -log10(marker_genes$p_val)*sign(marker_genes$avg_log2FC)
marker_genes[which(is.na(marker_genes$rank_metric)), "rank_metric"] <- 0
marker_genes <- marker_genes %>% arrange(rank_metric)
ranks <- marker_genes$rank_metric
names(ranks) <- marker_genes$X
ranks[1] <- -300 # findmarkers gives -inf, so set to very low value manually

set.seed(87)
fgseaRes <- fgsea(pathways = mouse_cp, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval),])

running_sum <- plotEnrichment(mouse_cp$REACTOME_FORMATION_OF_A_POOL_OF_FREE_40S_SUBUNITS, ranks) + 
  labs(title="Pathway: REACTOME_FORMATION_OF_A_POOL_OF_FREE_40S_SUBUNITS")
ggsave("figures/running_sum_reactome.png", running_sum, width = 8, height = 4)

plotEnrichment(mouse_cp$REACTOME_NONSENSE_MEDIATED_DECAY_NMD, ranks) + 
  labs(title="Reactome nonsense mediated decay nmd")

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                      mouse_cp, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(NES), pathway]
gsea_table <- plotGseaTable(mouse_cp[mainPathways], ranks, fgseaRes, 
                            gseaParam = 0.5, pathwayLabelStyle = list(size=10))
ggsave("figures/gsea_table.png", gsea_table, width = 15, height = 10)

# NES: normalized enrichment score

mouse_cp$REACTOME_INTERLEUKIN_RECEPTOR_SHC_SIGNALING
mouse_cp$REACTOME_FORMATION_OF_A_POOL_OF_FREE_40S_SUBUNITS
mouse_cp$WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS

