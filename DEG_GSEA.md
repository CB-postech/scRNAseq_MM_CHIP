Obtaining differentially expressed genes (DEGs) with and without CHIP for each cell type and using them to perform gene set enrichment analysis (GSEA)

```R
# Extracting DEGs -----
Idents(MM_hs_seurat3.final) <- MM_hs_seurat3.final$celltype
celltype_list <- levels(MM_hs_seurat3.final$celltype)

for(i in 1:length(celltype_list)) {
    DEG <- FindMarkers(MM_hs_seurat3.final, ident.1 = "TET2", ident.2 = "Control", group.by = 'sample_label_total', subset.ident = celltype_list[i], test.use = "MAST", latent.vars = 'sample_label', logfc.threshold = 0)
    DEG$ensembl_gene_id <- rownames(DEG)
    DEG$Symbol = ensemblGenes[rownames(DEG),]$external_gene_name
    assign(paste0('DEG_MAST_final_', celltype_list[i]), DEG)
}

# GSEA -----
for (i in 1:length(celltype_list)) {
    DEG = get(paste0('DEG_MAST_final_', celltype_list[i]))
    DEG <- DEG[c(order(-DEG$avg_log2FC)),]
    ranks <- DEG$avg_log2FC
    names(ranks) <- rownames(DEG)
    assign(paste0("ranks_DEG_MAST_final_", celltype_list[i]), ranks)
}

library(fgsea)
msigdbr_df_H <- msigdbr(species = "human", category = "H")
msigdbr_df_C5 <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
pathways_H = split(x = msigdbr_df_H$ensembl_gene, f = msigdbr_df_H$gs_name)
pathways_C5 = split(x = msigdbr_df_C5$ensembl_gene, f = msigdbr_df_C5$gs_name)

for (i in 1:length(celltype_list)) {
  ranks = get(paste0("ranks_DEG_MAST_final_", celltype_list[i]))
  fgsea <- fgsea(pathways = pathways_H,
                 stats = ranks,
                 nproc = 1)
  assign(paste0("fgsea_DEG_MAST_final_", celltype_list[i], "_H"), fgsea)
}

for (i in 1:length(celltype_list)) {
  ranks = get(paste0("ranks_DEG_MAST_final_", celltype_list[i]))
  fgsea <- fgsea(pathways = pathways_C5,
                 stats = ranks,
                 nproc = 1)
  assign(paste0("fgsea_DEG_MAST_final_", celltype_list[i], "_C5"), fgsea)
}
```
