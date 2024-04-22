To infer the activity of transcription factors, we utilized Dorothea. (Figure 4)

```R
### Loading packages -----
library(Seurat)
library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

### Sampling -----
clusters = as.vector(MM_hs_seurat3.final$celltype)
names(clusters) = names(MM_hs_seurat3.final$celltype)

cells = c()
for (i in names(table(clusters))) {
    if (length(names(clusters)[clusters == i]) > 1000) {
        set.seed(507)
        cells_sub = sample(names(clusters)[clusters == i], 1000)
        } else {
        cells_sub = names(clusters)[clusters == i]
        }
    cells = c(cells, cells_sub)
}

expr = MM_hs_seurat3.final@assays$originalexp@data
exprs_sample = expr[,cells]

MM_hs_seurat3.final$sampling = colnames(MM_hs_seurat3.final) %in% colnames(exprs_sample)
Idents(MM_hs_seurat3.final) <- MM_hs_seurat3.final$sampling
MM_hs_seurat3.final_sampling <- subset(MM_hs_seurat3.final, idents = TRUE)

### Dorothea -----
Idents(MM_hs_seurat3.final_sampling) <- MM_hs_seurat3.final_sampling$celltype

dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))
MM_hs_seurat3.final_sampling <- run_viper(MM_hs_seurat3.final_sampling, assay_key = 'originalexp', regulons = regulon,
                                          options = list(method = 'scale', minsize = 4, eset.filter = FALSE, verbose = FALSE))
viper_scores_df_sampling <- GetAssayData(MM_hs_seurat3.final_sampling, slot = "data", assay = "dorothea") %>%
  data.frame() %>%
  t()
CellsClusters <- data.frame(cell = names(Idents(MM_hs_seurat3.final_sampling)), 
                            cell_type = as.character(Idents(MM_hs_seurat3.final_sampling)),
                            stringsAsFactors = FALSE)
CellsClusters$cell <- gsub("-", ".", CellsClusters$cell)
viper_scores_clusters_sampling <- viper_scores_df_sampling  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

### TF activity -----
dorothea_NFKB1 <- viper_scores_clusters_sampling[viper_scores_clusters_sampling$tf == 'NFKB1',]
dorothea_NFKB1$cell <- gsub("[.]", "-", dorothea_NFKB1$cell)
rownames(dorothea_NFKB1) <- dorothea_NFKB1$cell
MM_hs_seurat3.final_sampling$dorothea_NFKB1 <- dorothea_NFKB1[colnames(MM_hs_seurat3.final_sampling),]$activity

dorothea_RELA <- viper_scores_clusters_sampling[viper_scores_clusters_sampling$tf == 'RELA',]
dorothea_RELA$cell <- gsub("[.]", "-", dorothea_RELA$cell)
rownames(dorothea_RELA) <- dorothea_RELA$cell
MM_hs_seurat3.final_sampling$dorothea_RELA <- dorothea_RELA[colnames(MM_hs_seurat3.final_sampling),]$activity
```
