After removing putative doublets and low quality cells (as described in manuscript), remaining cells were annotated as below (Figure 3)

```R
### Loading packages -----
library(Seurat)
library(scater)
library(scran)
library(SingleCellExperiment)
library(harmony)
library(ggplot2)

### Visulization functions -----
featureplot_gy <- function(seurat, feature, type = "ensembl") {
  feature_id <- ensemblGenes[ensemblGenes$external_gene_name == feature,]$ensembl_gene_id[1]
  
  if(sum(feature_id %in% rownames(seurat)) == 0) {
    print(paste0("no ", feature))
  } else {
    FeaturePlot(seurat, features = feature_id, pt.size = 0.5, reduction = "umap") + ggtitle(feature) +
      theme(title = element_text(size=20), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  }
}

### Annotation -----
 # Plasma
featureplot_gy(MM_hs_seurat3.final, "JCHAIN")
featureplot_gy(MM_hs_seurat3.final, "PRDM1")

# proB
featureplot_gy(MM_hs_seurat3.final, "MME") # CD10

# preB
featureplot_gy(MM_hs_seurat3.final, "MS4A1") # CD20

# pDC
featureplot_gy(MM_hs_seurat3.final, "CLEC4C")
featureplot_gy(MM_hs_seurat3.final, "IL3RA")
featureplot_gy(MM_hs_seurat3.final, "LILRA4")

# Macrophage
featureplot_gy(MM_hs_seurat3.final, "CD68")
featureplot_gy(MM_hs_seurat3.final, "CD163")
featureplot_gy(MM_hs_seurat3.final, "MRC1")
featureplot_gy(MM_hs_seurat3.final, "C1QA")
featureplot_gy(MM_hs_seurat3.final, "C1QB")
featureplot_gy(MM_hs_seurat3.final, "C1QC")
featureplot_gy(MM_hs_seurat3.final, "VCAM1")

# HSC
featureplot_gy(MM_hs_seurat3.final, "CD34")

# Erythrocyte
featureplot_gy(MM_hs_seurat3.final, "HBB")
featureplot_gy(MM_hs_seurat3.final, "HBA1")
featureplot_gy(MM_hs_seurat3.final, "HBA2")

# cDC
featureplot_gy(MM_hs_seurat3.final, "CLEC9A")
featureplot_gy(MM_hs_seurat3.final, "CD1C")
featureplot_gy(MM_hs_seurat3.final, "CLEC10A")

# GMP
featureplot_gy(MM_hs_seurat3.final, "MPO")

# Monocyte
featureplot_gy(MM_hs_seurat3.final, "CD14")
featureplot_gy(MM_hs_seurat3.final, "FCGR3A")

# T
featureplot_gy(MM_hs_seurat3.final, "CD3D")
featureplot_gy(MM_hs_seurat3.final, "CD3G")
featureplot_gy(MM_hs_seurat3.final, "CD3E")

featureplot_gy(MM_hs_seurat3.final, "TRAC")
featureplot_gy(MM_hs_seurat3.final, "TRDC")

MM_hs_seurat3.final$celltype <- MM_hs_seurat3.final$seurat_clusters

annots = c('Myeloid', 'Plasma', 'T.NK', 'Plasma', 'Plasma', 'T.NK',
           'T.NK', 'Plasma', 'Plasma', 'Erythroid', 'Myeloid',
           'Plasma', 'Plasma', 'B.pre', 'Plasma', 'Myeloid',
           'Myeloid', 'Erythroid', 'pDC', 'Myeloid', 'B.pro',
           'Mac')
length(annots)

for(i in rev(1:length(annots))) {
  specify = paste0("^", i-1, "$")
  MM_hs_seurat3.final$celltype = gsub(specify, annots[i], MM_hs_seurat3.final$celltype)
}

UMAPPlot(MM_hs_seurat3.final, group.by = "celltype", label = T)

```
