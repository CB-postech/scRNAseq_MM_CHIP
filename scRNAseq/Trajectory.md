Trajectory analysis for progenitor cells (Figure 3)

# Running Palantir
```R
import palantir
import scanpy as sc
import numpy as np
import pandas as pd
import os
import pickle
import rpy2
import harmonypy as hm

# Plotting
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Inline plotting
%matplotlib inline

# Reset random seed
np.random.seed(5)

# Load data
metadata = pd.read_csv("rdata/meta_Pro.csv", index_col = 0)
hvg_df=pd.read_csv("rdata/pro_hvg_normcounts_t.csv", index_col = 0)

# Running Palantir
start_cell = "BM5982_CD138Neg_AGTACTGCAGGTGAGT-1"

pcadims = [10, 30, 50, 70, 100]
dm = 10
tsne_perplexity = [100, 300, 500, 700, 1000]
num_wp = 500

filename = 'Pro'
savedir = 'R_analysis/Palantir_Pro/'
vars_use = ['sample_label']

for i in pcadims:
    
    pca,_ = palantir.utils.run_pca(hvg_df,n_components= i, use_hvg = False)
    ### Harmony
    ho = hm.run_harmony(pca, metadata, vars_use, max_iter_kmeans=50) # , sigma=50.0 defualt=0.1
    harmonyPCA=pd.DataFrame(ho.Z_corr)
    harmonyPCA=harmonyPCA.T
    harmonyPCA.index=hvg_df.index
    harmonyPCA.head()
    ####################
    dm_res = palantir.utils.run_diffusion_maps(harmonyPCA, n_components=dm)
    ms_data = palantir.utils.determine_multiscale_space(dm_res)
    ms_data.index = hvg_df.index
    imp_df=palantir.utils.run_magic_imputation(hvg_df, dm_res)
    with open(savedir + filename + '_ms_data_' + str(i) + '.pickle', 'wb') as f:
        pickle.dump(ms_data, f, pickle.HIGHEST_PROTOCOL)
    
    pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=num_wp)
    with open(savedir + filename + '_pr_res_' + str(i) + '.pickle', 'wb') as f:
        pickle.dump(pr_res, f, pickle.HIGHEST_PROTOCOL)
        
    pr_res.branch_probs.to_csv(savedir + filename + '_branch_probs_' + str(i) + '.csv')
    pr_res.pseudotime.to_csv(savedir + filename + '_pseudotime_' + str(i) + '.csv')
    
    for j in tsne_perplexity:
        
        tsne = palantir.utils.run_tsne(ms_data, perplexity = j)
        tsne.to_pickle(savedir + filename + '_tsne_' + str(i) + '_' + str(j) + '.p')
        
        tsne.to_csv(savedir + filename + '_tsne_' + str(i) + '_' + str(j) + '.csv')
        
        palantir.plot.highlight_cells_on_tsne(tsne, start_cell)
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_start_cell.png')

        palantir.plot.plot_palantir_results(pr_res, tsne)
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_pr_res.png')
        
        palantir.plot.plot_cell_clusters(tsne, metadata['sample_label_total'])
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_Condition.png')
        
        palantir.plot.plot_cell_clusters(tsne, metadata['celltype'])
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_celltype.png')
        
        cells = pr_res.branch_probs.columns
        palantir.plot.highlight_cells_on_tsne(tsne, cells)
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_path.png')

# Fine tuning
pcadims = [100]
dm = 10
tsne_perplexity = [700]
num_wp = [300, 400, 500, 600, 700, 800, 900, 1000]

filename = 'Pro'
savedir = 'R_analysis/Palantir_Pro/Fine_tuning/'

pca,_ = palantir.utils.run_pca(hvg_df,n_components= 100, use_hvg = False)

### Harmony
ho = hm.run_harmony(pca, metadata, vars_use, max_iter_kmeans=50) # , sigma=50.0 defualt=0.1
harmonyPCA=pd.DataFrame(ho.Z_corr)
harmonyPCA=harmonyPCA.T
harmonyPCA.index=hvg_df.index
harmonyPCA.head()
####################

dm_res = palantir.utils.run_diffusion_maps(harmonyPCA, n_components=10)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
imp_df=palantir.utils.run_magic_imputation(hvg_df, dm_res)

with open(savedir + filename + '_ms_data_' + str(i) + '.pickle', 'wb') as f:
    pickle.dump(ms_data, f, pickle.HIGHEST_PROTOCOL)

for i in num_wp:
    pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=i)
    with open(savedir + filename + '_pr_res_' + str(i) + '.pickle', 'wb') as f:
        pickle.dump(pr_res, f, pickle.HIGHEST_PROTOCOL)
        
    pr_res.branch_probs.to_csv(savedir + filename + '_branch_probs_' + str(i) + '.csv')
    pr_res.pseudotime.to_csv(savedir + filename + '_pseudotime_' + str(i) + '.csv')
    
    for j in tsne_perplexity:
        
        tsne = palantir.utils.run_tsne(ms_data, perplexity = j)
        tsne.to_pickle(savedir + filename + '_tsne_' + str(i) + '_' + str(j) + '.p')
        
        tsne.to_csv(savedir + filename + '_tsne_' + str(i) + '_' + str(j) + '.csv')
        
        palantir.plot.highlight_cells_on_tsne(tsne, start_cell)
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_start_cell.png')

        palantir.plot.plot_palantir_results(pr_res, tsne)
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_pr_res.png')
        
        palantir.plot.plot_cell_clusters(tsne, metadata['sample_label_total'])
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_Condition.png')
        
        palantir.plot.plot_cell_clusters(tsne, metadata['celltype'])
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_celltype.png')
        
        cells = pr_res.branch_probs.columns
        palantir.plot.highlight_cells_on_tsne(tsne, cells)
        plt.savefig(savedir + filename + '_pca_' + str(i) + '_' + str(j) + '_tsne_path.png')
```

# Trajectory analysis in R
```R
# pc 100 perplexity 700 num_wp 600
tsne_600_700 <- read.csv(paste0(ranaldir, 'Palantir_Pro/Fine_tuning/Pro_tsne_600_700.csv'), row.names = 'X')
colnames(tsne_600_700) <- c('Palantir_1', 'Palantir_2')
tsne_600_700 <- as.matrix(tsne_600_700)
MM_hs_seurat3.final_pro[["palantir"]] <- CreateDimReducObject(embeddings = tsne_600_700, key = 'Palantir_', assay = 'originalexp')

# Fig.3f
DimPlot(MM_hs_seurat3.final_pro,
             reduction = "palantir",
             group.by = "celltype",
             label = FALSE,
             raster = FALSE,
             cols = c('B.pro' = '#ffba08', 'HSC' = '#ffeeaa', 'M.pro' = '#d9ed92', 'E.pro' = '#bbdefb')) +
  theme(axis.line = element_blank(), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')

palantir_df <- Embeddings(object = MM_hs_seurat3.final_pro[["palantir"]])[(colnames(MM_hs_seurat3.final_pro)), c(1,2)]
palantir_df <- as.data.frame(palantir_df)
palantir_df$sample_label_total <- MM_hs_seurat3.final_pro$sample_label_total
palantir_df$celltype <- MM_hs_seurat3.final_pro$celltype

# Fig.3g
ggplot(palantir_df, aes(x = Palantir_1, y = Palantir_2)) +
  geom_point() +
  geom_density_2d_filled(alpha = 0.75) +
  facet_wrap(vars(sample_label_total)) +
  theme_classic() +
  theme(legend.position = 'none', axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())

branch_probs_600 <- read.csv(paste0(ranaldir, 'Palantir_Pro/Fine_tuning/Pro_branch_probs_600.csv'), row.names = 'X')
MM_hs_seurat3.final_pro$palantir_branch_E.pro <- branch_probs_600[colnames(MM_hs_seurat3.final_pro),]$BM5982_CD138Pos_GCCGTGATCGGCATAT
MM_hs_seurat3.final_pro$palantir_branch_B.pro <- branch_probs_600[colnames(MM_hs_seurat3.final_pro),]$BM5991_CD138Neg_AGCGTATGTGAATGTA
MM_hs_seurat3.final_pro$palantir_branch_M.pro <- branch_probs_600[colnames(MM_hs_seurat3.final_pro),]$BM6505_CD138Neg_CTACCTGCACTGGAAG

# Fig.3h
ggplot(metadata, aes(x = sample_label_total, y = palantir_branch_M.pro, fill = sample_label_total)) +
  scale_fill_manual(name = 'sample_label_total', guide = 'legend', values = c('#fed9b7', '#f07167')) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  theme_classic()
```
