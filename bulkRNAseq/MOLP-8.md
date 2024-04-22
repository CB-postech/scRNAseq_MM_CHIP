Bulk RNA sequencing analysis of the MOLP-8 cell line.

```R
### Load packages -----
library(DESeq2)
library(calibrate)
library(ggrepel)
library(topGO)
library(plyr)

### Load data -----
files = list.files("/home/gycho/Project/MM/bulk_MOLP8/STAR/")
files = files[grep("count.txt", files)]

sampleTable <- data.frame(sampleName = files,
                          fileName = files,
                          condition = c('MOLP8_NT', 'MOLP8_NT', 'MOLP8_TET2_KD', 'MOLP8_TET2_KD'),
                          sample = c('MOLP8_NT1', 'MOLP8_NT2', 'MOLP8_TET2_KD1', 'MOLP8_TET2_KD2'))

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "/home/gycho/Project/MM/bulk_MOLP8/STAR/",
                                       design= ~ condition)

keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
ddsHTSeq <- DESeq(ddsHTSeq)

### DEG -----
res <- results(ddsHTSeq, name=resultsNames(ddsHTSeq)[2])
res = res[!is.na(res$pvalue),]
resOrdered <- res[order(res$pvalue),] # ordering based on p-value
resSig <- subset(resOrdered, pvalue < 0.05) # significant genes
resOrdered$gene = ensemblGenes[match(rownames(resOrdered), ensemblGenes[,1],nomatch = 0),2]

### Gene Ontology (GO) analysis -----
backgroundGenes = rownames(res)
targetUPGenes = rownames(res[res$pvalue < 0.1 & !is.na(res$pvalue) & res$log2FoldChange > 0,])
targetDOWNGenes = rownames(res[res$pvalue < 0.1 & !is.na(res$pvalue) & res$log2FoldChange < 0,])

allGeneUP = factor(as.integer(backgroundGenes %in% targetUPGenes))
allGeneDOWN = factor(as.integer(backgroundGenes %in% targetDOWNGenes))

names(allGeneUP) = backgroundGenes
names(allGeneDOWN) = backgroundGenes

onts = c('MF', 'BP', 'CC')
tab = as.list(onts)
names(tab) = onts

tab_UP <- tab
tab_DOWN <- tab

for (i in 1:length(onts)) {
    tgd = new('topGOdata', ontology = onts[i], allGenes = allGeneUP, nodeSize = 5, annot = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'ensembl')
    resultTopGO.elim = runTest(tgd, algorithm = 'elim', statistic = 'Fisher')
    tab_UP[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim, orderBy = 'Fisher.classic', topNodes = 200)
}

for (i in 1:length(onts)) {
    tgd = new('topGOdata', ontology = onts[i], allGenes = allGeneDOWN, nodeSize = 5, annot = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'ensembl')
    resultTopGO.elim = runTest(tgd, algorithm = 'elim', statistic = 'Fisher')
    tab_DOWN[[i]] = GenTable(tgd, Fisher.elim = resultTopGO.elim, orderBy = 'Fisher.classic', topNodes = 200)
}
```
