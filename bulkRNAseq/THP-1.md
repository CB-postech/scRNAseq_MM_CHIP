Bulk RNA sequencing analysis of the THP-1 cell line.

```R
### Load packages -----
library(DESeq2)
library(calibrate)
library(ggrepel)

### Load data -----
files = list.files("./")
files = files[grep("count.txt", files)]

sampleTable <- data.frame(sampleName = strsplit(files, "_") %>% sapply(extract2, 1),
                          fileName = files,
                          condition = c("NT1", "NT1", "TET2_KD","TET2_KD","TET2_KD","TET2_KD"),
                          sample = c("NT1", "NT1", "TET2-1","TET2-1","TET2-2","TET2-2"))

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "./",
                                       design= ~ condition)

keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
ddsHTSeq <- DESeq(ddsHTSeq)

### DEG -----
res <- results(ddsHTSeq, name=resultsNames(ddsHTSeq)[2])
res = res[!is.na(res$padj),]
resOrdered <- res[order(res$pvalue),] # ordering based on p-value
resSig <- subset(resOrdered, padj < 0.05) # significant genes
resOrdered$gene = ensemblGenes[match(rownames(resOrdered), ensemblGenes[,1],nomatch = 0),2]
```
