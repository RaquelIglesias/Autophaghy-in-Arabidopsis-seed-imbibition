#KOATG7

#Dataset -? 
#open counts file and metadata
counts <- read.csv("gene_count_6.tsv", sep = "\t", header = TRUE)
samples <- read.csv("samples_6.tsv", sep = "\t", header = TRUE)

rownames(counts) <- counts$gene_id
counts <- counts[,2:dim(counts)[2]]

#DESeq2 dataseq
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData = samples,
                              design = ~ Mutant)

#run DGE analysis
deseq2_result <- DESeq(dds)

results_clean <- results(deseq2_result, contrast = c("Mutant", "Yes", "No"))

#clean rowns with NA's 
results_noNA <- results_clean[complete.cases(results_clean), ]



library(org.At.tair.db)
res.df <- as.data.frame(results_noNA)
res.df$symbol <- mapIds(org.At.tair.db, keys = rownames(res.df),keytype = "TAIR",
                        column = "SYMBOL")
res.df

install.packages("ggplot2")
library(ggplot2)

library(EnhancedVolcano)

EnhancedVolcano(res.df,lab = NA, xlim = c(-5,5), ylim = c(0,200), x = "log2FoldChange", 
                y = "padj", pCutoff = 10e-16, FCcutoff = 2)



EnhancedVolcano(res.df,lab = res.df$symbol, selectLab = c('CRC','AtCLO1', 'ALPHA-TIP',
                'SEIPIN1','AtCRT1b','ATPDI10'), 
                xlim = c(-5,5), ylim = c(0,200), x = "log2FoldChange", 
                y = "pvalue", pCutoff = 10e-16, FCcutoff = 2, drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black', labSize = 6.0,
                labCol = 'black',
                labFace = 'italic',
                boxedLabels = TRUE)

write.table(res.df, file = "Figure 4A_dataframe_.tsv", 
            sep = "\t", quote =FALSE, 
            col.name = TRUE, row.names = TRUE)


