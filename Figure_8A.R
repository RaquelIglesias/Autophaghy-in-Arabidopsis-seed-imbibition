counts <- read.csv("gene_count_Fig_8.tsv", sep = "\t", header = TRUE)
samples <- read.csv("samples_Fig_8.tsv", sep = "\t", header = TRUE)
samples$Treatment <- factor(samples$Treatment, levels = c("No", "ABA"))

row.names(counts) <- counts$gene_id
counts <- counts[,2:dim(counts)[2]]


#DESeq2 dataseq
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData = samples,
                              design = ~ Treatment + Mutant + Treatment:Mutant)

deseq2_result <- DESeq(dds, test = "LRT", reduce = ~ Treatment + Mutant)


resultsNames(deseq2_result)

results_clean <- results(deseq2_result, name = "TreatmentABA.MutantYes")
#clean rowns with NA's 
results_noNA <- results_clean[complete.cases(results_clean), ]

write.table(results_noNA, file = "results_Figure_8.tsv", 
            sep = "\t", quote =FALSE, 
            col.name = TRUE, row.names = TRUE)

#Vulcano Plot 
library(org.At.tair.db)
library(ggplot2)
library(EnhancedVolcano)

res.df <- as.data.frame(results_noNA)
res.df$symbol <- mapIds(org.At.tair.db, keys = rownames(res.df),keytype = "TAIR",
                        column = "SYMBOL")
res.df

write.table(res.df, file = "res.df_.tsv", 
            sep = "\t", quote =FALSE, 
            col.name = TRUE, row.names = TRUE)

EnhancedVolcano(results_noNA, lab = row.names(results_noNA), xlim = c(-10,10), x = "log2FoldChange", 
                selectLab = c('AT3G51810', 'AT1G69260', 'AT5G35660', 'AT4G25140',
                              'AT5G54270', 'AT3G11050', 'AT1G70670', 'AT2G37770', 'AT2G38170', 'AT2G21490'),
                y = "padj", pCutoff = 0.05, FCcutoff = 1, drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black', labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE)

#Some Symbols modified in res.df.tsv

res.df_1 <- read.csv("Figure_8A_ Vulgano_Symbol.tsv", sep = "\t", header = TRUE)



EnhancedVolcano(res.df_1,
                lab = res.df_1$symbol, 
                selectLab = c('EM1','AtCLO4', 'AT5G35660','ATFER2','AFP1',
                              'ATCAX1', 'OLE1', 'LEA', 'AKR4C9', 'LHCB3'), 
                xlim = c(-5,5), ylim = c(0,150), 
                x = "log2FoldChange", 
                y = "padj", 
                pCutoff = 0.05, 
                FCcutoff = 1, 
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black', 
                labSize = 8.0,       
                labCol = 'black', 
                labFace = 'italic',  
                boxedLabels = TRUE) 


