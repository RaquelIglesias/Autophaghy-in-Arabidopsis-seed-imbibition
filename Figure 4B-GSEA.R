#KOATG7

#Dataset -? 
#open counts file and metadata
counts <- read.csv("gene_count_6.tsv", sep = "\t", header = TRUE)
samples <- read.csv("samples_6.tsv", sep = "\t", header = TRUE)

row.names(counts) <- counts$gene_id
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

write.table(results_noNA, file = "results_Figure_4.tsv", 
            sep = "\t", quote =FALSE, 
            col.name = TRUE, row.names = TRUE)


BiocManager::install(version = "3.17")
BiocManager::install("org.At.tair.db")
BiocManager::install("clusterProfiler")


#GSEA analysis
library(clusterProfiler)
library(org.At.tair.db)

sorted_vector_L2FC <- sort(results_noNA$log2FoldChange, na.last = TRUE, decreasing = TRUE)

order_vector <- order(results_noNA$log2FoldChange, na.last = TRUE, decreasing = TRUE)

results_clean_ordered <- results_noNA[order_vector,]

geneIDs <- rownames(results_clean_ordered)

#names(x) <- vector to assign names to the object
names(sorted_vector_L2FC) <- geneIDs

head(sorted_vector_L2FC)
gseGO_result <- gseGO(sorted_vector_L2FC, org.At.tair.db, keyType = "TAIR", 
               ont = "BP", eps = 0)
gseGO_result <- simplify(gseGO_result)
gseGO_df <- as.data.frame(gseGO_result)
head(gseGO_df)
write.table(gseGO_df, file = "Figure_4B_GSEA_GO_df.tsv", 
            sep = "\t", quote =FALSE, 
            col.name = TRUE, row.names = FALSE)

# check top and bottom NES scores
sig_gseGO <- gseGO_df[gseGO_df$p.adjust < 0.05,]
View(sig_gseGO[order(sig_gseGO$NES),])

install.packages("ggridges")
ridgeplot(gseGO_result, showCategory = 20, label_format = 30)
?ridgeplot


ridgeplot(gseGO_result, showCategory = 15, label_format = 30)




