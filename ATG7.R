#KOATG7

#Dataset -? 
#open counts file and metadata
counts <- read.csv("gene_count.tsv", sep = "\t", header = TRUE)
samples <- read.csv("samples_2.tsv", sep = "\t", header = TRUE)

row.names(counts) <- counts$gene_id
counts <- counts[,2:dim(counts)[2]]

samples$Time <- factor(samples$Time)

#DESeq2 dataseq
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData = samples,
                              design = ~ Treatment + Mutant + Time)

#Clustering

# apply VST transformation
dds_vst <- vst(dds)
vst_clean <- assay(dds_vst)
# make a new version of the VST values 
# where you only have the 500 most variable genes
vars <- apply(vst_clean, 1, var)
sorted_vars <- sort(vars, decreasing = TRUE)[1:500]
topgenes <- names(sorted_vars)
top_vars <- vst_clean[topgenes,]
top_vars_scaled <- t(scale(t(top_vars)))


# make heatmap
install.packages("gplots")
library(gplots)
palette <- hcl.colors(100, palette = "RdBu")
heatmap.2(top_vars_scaled, trace = "none",
          scale = "none", col = palette)

dev.off()
mutant_colours <- rep("steelblue", 30)
mutant_colours[samples$Mutant == "Yes"] <- "coral1" 
  heatmap.2(top_vars_scaled, trace = "none",
            scale = "row", col = palette,
            ColSideColors = mutant_colours, cexRow = 1.2,cexCol = 1.2, 
            margins = c(10,10)) 
  

# PCA
pca <- prcomp(t(top_vars_scaled)) 
scores <- pca$x
dim(scores)
head(scores)
           
install.packages("ggplot2")
library(ggplot2) 

# PCA Mutant
library(factoextra)
fviz_pca_ind(pca, axes = c(1,2),
col.ind = samples$Mutant,
label = "none", pointsize = 8) +
theme(text = element_text(size = 24),
axis.title = element_text(size = 24),
axis.text = element_text(size = 24))

         
?fviz_pca_ind    
fviz_eig(pca) + theme(text = element_text(size = 25),
axis.title = element_text(size = 25),
axis.text = element_text(size = 25))
            
loadings <- pca$rotation
dim(loadings) 

# 500 genes, 30 PCs
head(loadings)
            
            
            