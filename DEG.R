## DESeq 2

library(DESeq2)
library(dplyr)

setwd("~/Documents/practiques/data")


countData <- read.csv("unproc_rnaseq.csv")

metaData <- read.csv('donor.csv', header = TRUE, sep = ",")
metaData


filtered_df <- subset(metaData, record_id %in% colnames(countData))
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=filtered_df, 
                              design= ~ diagnosis, tidy = TRUE)

dds
dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE))

#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)

#Filtrar Padj mes petit de 0,05
res_padj<- subset(res, padj<0.05)
res_padj



#Annotate genes
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(annotate)
library(clusterProfiler)

rownames(res_padj)

Change <- clusterProfiler:: bitr(rownames(res_padj),
     fromType = "ENTREZID",
     toType = "SYMBOL",
     OrgDb = org.Hs.eg.db,
     drop = TRUE)

rownames(res_padj) <- Change$SYMBOL
rownames(res_padj) 




#Volcanoplot 

library(ggpubr)

ggplot(data = res_padj, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point()


  #res_padj$diffexpressed <- "NO"
  #res_padj$diffexpressed[res_padj$log2FoldChange > 0.6 & res_padj$pvalue < 0.05] <- "UP"
  #res_padj$diffexpressed[res_padj$log2FoldChange < -0.6 & res_padj$pvalue < 0.05] <- "DOWN"
  #head(res_padj[order(res_padj$padj) & res_padj$diffexpressed == 'DOWN', ])
    
  
  #ggplot(data = res_padj, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) 
    #geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') 
    #geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') 
    #geom_point(size = 2) 
    #scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), # to set the colours of our variable<br /><br /><br />
                       #labels = c("Downregulated", "Not significant", "Upregulated")) # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)</p><br /><br />
    

# Mirar si coincideixen amb els de HI
hi_results <- read.csv("dea_results_2025-03-10T11-43-06.csv")
hi_results

intersect(rownames(res_padj), hi_results$Feature)
#577 intersect
#619 res_padj

hi_results_adj<- subset(hi_results, Adjusted.p_value<0.05)

intersect(rownames(res_padj), hi_results_adj$Feature)
#187