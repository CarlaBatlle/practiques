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


# Assignació correcta de l'estat d'expressió diferencial
res_padj$diffexpressed <- "NO"
res_padj$diffexpressed[res_padj$log2FoldChange > 0.6 & res_padj$pvalue < 0.05] <- "UP"
res_padj$diffexpressed[res_padj$log2FoldChange < -0.6 & res_padj$pvalue < 0.05] <- "DOWN"

# Comprovació de valors
table(res_padj$diffexpressed)  # Comprovar si hi ha valors "UP" i "DOWN"

# Eliminació de valors NA per evitar errors en la visualització
res_padj <- res_padj[!is.na(res_padj$log2FoldChange) & !is.na(res_padj$pvalue), ]

# Selecció dels gens més diferencialment expressats
head(res_padj[order(res_padj$padj), ][res_padj[order(res_padj$padj), ]$diffexpressed == 'DOWN', ])

# Filtrar només els gens amb log2FC < -0.5
res_padj_down <- subset(res_padj, log2FoldChange < -5)

# Gràfic de volcà amb anotació
ggplot(data = res_padj, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  # Afegir noms als gens amb log2FC < -0.5
  geom_text(data = res_padj_down, 
            aes(label = rownames(res_padj_down)), 
            vjust = 1.5, hjust = 0.5, size = 3, check_overlap = TRUE) +
  theme_minimal() 


# Mirar si coincideixen amb els de HI
hi_results <- read.csv("dea_results_2025-03-10T11-43-06.csv")
hi_results

intersect(rownames(res_padj), hi_results$Feature)

hi_results_adj<- subset(hi_results, Adjusted.p_value<0.05)

intersect<-intersect(rownames(res_padj), hi_results_adj$Feature)
summary(intersect)

#188


## FUNCTIONAL ANALYSIS

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(forcats)
library(ggupset)

signif_genes <- rownames(res_padj)
head(signif_genes)

ego <- enrichGO(gene = signif_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db,
                ont = "BP", # BP: Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = F
)
# use simplify() to remove redundant terms
ego2 <- simplify(ego, 
                 cutoff = 0.7,  #llindar del minim redundant
                 by = "p.adjust", 
                 select_fun = min)
# GRÀFICS
barplot(ego2, showCategory = 10) +
  theme(axis.text.y = element_text(size = 8))

dotplot(ego2, showCategory = 10) +
  theme(axis.text.y = element_text(size = 8))

#RICH FACTOR
ego3 <- mutate(ego2, 
               richFactor = Count/as.numeric(
                 sub("/\\d+", "", BgRatio)
               )
)

ggplot(ego3, showCategory = 10,
       aes(richFactor, fct_reorder(Description, richFactor))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"), 
                        trans = "log10",
                        guide = guide_colorbar(reverse = T, order = 1)) +
  scale_size_continuous(range = c(2,10)) +
  #theme_classic(text = element_text(size = 14)) +
  xlab("Rich Factor") +
  labs(y = NULL) +
  ggtitle("Biological Processes")


#Gene Set Enrichment Analysis
geneList <- res_padj$log2FoldChange
names(geneList) <- rownames(res_padj)
head(geneList)
geneList <- geneList[order(geneList, decreasing = T)]
ggo <- gseGO(geneList, keyType = "SYMBOL",
             OrgDb = org.Hs.eg.db, pvalueCutoff = 1)
ggo2 <- simplify(ggo, cutoff = 0.7, by = "p.adjust", select_fun = min)

ggo2 <- pairwise_termsim(ggo2)

#Gene-Concept Network
cnetplot(ego2, foldChange = geneList, showCategory = 2,
         cex_label_gene = 0.2,
         cex_label_category = 0.5)

#Heatplot
heatplot(ggo2, foldChange = geneList, showCategory = 10) +
  theme(axis.text.x = element_text(angle = 90, size = 5))

#Upseat plot
upsetplot(ego, n = 5)



