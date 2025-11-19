library(DESeq2)
library(dplyr)
library(plyr)
library(tidyverse)
library(ggpubr)
library(EnhancedVolcano)
library(ggplot2)
library(tidyverse)
library(CorLevelPlot)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(pathview)
library(AnnotationDbi)
library(AnnotationHub)
library(ggridges)
library(enrichplot)
library(KEGGREST)
library(EnrichmentBrowser)
library(pheatmap)


# load data
metaData = read.csv("All_phenotypes.csv")
countData<-read.csv("counts.csv",sep=",", header=T, check.names=F)

orig_names <- names(countData) # keep a back-up copy of the original names
geneID <- countData[,1] # Convert count data to a matrix of appropriate form that DEseq2 can read
countData <- as.matrix(countData[ , -1]) 
sampleIndex <- colnames(countData)
countData <- as.matrix(countData[,sampleIndex])
rownames(countData) <- geneID
countData2<-countData


# Calculate total number of columns
total_columns <- ncol(countData2)
zero_counts <- rowSums(countData2 == 0)
length(zero_counts)
table(zero_counts)
countData2 <- as.data.frame(countData2) %>% dplyr::filter(zero_counts < 0.7 * total_columns)
mycounts <- countData2

rownames(metaData) <- metaData$LambID 
metaData$LambID  <- factor(metaData$LambID)
colnames(countData2) == metaData$LambID 

# Select the columns from countData2 that are in metadata
countData2 <- countData2[, colnames(countData2) %in% metaData$LambID]

# construct PCA to see if any environmental effects or not
deseq2Data <- DESeqDataSetFromMatrix(countData=countData2, colData=metaData, design= ~ 1)# no model only for PCA
vsd <- vst(deseq2Data, blind = TRUE)
rlog <- rlog(deseq2Data, blind = TRUE)

# pca
colnames(metaData)
plotPCA(rlog, intgroup = "Treatment", ntop=15000)

pca <- prcomp(t(assay(vsd)))
pca_df <- data.frame(Sample = rownames(pca$x),
                     PC1 = pca$x[,1],
                     PC2 = pca$x[,2]) %>% left_join(metaData, by = c("Sample" = "LambID"))

ggplot(pca_df, aes(x=PC1, y=PC2, color=Ave_MO)) +
  geom_point(size = 4) +
  theme_bw()+
  labs(title = "PCA plot of the RNASeq samples")

# Differential gene expression analysis
deseq2Data <- DESeqDataSetFromMatrix(countData=countData2, colData=metaData, design= ~ DMI)
deseq2Data <- DESeq(deseq2Data)

resultsNames(deseq2Data)
pval = 0.05
lfc = 2
results = resultsNames(deseq2Data)
upresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"upDEGs"))
downresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"downDEGs"))

for(i in 2:length(results)){
  
  res = results(deseq2Data, 
                name = results[i])# independent filtering occurs in this step to save you from multiple test correction on genes with no power
  resorder <- res[order(res$padj),]
  upDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange > lfc))))
  downDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange < -lfc))))
  resSig = subset(resorder, padj < pval & log2FoldChange > lfc | padj < pval & log2FoldChange < -lfc)
  write.csv(resSig , file=paste0(out_dir,results[i],".0.05P.2LFC.updownDEGs.csv"), row.names = T)
  res_pval = subset(resorder, padj < pval )
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

summary(res)

# Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Create the volcano plot
EnhancedVolcano(res_df,
    lab = res_df$gene,
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 4,
    ylim=c(0,2),
    pointSize = 1.5,
    labSize = 2.5,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.7,
    legendLabels = c('NS','Log2FC','p-value','p-value & Log2FC'),
    legendPosition = 'right',
    title = 'Volcano Plot of DESeq2 Results',
    subtitle = 'Differential expression analysis on dry matter intake',
    caption = 'Thresholds: p < 0.05, |log2FC| > 2',
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    max.overlaps = Inf
)
 
# Over representation analysis
# extract annotation DB of sheep

ah <- AnnotationHub()
AnnotationHub::query(ah, c("Ovis", "aries"))
Oaries <- ah[["AH114633"]] # old one which worked with Biocversion 3.18
columns(Oaries)
res_original = res
resSig_original = resSig
# Taking the resSig variable from the above deseq2 as it has all the genes with corresponding pvalue and log2fc values
resSig$GeneID = rownames(resSig)
resSigEntrez <- AnnotationDbi::select(Oaries, keys =  rownames(resSig),
                                columns = c('ENTREZID','GENENAME'), keytype = 'SYMBOL') # Oaries is the annotation db created from the geneset_enrichment_analysis.R code
# Now replace the NA values in entrezid column with the values from first column. 
resSigEntrez$ENTREZID <- ifelse(is.na(resSigEntrez$ENTREZID), resSigEntrez$SYMBOL, resSigEntrez$ENTREZID)
colnames(resSigEntrez) = c("GeneID", "ENTREZID")
resSig_with_entrez = merge(as.data.frame(resSig), resSigEntrez, by = "GeneID") #All entrezIDs have been retrieved for the 36 genes with p<0.1 and lfc=0 threshold
# Extract entrezIDs for the res variable
res_allgenes<-res
res_allgenes$GeneID<-rownames(res)
# removing "LOC" info
res_allgenes$GeneID <- stringr::str_remove(res_allgenes$GeneID, "LOC")
rownames(res_allgenes) = stringr::str_remove(rownames(res_allgenes), "LOC")

res_allgenesEntrez <- AnnotationDbi::select(Oaries, keys =  res_allgenes$GeneID,
                                            columns = c('ENTREZID','GENENAME'), keytype = 'SYMBOL')
# Now replace the NA values in entrezid column with the values from first column.
res_allgenesEntrez$ENTREZID <- ifelse(is.na(res_allgenesEntrez$ENTREZID), res_allgenesEntrez$SYMBOL, res_allgenesEntrez$ENTREZID)
colnames(res_allgenesEntrez) = c("GeneID", "ENTREZID")
res_allgenes_with_entrez = merge(as.data.frame(res_allgenes), res_allgenesEntrez, by = "GeneID") #All entrezIDs have been retrieved for the DEGs
dim(res_allgenes_with_entrez)
#we want the log2 fold change 
original_gene_list <- res_allgenes_with_entrez$log2FoldChange
# name the vector
names(original_gene_list) <- res_allgenes_with_entrez$ENTREZID
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
# Exctract significant results (padj < 0.05)
sig_genes_df = subset(res_allgenes_with_entrez, padj < 0.05)
# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange
# Name the vector
names(genes) <- sig_genes_df$ENTREZID
# omit NA values
genes <- na.omit(genes)
# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]
length(genes)
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = Oaries, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

as.data.frame(go_enrich)
Ora_go_df = file.path(out_dir, "DMI_ORA_GO_enrichments.csv")
write.csv(as.data.frame(go_enrich),Ora_go_df, row.names=T)

# Exctract significant results from df2
kegg_sig_genes_df = subset(resSig_with_entrez, padj < 0.05)
# From significant results, we want to filter on log2fold change
kegg_genes <- kegg_sig_genes_df$log2FoldChange
# Name the vector with the CONVERTED ID!
names(kegg_genes) <- kegg_sig_genes_df$ENTREZID
# omit NA values
kegg_genes <- na.omit(kegg_genes)
# filter on log2fold change (PARAMETER)
kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]
Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE)
kegg_organism = "oas"
kk <- enrichKEGG(gene=kegg_genes, universe=names(gene_list),organism=kegg_organism, pvalueCutoff = 0.05)
as.data.frame(kk)
Ora_kegg_df = file.path(out_dir, "DMI_ORA_KEGG_enrichments.csv")
write.csv(as.data.frame(kk),Ora_kegg_df, row.names=T)
