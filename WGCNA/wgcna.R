library(WGCNA)
library(flashClust)
library(curl)
library(DESeq2)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(CorLevelPlot)
library(clusterProfiler)
library(pathview)
library(stringr)
library(AnnotationHub)
library(ggridges)
library(enrichplot)



setwd("C:/Users/afbi-roses/Steffi_sheep_transcriptomics/WGCNA_DMI_ADG/")
out_dir="C:/Users/afbi-roses/Steffi_sheep_transcriptomics/WGCNA_DMI_ADG/"

countData=read.csv("33lambs_featureCounts.csv",header=T,row.names=1,sep=",", check.names = FALSE)
head(countData)
# Calculate total number of columns
# total_columns <- ncol(countData)
# zero_counts <- rowSums(countData == 0)
# length(zero_counts)
# table(zero_counts)
# data <- as.data.frame(countData) %>% filter(zero_counts < 0.7 * total_columns)
# dim(data)
keep_genes <- rowSums( countData > 5 ) >= 7
countData2 <- countData[keep_genes,]
dim(countData)
dim(countData2)

# Read the metadata
sample_metadata = read.csv(file = "33lambs_phenotypes.csv", sep=",", row.names = 1)
head(sample_metadata)


#####################################################################################
# Normalization
###################################################################################
# create a deseq2 dataset
data = countData2
data <- data[,unique(rownames(sample_metadata))]
all(colnames(data) == rownames(sample_metadata))
# create dds
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = sample_metadata,
                              design = ~ 1) # not spcifying model
dim(dds)
# perform variance stabilization
dds_norm <- vst(dds)
# get normalized counts
norm.counts <- assay(dds_norm) %>% t()

###########################################################################################
# QC - outlier detection
###########################################################################################

# detect outlier genes

gsg <- goodSamplesGenes(t(norm.counts))
summary(gsg)
gsg$allOK
table(gsg$goodGenes) # 25030 genes passed
table(gsg$goodSamples)
# if allOK returen false, remove genes that are detectd as outliers
data <- norm.counts[gsg$goodGenes == TRUE,]
# detect outlier samples - hierarchical clustering - method 1

sampleTree <- hclust(dist(data), method = "average") #Clustering samples based on distance 
sampleConn <- adjacency(t(data), type ="unsigned")
Z.k <- scale(apply(sampleConn, 2, sum))
Z.k

# Outlier sample removal
outlier_samples <- c("7094", "6977", "6764", "6976", "6856", "6968", "7153", "7196") 
datExpr_clean <- data[!rownames(data) %in% outlier_samples, ]

# rerun sample clustering to see any change
sampleTree <- hclust(dist(datExpr_clean), method = "average") #Clustering samples based on distance 
hclust = pdf("1.hclust-samples_on_normalized_counts_7094_6967_6764_removed.pdf")
par(cex = 0.6);
par(mar = c(0,4,2,0))
# Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
             cex.axis = 1.5, cex.main = 2)
dev.off()



##################################################################
# PCA
#################################################################

data <- datExpr_clean
# detect outlier samples - pca - method 2
pca <- prcomp(data)
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
pcaplot = file.path(out_dir, "2.PCA_norm_counts.pdf")
plot2=ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
ggsave(filename = pcaplot, plot = plot2) 



###################################################################################################
# Choosing the soft-thresholding power: analysis of network topology
###################################################################################################
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts, powerVector = powers, verbose = 5)
# Plot the results:
pdf(paste0(out_dir,"3.soft-power.pdf"),height=7,width=7.5)
par(mfrow = c(1,2))
# Set some parameters
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot3 = plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
             xlab="Soft Threshold (powers)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
             main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot3 = plot(sft$fitIndices[,1], sft$fitIndices[,5], 
             xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
             main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
cex1 = 0.9
print(plot3)
dev.off()

sft

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)
softPower <- sft$powerEstimate #7
# calling adjacency function
adjacency <- adjacency(norm.counts, power = softPower, type="signed")
# TOM
TOM <- TOMsimilarity(adjacency, TOMType = "signed")#This gives similarity between genes
TOM.dissimilarity <- 1-TOM # get dissimilarity matrix

# Hierarchical Clustering Analysis
#The dissimilarity/distance measures are then clustered using linkage hierarchical clustering and a dendrogram (cluster tree) of genes is constructed.

# creating the dendrogram 
pdf(paste0(out_dir, "4.dendrogram.gene.clustering.TOM.dissimilarity.pdf"))
hclustGeneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

# Plot the resulting clustering tree (dendogram)

plot4 = plot(hclustGeneTree, xlab = "", sub = "", 
             main = "Gene Clustering on TOM-based disssimilarity", 
             labels = FALSE, hang = 0.04)
print(plot4)
dev.off()

##############################################################################
# identify modules
##############################################################################
# Make the modules larger, so set the minimum higher
minModuleSize <- 30
# Module ID using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = hclustGeneTree, cutHeight=0.975,method="hybrid",
                             distM = TOM.dissimilarity,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)#returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module.
# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)#returns the counts for each color (aka the number of genes within each module) 
pdf(paste0(out_dir, "5.dendrogram.with.module.colors.pdf"))
# Plot the dendrogram and colors underneath
plot5 = plotDendroAndColors(hclustGeneTree, dynamicColors, "Dynamic Tree Cut", 
                            dendroLabels = FALSE, hang = 0.03, 
                            addGuide = TRUE, guideHang = 0.05, 
                            main = "Gene dendrogram and module colors")
print(plot5)
dev.off()

##############################################################
# Module Eigengenes
#############################################################
# Calculate the module eigengenes
dynamic_MEList <- moduleEigengenes(norm.counts, colors = dynamicColors)
dynamic_MEs <- dynamic_MEList$eigengenes
head(dynamic_MEs)
#To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations 
#and merge the modules that have similar expression profiles.
# Calculate dissimilarity of module eigengenes
dynamic_MEDiss <- 1-cor(dynamic_MEs) #Calculate eigengene dissimilarity
dynamic_METree <- hclust(as.dist(dynamic_MEDiss), method = "average")#Clustering eigengenes 
pdf(paste0(out_dir, "6.clustering of module eigen genes.pdf"))
# Plot the hclust
plot6 = plot(dynamic_METree, main = "Dynamic Clustering of module eigengenes",
             xlab = "", sub = "",)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
print(plot6)
dev.off()

########################################################################
# Merge similar modules
########################################################################
dynamic_MEDissThres <- 0.25
# Call an automatic merging function
merge_dynamic_MEDs <- mergeCloseModules(norm.counts, dynamicColors, cutHeight = dynamic_MEDissThres, verbose = 3)
# The Merged Colors
dynamic_mergedColors <- merge_dynamic_MEDs$colors
# Eigen genes of the new merged modules
mergedMEs <- merge_dynamic_MEDs$newMEs
table(dynamic_mergedColors)
# dendrogram with original and merged modules
pdf(paste0(out_dir, "7.Dendrogram-with-original_and-merged-modules.pdf"))
plot7 = plotDendroAndColors(hclustGeneTree, cbind(dynamicColors, dynamic_mergedColors),
                            c("Dynamic Tree Cut", "Merged dynamic"),
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors for original and merged modules")
print(plot7)
dev.off()
# Rename Module Colors 
moduleColors <- dynamic_mergedColors

# Construct numerical labels corresponding to the colors 
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs
#save(MEs, moduleLabels, moduleColors, hclustGeneTree, file = "wgcna-networkConstruction.RData")


#####################################################################################
# Relating modules to external traits
#####################################################################################
# pull out all continuous traits. Before that convertiing ID column which is a factor to numeric type
sample_metadata_original = sample_metadata
#str(sample_metadata)
allTraits <- sample_metadata[1:3]
table(rownames(MEs) == rownames(allTraits))

# define numbers of genes and samples
nGenes <- ncol(norm.counts)
nSamples <- nrow(norm.counts)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(norm.counts, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
names(MEs) <- substring(names(MEs), 3)
moduleTraitCor <- cor(MEs, allTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# create module-trait heatmap
# Will display correlations and their p-values
pdf(paste0(out_dir, "8.Module-trait-heatmap.pdf"), width = 17, height=10)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap
plot8 = labeledHeatmap(Matrix = moduleTraitCor, 
                       xLabels = names(allTraits),
                       yLabels = names(MEs), 
                       ySymbols = names(MEs), 
                       colorLabels = FALSE, 
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix, 
                       setStdMargins = FALSE,
                       cex.text = 0.8,
                       zlim = c(-1,1),
                       main = paste("Module-trait Relationships"))
print(plot8)	
dev.off()
#Each row corresponds to a module eigengene, and the columns correspond to a trait. 
#Each cell contains a p-value and correlation. Those with strong positive correlations are shaded a darker red while those with stronger negative correlations become more blue. 
# heatmap with significance as stars (*)
heatmap.data <- merge(MEs , allTraits, by = 'row.names')
head(heatmap.data)
heatmap.data <- heatmap.data %>%   column_to_rownames(var = 'Row.names')
pdf(paste0(out_dir, "9.module trait relationships with significance.pdf"), width = 17, height = 10)
plot9 = CorLevelPlot(heatmap.data,
                     x = names(heatmap.data)[19:21],
                     y = names(heatmap.data)[1:18],
                     col = c("blue1", "skyblue", "white", "pink", "red"),
                     signifSymbols = c("***", "**", "*", ""),
                     signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     rotLabX = 30, rotLabY = 30)
print(plot9)
dev.off()
table(dynamic_mergedColors)

###################################
#LIST OF GEnes in each modules
#####################################
# LIST OF GENES IN EACH MODUALS
topGOgenes <- names(norm.counts)[moduleColors =="midnightblue"]
table(topGOgenes) 
write.csv(topGOgenes,"midnightblue.txt")

##############################################
#13. correlation values within a heatmap plot #
##############################################
dissTOM = 1-TOMsimilarityFromExpr(norm.counts, power = 6);
plotTOM = dissTOM^7;
diag(plotTOM) = NA;
# Call the plot function
memory.limit(40000)
sizeGrWindow(9,9)
pdf(paste0(out_dir, "11.Network heatmap plot.pdf"), width = 17, height = 10)
TOMplot(plotTOM, hclustGeneTree, dynamicColors, main = "Network heatmap plot")
dev.off()

#####################################################################################################
# Target gene identification
# Gene relationship to trait and important modules
#####################################################################################################
# To Quantify associations of individual gene with trait of interest (methane production)
# For this first define gene significance (GS) as the absolute value of the correlation between the gene and the trait.
# For each module, also define a quantitative measure of module membership (MM) as the correlation of the module eigngene and the gene expression profile. Allows us to quantify the similarity of all genes on the array to every module.
# Define variable methane production 
metpro <- as.data.frame(allTraits$ADG)
names(metpro) <- "ADG" # rename
# Calculate the correlations between modules
geneModuleMembership <- as.data.frame(WGCNA::cor(norm.counts, MEs, use = "p"))
# What are the p-values for each correlation?
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# What's the correlation for the trait: methane production?
geneTraitSignificance <- as.data.frame(cor(norm.counts, metpro, use = "p"))
# What are the p-values for each correlation?
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nSamples))
names(geneTraitSignificance) <- paste("GS.", names(metpro), sep = "")
names(GSPvalue) <- paste("p.GS.", names(metpro), sep = "")

####################################################
# IIntramodular analysis: identify genes with high GS and MM
####################################################
modNames <- names(geneModuleMembership)
par(mfrow = c(2,3)) 
# Initialize for loop
# NEED: modNames
pdf(paste0(out_dir, "10.Gene_signicance_Module_membership.pdf"))
for (i in names(geneModuleMembership)) {
  # Pull out the module we're working on
  module <- i
  print(module) 
  # Find the index in the column of the dataframe 
  column <- match(module, modNames) #print(column)
  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors == module
  genenames = rownames(geneTraitSignificance)
  #print(paste("There are ", length(genenames[moduleGenes]), " genes in the ", module, " module.", sep = ""))
  #print(genenames[moduleGenes])
  # Make the plot
  plot10 = verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                              abs(geneTraitSignificance[moduleGenes, 1]),
                              xlab = paste("Module Membership in", module, "module"),
                              ylab = "Gene significance for Total VFA",
                              main = paste("Module membership vs. gene significnace \n"),
                              cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}    
print(plot10)
dev.off()

# # Print the number of genes in each module
for (i in names(geneModuleMembership)) {
  # Pull out the module we're working on
  module <- i
  # Find the index in the column of the dataframe
  column <- match(module, modNames)
  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors == module
  genenames = rownames(geneTraitSignificance)
  print(paste("There are ", length(genenames[moduleGenes]), " genes in the ", module, " module.", sep = ""))

  # NOTE: This makes hidden variables with the gene names
  assign(paste(module, "_genes", sep = ""), genenames[moduleGenes])
}

# merge this important module information with gene annotation information and write out a file that summarizes the results.
# Combine pval, module membership, and gene significance into one dataframe
# Prepare pvalue df
# GSpval <- GSPvalue %>%
#   tibble::rownames_to_column(var = "gene")
# # Prepare module membership df
# gMM_df <- geneModuleMembership %>%
#   tibble::rownames_to_column(var = "gene") %>%
#   gather(key = "moduleColor", value = "moduleMemberCorr", -gene) 
# # Prepare gene significance df
# GS_metprod_df <- geneTraitSignificance %>%
#   data.frame() %>%
#   tibble::rownames_to_column(var = "gene")
# # Put everything together 
# allData_df <- gMM_df %>%
#   left_join(GS_metprod_df, by = "gene") %>%
#   left_join(GSpval, by = "gene") 

# Write a file 
#write.csv(allData_df, file = "allData_df_WGCNA.csv")


#################################################
# Hub gene in each module
#################################################
hub = chooseTopHubInEachModule(norm.counts, moduleColors)
hubgene = file.path(out_dir,"hub_genes_in_each_module.csv" )
write.csv(hub, file = hubgene, row.names = T) 

#########################################################
# Writing all genes into a file
########################################################

#Identifying most important genes for one determined characteristic inside of the cluster
probes = colnames(norm.counts)
geneInfo = data.frame(Genes = probes,
                      moduleColor = moduleColors,
                      geneTraitSignificance,
                      GSPvalue)
#Order modules by their significance for trait
modOrder = order(-abs(cor(MEs, metpro, use = "p")))

for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo)
  geneInfo0 = data.frame(geneInfo, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.ADG))
  geneInfo = geneInfo0[geneOrder, ]
  
}
geneInfowrite = file.path(out_dir, "geneInfo_ADG.csv")
write.csv(geneInfo, file = geneInfowrite, row.names=T)
save.image(file = "WGCNA_Norm_counts_workspace.RData")
