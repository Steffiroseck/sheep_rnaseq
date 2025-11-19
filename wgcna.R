library(WGCNA)
library(flashClust)
library(curl)
library(DESeq2)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyverse)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads()          # allow multi-threading (optional)

data = counts #(from deseq2)

# detect outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
  
table(gsg$goodGenes) # 25030 genes passed
table(gsg$goodSamples)
  
# if allOK returen false, remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]
  
# detect outlier samples - hierarchical clustering - method 1
sampleTree <- hclust(dist(t(data)), method = "average") #Clustering samples based on distance 
hclust = file.path(out_dir, "1.hclust-samples.pdf")
par(cex = 0.6);
par(mar = c(0,4,2,0))
# Plotting the cluster dendrogram
plot1 = plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
ggsave(filename = hclust, plot = plot1)

# detect outlier samples - pca - method 2
pca <- prcomp(t(data))
pca.dat <- pca$x
  
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
  
pca.dat <- as.data.frame(pca.dat)
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
        y = paste0('PC2: ', pca.var.percent[2], ' %'))

# Normalization
# create a deseq2 dataset
# making the rownames and column names identical
# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
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
  norm.counts <- assay(dds_norm) %>% 
  t()

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts, powerVector = powers, verbose = 5)

# Set some parameters
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (powers)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

softPower <- 6
# calling adjacency function
adjacency <- adjacency(norm.counts, power = softPower, type="signed")
# TOM

TOM <- TOMsimilarity(adjacency, TOMType = "signed")#This gives similarity between genes
TOM.dissimilarity <- 1-TOM # get dissimilarity matrix

# Hierarchical Clustering Analysis
#The dissimilarity/distance measures are then clustered using linkage hierarchical clustering and a dendrogram (cluster tree) of genes is constructed.

# creating the dendrogram 
hclustGeneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
plot(hclustGeneTree, xlab = "", sub = "", 
     main = "Gene Clustering on TOM-based disssimilarity", 
     labels = FALSE, hang = 0.04)
print(plot4)

# Make the modules larger, so set the minimum higher

minModuleSize <- 30

# Module ID using dynamic tree cut

dynamicMods <- cutreeDynamic(dendro = hclustGeneTree, 
                             distM = TOM.dissimilarity,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
							 
table(dynamicMods)#returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module.

# Convert numeric labels into colors

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)#returns the counts for each color (aka the number of genes within each module) 

# Plot the dendrogram and colors underneath
plotDendroAndColors(hclustGeneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

# Calculate the module eigengenes

dynamic_MEList <- moduleEigengenes(norm.counts, colors = dynamicColors)
dynamic_MEs <- dynamic_MEList$eigengenes
head(dynamic_MEs)

#To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations 
#and merge the modules that have similar expression profiles.
# Calculate dissimilarity of module eigengenes

dynamic_MEDiss <- 1-cor(dynamic_MEs) #Calculate eigengene dissimilarity
dynamic_METree <- hclust(as.dist(dynamic_MEDiss), method = "average")#Clustering eigengenes 

# Plot the hclust
plot(dynamic_METree, main = "Dynamic Clustering of module eigengenes",
     xlab = "", sub = "",)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75

dynamic_MEDissThres <- 0.25

# Call an automatic merging function
merge_dynamic_MEDs <- mergeCloseModules(norm.counts, dynamicColors, cutHeight = dynamic_MEDissThres, verbose = 3)

# The Merged Colors
dynamic_mergedColors <- merge_dynamic_MEDs$colors

# Eigen genes of the new merged modules
mergedMEs <- merge_dynamic_MEDs$newMEs

table(dynamic_mergedColors)

# dendrogram with original and merged modules
plotDendroAndColors(hclustGeneTree, cbind(dynamicColors, dynamic_mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors for original and merged modules")

# Rename Module Colors 
moduleColors <- dynamic_mergedColors

# Construct numerical labels corresponding to the colors 
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

# Trait matching
# pull out all continuous traits. Before that convertiing ID column which is a factor to numeric type
sample_metadata = metaData
sample_metadata$ID = as.numeric(levels(sample_metadata$ID))[sample_metadata$ID]
str(sample_metadata)
allTraits <- sample_metadata[2:3]
rownames(allTraits)
#some of the columns have 0 values. better to remove them prior to plotting heatmap

# sample names should be consistent in eigen genes and traits !!!!
#allTraits = allTraits[match(rownames(MEs), rownames(allTraits)), ]
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
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Display the correlation values within a heatmap
labeledHeatmap(Matrix = moduleTraitCor, 
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

heatmap.data <- merge(MEs , allTraits, by = 'row.names')
head(heatmap.data)
heatmap.data <- heatmap.data %>% 
column_to_rownames(var = 'Row.names')
  
 
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[16:45],
             y = names(heatmap.data)[1:15],
             col = c("blue1", "skyblue", "white", "pink", "red"),
			       signifSymbols = c("***", "**", "*", ""),
             signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
             rotLabX = 30, rotLabY = 30)


# Target gene identification
# Define variable methane production 
metpro <- as.data.frame(allTraits$DMI)
names(metpro) <- "DMI" # rename

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


modNames <- names(geneModuleMembership)


# Initialize for loop
# NEED: modNames
for (i in names(geneModuleMembership)) {
  
  # Pull out the module we're working on
  module <- i
  print(module)   
  
  # Find the index in the column of the dataframe 
  column <- match(module, modNames)
  #print(column)
  
  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors == module
  genenames = rownames(geneTraitSignificance)
  #print(paste("There are ", length(genenames[moduleGenes]), " genes in the ", module, " module.", sep = ""))
  #print(genenames[moduleGenes])
  
  # Make the plot
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                 abs(geneTraitSignificance[moduleGenes, 1]),
                 xlab = paste("Module Membership in", module, "module"),
                 ylab = "Gene significance for Total VFA",
                 main = paste("Module membership vs. gene significnace \n"),
                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}    


# Print the number of genes in each module
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
GSpval <- GSPvalue %>%
  tibble::rownames_to_column(var = "gene")

# Prepare module membership df
gMM_df <- geneModuleMembership %>%
  tibble::rownames_to_column(var = "gene") %>%
  gather(key = "moduleColor", value = "moduleMemberCorr", -gene) 

# Prepare gene significance df
GS_metprod_df <- geneTraitSignificance %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "gene")

# Put everything together 
allData_df <- gMM_df %>%
  left_join(GS_metprod_df, by = "gene") %>%
  left_join(GSpval, by = "gene") 
