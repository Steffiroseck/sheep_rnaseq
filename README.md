# 🧬 Multi-Omics Analysis in R: DESeq2, ClusterProfiler, and WGCNA

This repository contains bash and R scripts and workflows for transcriptomic analysis using **DESeq2**, **functional enrichment via ClusterProfiler**, and **network construction with WGCNA**. The analyses were performed on RNA-seq data from lambs samples, with a focus on differential expression, biological interpretation, and co-expression module detection.

---


## 🧪 Analyses Overview

### 1. Sequence pre-processing, alignment, counts generation
- Initial QC
- Alignment to reference genome
- Generation of read counts

### 2. Differential Expression (DESeq2)
- Normalization and dispersion estimation
- Wald test for differential expression
- PCA plot, volcano plot
- Output: `deseq2_results`, `plotPCA.png`

### 3. Functional Enrichment (ClusterProfiler)
- Overrepresentation analysis (ORA) for GO and KEGG
- Visualization: dotplot, barplot, enrichment map
- Output: `go_enrichment.csv`, `kegg_enrichment.csv`, `dotplot.png`

### 4. Co-expression Network (WGCNA)
- Sample clustering and outlier detection
- Soft-threshold power selection
- Module detection and trait correlation
- Output: `module_colors.csv`, `eigengene_heatmap.png`, `trait_correlations.csv`

---

## 🚀 Getting Started

### Prerequisites for R
Install required packages:
```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "WGCNA", "org.Hs.eg.db"))

# DESeq2 and ClusterProfiler
source("Deseq2/2.DGE_and_ORA.R")

# WGCNA
source("WGCNA/wgcna.R")


---

## 🔄 Workflow Overview

```mermaid
flowchart TD
    A[Raw FASTQ files] --> B[FastQC: Quality check]
    B --> C[Trimmomatic: Adapter trimming & filtering]
    C --> D[Hisat2: Alignment to reference genome]
    D --> E[FeatureCounts: Gene-level quantification]
    E --> F[DESeq2: Differential expression analysis]
    F --> G[ClusterProfiler: Functional enrichment (GO/KEGG)]
    G --> H[WGCNA: Co-expression network & trait correlation]
