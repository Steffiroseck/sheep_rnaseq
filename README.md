# 🧬 RNA Sequencing Analysis

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
- Output: `deseq2_results`, `plotPCA.png`, `Volcano.png`

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

### Prerequisites for bash
install via conda por apt.
- fatsqc
- trimmomatic
- hisat2
- subread(featureCounts)

Make sure you have installed the necessary softwares required to run the pipeline. These are mentioned in detail below.

**1. FastQC:** To generate quality reports for the raw fastq files coming from the sequencing platforms.
  
  ```
  conda install bioconda::fastqc
  ```
   or
  
   ```
   sudo apt-install fastqc
  ```
   or install from source
   ```
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
  unzip fastqc_v0.12.1.zip
  cd FastQC/
  chmod 755 fastqc
  ./fastqc --help

  ```
   
**2. Trimmomatic:** 
  ```
  conda install bioconda::trimmomatic
  ```
  or install from source. Refer to the link below.
  ```
  http://www.usadellab.org/cms/?page=trimmomatic
  ```

**3. Hisat2:** 
  ```
  sudo apt-get update
  sudo apt-get -y install hisat2
  ```
  or you can simply download from the source
  ```
  curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
  unzip hisat2-2.2.1-Linux_x86_64.zip
  cd hisat2-2.2.1
  ./hisat2 -h
  ```
  **5. featureCounts:**
  The installation instructions can be found at : https://subread.sourceforge.net/featureCounts.html


  ```
  
  bash NGS_data_analysis/1.ngs_data_analysis.sh
  
  ```

### Prerequisites for R
Install required packages:
```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "WGCNA"))

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

```

