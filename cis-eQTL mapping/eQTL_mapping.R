# ---------------------------
# 1. Load data
# ---------------------------
counts_raw <- read.table("28_lambs_adg_dmi_twas.csv", header = TRUE, sep=",", check.names = FALSE, row.names = 1)
colnames(counts_raw) <- paste0("UK17906910", colnames(counts_raw))
head(counts_raw)



sample_info <- fread("28lambs_metadata.csv", sep=",", header=TRUE)
# Expect columns: SampleID, batch, sex, age, RIN, etc.

# Reorder counts to match sample_info
stopifnot(all(sample_info$Genotype %in% colnames(counts_raw)))
counts_raw <- counts_raw[, sample_info$Genotype]

# ---------------------------
# 2. Filter lowly expressed genes
# ---------------------------
# Keep genes with CPM > 1 in at least 50% samples (>=14 for n=28)
# dge <- DGEList(counts = counts)
# dge <- calcNormFactors(dge, method = "TMM")
# 
# keep <- rowSums(cpm(dge) > 1) >= ceiling(ncol(counts) * 0.5)
# dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
# dge_filt <- calcNormFactors(dge_filt, method = "TMM")

# Filter genes
keep_genes <- rowSums( counts_raw > 5 ) >= 7
counts <- counts_raw[keep_genes,]
dim(counts)

# ---------------------------
# 3. Normalize and transform expression
# ---------------------------
library(DESeq2)

# Assuming 'countData' is your filtered raw count matrix (genes × samples)
# and 'colData' is a data frame with sample metadata (e.g., Treatment)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ 1)  # No design needed for normalization only


dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

vst_dds <- vst(dds, blind = TRUE)         # or blind = FALSE if you want to preserve design signal
vst_mat <- assay(vst_dds)

# Use voom (log2 CPM with precision weights) or simple log2 CPM.
# v <- voom(dge_filt, design = NULL, plot = FALSE)
# expr_logcpm <- v$E  # genes x samples (log2-CPM)

# compute PCs and keep top PCs as covariates
pca <- prcomp(t(norm_counts), center = TRUE, scale. = FALSE)
K <- 10
peer_like_factors <- pca$x[, 1:K]    # samples x K matrix
# inspect:
cor(peer_like_factors[,1:10], sample_info$Ave_MO, use = "pairwise.complete.obs")

# ---------------------------
# 4. Covariates from genotypes
# ---------------------------

# Load eigenmatrix from PCA
cov_geno <- read.table("genotype_PCA/Twas_afterQC_PCA.eigenvec")
colnames(cov_geno) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
cov_geno$FID<- NULL
rownames(cov_geno) <- cov_geno$IID
cov_geno$IID <- NULL

# use only first three PCs
covrt <- cov_geno[,c("PC1", "PC2", "PC3")]
rownames(covrt) <- rownames(cov_geno)

# Covariates addition did not improved results in this study, hence was not used in subsequent codes.

# ---------------------------
# 5. Load genotype matrix and positions
# ---------------------------
# geno_matrix.txt: rows=SNPs, cols=samples, values 0/1/2; first column "SNP"
geno <- fread("Twas_afterQC_dosage.raw")
geno_matrix <- as.data.frame(geno[,-c(1:6)])
rownames(geno_matrix) <- geno$IID
geno_matrix <- t(geno_matrix) #snps as rows and samples as columns
rownames(geno_matrix) <- gsub("_[A-Z]$", "", rownames(geno_matrix))

snp_ids <- rownames(geno_matrix)
geno_mat <- as.matrix(geno_matrix[, sample_info$Genotype])  # SNPs x samples
rownames(geno_mat) <- snp_ids

# Load SNP and gene positions
snp_pos <- fread("SNP_annotation.csv", sep=",", header=T)   # columns: SNP, chr, pos
gene_pos <- fread("Gene_annotation.csv", sep=",", header=T) # columns: gene, chr, start, end, TSS

common_snps <- intersect(rownames(geno_mat), snp_pos$SNPID)
geno_mat <- geno_mat[common_snps, ]
snp_pos  <- snp_pos[snp_pos$SNPID %in% common_snps, ]

stopifnot(all(rownames(norm_counts) %in% gene_pos$GeneID))
stopifnot(all(rownames(geno_mat) %in% snp_pos$SNPID))

# Align orders
gene_pos <- gene_pos[match(rownames(norm_counts), gene_pos$GeneID)]
snp_pos  <- snp_pos[match(rownames(geno_mat), snp_pos$SNPID)]

str(snp_pos)
str(gene_pos)

gene_pos$CHR <- as.integer(gene_pos$CHR)

common_genes <- intersect(rownames(vst_mat), gene_pos$GeneID)
vst_mat <- vst_mat[common_genes, , drop = FALSE]
gene_pos_sub <- gene_pos[match(common_genes, gene_pos$GeneID), ]


#snp_gene_pos <- merge(as.data.frame(gene_pos), as.data.frame(snp_pos), by="CHR")
#write.csv(snp_gene_pos, "SNP_GENE_ANNOTATION.csv")

# ---------------------------
# 6. Prepare SlicedData objects for MatrixEQTL
# ---------------------------
# Expression
geneSD <- SlicedData$new()
geneSD$CreateFromMatrix(as.matrix(vst_mat))    # genes x samples

# Genotypes
snpSD <- SlicedData$new()
snpSD$CreateFromMatrix(geno_mat)                  # SNPs x samples

# Covariates for MatrixEQTL (if using non-residualized expression)
# Here we residualized already; we pass an empty covariate matrix.
cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(t(peer_like_factors))   # transpose: PCs as rows

# ---------------------------
# 7. Run SNP-level cis-eQTL mapping
# ---------------------------
useModel <- modelLINEAR
cisDist <-  1e6  # 1 Mb (1e6) window; adjust as needed

# Position data frames must be: snps: SNP, chr, pos; genes: gene, chr, TSS
snpspos  <- data.frame(snp = snp_pos$SNPID, chr = snp_pos$CHR, pos = snp_pos$POSITION)
genepos  <- data.frame(gene = gene_pos$GeneID, chr = gene_pos$CHR, s1 = gene_pos$START, s2 = gene_pos$END)

pvThreshold_cis <- 1  # liberal output threshold to collect candidates
pvThreshold_trans <- 0   # no trans output in this run. If 1 then # keep all cis results, filter later

cisOutFile <- "cis_eqtls.raw.csv" 

cis_res <- Matrix_eQTL_main(
  snps = snpSD,
  gene = geneSD,
  #cvrt = cvrt,
  output_file_name     = NULL,
  pvOutputThreshold    = pvThreshold_trans,
  useModel             = useModel,
  errorCovariance      = numeric(),
  verbose              = TRUE,
  output_file_name.cis = cisOutFile,
  pvOutputThreshold.cis= pvThreshold_cis,
  snpspos              = snpspos,
  genepos              = genepos,
  cisDist              = cisDist,
  pvalue.hist          = "qqplot",
  min.pv.by.genesnp    = TRUE,
  noFDRsaveMemory      = FALSE
)

plot(cis_res)
cis_tab <- cis_res$cis$eqtls #ALL cis tests + FDR
cat(sprintf("cis-eQTL tests performed: %d; reported (p<=%g): %d\n",
            cis_res$cis$number_of_tests, pvThreshold_cis, nrow(cis_tab)))

cat(sprintf("Genes retained: %d of %d\n", nrow(counts), nrow(counts_raw)))

# ---------------------------
# 8. Multiple testing correction (cis)
# ---------------------------
if (nrow(cis_tab) > 0) {
  significant_eqtls <- cis_res$cis$eqtls[cis_res$cis$eqtls$FDR < 0.05, ]
  cis_tab$FDR <- as.numeric(as.character(cis_tab$FDR))
  fwrite(significant_eqtls, "Original_cis_eqtls.fdr0.05.tsv", sep = "\t")
  cat(sprintf("cis-eQTLs at FDR<=0.05: %d\n", nrow(significant_eqtls)))
} else {
  cat("No cis-eQTLs passed the reporting threshold.\n")
}


# ---------------------------
# 9. Manhattan Ploting
# ---------------------------

library(qqman)
# Rename p-value column if needed
colnames(cis_eqtl)[colnames(cis_eqtl) == "pvalue"] <- "P"
colnames(cis_eqtl)[colnames(cis_eqtl) == "snps"] <- "SNPID"

# get CHR info from snppos
cis_eqtl_chr <- merge(cis_eqtl, snp_pos, by="SNPID")

# Convert chromosome to numeric if needed
cis_eqtl_chr$CHR <- as.numeric(as.character(cis_eqtl_chr$CHR))

# Label genes
cis_eqtl_chr$logP <- -log10(cis_eqtl_chr$P)
genome_threshold <- 5e-8
cis_eqtl_chr$significant <- cis_eqtl_chr$P < genome_threshold


# Manhattan plot
pdf("Original_cis-eQTL_ManhattanPlot.pdf", width=6, height=5)
manhattan(
  cis_eqtl_chr,
  chr = "CHR",
  bp = "POSITION",
  p = "P",
  snp = "gene",
  col = c("black", "darkgrey"),
  #genomewideline = -log10(5e-8),
  #suggestiveline = -log10(1e-5),
  genomewideline = FALSE,
  suggestiveline = FALSE,
  #cex = 0.9,
  ylim=c(0,11),
  cex.axis = 0.9,
  main = "Manhattan Plot of cis-eQTLs",
  annotateTop = 8, 
  annotatePval =5e-7
)
dev.off()

# ---------------------------
# 10. Gene-level analysis (FastQTL style)
# ---------------------------
set.seed(123)
n_perm <- 1000   # increase to 1000+ for publication
 
perm_results <- vector("list", n_perm)
 
for (i in 1:n_perm) {
 
  cat("Permutation", i, "\n")
  
  # Permute expression across samples 
  perm_idx <- sample(ncol(vst_mat)) # random permutation of sample indices
 
  # Permute expression across samples
  perm_expr <- SlicedData$new()
  perm_expr$CreateFromMatrix(as.matrix(vst_mat[, perm_idx]))
 
  me_perm <- Matrix_eQTL_main(
    snps = snpSD,
    gene = perm_expr,
    cvrt = cvrt,
    pvOutputThreshold.cis = 1,
    useModel = modelLINEAR,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = 1e6,
    verbose = FALSE
  )
 
  perm_results[[i]] <- me_perm$cis$eqtls
}


# original cis eqtl
nrow(cis_res$cis$eqtls)
# perm based
nrow(perm_results[[1]])
summary(cis_res$cis$eqtls$pvalue)
summary(perm_results[[1]]$pvalue)
summary(perm_results[[2]]$pvalue)

#cis_res            # observed cis-eQTLs
#perm_results        # list of permutation cis-eQTL tables


#Gene level min p
eqtl_obs <- as.data.table(cis_res$cis$eqtls)
 
obs_minp <- eqtl_obs[
  , .(obs_min_p = min(pvalue, na.rm = TRUE)),
  by = gene
]

#gene level min p from permutations
#length(perm_results) = 1000
 
genes <- obs_minp$gene
n_perm <- length(perm_results)
 
perm_minp <- matrix(NA, nrow = length(genes), ncol = n_perm)
rownames(perm_minp) <- genes


for (i in seq_len(n_perm)) {

  perm_dt <- as.data.table(perm_results[[i]])

  perm_gene_minp <- perm_dt[
    , .(perm_min_p = min(pvalue, na.rm = TRUE)),
    by = gene
  ]

  # match gene names to perm_minp rows
  idx <- match(perm_gene_minp$gene, rownames(perm_minp))

  # only update rows that exist
  valid <- !is.na(idx)

  perm_minp[idx[valid], i] <- perm_gene_minp$perm_min_p[valid]
}


#empirical p value for fastqtl style gene level
empirical_p <- sapply(genes, function(g) {
 
obs_p <- obs_minp[gene == g, obs_min_p]
perm_p <- perm_minp[g, ]
 
if (is.na(obs_p)) return(NA)
 
  (sum(perm_p <= obs_p, na.rm = TRUE) + 1) /
  (sum(!is.na(perm_p)) + 1)
})

gene_eqtl <- data.table(
  gene = genes,
  obs_min_p = obs_minp$obs_min_p,
  empirical_p = empirical_p
)

#gene level BH FDR
gene_eqtl$gene_FDR = p.adjust(gene_eqtl$empirical_p, method = "BH")
#gene-level BH on obs p
gene_eqtl$gene_obs_p_FDR = p.adjust(gene_eqtl$obs_min_p, method = "BH")
write.csv(gene_eqtl, "Fastqtl_style_all_genes_and_snps.csv")
#write.csv(empirical_genes_FDR05, "Candidate_genes_with_emp_p_fdr_lessthan_0.05_FastqlStyle.csv")

sig_genes_FDR05 <- gene_eqtl[gene_FDR < 0.05]
sig_genes_FDR10 <- gene_eqtl[gene_FDR < 0.20]

# get snp info for these eGenes
sig_genes_FDR10_snp <- merge(sig_genes_FDR10,cis_eqtl_chr, by="gene")
colnames(sig_genes_FDR10_snp)[colnames(sig_genes_FDR10_snp) == "SNPID"] <- "snp"
write.csv(sig_genes_FDR10, "Significant_Genes_FDR_0.20_Fastql_style_permutation.csv")
 
nrow(sig_genes_FDR05)
nrow(sig_genes_FDR10)
 
summary(gene_eqtl$obs_min_p)
summary(gene_eqtl$empirical_p)
summary(gene_eqtl$gene_FDR)

hist(perm_results[[1]]$pvalue, breaks = 50) # should be flat and mean ~ 0.5. If skewed left, permutation is wrong
hist(gene_eqtl$empirical_p, breaks = 50)

# ---------------------------
# 11. GTEx style permutation (SNP-level)
# ---------------------------
#collect real p from initial run
real_p <- cis_res$cis$eqtls$pvalue
# permuted p 
null_p <- unlist(lapply(perm_results, function(dt) dt$pvalue))

length(real_p)
length(null_p)

summary(real_p)
summary(null_p)

hist(real_p, 100, main = "Real p-values")
hist(null_p, 100, main = "Permuted p-values")


#sort them
real_sorted <- sort(real_p)
null_sorted <- sort(null_p)

#compute empirical FDR curve
# cumulative fraction of real p-values
real_cum <- seq_along(real_sorted) / length(real_sorted)
# cumulative fraction of null p-values
null_counts <- findInterval(real_sorted, null_sorted)
null_cum <- null_counts / length(null_sorted)
# empirical FDR
emp_fdr_curve <- null_cum / real_cum

# map each real p-value to its empirical FDR
emp_fdr <- emp_fdr_curve[findInterval(real_p, real_sorted)] #SNP‑level empirical FDR for each cis‑eQTL.

# attaching it to real dataset cis_tab
cis_tab$emp_fdr_gtex <- emp_fdr
sig_eqtls_gtex <- cis_tab[cis_tab$emp_fdr_gtex < 0.05, ] #100
write.csv(sig_eqtls_gtex, "Significant_SNP-Gene_pairs_FDR_0.05_GTEx_style_permutation.csv")

                        
