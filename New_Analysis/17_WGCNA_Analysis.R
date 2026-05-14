# ==============================================================================
# 文件名称: 17_WGCNA_Analysis.R
# 功能描述: 加权基因共表达网络分析 (WGCNA)
# 联动说明: 作为 GRN 的补充，在 R 中运行 WGCNA 提取高度协同表达的基因模块，识别与癫痫/病理特征高度相关的 Hub 基因。
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath("New_Analysis", winslash = "/", mustWork = FALSE)
}
source(file.path(script_dir, "analysis_paths.R"))

# ================= Configuration =================
# Set paths based on analysis_config.py
BASE_DIR <- base_dir

# Configure reticulate to use the current project .venv Python interpreter
venv_python_path <- file.path(BASE_DIR, ".venv", "Scripts", "python.exe")
Sys.setenv(RETICULATE_PYTHON = venv_python_path)
Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")

suppressPackageStartupMessages(library(reticulate))
use_python(venv_python_path, required = TRUE)

INPUT_FILE <- pipeline_files$annotated_h5ad
OUTPUT_DIR <- pipeline_files$wgcna_dir

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# WGCNA settings
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Ensure PDFs generate editable text (Type 42 font equivalent in R)
pdf.options(useDingbats = FALSE)

# ================= 1. Load Data =================
cat("Loading h5ad data...\n")
# Instead of using the `anndata` R package which tries to find conda,
# we use reticulate to directly import the python `anndata` module.
ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad(INPUT_FILE)

# Extract expression matrix (genes x cells/spots)
# Using normalized data in .X or .raw
if (!is.null(adata$raw)) {
  expr_mat <- py_to_r(adata$raw$X)
  gene_names <- py_to_r(adata$raw$var_names$to_list())
} else {
  expr_mat <- py_to_r(adata$X)
  gene_names <- py_to_r(adata$var_names$to_list())
}

# Convert sparse matrix to dense if necessary, or sample to avoid memory issues
# For WGCNA, we transpose to (samples x genes)
# To save memory and time, we filter for top highly variable genes or top expressed genes
cat("Filtering genes for WGCNA...\n")
# For spatial data, we might have 50k spots. WGCNA on 50k spots is very slow.
# We will aggregate by spatial clusters or sample a subset of spots.
# Let's pseudo-bulk by cell type and sample, or randomly sample 5000 spots.
metadata <- as.data.frame(adata$obs)
celltype_col <- if ("celltype" %in% colnames(metadata)) "celltype" else "cell_type"
cell_types <- as.character(metadata[[celltype_col]])

# We will use top 3000 highly variable genes
if (!is.null(adata$var$highly_variable)) {
  hvg_idx <- which(adata$var$highly_variable == TRUE)
  if (length(hvg_idx) > 3000) {
    hvg_idx <- hvg_idx[1:3000]
  }
} else {
  # If no HVG, just take top 3000 by variance
  gene_vars <- apply(expr_mat, 2, var)
  hvg_idx <- order(gene_vars, decreasing = TRUE)[1:3000]
}

datExpr <- as.matrix(expr_mat[, hvg_idx])
colnames(datExpr) <- gene_names[hvg_idx]

# If too many spots (>5000), sample to save memory
if (nrow(datExpr) > 5000) {
  set.seed(42)
  sample_idx <- sample(1:nrow(datExpr), 5000)
  datExpr <- datExpr[sample_idx, ]
  metadata <- metadata[sample_idx, ]
}

# Check for missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# ================= 2. Soft Thresholding =================
cat("Calculating soft threshold power...\n")
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot soft threshold
pdf(file.path(OUTPUT_DIR, "1_Soft_Thresholding.pdf"), width = 9, height = 5)
par(mfrow = c(1,2))
cex1 <- 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6

# ================= 3. Network Construction & Module Detection =================
cat("Constructing co-expression network (Power =", softPower, ")...\n")
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# Plot Dendrogram
pdf(file.path(OUTPUT_DIR, "2_Module_Dendrogram.pdf"), width = 12, height = 9)
plotDendroAndColors(geneTree, moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# ================= 4. Relate Modules to Traits (Cell Types) =================
cat("Relating modules to cell types...\n")
# Create dummy variables for cell types
trait_formula <- as.formula(paste0("~ 0 + ", celltype_col))
traits <- model.matrix(trait_formula, data = metadata)
colnames(traits) <- str_replace(colnames(traits), paste0("^", celltype_col), "")

# Calculate correlation
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Plot Heatmap
pdf(file.path(OUTPUT_DIR, "3_Module_Trait_Relationships.pdf"), width = 10, height = 8)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(8, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# ================= 5. Identify Hub Genes =================
cat("Identifying Hub Genes...\n")
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

# Export Top Hub Genes for each module
hub_genes <- data.frame()
for (mod in unique(moduleColors)) {
  if (mod == "grey") next
  
  modGenes <- (moduleColors == mod)
  mod_genes_names <- colnames(datExpr)[modGenes]
  
  # Get kME for this module
  kME_col <- paste0("ME", mod)
  if (kME_col %in% colnames(geneModuleMembership)) {
    mod_kME <- geneModuleMembership[modGenes, kME_col]
    
    # Sort by kME
    sorted_idx <- order(mod_kME, decreasing = TRUE)
    top_genes <- mod_genes_names[sorted_idx][1:min(20, length(sorted_idx))]
    
    hub_genes <- rbind(hub_genes, data.frame(
      Module = mod,
      Gene = top_genes,
      kME = mod_kME[sorted_idx][1:min(20, length(sorted_idx))]
    ))
  }
}

write.csv(hub_genes, file.path(OUTPUT_DIR, "Hub_Genes_Per_Module.csv"), row.names = FALSE)
cat("WGCNA analysis completed. Results saved to:", OUTPUT_DIR, "\n")


