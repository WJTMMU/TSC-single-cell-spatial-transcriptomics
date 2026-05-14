# 14_CellChat_Analysis.R
# This script reads the CellPhoneDB formatted data exported from Python
# and runs a comprehensive CellChat analysis to uncover communication mechanisms.

# Set working directory to the project root (adjust if needed)
# setwd("d:/Data_Analysis/空间转录组分析/ST_analysis")

library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(Matrix)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else normalizePath("New_Analysis", winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "analysis_paths.R"))

print("Step 12: Running CellChat Analysis on GC/DN Neighborhoods")

# 1. Define paths
input_dir <- pipeline_files$communication_input_dir
out_dir <- pipeline_files$cellchat_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load the exported data
print("Loading data...")
# Read counts and meta
counts_mtx_file <- file.path(input_dir, "counts.mtx")
counts_txt_file <- file.path(input_dir, "counts.txt")
genes_file <- file.path(input_dir, "genes.tsv")
barcodes_file <- file.path(input_dir, "barcodes.tsv")
meta_file <- file.path(input_dir, "meta.txt")

# Read metadata
meta <- read.table(meta_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Force row names and column names to be unique and perfectly matched
new_cell_names <- meta$Cell
rownames(meta) <- new_cell_names

if (file.exists(counts_mtx_file)) {
    # Read counts from mtx
    print("Reading mtx file...")
    data.input <- readMM(counts_mtx_file)
    data.input <- as(data.input, "CsparseMatrix")
    
    # Read genes and barcodes
    genes <- read.table(genes_file, header=FALSE, stringsAsFactors=FALSE)$V1
    barcodes <- read.table(barcodes_file, header=FALSE, stringsAsFactors=FALSE)$V1
    
    rownames(data.input) <- genes
    colnames(data.input) <- barcodes
} else if (file.exists(counts_txt_file)) {
    print("mtx file not found. Falling back to reading counts.txt using data.table::fread...")
    if (!requireNamespace("data.table", quietly = TRUE)) {
        install.packages("data.table", repos="https://cloud.r-project.org")
    }
    library(data.table)
    
    counts_df <- fread(counts_txt_file, data.table=FALSE)
    
    genes <- counts_df[, 1]
    counts_df <- counts_df[, -1]
    
    # Fix column mismatch if any
    if (ncol(counts_df) > nrow(meta)) {
        counts_df <- counts_df[, 1:nrow(meta)]
    }
    
    colnames(counts_df) <- new_cell_names
    rownames(counts_df) <- make.unique(as.character(genes))
    
    print("Converting to sparse matrix...")
    data.input <- as(as.matrix(counts_df), "sparseMatrix")
} else {
    stop("Neither counts.mtx nor counts.txt was found!")
}

# Ensure metadata rownames match data.input colnames
meta <- meta[colnames(data.input), , drop=FALSE]

# Format cell types
meta$labels <- as.factor(meta$cell_type)

print(paste("Loaded data with", nrow(data.input), "genes and", ncol(data.input), "cells."))

# 3. Initialize CellChat object
print("Initializing CellChat object...")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# Set the database to Human
CellChatDB <- CellChatDB.human 
# Use all CellChatDB (including Secreted Signaling, ECM-Receptor, Cell-Cell Contact)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use

# 4. Preprocessing and Computing
print("Preprocessing and computing communications...")
cellchat <- subsetData(cellchat) # Subset the expression data of signaling genes
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project gene expression data onto PPI (Protein-Protein Interaction) network to mitigate dropout effects
print("Projecting data onto PPI network...")
cellchat <- projectData(cellchat, PPI.human)

# Compute communication probability
# Default is triMean. If data is NOT log-normalized, CellChat recommends normalizeData() first, 
# but if the Python script exported raw counts, CellChat needs normalization.
# Let's explicitly normalize to be safe.
print("Normalizing data...")
# Use Seurat's NormalizeData as a more robust alternative if CellChat's fails on sparse matrices
if (!requireNamespace("Seurat", quietly = TRUE)) {
    install.packages("Seurat")
}
library(Seurat)

# ==============================================================================
# 全局随机种子设置 (Global Seed Setting)
# 确保分析结果的完全可重复性
# ==============================================================================
set.seed(42)
# ==============================================================================

cellchat@data.signaling <- NormalizeData(cellchat@data.signaling)

cellchat <- computeCommunProb(cellchat, type = "triMean")

# Filter out communications with few cells in a group
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer cell-cell communication at signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Extract and save the inferred communication network dataframe
print("Saving subset communication network dataframe...")
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file.path(out_dir, "cellchat_inferred_network.csv"), row.names = FALSE)

# Save the object
saveRDS(cellchat, file = file.path(out_dir, "cellchat_neighborhood.rds"))

# ==============================================================================
# 5. Visualization (Nature-level Plots)
# ==============================================================================
print("Generating visualizations...")

# Clear any leftover graphics devices to avoid corrupted PDFs
graphics.off()

# Get group sizes for the vertex weights in the circle plot
groupSize <- as.numeric(table(cellchat@idents))

# Color palette
cell_types <- levels(cellchat@idents)
num_clusters <- length(cell_types)
my_palette <- scales::hue_pal()(num_clusters)
names(my_palette) <- cell_types

# Create a dedicated directory for plots to keep things clean
plots_dir <- file.path(out_dir, "Plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# (A) Circle Plot (Aggregated network)
pdf(file.path(plots_dir, "1_Aggregated_Network_Circle.pdf"), width = 8, height = 8)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
dev.off()

png(file.path(plots_dir, "1_Aggregated_Network_Circle.png"), width = 2400, height = 2400, res=300)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
dev.off()

# (A-2) Individual Circle Plots for each cell type (Count and Weight)
mat_count <- cellchat@net$count
pdf(file.path(plots_dir, "1_2_NetVisual_Circle_Individual_Count.pdf"), width = 12, height = 10)
par(mfrow = c(ceiling(nrow(mat_count)/4), 4), xpd=TRUE)
for (i in 1:nrow(mat_count)) {
    mat2 <- matrix(0, nrow = nrow(mat_count), ncol = ncol(mat_count), dimnames = dimnames(mat_count))
    mat2[i, ] <- mat_count[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(mat_count), title.name = rownames(mat_count)[i])
}
dev.off()

mat_weight <- cellchat@net$weight
pdf(file.path(plots_dir, "1_3_NetVisual_Circle_Individual_Weight.pdf"), width = 12, height = 10)
par(mfrow = c(ceiling(nrow(mat_weight)/4), 4), xpd=TRUE)
for (i in 1:nrow(mat_weight)) {
    mat2 <- matrix(0, nrow = nrow(mat_weight), ncol = ncol(mat_weight), dimnames = dimnames(mat_weight))
    mat2[i, ] <- mat_weight[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(mat_weight), title.name = rownames(mat_weight)[i])
}
dev.off()

# (A-3) Interaction Heatmaps
pdf(file.path(plots_dir, "1_4_NetVisual_Heatmap.pdf"), width = 12, height = 8)
h1 <- netVisual_heatmap(cellchat, measure = "count", title.name = "Number of Interactions (Heatmap)")
h2 <- netVisual_heatmap(cellchat, measure = "weight", title.name = "Interaction Weight (Heatmap)")
print(h1 + h2)
dev.off()

# (B) Identify major signaling roles (Sender/Receiver)
# Compute network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf(file.path(plots_dir, "2_Signaling_Roles_Scatter.pdf"), width = 8, height = 6)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

png(file.path(plots_dir, "2_Signaling_Roles_Scatter.png"), width = 2400, height = 1800, res=300)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

# (C) DotPlot (Bubble Plot) for Interactions
# First, a comprehensive Bubble Plot for all cell types
pdf(file.path(plots_dir, "3_Bubble_Plot_All.pdf"), width = 14, height = 10)
print(netVisual_bubble(cellchat, remove.isolate = FALSE) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

# Then focus on interactions where GC or DN are senders or receivers
targets <- grep("Giant Cells|Dysmorphic Neurons", cell_types, value=TRUE)
others <- setdiff(cell_types, targets)

# Clean up cell type names to make valid filenames
make_safe_name <- function(name) {
    gsub("[^A-Za-z0-9]", "_", name)
}

if(length(targets) > 0) {
    
    # ---------------------------------------------------------
    # Target as Sender (Target -> Others) - Separated by cell type
    # ---------------------------------------------------------
    for (target in targets) {
        safe_target <- make_safe_name(target)
        
        pdf(file.path(plots_dir, paste0("3_2_DotPlot_", safe_target, "_as_Sender.pdf")), width = 10, height = 8)
        # remove.isolate=TRUE removes columns with no communication
        tryCatch({
            print(netVisual_bubble(cellchat, sources.use = target, targets.use = others, remove.isolate = TRUE) + 
                  ggtitle(paste(target, "as Sender")))
        }, error = function(e) {
            plot.new(); text(0.5, 0.5, paste("No significant interactions found for", target, "as Sender"))
        })
        dev.off()
        
        png(file.path(plots_dir, paste0("3_2_DotPlot_", safe_target, "_as_Sender.png")), width = 3000, height = 2400, res=300)
        tryCatch({
            print(netVisual_bubble(cellchat, sources.use = target, targets.use = others, remove.isolate = TRUE) + 
                  ggtitle(paste(target, "as Sender")))
        }, error = function(e) {
            plot.new(); text(0.5, 0.5, paste("No significant interactions found for", target, "as Sender"))
        })
        dev.off()
    }
    
    # ---------------------------------------------------------
    # Target as Receiver (Others -> Target) - Separated by cell type
    # ---------------------------------------------------------
    for (target in targets) {
        safe_target <- make_safe_name(target)
        
        pdf(file.path(plots_dir, paste0("4_DotPlot_", safe_target, "_as_Receiver.pdf")), width = 10, height = 8)
        tryCatch({
            print(netVisual_bubble(cellchat, sources.use = others, targets.use = target, remove.isolate = TRUE) + 
                  ggtitle(paste(target, "as Receiver")))
        }, error = function(e) {
            plot.new(); text(0.5, 0.5, paste("No significant interactions found for", target, "as Receiver"))
        })
        dev.off()
        
        png(file.path(plots_dir, paste0("4_DotPlot_", safe_target, "_as_Receiver.png")), width = 3000, height = 2400, res=300)
        tryCatch({
            print(netVisual_bubble(cellchat, sources.use = others, targets.use = target, remove.isolate = TRUE) + 
                  ggtitle(paste(target, "as Receiver")))
        }, error = function(e) {
            plot.new(); text(0.5, 0.5, paste("No significant interactions found for", target, "as Receiver"))
        })
        dev.off()
    }
}

# (D) Detailed Pathway Visualizations
# Find the top 10 most significant pathways
pathways_dir <- file.path(plots_dir, "Top10_Pathways")
dir.create(pathways_dir, recursive = TRUE, showWarnings = FALSE)

pathways.show.all <- cellchat@netP$pathways
if(length(pathways.show.all) > 0) {
    # Expand to top 10 pathways
    top_pathways <- head(pathways.show.all, 10)
    
    for (pw in top_pathways) {
        # Circle plot for specific pathway
        pdf(file.path(pathways_dir, paste0("Pathway_", pw, "_Circle.pdf")), width = 8, height = 8)
        netVisual_aggregate(cellchat, signaling = pw, layout = "circle")
        dev.off()
        
        png(file.path(pathways_dir, paste0("Pathway_", pw, "_Circle.png")), width = 2400, height = 2400, res=300)
        netVisual_aggregate(cellchat, signaling = pw, layout = "circle")
        dev.off()
        
        # Chord diagram for specific pathway
        pdf(file.path(pathways_dir, paste0("Pathway_", pw, "_Chord.pdf")), width = 8, height = 8)
        tryCatch({
            netVisual_aggregate(cellchat, signaling = pw, layout = "chord")
            title(paste0(pw, " (Chord)"))
        }, error = function(e) {
            plot.new(); text(0.5, 0.5, paste("Chord plot error:", e$message))
        })
        dev.off()
        
        # Contribution of L-R pairs to the pathway
        pdf(file.path(pathways_dir, paste0("Pathway_", pw, "_LR_Contribution.pdf")), width = 6, height = 4)
        print(netAnalysis_contribution(cellchat, signaling = pw))
        dev.off()
        
        png(file.path(pathways_dir, paste0("Pathway_", pw, "_LR_Contribution.png")), width = 1800, height = 1200, res=300)
        print(netAnalysis_contribution(cellchat, signaling = pw))
        dev.off()
    }
}

print("CellChat analysis complete! Visualizations saved in CellChat_Neighborhood directory.")
