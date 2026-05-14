# 14_NicheNet_Analysis.R
# This script performs NicheNet analysis to predict ligand-receptor interactions
# and their downstream target genes, focusing on GC/DN as receivers in the
# neighborhood subset exported by `11_GC_DN_Neighborhood_Communication.py`.

library(nichenetr)
library(Seurat)
library(tidyverse)
library(Matrix)
library(data.table)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else normalizePath("New_Analysis", winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "analysis_paths.R"))

print("Step 14: Running NicheNet Analysis on GC/DN Neighborhoods")

# 1. Define paths
input_dir <- pipeline_files$communication_input_dir
out_dir <- pipeline_files$nichenet_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load the exported data (Matrix Market format from Step 05)
print("Loading data...")
counts_mtx_file <- file.path(input_dir, "counts.mtx")
genes_file <- file.path(input_dir, "genes.tsv")
barcodes_file <- file.path(input_dir, "barcodes.tsv")
meta_file <- file.path(input_dir, "meta.txt")

meta <- read.table(meta_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(meta) <- meta$Cell

# ==============================================================================
# 全局随机种子设置 (Global Seed Setting)
# 确保分析结果的完全可重复性
# ==============================================================================
set.seed(42)
# ==============================================================================

if (!file.exists(counts_mtx_file) || !file.exists(genes_file) || !file.exists(barcodes_file)) {
    stop("counts.mtx / genes.tsv / barcodes.tsv was not found. Please rerun Step 11.")
}

data.input <- readMM(counts_mtx_file)
data.input <- as(data.input, "CsparseMatrix")
genes <- read.table(genes_file, header=FALSE, stringsAsFactors=FALSE)$V1
barcodes <- read.table(barcodes_file, header=FALSE, stringsAsFactors=FALSE)$V1
rownames(data.input) <- make.unique(as.character(genes))
colnames(data.input) <- barcodes

meta <- meta[colnames(data.input), , drop=FALSE]

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = data.input, meta.data = meta)
Idents(seurat_obj) <- "cell_type"
seurat_obj <- NormalizeData(seurat_obj)

# 3. Load NicheNet Prior Models
# Note: NicheNet requires downloading background models (~100MB total). 
print("Loading NicheNet prior models from local directory...")
models_dir <- "NicheNet_Models"
if (!file.exists(file.path(models_dir, "ligand_target_matrix.rds")) ||
    !file.exists(file.path(models_dir, "lr_network.rds")) ||
    !file.exists(file.path(models_dir, "weighted_networks.rds"))) {
    stop("NicheNet prior models not found in NicheNet_Models/. Please prepare the models locally before running.")
}

ligand_target_matrix <- readRDS(file.path(models_dir, "ligand_target_matrix.rds"))
lr_network <- readRDS(file.path(models_dir, "lr_network.rds"))
weighted_networks <- readRDS(file.path(models_dir, "weighted_networks.rds"))

# 4. Define NicheNet Analysis Parameters
# We will treat GC and DN as "Receivers" and all other cells as "Senders"
cell_types <- levels(Idents(seurat_obj))
receivers <- grep("Giant Cells|Dysmorphic Neurons", cell_types, value=TRUE)
senders <- setdiff(cell_types, receivers)

if(length(receivers) == 0) {
    stop("No Giant Cells or Dysmorphic Neurons found in the dataset.")
}

print(paste("Receivers:", paste(receivers, collapse=", ")))
print(paste("Senders:", paste(senders, collapse=", ")))

# 5. Define genes of interest
# Find marker genes of the receivers (compared to senders) to see what ligands are driving their specific expression
print("Finding marker genes for receivers to define genes of interest...")
receiver_markers <- FindMarkers(seurat_obj, ident.1 = receivers, ident.2 = senders, only.pos = TRUE, logfc.threshold = 0.25)
genes_of_interest <- rownames(receiver_markers)[receiver_markers$p_val_adj < 0.05]

if(length(genes_of_interest) == 0) {
    print("Warning: No significant DEGs found. Using top 100 genes by fold change.")
    genes_of_interest <- rownames(head(receiver_markers, 100))
}

# Define background genes (all expressed genes in receivers)
expressed_genes_receiver <- get_expressed_genes(receivers, seurat_obj, pct = 0.10)
background_genes <- expressed_genes_receiver

# 6. Define potential ligands from senders
expressed_genes_sender <- get_expressed_genes(senders, seurat_obj, pct = 0.10)
expressed_ligands <- intersect(lr_network$from, expressed_genes_sender)
expressed_receptors <- intersect(lr_network$to, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# 7. Perform NicheNet ligand activity prediction
print("Predicting ligand activity...")
ligand_activities <- predict_ligand_activities(genes = genes_of_interest, background_express = background_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
top_ligands <- ligand_activities %>% top_n(20, pearson) %>% pull(test_ligand) %>% unique()

# Save ligand activity results
write.csv(ligand_activities, file.path(out_dir, "Ligand_Activities.csv"), row.names=FALSE)

# 8. Visualization
print("Generating NicheNet visualizations...")
graphics.off()

# Create a dedicated directory for plots to keep things clean
plots_dir <- file.path(out_dir, "Plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# (A) Ligand activity bar plot
p_ligand_activity <- ligand_activities %>% top_n(20, pearson) %>% 
    mutate(test_ligand = fct_reorder(test_ligand, pearson)) %>% 
    ggplot(aes(x = pearson, y = test_ligand)) + 
    geom_bar(stat = "identity", fill="steelblue") + 
    theme_minimal() + 
    labs(x = "Pearson correlation (Ligand Activity)", y = "Top Ligands")

pdf(file.path(plots_dir, "1_Top_Ligands_Activity.pdf"), width = 6, height = 6)
print(p_ligand_activity)
dev.off()

png(file.path(plots_dir, "1_Top_Ligands_Activity.png"), width = 1800, height = 1800, res=300)
print(p_ligand_activity)
dev.off()

# NEW: Determine which ligands are highly expressed in GC vs DN vs Others
print("Determining ligand expression specificity (GC vs DN)...")
avg_exp <- AverageExpression(seurat_obj, features = top_ligands, group.by = "cell_type")$RNA
gc_col <- grep("Giant Cells", colnames(avg_exp), value=TRUE)
dn_cols <- grep("Dysmorphic Neurons", colnames(avg_exp), value=TRUE)

source_df <- data.frame(Ligand = top_ligands, Source = rep("Unknown", length(top_ligands)))

if(length(gc_col) > 0 && length(dn_cols) > 0) {
    # Simple logic: for each ligand, find which cell type has max expression among GC, DN(avg), and others
    dn_avg <- rowMeans(avg_exp[, dn_cols, drop=FALSE])
    gc_exp <- avg_exp[, gc_col]
    
    ligand_source <- sapply(top_ligands, function(lig) {
        gc_val <- gc_exp[lig]
        dn_val <- dn_avg[lig]
        if(gc_val > dn_val * 1.5) return("GC Enriched")
        if(dn_val > gc_val * 1.5) return("DN Enriched")
        return("Shared / Other")
    })
    
    source_df$Source <- ligand_source
    write.csv(source_df, file.path(out_dir, "Top_Ligands_Source_GC_vs_DN.csv"), row.names=FALSE)
    print("Saved Ligand Source classification to Top_Ligands_Source_GC_vs_DN.csv")
}

# (B) Ligand-Target matrix heatmap
# Ensure we only use genes that actually exist in the ligand_target_matrix
valid_genes_of_interest <- intersect(genes_of_interest, rownames(ligand_target_matrix))
valid_top_ligands <- intersect(top_ligands, colnames(ligand_target_matrix))

if (length(valid_genes_of_interest) > 0 && length(valid_top_ligands) > 0) {
    active_ligand_target_links_df <- ligand_target_matrix[valid_genes_of_interest, valid_top_ligands, drop=FALSE]
    active_ligand_target_links_df[active_ligand_target_links_df < 0.05] <- 0
    # Keep only rows (genes) and columns (ligands) that have at least one non-zero interaction
    active_ligand_target_links_df <- active_ligand_target_links_df[rowSums(active_ligand_target_links_df) > 0, , drop=FALSE]
    active_ligand_target_links_df <- active_ligand_target_links_df[, colSums(active_ligand_target_links_df) > 0, drop=FALSE]
    
    if(nrow(active_ligand_target_links_df) > 0 && ncol(active_ligand_target_links_df) > 0) {
        library(pheatmap)
        
        pdf(file.path(plots_dir, "2_Ligand_Target_Heatmap.pdf"), width = max(6, ncol(active_ligand_target_links_df)*0.5), height = max(4, nrow(active_ligand_target_links_df)*0.2))
        pheatmap(t(active_ligand_target_links_df), main="Ligand-Target Regulatory Potential")
        dev.off()
        
        png(file.path(plots_dir, "2_Ligand_Target_Heatmap.png"), width = max(1800, ncol(active_ligand_target_links_df)*150), height = max(1200, nrow(active_ligand_target_links_df)*60), res=300)
        pheatmap(t(active_ligand_target_links_df), main="Ligand-Target Regulatory Potential")
        dev.off()
    } else {
        print("Warning: No strong regulatory links (>0.05) found to plot in heatmap.")
    }
} else {
    print("Warning: No valid intersecting genes/ligands found for heatmap.")
}

# NEW: Downstream target gene visualization per ligand
print("Visualizing downstream target genes for each ligand...")

# Create a subdirectory for downstream targets
targets_dir <- file.path(plots_dir, "Downstream_Targets")
dir.create(targets_dir, recursive = TRUE, showWarnings = FALSE)

# Iterate over each ligand in top_ligands
for(lig in top_ligands) {
    if(lig %in% colnames(ligand_target_matrix)) {
        # Get target scores for this ligand
        ligand_targets <- ligand_target_matrix[, lig]
        
        # Filter for genes of interest (DEGs in receivers)
        target_genes <- intersect(names(ligand_targets), valid_genes_of_interest)
        
        if(length(target_genes) > 0) {
            ligand_targets <- ligand_targets[target_genes]
            
            # Sort scores from high to low
            ligand_targets <- sort(ligand_targets, decreasing = TRUE)
            
            # Keep top 20 targets for visualization
            top_targets <- head(ligand_targets, 20)
            
            # Only plot if max score is > 0
            if(max(top_targets) > 0) {
                target_df <- data.frame(
                    Target = factor(names(top_targets), levels = rev(names(top_targets))),
                    Score = top_targets
                )
                
                # Get the source annotation for the plot title
                lig_source <- source_df$Source[source_df$Ligand == lig]
                if(length(lig_source) == 0) lig_source <- "Unknown"
                
                p_target <- ggplot(target_df, aes(x = Score, y = Target)) +
                    geom_bar(stat = "identity", fill = ifelse(grepl("GC", lig_source), "darkred", 
                                                       ifelse(grepl("DN", lig_source), "darkblue", "grey"))) +
                    theme_minimal() +
                    labs(title = paste0("Targets of ", lig, " (", lig_source, ")"),
                         x = "Regulatory Potential Score", y = "Target Gene")
                
                # Sanitize filename
                safe_lig <- gsub("[^A-Za-z0-9_]", "_", lig)
                
                pdf(file.path(targets_dir, paste0("Targets_", safe_lig, ".pdf")), width = 5, height = 5)
                print(p_target)
                dev.off()
                
                png(file.path(targets_dir, paste0("Targets_", safe_lig, ".png")), width = 1500, height = 1500, res=300)
                print(p_target)
                dev.off()
            }
        }
    }
}
print("Downstream target visualizations completed.")

# (C) Ligand expression dotplot in sender cells
p_dotplot <- DotPlot(seurat_obj, features = top_ligands %>% rev(), group.by = "cell_type") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(plots_dir, "3_Ligand_Expression_DotPlot.pdf"), width = 8, height = 6)
print(p_dotplot)
dev.off()

png(file.path(plots_dir, "3_Ligand_Expression_DotPlot.png"), width = 2400, height = 1800, res=300)
print(p_dotplot)
dev.off()

# (D) Receptors and Intracellular Signaling network
# Predict receptors for the top ligands
print("Predicting Receptors and mapping intracellular networks...")
lr_network_top <- lr_network %>% filter(from %in% top_ligands & to %in% expressed_receptors) %>% filter(from %in% valid_top_ligands)
top_receptors <- lr_network_top %>% pull(to) %>% unique()

if(length(top_receptors) > 0) {
    # Dotplot for Receptor expression in Receiver cells
    p_rec_dotplot <- DotPlot(seurat_obj, features = top_receptors %>% rev(), group.by = "cell_type") + 
        coord_flip() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    pdf(file.path(plots_dir, "4_Receptor_Expression_DotPlot.pdf"), width = 8, height = max(4, length(top_receptors)*0.3))
    print(p_rec_dotplot)
    dev.off()
    
    png(file.path(plots_dir, "4_Receptor_Expression_DotPlot.png"), width = 2400, height = max(1200, length(top_receptors)*90), res=300)
    print(p_rec_dotplot)
    dev.off()
    
    # Extract Ligand-Receptor-Target paths
    # For a high level visualization, we can use the weighted networks
    active_receptors <- lr_network_top$to
    
    # Save the L-R network links
    write.csv(lr_network_top, file.path(out_dir, "Ligand_Receptor_Links.csv"), row.names=FALSE)
    
    # Optional: Ligand-Receptor heatmap
    lr_matrix <- matrix(0, nrow=length(top_ligands), ncol=length(top_receptors))
    rownames(lr_matrix) <- top_ligands
    colnames(lr_matrix) <- top_receptors
    for(i in 1:nrow(lr_network_top)) {
        lr_matrix[lr_network_top$from[i], lr_network_top$to[i]] <- 1
    }
    
    # Filter empty rows/cols
    lr_matrix <- lr_matrix[rowSums(lr_matrix) > 0, colSums(lr_matrix) > 0, drop=FALSE]
    if(nrow(lr_matrix) > 0 && ncol(lr_matrix) > 0) {
        pdf(file.path(plots_dir, "5_Ligand_Receptor_Network.pdf"), width = max(6, ncol(lr_matrix)*0.4), height = max(4, nrow(lr_matrix)*0.4))
        pheatmap(lr_matrix, cluster_rows=TRUE, cluster_cols=TRUE, main="Ligand-Receptor Network (Binary)", color = colorRampPalette(c("white", "darkred"))(2), legend=FALSE)
        dev.off()
        
        png(file.path(plots_dir, "5_Ligand_Receptor_Network.png"), width = max(1800, ncol(lr_matrix)*120), height = max(1200, nrow(lr_matrix)*120), res=300)
        pheatmap(lr_matrix, cluster_rows=TRUE, cluster_cols=TRUE, main="Ligand-Receptor Network (Binary)", color = colorRampPalette(c("white", "darkred"))(2), legend=FALSE)
        dev.off()
    }
}

print("NicheNet analysis complete! Results and plots saved.")


