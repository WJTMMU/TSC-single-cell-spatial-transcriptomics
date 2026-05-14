# 21_Targeted_NicheNet_Downstream.R
# This script performs a TARGETED NicheNet analysis to predict downstream 
# target genes of specific key ligands acting on GC/DN.

# Install nichenetr if not present
if (!requireNamespace("nichenetr", quietly = TRUE)) {
    devtools::install_github("saeyslab/nichenetr")
}

library(nichenetr)
library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(data.table)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else normalizePath("New_Analysis", winslash = "/", mustWork = FALSE)
source(file.path(script_dir, "analysis_paths.R"))

print("Step 21: Running Targeted NicheNet Analysis on Key Molecules")

# ==============================================================================
# 全局绘图风格设置 (Nature NPG 色系, 高DPI, 无插值)
# 与 analysis_config.py 中的 NPG_COLORS 保持严格一致
# ==============================================================================
NPG_COLORS <- c(
    '#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', 
    '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85',
    '#2E2A2B', '#E5A5C2', '#586E6E', '#C0C0C0', '#D4E2E0',
    '#8A2BE2', '#D2691E', '#20B2AA', '#FF1493', '#32CD32'
)
theme_set(theme_classic() + theme(
    text = element_text(family = "sans", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
))

# 全局随机种子设置 (Global Seed Setting)
set.seed(42)

# ==============================================================================
# 1. 定义路径与输入数据
# ==============================================================================
input_dir <- pipeline_files$communication_input_dir
out_dir <- pipeline_files$targeted_nichenet_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 核心靶点中的配体 (从之前筛选的 28 个双重验证核心基因中挑出的外分泌/膜表面蛋白)
core_ligands <- c("EFNA5", "LAMB1", "TNR", "CDH9", "CDH12", "CDH18", "CDH22", "ROBO1", "ROBO2")

# 在 Figure 6 的逻辑中，DN 是 Sender，我们想看这些配体对“微环境细胞 (Receivers)”的影响！
# 所以接收细胞列表应该改为微环境中的胶质细胞或免疫细胞，而不是 DN 本身。
receiver_cells_list <- list(
    "Astrocytes" = "Astrocytes",
    "Microglia" = "Microglia",
    "Endothelial Cells" = "Endothelial Cells"
)


# ==============================================================================
# 2. 加载表达数据与 Metadata
# ==============================================================================
print("Loading expression data and metadata...")
counts_mtx_file <- file.path(input_dir, "counts.mtx")
genes_file <- file.path(input_dir, "genes.tsv")
barcodes_file <- file.path(input_dir, "barcodes.tsv")
meta_file <- file.path(input_dir, "meta.txt")

if (!file.exists(counts_mtx_file) || !file.exists(genes_file) || !file.exists(barcodes_file) || !file.exists(meta_file)) {
    stop("Input expression data not found. Please run earlier communication steps first.")
}

meta <- read.table(meta_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(meta) <- meta$Cell

data.input <- readMM(counts_mtx_file)
data.input <- as(data.input, "CsparseMatrix")
genes <- read.table(genes_file, header=FALSE, stringsAsFactors=FALSE)$V1
barcodes <- read.table(barcodes_file, header=FALSE, stringsAsFactors=FALSE)$V1
rownames(data.input) <- make.unique(as.character(genes))
colnames(data.input) <- barcodes
meta <- meta[colnames(data.input), , drop=FALSE]
rownames(counts_df) <- make.unique(as.character(genes))

data.input <- as(as.matrix(counts_df), "sparseMatrix")
meta <- meta[colnames(data.input), , drop=FALSE]

seurat_obj <- CreateSeuratObject(counts = data.input, meta.data = meta)
Idents(seurat_obj) <- "cell_type"
seurat_obj <- NormalizeData(seurat_obj)

# ==============================================================================
# 3. 加载 NicheNet 先验模型 (Prior Models)
# ==============================================================================
print("Loading NicheNet prior models from local directory...")
models_dir <- "NicheNet_Models"
if (!file.exists(file.path(models_dir, "ligand_target_matrix.rds"))) {
    stop("NicheNet models not found. Please ensure they are downloaded in NicheNet_Models/.")
}
ligand_target_matrix <- readRDS(file.path(models_dir, "ligand_target_matrix.rds"))
lr_network <- readRDS(file.path(models_dir, "lr_network.rds"))
weighted_networks <- readRDS(file.path(models_dir, "weighted_networks.rds"))

# ==============================================================================
# 4. 循环对 GC 和 DN 分别进行分析
# ==============================================================================
for (cell_label in names(receiver_cells_list)) {
    receiver_cells <- receiver_cells_list[[cell_label]]
    current_key_ligands <- core_ligands
    
    print(paste0("=================================================="))
    print(paste0("Starting analysis for ", cell_label, " (", receiver_cells, ")"))
    print(paste0("=================================================="))
    
    print("Defining background genes and genes of interest (DEGs of receivers)...")

    # 背景基因：在 Receiver 细胞中至少表达在 10% 的细胞中的基因
    receiver_cells_idx <- WhichCells(seurat_obj, idents = receiver_cells)
    
    if(length(receiver_cells_idx) == 0) {
        print(paste("Skipping", cell_label, "- No cells found."))
        next
    }
    
    expressed_genes_receiver <- seurat_obj@assays$RNA@counts[, receiver_cells_idx, drop=FALSE] %>% rowSums()
    background_expressed_genes <- names(expressed_genes_receiver)[expressed_genes_receiver > (length(receiver_cells_idx) * 0.1)]

    # 感兴趣基因集 (geneset_oi)：计算 DEGs
    all_cell_types <- levels(Idents(seurat_obj))
    reference_cells <- setdiff(all_cell_types, receiver_cells)
    
    print(paste("Calculating DEGs for", cell_label, "vs Reference..."))
    degs <- FindMarkers(seurat_obj, ident.1 = receiver_cells, ident.2 = reference_cells, logfc.threshold = 0.25, min.pct = 0.1)
    geneset_oi <- rownames(degs)[degs$p_val_adj < 0.05 & degs$avg_log2FC > 0.5]

    if (length(geneset_oi) < 5) {
        print("Warning: Very few DEGs found. Expanding the geneset_oi by lowering the threshold.")
        geneset_oi <- rownames(degs)[degs$p_val_adj < 0.1 & degs$avg_log2FC > 0.25]
    }
    print(paste("Number of genes in geneset_oi:", length(geneset_oi)))

    # 确保 key_ligands 存在于配体模型中
    valid_key_ligands <- intersect(current_key_ligands, colnames(ligand_target_matrix))
    if(length(valid_key_ligands) == 0) {
        print(paste("Skipping", cell_label, "- None of the key ligands are present in the NicheNet prior model."))
        next
    }
    print(paste("Targeted Ligands for", cell_label, ":", paste(valid_key_ligands, collapse=", ")))

    # ==============================================================================
    # 5. 提取靶向配体-靶基因网络 (Ligand-Target Links)
    # ==============================================================================
    print("Extracting ligand-target regulatory links...")

    active_ligand_target_links_df <- valid_key_ligands %>% 
        lapply(get_weighted_ligand_target_links, 
               geneset = geneset_oi, 
               ligand_target_matrix = ligand_target_matrix, 
               n = 200) %>% 
        bind_rows() %>% 
        drop_na()

    # 过滤低权重的连接
    active_ligand_target_links_df <- active_ligand_target_links_df[active_ligand_target_links_df$weight > 0.05, ]
    
    if(nrow(active_ligand_target_links_df) == 0) {
        print(paste("No strong ligand-target links found for", cell_label))
        next
    }

    # 准备绘图用的 Heatmap 矩阵
    active_ligand_target_links <- prepare_ligand_target_matrix(
        ligand_target_df = active_ligand_target_links_df, 
        ligands = valid_key_ligands, 
        targets = active_ligand_target_links_df$target
    )

    # 选取得分最高的前 50 个 Target 进行展示
    top_targets <- active_ligand_target_links_df %>% 
        group_by(target) %>% 
        summarise(max_weight = max(weight)) %>% 
        arrange(desc(max_weight)) %>% 
        head(50) %>% 
        pull(target)

    active_ligand_target_links_subset <- active_ligand_target_links[top_targets, valid_key_ligands, drop=FALSE]

    # ==============================================================================
    # 6. 可视化与保存 (Nature 风格)
    # ==============================================================================
    print(paste("Visualizing Network for", cell_label, "..."))

    heatmap_data <- as.data.frame(as.table(active_ligand_target_links_subset))
    colnames(heatmap_data) <- c("Target", "Ligand", "Weight")

    p_heatmap <- ggplot(heatmap_data, aes(x = Ligand, y = Target, fill = Weight)) +
        geom_tile(color = "white") +
        scale_fill_gradientn(colors = c("white", NPG_COLORS[1]), name = "Regulatory\nPotential") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(face = "italic", size = 8),
            axis.title = element_blank()
        ) +
        ggtitle(paste("Targeted Ligand-Target Network in", cell_label))

    pdf_file <- file.path(out_dir, paste0("Targeted_NicheNet_Downstream_", cell_label, "_Heatmap.pdf"))
    png_file <- file.path(out_dir, paste0("Targeted_NicheNet_Downstream_", cell_label, "_Heatmap.png"))

    ggsave(pdf_file, plot = p_heatmap, width = 6, height = 10, dpi = 600, device = cairo_pdf)
    ggsave(png_file, plot = p_heatmap, width = 6, height = 10, dpi = 600)

    write.csv(active_ligand_target_links_df, file.path(out_dir, paste0("Targeted_Downstream_Targets_", cell_label, ".csv")), row.names = FALSE)

    print(paste(cell_label, "分析完成！"))
}

print("所有靶向分析流程结束！")


