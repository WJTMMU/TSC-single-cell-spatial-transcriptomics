# ==============================================================================
# 文件名称: Tool_Custom_Enrichment_Plot.R
# 功能描述: 工具脚本: 自定义富集结果绘图 (R)
# 联动说明: 独立工具，用于将 Python 差异分析或 NMF 分析生成的 CSV 结果在 R 中绘制高质量的柱状图。
# ==============================================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(readr)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Set your working directory to where your CSV files are located
# setwd("D:/Data_Analysis/空间转录组分析/ST_analysis/Results/01_Basis/05_SingleCell_DEGs")
# OR
# setwd("D:/Data_Analysis/空间转录组分析/ST_analysis/Results/02_Niche/NMF_Interpretation")

# Provide the path to your CSV file here
# Example for 11 DEGs:
# csv_file_path <- "Enrichment_DN_vs_Normal_Up.csv"
# Example for 13 NMF:
# csv_file_path <- "Combined_AllSamples_NMF_Factors_Enrichment.csv"

# For demonstration, we will create a function that takes the file path as input
# ==============================================================================

plot_custom_enrichment <- function(csv_file_path, output_prefix = "Custom_Enrichment", top_n = 15, is_nmf = FALSE, target_factor = NULL, custom_terms = NULL) {
  
  # Ensure PDFs generate editable text (Type 42 font equivalent in R)
  pdf.options(useDingbats = FALSE)
  
  if (!file.exists(csv_file_path)) {
    stop(paste("File not found:", csv_file_path))
  }
  
  # Read the CSV file
  df <- read_csv(csv_file_path, show_col_types = FALSE)
  
  # If it's NMF results, filter by the specific factor if provided
  if (is_nmf && !is.null(target_factor)) {
    if ("Factor" %in% colnames(df)) {
      df <- df %>% filter(Factor == target_factor)
      output_prefix <- paste0(output_prefix, "_", target_factor)
    } else {
      warning("Column 'Factor' not found in NMF results. Proceeding with all data.")
    }
  }
  
  # Check required columns (handles both Python gseapy output formats)
  # Standard gseapy output has: Gene_set, Term, P-value, Adjusted P-value
  # NMF output has: Database, Term, P-value, Adjusted P-value
  
  db_col <- if("Gene_set" %in% colnames(df)) "Gene_set" else if("Database" %in% colnames(df)) "Database" else NULL
  pval_col <- if("Adjusted P-value" %in% colnames(df)) "Adjusted P-value" else if("pvals_adj" %in% colnames(df)) "pvals_adj" else NULL
  term_col <- "Term"
  
  if (is.null(db_col) || is.null(pval_col) || !(term_col %in% colnames(df))) {
    stop("Required columns not found. Ensure CSV has Term, P-value/Adjusted P-value, and Gene_set/Database columns.")
  }
  
  # Rename columns for consistency
  df <- df %>%
    rename(
      Database = !!sym(db_col),
      Adj_P_Value = !!sym(pval_col),
      Term_Name = !!sym(term_col)
    )
  
  # Calculate -log10(P-value)
  df <- df %>%
    mutate(Log_P = -log10(Adj_P_Value + 1e-300)) # Prevent log(0)
  
  # Clean up Term names (remove GO IDs for cleaner plots)
  # e.g., "synaptic transmission (GO:0007268)" -> "synaptic transmission"
  df <- df %>%
    mutate(Clean_Term = str_replace(Term_Name, "\\s*\\(GO:[0-9]+\\)", "")) %>%
    mutate(Clean_Term = str_replace(Clean_Term, "\\s*Homo sapiens\\s*hsa[0-9]+", "")) # Clean KEGG if needed
  
  # Truncate very long terms
  df <- df %>%
    mutate(Display_Term = ifelse(nchar(Clean_Term) > 50, 
                                 paste0(substr(Clean_Term, 1, 47), "..."), 
                                 Clean_Term))
  
  # Simplify Database names
  df <- df %>%
    mutate(Database_Type = case_when(
      str_detect(Database, "(?i)GO") ~ "GO",
      str_detect(Database, "(?i)KEGG") ~ "KEGG",
      TRUE ~ "Other"
    ))
  
  # Select Terms based on custom list OR Top N
  if (!is.null(custom_terms) && length(custom_terms) > 0) {
    # Filter by user-provided custom terms (matching raw or cleaned names)
    top_df <- df %>%
      filter(Term_Name %in% custom_terms | Clean_Term %in% custom_terms) %>%
      arrange(Log_P) # Re-arrange for plotting bottom-to-top
      
    if (nrow(top_df) == 0) {
      warning("None of the provided custom_terms were found in the data.")
      return(NULL)
    }
  } else {
    # Default behavior: Select Top N terms based on Adjusted P-value
    top_df <- df %>%
      arrange(Adj_P_Value) %>%
      head(top_n) %>%
      arrange(Log_P) # Re-arrange for plotting bottom-to-top
  }
  
  # Ensure factors are ordered correctly in the plot
  top_df$Display_Term <- factor(top_df$Display_Term, levels = top_df$Display_Term)
  
  # Define distinct colors
  my_colors <- c("GO" = "#1f77b4", "KEGG" = "#ff7f0e", "Other" = "#2ca02c")
  
  # Create the plot
  p <- ggplot(top_df, aes(x = Log_P, y = Display_Term, fill = Database_Type)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
    scale_fill_manual(values = my_colors) +
    theme_classic(base_size = 14) +
    labs(
      title = paste("Enrichment Analysis:", output_prefix),
      x = bquote("-Log"[10]~"(Adjusted P-value)"),
      y = NULL,
      fill = "Database"
    ) +
    theme(
      axis.text.y = element_text(color = "black", size = 11),
      axis.text.x = element_text(color = "black", size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank()
    )
  
  # Save the plot
  output_pdf <- paste0(output_prefix, "_Custom_Barplot.pdf")
  output_png <- paste0(output_prefix, "_Custom_Barplot.png")
  
  ggsave(output_pdf, plot = p, width = 8, height = max(5, length(top_df$Display_Term) * 0.3), dpi = 300)
  ggsave(output_png, plot = p, width = 8, height = max(5, length(top_df$Display_Term) * 0.3), dpi = 300)
  
  message(paste("Successfully saved custom plots to:", output_pdf, "and", output_png))
  
  return(p)
}

# ==============================================================================
# EXAMPLES OF HOW TO RUN (可以直接取消注释运行)
# ==============================================================================
cat("\n# === Instructions for running in RStudio (使用说明) === #\n")
cat("# 第一步：选中上面的所有代码并运行 (Run)，将函数加载到 R 的内存中。\n")
cat("# 第二步：修改下面的路径和参数，然后运行对应的例子。\n\n")

# ---------------------------------------------------------
# 例子 1：常规使用（自动取 P-value 最小的前 N 个通路画图）
# ---------------------------------------------------------
# # 1. 设置工作目录到你存放 CSV 文件的文件夹
# setwd("D:/Data_Analysis/空间转录组分析/ST_analysis/Results/01_Basis/05_SingleCell_DEGs")
# 
# # 2. 调用函数，画出最显著的 Top 15 个通路
# plot_custom_enrichment(
#   csv_file_path = "Enrichment_DN_vs_Normal_Up.csv", 
#   output_prefix = "DN_vs_Normal_Up", 
#   top_n = 15
# )

# ---------------------------------------------------------
# 例子 2：【高级】自定义挑选特定的通路画图（解决 Top N 无法完全展示的问题）
# ---------------------------------------------------------
# # 1. 设置工作目录
# setwd("D:/Data_Analysis/空间转录组分析/ST_analysis/Results/01_Basis/05_SingleCell_DEGs")
#
# # 2. 定义你想要重点展示的通路名称列表 (注意：名字要和 CSV 里 Term 列中的名称匹配，可以带 GO 编号也可以不带)
# my_picked_terms <- c(
#   "synaptic transmission", 
#   "axon guidance", 
#   "calcium ion binding",
#   "Pathways in cancer"
# )
# 
# # 3. 调用函数并传入 custom_terms 参数
# plot_custom_enrichment(
#   csv_file_path = "Enrichment_DN_vs_Normal_Up.csv", 
#   output_prefix = "DN_Custom_Selected", 
#   custom_terms = my_picked_terms
# )

# ---------------------------------------------------------
# 例子 3：处理 NMF 的结果 (由于 NMF 所有 Factor 都在一个表里，需要指定 target_factor)
# ---------------------------------------------------------
# # 1. 设置工作目录到 NMF 结果文件夹
# setwd("D:/Data_Analysis/空间转录组分析/ST_analysis/Results/02_Niche/NMF_Interpretation")
# 
# # 2. 调用函数，指定 is_nmf=TRUE 并告诉它画哪一个 Factor
# plot_custom_enrichment(
#   csv_file_path = "Combined_AllSamples_NMF_Factors_Enrichment.csv", 
#   output_prefix = "NMF", 
#   is_nmf = TRUE, 
#   target_factor = "Factor_0", 
#   top_n = 10
# )


