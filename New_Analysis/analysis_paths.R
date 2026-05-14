args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (basename(script_dir) != "New_Analysis") {
    script_dir <- file.path(script_dir, "New_Analysis")
  }
}

base_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
results_dir <- file.path(base_dir, "New_Analysis", "Results")

module_dirs <- list(
  basis = file.path(results_dir, "01_Basis"),
  annotation = file.path(results_dir, "02_Annotation"),
  pathology = file.path(results_dir, "03_Pathology_DEG_Enrichment"),
  niche = file.path(results_dir, "04_Niche_Neighborhood"),
  communication = file.path(results_dir, "05_Communication"),
  regulation = file.path(results_dir, "06_Regulation_Coexp"),
  trajectory = file.path(results_dir, "07_Trajectory"),
  cnv = file.path(results_dir, "08_CNV"),
  targeted = file.path(results_dir, "09_Targeted_Networks")
)

step_layout <- list(
  "01" = list(module = "basis", folder = "01_Data_Preprocessing"),
  "02" = list(module = "annotation", folder = "02_Manual_Annotation"),
  "03" = list(module = "annotation", folder = "03_Manual_Annotation_Comb2"),
  "04" = list(module = "annotation", folder = "04_Check_Gene_Expression"),
  "05" = list(module = "pathology", folder = "05_SingleCell_DEGs"),
  "06" = list(module = "pathology", folder = "06_Functional_Enrichment"),
  "07" = list(module = "niche", folder = "07_Spatial_Niche_NMF"),
  "08" = list(module = "niche", folder = "08_NMF_Downstream"),
  "09" = list(module = "niche", folder = "09_Venn_Factor7_DNExcit"),
  "10" = list(module = "niche", folder = "10_Core_Genes_Visualization"),
  "11" = list(module = "communication", folder = "11_GC_DN_Neighborhood_Communication"),
  "12" = list(module = "communication", folder = "12_CellChat_Analysis"),
  "13" = list(module = "communication", folder = "13_Spatial_Communication_Visualization"),
  "14" = list(module = "communication", folder = "14_NicheNet_Analysis"),
  "15" = list(module = "regulation", folder = "15_Spatial_GRN_Squidpy"),
  "16" = list(module = "regulation", folder = "16_Spatial_GRN_Comparison"),
  "17" = list(module = "regulation", folder = "17_WGCNA_Analysis"),
  "18" = list(module = "trajectory", folder = "18_Trajectory_Pseudotime"),
  "19" = list(module = "cnv", folder = "19_Trajectory_CNV"),
  "20" = list(module = "targeted", folder = "20_Targeted_SCENIC_Upstream"),
  "21" = list(module = "targeted", folder = "21_Targeted_NicheNet_Downstream")
)

get_step_dir <- function(step_id) {
  step_id <- sprintf("%02d", as.integer(step_id))
  layout <- step_layout[[step_id]]
  if (is.null(layout)) {
    stop(paste("Unknown pipeline step:", step_id))
  }
  step_dir <- file.path(module_dirs[[layout$module]], layout$folder)
  dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
  step_dir
}

get_step_file <- function(step_id, filename) {
  file.path(get_step_dir(step_id), filename)
}

pipeline_files <- list(
  processed_h5ad = get_step_file("01", "adata_tsc_processed.h5ad"),
  cluster_markers_top50 = get_step_file("01", "cluster_markers_top50.csv"),
  annotated_h5ad = get_step_file("03", "adata_tsc_annotated.h5ad"),
  communication_input_dir = file.path(get_step_dir("11"), "CellPhoneDB_Input"),
  cellchat_dir = get_step_dir("12"),
  nichenet_dir = get_step_dir("14"),
  grn_comparison_dir = get_step_dir("16"),
  wgcna_dir = get_step_dir("17"),
  targeted_nichenet_dir = get_step_dir("21")
)

dir.create(pipeline_files$communication_input_dir, recursive = TRUE, showWarnings = FALSE)
