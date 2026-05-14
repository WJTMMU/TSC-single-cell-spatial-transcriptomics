
import os
import warnings

# 忽略警告
warnings.filterwarnings("ignore")

# ================= 数据路径配置 =================
BASE_DIR = r"d:\Data_Analysis\空间转录组分析\ST_analysis"
TSC_DATA_DIR = os.path.join(BASE_DIR, "结节性硬化数据")

# 原始方法代码路径 (供参考)
METHOD_DIR = os.path.join(BASE_DIR, "空转-2025")

# 样本列表
TSC_SAMPLES = [
    "25060668_JinJZ3_expression",
    "25060668_TSC_1_expression",
    "25060668_YanXR3_expression"
]

# ================= 输出路径配置 =================
# 主输出目录
OUTPUT_DIR = os.path.join(BASE_DIR, "New_Analysis", "Results")

# 各分析模块子目录
# 说明:
# - 01: 预处理与基础对象
# - 02: 注释与表达验证
# - 03: 病理差异与功能富集
# - 04: 空间生态位与核心基因
# - 05: 细胞通讯
# - 06: 调控网络
# - 07: 轨迹
# - 08: CNV
# - 09: 靶向整合
DIRS = {
    "basis": os.path.join(OUTPUT_DIR, "01_Basis"),
    "annotation": os.path.join(OUTPUT_DIR, "02_Annotation"),
    "pathology": os.path.join(OUTPUT_DIR, "03_Pathology_DEG_Enrichment"),
    "niche": os.path.join(OUTPUT_DIR, "04_Niche_Neighborhood"),
    "communication": os.path.join(OUTPUT_DIR, "05_Communication"),
    "regulation": os.path.join(OUTPUT_DIR, "06_Regulation_Coexp"),
    "trajectory": os.path.join(OUTPUT_DIR, "07_Trajectory"),
    "cnv": os.path.join(OUTPUT_DIR, "08_CNV"),
    "targeted": os.path.join(OUTPUT_DIR, "09_Targeted_Networks"),
}

# 自动创建目录
for d in DIRS.values():
    if not os.path.exists(d):
        os.makedirs(d)

# 统一的脚本级输出目录
STEP_LAYOUT = {
    "01": {"module": "basis", "folder": "01_Data_Preprocessing"},
    "02": {"module": "annotation", "folder": "02_Manual_Annotation"},
    "03": {"module": "annotation", "folder": "03_Manual_Annotation_Comb2"},
    "04": {"module": "annotation", "folder": "04_Check_Gene_Expression"},
    "05": {"module": "pathology", "folder": "05_SingleCell_DEGs"},
    "06": {"module": "pathology", "folder": "06_Functional_Enrichment"},
    "07": {"module": "niche", "folder": "07_Spatial_Niche_NMF"},
    "08": {"module": "niche", "folder": "08_NMF_Downstream"},
    "09": {"module": "niche", "folder": "09_Venn_Factor7_DNExcit"},
    "10": {"module": "niche", "folder": "10_Core_Genes_Visualization"},
    "11": {"module": "communication", "folder": "11_GC_DN_Neighborhood_Communication"},
    "12": {"module": "communication", "folder": "12_CellChat_Analysis"},
    "13": {"module": "communication", "folder": "13_Spatial_Communication_Visualization"},
    "14": {"module": "communication", "folder": "14_NicheNet_Analysis"},
    "15": {"module": "regulation", "folder": "15_Spatial_GRN_Squidpy"},
    "16": {"module": "regulation", "folder": "16_Spatial_GRN_Comparison"},
    "17": {"module": "regulation", "folder": "17_WGCNA_Analysis"},
    "18": {"module": "trajectory", "folder": "18_Trajectory_Pseudotime"},
    "19": {"module": "cnv", "folder": "19_Trajectory_CNV"},
    "20": {"module": "targeted", "folder": "20_Targeted_SCENIC_Upstream"},
    "21": {"module": "targeted", "folder": "21_Targeted_NicheNet_Downstream"},
}


def get_step_dir(step_id):
    """返回某个脚本步骤的标准输出目录，并自动创建。"""
    step_id = str(step_id).zfill(2)
    if step_id not in STEP_LAYOUT:
        raise KeyError(f"Unknown pipeline step: {step_id}")
    layout = STEP_LAYOUT[step_id]
    step_dir = os.path.join(DIRS[layout["module"]], layout["folder"])
    os.makedirs(step_dir, exist_ok=True)
    return step_dir


def get_step_file(step_id, filename):
    """返回某个脚本步骤下的标准文件路径。"""
    return os.path.join(get_step_dir(step_id), filename)


PIPELINE_FILES = {
    "processed_h5ad": get_step_file("01", "adata_tsc_processed.h5ad"),
    "cluster_markers_top50": get_step_file("01", "cluster_markers_top50.csv"),
    "annotated_h5ad": get_step_file("03", "adata_tsc_annotated.h5ad"),
    "deg_dir": get_step_dir("05"),
    "nmf_dir": get_step_dir("07"),
    "nmf_downstream_dir": get_step_dir("08"),
    "venn_dir": get_step_dir("09"),
    "core_genes_dir": get_step_dir("10"),
    "communication_dir": get_step_dir("11"),
    "communication_input_dir": os.path.join(get_step_dir("11"), "CellPhoneDB_Input"),
    "cellchat_dir": get_step_dir("12"),
    "spatial_comm_viz_dir": get_step_dir("13"),
    "nichenet_dir": get_step_dir("14"),
    "grn_dir": get_step_dir("15"),
    "grn_comparison_dir": get_step_dir("16"),
    "wgcna_dir": get_step_dir("17"),
    "trajectory_dir": get_step_dir("18"),
    "cnv_dir": get_step_dir("19"),
    "targeted_scenic_dir": get_step_dir("20"),
    "targeted_nichenet_dir": get_step_dir("21"),
}

os.makedirs(PIPELINE_FILES["communication_input_dir"], exist_ok=True)

# ================= Marker基因配置 =================
# 包含主要的脑细胞类型
MARKERS = {
    "Excitatory_Neurons": ["SLC17A7", "CAMK2A", "NRGN", "SATB2"],
    "Inhibitory_Neurons": ["GAD1", "GAD2", "SLC6A1", "PVALB", "SST"],
    "Astrocytes": ["GFAP", "AQP4", "ALDH1L1", "SLC1A2"],
    "Oligodendrocytes": ["MBP", "PLP1", "MOG", "CNP"],
    "OPC": ["PDGFRA", "CSPG4", "OLIG1"],
    "Microglia": ["AIF1", "TMEM119", "C1QA", "CX3CR1"],
    "Endothelial": ["CLDN5", "FLT1", "PECAM1", "VWF"],
    "Mural": ["ACTA2", "TAGLN", "PDGFRB"],
    "Ependymal": ["FOXJ1", "PIFO"]
}

# ================= Nature 级别全局配色配置 =================
# NPG (Nature Publishing Group) 经典调色板
NPG_COLORS = [
    '#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', 
    '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85',
    '#2E2A2B', '#E5A5C2', '#586E6E', '#C0C0C0', '#D4E2E0',
    '#8A2BE2', '#D2691E', '#20B2AA', '#FF1493', '#32CD32'
]

# 统一的连续色带 (连续值使用，如表达量)
# 推荐使用 viridis (双向) 或 plasma/magma (单向)

def set_nature_style():
    """设置全局 matplotlib 样式为 Nature 风格"""
    import matplotlib.pyplot as plt
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=NPG_COLORS)
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12

# ================= 辅助函数 =================
def check_paths():
    """检查数据路径是否存在"""
    print(f"Checking TSC Data Directory: {TSC_DATA_DIR}")
    if os.path.exists(TSC_DATA_DIR):
        print("  [OK] TSC Directory found.")
    else:
        print("  [ERROR] TSC Directory NOT found.")

    print("\nOutput Directories:")
    for name, path in DIRS.items():
        print(f"  {name}: {path}")

if __name__ == "__main__":
    check_paths()
