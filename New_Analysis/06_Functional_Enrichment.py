# ==============================================================================
# 文件名称: 06_Functional_Enrichment.py
# 功能描述: 功能富集分析
# 联动说明: 承接 05 步差异分析和注释后的表达对象，运行通路活性分析、手动基因集打分和空间自相关分析，揭示病理细胞及其空间区域的功能改变。
# ==============================================================================

import sys
import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# ==============================================================================
# 全局随机种子设置 (Global Seed Setting)
# 确保分析结果的完全可重复性
# ==============================================================================
import random
import numpy as np
try:
    import scanpy as sc
except ImportError:
    pass
try:
    import torch
except ImportError:
    pass

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
try:
    sc.settings.seed = SEED
except NameError:
    pass
try:
    torch.manual_seed(SEED)
    torch.cuda.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
except NameError:
    pass
# ==============================================================================


# 添加当前目录到 sys.path 以便导入配置
sys.path.append(os.getcwd())
try:
    from analysis_config import *
except ImportError:
    # Fallback if running from a different directory
    sys.path.append(r"d:\Data_Analysis\空间转录组分析\ST_analysis\New_Analysis")
    from analysis_config import *
from analysis_contract import compute_samplewise_spatial_neighbors, require_data_contract

# 忽略警告
warnings.filterwarnings("ignore")

# 绘图设置
sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (8, 8)
sc.set_figure_params(dpi=300, facecolor='white', vector_friendly=True)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# ================= 1. 定义关键通路基因集 (Fallback Knowledge Base) =================
# 如果 Decoupler 未安装，将使用这些预定义基因集进行打分
# 这些通路与 TSC/FCD 疾病机制高度相关 (mTOR, Inflammation, Hypoxia)
MANUAL_PATHWAYS = {
    "mTOR_Signaling": [
        "MTOR", "RPS6KB1", "RPS6", "EIF4EBP1", "AKT1", "PIK3CA", "TSC1", "TSC2", 
        "RHEB", "MLST8", "RPTOR", "RICTOR", "EIF4E"
    ],
    "Inflammation_Microglia": [
        "AIF1", "TMEM119", "C1QA", "C1QB", "HLA-DRA", "CD68", "TNF", "IL1B", 
        "IL6", "CCL2", "TLR4", "NFKB1"
    ],
    "Hypoxia_VEGF": [
        "HIF1A", "VEGFA", "SLC2A1", "LDHA", "PGK1", "CA9", "BNIP3", "DDIT4"
    ],
    "Neuronal_Excitability": [
        "SLC17A7", "CAMK2A", "GRIN1", "GRIN2B", "SCN1A", "SCN2A", "KCNQ2", "GABRA1"
    ],
    "Astrogliosis": [
        "GFAP", "VIM", "S100B", "AQP4", "STAT3", "OSMR", "LCN2", "TIMP1"
    ]
}

def save_plot(filename, fig=None):
    """同时保存 PDF 和 PNG，确保出版级质量"""
    if fig is None:
        fig = plt.gcf()
    
    base, ext = os.path.splitext(filename)
    # 强制保存 PDF 和 PNG
    try:
        fig.savefig(f"{base}.pdf", bbox_inches='tight', dpi=300)
        fig.savefig(f"{base}.png", bbox_inches='tight', dpi=300)
        print(f"Saved plots to {base}.pdf/png")
    except Exception as e:
        print(f"Error saving plots {filename}: {e}")

def run_pathway_analysis(adata, sample_type, output_dir):
    """
    运行通路活性分析。
    优先使用 Decoupler (PROGENy/CollecTRI)，如果未安装则使用 Scanpy score_genes。
    """
    print(f"\n--- Running Pathway Analysis for {sample_type} ---")
    require_data_contract(adata, "06_Functional_Enrichment.py")
    
    # 尝试导入 Decoupler
    has_decoupler = False
    try:
        import decoupler as dc
        has_decoupler = True
        print("Decoupler detected. Using PROGENy and CollecTRI for inference.")
    except ImportError:
        print("Decoupler NOT found. Using manual gene set scoring (sc.tl.score_genes).")

    if has_decoupler:
        try:
            # 1. PROGENy (Pathways)
            print("Estimating Pathway Activities (PROGENy)...")
            # 自动下载模型 (human) - 注意 TSC/FCD 通常也是 Human 数据
            # 如果是 Mouse，需改 organism='mouse'
            # 假设数据是 Human (基于基因名全大写判断)
            is_human = any(g.isupper() for g in adata.var_names[:10])
            organism = 'human' if is_human else 'mouse'
            
            progeny = dc.get_progeny(organism=organism, top=100)
            
            # 运行 MLM (Multivariate Linear Model)
            dc.run_mlm(mat=adata, net=progeny, source='source', target='target', weight='weight', verbose=False)
            
            # 提取结果并保存到 obsm
            acts = dc.get_acts(adata, obsm_key='mlm_estimate')
            adata.obsm['progeny_mlm'] = acts.copy()
            
            # 绘图：Top 3 Pathways
            # 选取方差最大的通路
            top_paths = acts.var.std().sort_values(ascending=False).head(3).index.tolist()
            
            # 空间绘图 (如果可用)
            if 'spatial' in adata.obsm:
                # 按样本循环绘图
                samples = adata.obs['sample'].unique()
                for sample_id in samples:
                    adata_sub = adata[adata.obs['sample'] == sample_id].copy()
                    # 优化可视化参数
                    sc.pl.embedding(adata_sub, basis="spatial", color=top_paths, cmap='vlag', vcenter=0, size=5, frameon=False, show=False)
                    save_plot(os.path.join(output_dir, f"{sample_type}_{sample_id}_progeny_spatial.pdf"))
                    plt.close()
                    del adata_sub
            
            # 2. CollecTRI (Transcription Factors)
            print("Estimating TF Activities (CollecTRI)...")
            collectri = dc.get_collectri(organism=organism)
            dc.run_ulm(mat=adata, net=collectri, source='source', target='target', verbose=False)
            
            tf_acts = dc.get_acts(adata, obsm_key='ulm_estimate')
            adata.obsm['collectri_ulm'] = tf_acts.copy()
            
            # 绘图：Top 3 TFs
            top_tfs = tf_acts.var.std().sort_values(ascending=False).head(3).index.tolist()
            if 'spatial' in adata.obsm:
                # 按样本循环绘图
                samples = adata.obs['sample'].unique()
                for sample_id in samples:
                    adata_sub = adata[adata.obs['sample'] == sample_id].copy()
                    sc.pl.embedding(adata_sub, basis="spatial", color=top_tfs, cmap='RdBu_r', vcenter=0, show=False)
                    save_plot(os.path.join(output_dir, f"{sample_type}_{sample_id}_collectri_tf_spatial.pdf"))
                    plt.close()
                    del adata_sub

        except Exception as e:
            print(f"Error running Decoupler: {e}")
            print("Falling back to manual scoring.")
            has_decoupler = False # Fallback to manual

    if not has_decoupler:
        # 使用 sc.tl.score_genes
        print("Calculating scores for manual gene sets...")
        
        # 检查是否存在 raw 数据以获取完整基因列表
        scoring_adata = adata
        use_raw = False
        if adata.raw is not None:
            print(f"Using adata.raw for scoring ({adata.raw.shape[1]} genes available).")
            # 创建一个临时的 adata 用于打分
            scoring_adata = adata.raw.to_adata()
            use_raw = True
        
        for name, genes in MANUAL_PATHWAYS.items():
            # 过滤存在的基因
            valid_genes = [g for g in genes if g in scoring_adata.var_names]
            if len(valid_genes) < 2:
                print(f"Skipping {name}: too few genes found ({len(valid_genes)}).")
                continue
            
            try:
                sc.tl.score_genes(scoring_adata, gene_list=valid_genes, score_name=name)
                # 如果使用了临时 adata，将结果拷回主 adata
                if use_raw:
                    adata.obs[name] = scoring_adata.obs[name]
                print(f"  Scored {name} with {len(valid_genes)} genes.")
            except Exception as e:
                print(f"Error scoring {name}: {e}")
        
        # 释放内存
        if use_raw:
            del scoring_adata
            import gc
            gc.collect()
            
        # 绘图
        pathways_to_plot = [p for p in MANUAL_PATHWAYS.keys() if p in adata.obs.columns]
        if pathways_to_plot and 'spatial' in adata.obsm:
            samples = adata.obs['sample'].unique()
            for sample_id in samples:
                adata_sub = adata[adata.obs['sample'] == sample_id].copy()
                # 优化：缩小点大小，去除边框，使用发散高对比度色板
                sc.pl.embedding(adata_sub, basis="spatial", color=pathways_to_plot, cmap='coolwarm', vcenter=0, size=5, frameon=False, show=False)
                save_plot(os.path.join(output_dir, f"{sample_type}_{sample_id}_manual_pathways_spatial.pdf"))
                plt.close()
                del adata_sub
            
            # 另外画一个小提琴图对比不同 Cluster
            cluster_key = 'leiden' if 'leiden' in adata.obs else 'celltype'
            if cluster_key in adata.obs:
                sc.pl.violin(adata, keys=pathways_to_plot, groupby=cluster_key, rotation=90, show=False)
                save_plot(os.path.join(output_dir, f"{sample_type}_manual_pathways_violin.pdf"))
                plt.close()

def analyze_spatial_modules(adata, sample_type, output_dir):
    """
    分析空间高变基因 (Spatially Variable Genes, SVGs)。
    使用 Squidpy Moran's I。
    """
    print(f"\n--- Running Spatial Module Analysis for {sample_type} ---")
    
    # 确保有空间邻居图
    if 'spatial' not in adata.obsm:
        print("No spatial coordinates found. Skipping.")
        return

    try:
        # 按样本循环计算 Moran's I
        samples = adata.obs['sample'].unique()
        for sample_id in samples:
            print(f"\n[{sample_id}] Calculating Moran's I for spatial autocorrelation...")
            adata_sub = adata[adata.obs['sample'] == sample_id].copy()
            
            # 计算空间邻居
            compute_samplewise_spatial_neighbors(adata_sub, spatial_key="spatial", sample_key="sample")
            
            # 计算 Moran's I
            sq.gr.spatial_autocorr(adata_sub, mode="moran")
            
            # 获取 Top 空间基因
            moran_df = adata_sub.uns["moranI"]
            moran_df.to_csv(os.path.join(output_dir, f"{sample_type}_{sample_id}_moranI.csv"))
            
            top_svgs = moran_df.index[:4].tolist()
            print(f"[{sample_id}] Top SVGs: {top_svgs}")
            
            # 绘图
            if top_svgs:
                sc.pl.embedding(adata_sub, basis="spatial", color=top_svgs, cmap='viridis', show=False)
                save_plot(os.path.join(output_dir, f"{sample_type}_{sample_id}_top_svgs.pdf"))
                plt.close()
            
            del adata_sub
            import gc
            gc.collect()
        
    except Exception as e:
        print(f"Error in spatial module analysis: {e}")

def main():
    print("Step 06: Functional Enrichment & Regulation...")
    
    output_dir = get_step_dir("06")

    # --- 1. TSC Analysis ---
    tsc_path = PIPELINE_FILES["annotated_h5ad"]
    if os.path.exists(tsc_path):
        print(f"Loading TSC data from {tsc_path}...")
        adata_tsc = sc.read_h5ad(tsc_path)
        run_pathway_analysis(adata_tsc, "TSC", output_dir)
        analyze_spatial_modules(adata_tsc, "TSC", output_dir)
        del adata_tsc
    else:
        print("TSC data not found. Skipping.")

    print("\nStep 06 Done.")

if __name__ == "__main__":
    main()

