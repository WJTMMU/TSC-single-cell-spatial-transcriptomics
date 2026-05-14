# ==============================================================================
# 文件名称: 19_Trajectory_CNV.py
# 功能描述: 结合拷贝数变异 (CNV) 的轨迹分析
# 联动说明: 在当前流程的轨迹分析基础上，引入 CNV 分析，验证沿着拟时序演进过程中是否存在基因组层面的异常改变。
# ==============================================================================

import sys
import os
import matplotlib
matplotlib.use('Agg') # Force non-interactive backend
import matplotlib.pyplot as plt
# Set font type for editable PDFs
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

import scanpy as sc
# Set n_jobs=1 globally to avoid Windows multiprocessing issues
sc.settings.n_jobs = 1
import numpy as np
import pandas as pd
import seaborn as sns
import infercnvpy as cnv
import mygene

# Add current directory to path to import config
sys.path.append(os.getcwd())
from analysis_config import *

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


def save_plot(filename, fig=None):
    """Save plot to PDF and PNG with high resolution."""
    if fig is None:
        fig = plt.gcf()
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    base, ext = os.path.splitext(filename)
    fig.savefig(f"{base}.pdf", bbox_inches='tight', dpi=300)
    fig.savefig(f"{base}.png", bbox_inches='tight', dpi=300)
    print(f"Saved plots to {base}.pdf/png")

def build_cnv_input_adata(adata):
    """
    尽可能恢复适合 CNV 的全基因原始计数矩阵。
    优先使用 raw 中的全基因 var，再使用 counts 层作为表达矩阵。
    """
    if 'counts' in adata.layers and adata.layers['counts'].shape == adata.shape:
        print("Building CNV input from full-gene counts layer...")
        adata_cnv = adata.copy()
        adata_cnv.X = adata.layers['counts'].copy()
        return adata_cnv

    if adata.raw is not None and hasattr(adata.raw, "var") and adata.raw.n_vars > adata.n_vars:
        print("Building CNV input from full-gene raw object...")
        adata_cnv = adata.raw.to_adata()
        adata_cnv.obs = adata.obs.copy()
        adata_cnv.obsm = {k: v.copy() for k, v in adata.obsm.items()}
        if 'counts' in adata.layers and adata.layers['counts'].shape == adata_cnv.shape:
            adata_cnv.X = adata.layers['counts'].copy()
        else:
            print("Warning: Full-gene counts layer is unavailable. Falling back to full-gene log-normalized raw matrix for CNV input.")
        return adata_cnv

    print("Warning: Full-gene raw.var is unavailable. CNV will use the current feature set only.")
    return adata.copy()

def fetch_gene_positions(adata):
    """
    Fetch genomic positions for genes in adata using MyGene.info.
    Adds 'chromosome', 'start', 'end' to adata.var.
    """
    print("Fetching genomic positions using MyGene.info...")
    mg = mygene.MyGeneInfo()
    genes = adata.var_names.tolist()
    
    # Batch query (species='human' for TSC/FCD)
    # scopes='symbol' assumes var_names are gene symbols
    # fields='genomic_pos,symbol'
    print(f"Querying {len(genes)} genes...")
    results = mg.querymany(genes, scopes='symbol', fields='genomic_pos', species='human', verbose=False)
    
    # Parse results
    gene_pos = {}
    for res in results:
        if 'notfound' in res:
            continue
        query = res['query']
        if 'genomic_pos' in res:
            pos = res['genomic_pos']
            # genomic_pos can be a list (multiple loci) or dict (single locus)
            # We take the first one if list, or the dict itself
            if isinstance(pos, list):
                pos = pos[0]
            
            if isinstance(pos, dict):
                chrom = pos.get('chr')
                start = pos.get('start')
                end = pos.get('end')
                
                # Check validity
                if chrom and start and end:
                    # Normalize chromosome name (e.g., '1' -> 'chr1')
                    # infercnvpy usually expects just '1', '2', 'X' or 'chr1' depending on config.
                    # Standardize to 'chrX' format or just 'X'?
                    # Let's check infercnvpy defaults. It handles both usually.
                    gene_pos[query] = {
                        'chromosome': str(chrom),
                        'start': int(start),
                        'end': int(end)
                    }
    
    # Update adata.var
    print(f"Found positions for {len(gene_pos)} genes.")
    sys.stdout.flush()
    
    # Create DataFrame to merge
    try:
        df_pos = pd.DataFrame.from_dict(gene_pos, orient='index')
        print("Created DataFrame from positions.")
        sys.stdout.flush()
    except Exception as e:
        print(f"Error creating DataFrame: {e}")
        return adata

    # Merge with adata.var
    # We only keep genes with positions for CNV
    original_gene_count = adata.n_vars
    try:
        mask = adata.var_names.isin(df_pos.index)
        print(f"Slicing adata with mask (keeping {sum(mask)} genes)...")
        sys.stdout.flush()
        adata_cnv = adata[:, mask].copy()
        print("Created adata_cnv copy.")
        sys.stdout.flush()
    except Exception as e:
        print(f"Error slicing adata: {e}")
        return adata
    
    # Add columns
    adata_cnv.var['chromosome'] = df_pos.loc[adata_cnv.var_names, 'chromosome']
    adata_cnv.var['start'] = df_pos.loc[adata_cnv.var_names, 'start']
    adata_cnv.var['end'] = df_pos.loc[adata_cnv.var_names, 'end']
    
    # Filter chromosomes
    # Keep only standard chromosomes 1-22, X, Y, MT
    # Normalize to "chrX" format
    def normalize_chrom(c):
        c = str(c)
        if c.startswith('chr'):
            return c
        if c in ['X', 'Y', 'MT', 'M']:
            return f"chr{c}"
        if c.isdigit() and 1 <= int(c) <= 22:
            return f"chr{c}"
        return None # Filter out others
    
    adata_cnv.var['chromosome'] = adata_cnv.var['chromosome'].apply(normalize_chrom)
    
    # Drop genes with invalid chromosomes
    if adata_cnv.var['chromosome'].isna().any():
        n_dropped = adata_cnv.var['chromosome'].isna().sum()
        print(f"Dropping {n_dropped} genes with non-standard chromosomes.")
        adata_cnv = adata_cnv[:, ~adata_cnv.var['chromosome'].isna()].copy()
        
    print(f"Filtered adata for CNV: {adata_cnv.n_vars} genes (from {original_gene_count}).")
    print("Chromosome distribution:")
    print(adata_cnv.var['chromosome'].value_counts())
    sys.stdout.flush()
    return adata_cnv

def analyze_cnv(adata, output_dir):
    """
    Perform CNV analysis using infercnvpy.
    """
    print("Running CNV Analysis (InferCNVpy)...")
    
    # 0. Build a dedicated CNV input object
    adata = build_cnv_input_adata(adata)

    # 1. Prepare Genomic Positions
    # Check if positions already exist
    required_cols = ['chromosome', 'start', 'end']
    if not all(col in adata.var.columns for col in required_cols):
        adata = fetch_gene_positions(adata)
        if adata.n_vars < 100:
            print("Error: Too few genes with genomic positions found. Skipping CNV.")
            return

    # 2. Setup Reference
    # We need a reference group (normal cells). 
    # In TSC/FCD, maybe 'Microglia' or 'Endothelial' or 'Oligodendrocytes' can be reference if assumed normal?
    # Or use a generic average.
    # Let's try to find a suitable reference from 'celltype' or 'leiden'.
    
    cluster_key = "celltype"
    if cluster_key not in adata.obs.columns:
        cluster_key = "leiden"
    
    print(f"Using '{cluster_key}' for grouping.")
    groups = adata.obs[cluster_key].unique().tolist()
    print(f"Available groups: {groups}")
    
    # Heuristic for reference: Microglia or Endothelial often used as reference in brain
    ref_group = []
    for g in groups:
        g_str = str(g).lower()
        if "microglia" in g_str or "endothelial" in g_str or "immune" in g_str:
            ref_group.append(g)
    
    if not ref_group:
        # If no specific reference found, pick the largest group or just skip reference (all vs average)
        # infercnvpy requires reference_key usually.
        # If no reference, we can use the rest as reference for each cluster?
        # Let's pick the first cluster as reference if nothing better found, or user defined.
        print("No standard reference cell type found (Microglia/Endothelial). Using all cells as reference (average).")
        ref_group = None 
    else:
        print(f"Using reference groups: {ref_group}")

    # 3. Run InferCNV
    try:
        # infercnvpy.tl.infercnv
        # window_size: 100-200 for single cell
        print("Running infercnvpy.tl.infercnv...")
        # n_jobs=1 to avoid multiprocessing issues on Windows
        # Use a smaller window size if gene count is low
        # But for CNV, usually window=100.
        # Filter out chromosomes with < 20 genes?
        # infercnvpy might handle it, but let's be safe.
        
        counts = adata.var['chromosome'].value_counts()
        small_chroms = counts[counts < 10].index.tolist()
        if small_chroms:
            print(f"Dropping chromosomes with < 10 genes: {small_chroms}")
            adata = adata[:, ~adata.var['chromosome'].isin(small_chroms)].copy()
        
        cnv.tl.infercnv(
            adata, 
            reference_key=cluster_key, 
            reference_cat=ref_group, 
            window_size=100,
            step=10,
            exclude_chromosomes=['chrX', 'chrY', 'chrM'], # Standardized names
            n_jobs=1
        )
        
        # 4. Clustering on CNV profiles
        print("Clustering CNV profiles...")
        cnv.tl.pca(adata)
        cnv.pp.neighbors(adata)
        cnv.tl.leiden(adata, key_added="cnv_leiden")
        
        # 5. Visualization
        print("Plotting CNV Heatmap...")
        plt.figure(figsize=(12, 10))
        # heatmap requires 'cnv' in obsm or X_cnv? 
        # infercnvpy stores in adata.obsm['X_cnv'] usually.
        # Plotting
        # Note: the `save` parameter in scanpy/infercnvpy plotting functions 
        # is often just appended to the default prefix, not treated as an absolute path.
        # So we pass True or a string suffix and move it later, or just save the figure directly.
        cnv.pl.chromosome_heatmap(
            adata, 
            groupby=cluster_key, 
            show=False
        )
        save_plot(os.path.join(output_dir, "TSC_CNV_heatmap_by_celltype.pdf"))
        plt.close()

        cnv.pl.chromosome_heatmap(
            adata, 
            groupby="cnv_leiden", 
            show=False
        )
        save_plot(os.path.join(output_dir, "TSC_CNV_heatmap_by_cnv_cluster.pdf"))
        plt.close()
        
        # Plot CNV score on UMAP
        # Calculate CNV score (e.g. mean absolute deviation from baseline)
        # infercnvpy has 'cnv_score' usually? Or we calculate it.
        # simple score: mean of squared differences from 0 (if normalized to 0)
        if 'X_cnv' in adata.obsm:
            cnv_matrix = adata.obsm['X_cnv']
            # assuming X_cnv is log-ratio or similar, centered at 0?
            # infercnvpy result is usually smoothed expression.
            # We can calculate a 'CNV Load'
            cnv_score = np.mean(np.abs(cnv_matrix), axis=1)
            adata.obs['cnv_score'] = cnv_score
            
            plt.figure(figsize=(8, 8))
            sc.pl.umap(adata, color=['cnv_score', 'cnv_leiden'], show=False)
            save_plot(os.path.join(output_dir, "TSC_CNV_score_UMAP.pdf"))
            plt.close()
            
            # 核心图表：双谱系收敛图 (Lineage Convergence Hypothesis)
            # 将兴奋性和抑制性两条独立的轨迹绘制在同一张图上
            print("Plotting Converging CNV Burden vs Pseudotime...")
            
            lineages = {
                'Excitatory': ['Excitatory Neurons', 'Dysmorphic Neurons(Excit)', 'Giant Cells'],
                'Inhibitory': ['Inhibitory Neurons', 'Dysmorphic Neurons(Inhib)', 'Giant Cells']
            }
            
            df_list = []
            for lin_name, lin_cats in lineages.items():
                dpt_col = f'dpt_pseudotime_{lin_name}'
                if dpt_col in adata.obs.columns and cluster_key in adata.obs.columns:
                    mask = adata.obs[cluster_key].isin(lin_cats) & adata.obs[dpt_col].notna()
                    df_sub = adata.obs[mask].copy()
                    df_sub['Plot_Pseudotime'] = df_sub[dpt_col]
                    df_sub['Lineage_Path'] = lin_name
                    df_list.append(df_sub)
            
            if df_list:
                df_plot = pd.concat(df_list)
                
                # 移除未使用的分类因子
                if hasattr(df_plot[cluster_key], 'cat'):
                    df_plot[cluster_key] = df_plot[cluster_key].cat.remove_unused_categories()
                
                # 拉长图片比例
                plt.figure(figsize=(14, 6), dpi=600)
                
                # 绘制散点（基于索引去重，防止共同终点 GC 的点被画两次重叠变深）
                df_scatter = df_plot[~df_plot.index.duplicated(keep='first')]
                
                sns.scatterplot(
                    data=df_scatter, 
                    x='Plot_Pseudotime', 
                    y='cnv_score', 
                    hue=cluster_key, 
                    palette='Set1',
                    alpha=0.6, 
                    s=50,
                    edgecolor='black',
                    linewidth=0.2
                )
                
                # 分别绘制两条趋势线
                # 1. 兴奋性路径
                df_exc = df_plot[df_plot['Lineage_Path'] == 'Excitatory']
                if not df_exc.empty:
                    sns.regplot(
                        data=df_exc, 
                        x='Plot_Pseudotime', 
                        y='cnv_score', 
                        scatter=False, 
                        color='black', 
                        lowess=True, 
                        line_kws={'linestyle': '-', 'linewidth': 2.5, 'label': 'Excitatory Path'}
                    )
                
                # 2. 抑制性路径
                df_inh = df_plot[df_plot['Lineage_Path'] == 'Inhibitory']
                if not df_inh.empty:
                    sns.regplot(
                        data=df_inh, 
                        x='Plot_Pseudotime', 
                        y='cnv_score', 
                        scatter=False, 
                        color='black', 
                        lowess=True, 
                        line_kws={'linestyle': '--', 'linewidth': 2.5, 'label': 'Inhibitory Path'}
                    )
                
                plt.title("Converging Pathological Trajectories Driven by Genomic Instability", fontweight='bold', fontsize=15)
                plt.xlabel("Pathological Deviation (Pseudotime)", fontweight='bold', fontsize=12)
                plt.ylabel("CNV Burden Score", fontweight='bold', fontsize=12)
                
                # 调整图例
                handles, labels = plt.gca().get_legend_handles_labels()
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Cell Lineage & Paths")
                sns.despine()
                
                save_plot(os.path.join(output_dir, "Hypothesis_CNV_vs_Pseudotime_Converging.pdf"))
                plt.close()
                print("Saved Hypothesis_CNV_vs_Pseudotime_Converging.pdf")
            else:
                print("Warning: Missing lineage pseudotime data. Make sure analyze_trajectory ran successfully.")
            
    except Exception as e:
        print(f"Error in CNV workflow: {e}")
        import traceback
        traceback.print_exc()

def analyze_trajectory(adata, output_dir):
    """
    Perform trajectory analysis using PAGA and DPT.
    """
    print("Running Trajectory Analysis (PAGA & DPT)...")
    
    # Save a reference to the original adata so we can map pseudotime back
    original_adata = adata
    
    # 1. Choose clustering key
    cluster_key = "celltype"
    if cluster_key not in adata.obs.columns:
        print(f"'{cluster_key}' not found. Falling back to 'leiden'.")
        cluster_key = "leiden"
    
    # Optional: Filter out non-neural cell types (Endothelial, Microglia, T_Cells, etc.)
    # because they don't share a developmental lineage with GC/DN/Neurons
    # and their inclusion can cause spurious connections in PAGA.
    if cluster_key in adata.obs.columns:
        cats = adata.obs[cluster_key].unique().tolist()
        exclude_keywords = ['endothelial', 'microglia', 'mural', 'immune', 'macrophage', 'fibroblast', 'blood', 'vascular', 'ependymal', 't_cell', 'b_cell', 't_cells', 't cells']
        
        keep_cats = []
        for cat in cats:
            cat_str = str(cat).lower()
            # 确保巨细胞 (Giant Cells) 不会被意外过滤
            if 'giant' in cat_str or 'gc' in cat_str:
                keep_cats.append(cat)
            elif not any(kw in cat_str for kw in exclude_keywords):
                keep_cats.append(cat)
                
        if len(keep_cats) < len(cats):
            print(f"Filtering out non-neural lineages to prevent spurious connections.")
            print(f"Keeping categories: {keep_cats}")
            adata = adata[adata.obs[cluster_key].isin(keep_cats)].copy()
            # Remove unused categories to prevent plotting errors
            if hasattr(adata.obs[cluster_key], 'cat'):
                adata.obs[cluster_key] = adata.obs[cluster_key].cat.remove_unused_categories()
        else:
            print("No non-neural cell types found to filter.")

    if cluster_key not in adata.obs.columns:
        print("No clustering found (leiden/celltype). Running Leiden clustering...")
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
        cluster_key = "leiden"

    # 2. PAGA (Partition-based Graph Abstraction)
    print(f"Calculating PAGA using '{cluster_key}'...")
    sys.stdout.flush()
    # Always recompute neighbors if we subsetted the data
    print("Computing neighbors for the trajectory lineage...")
    try:
        # It's safer to recompute PCA and neighbors on the subsetted neuro-lineage data
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, use_rep='X_pca')
    except Exception as e:
        print(f"Error computing neighbors: {e}")
        
    print("Running PAGA...")
    sys.stdout.flush()
    sc.tl.paga(adata, groups=cluster_key)
    
    # Plot PAGA
    plt.figure(figsize=(10, 8))
    sc.pl.paga(adata, color=cluster_key, show=False)
    save_plot(os.path.join(output_dir, "TSC_PAGA_graph.pdf"))
    plt.close() # sc.pl.paga might not create a figure object returned by gcf() cleanly if used with frame, but usually fine.
    
    # 3. Diffusion Pseudotime (DPT) - 计算两条独立轨迹
    print("Calculating Diffusion Pseudotime (DPT) for Split Lineages...")
    
    lineages = {
        'Excitatory': ['Excitatory Neurons', 'Dysmorphic Neurons(Excit)', 'Giant Cells'],
        'Inhibitory': ['Inhibitory Neurons', 'Dysmorphic Neurons(Inhib)', 'Giant Cells']
    }
    
    for lin_name, lin_cats in lineages.items():
        print(f"\n--- Processing {lin_name} Lineage ---")
        # 为每条轨迹提取对应的细胞子集
        lin_adata = original_adata[original_adata.obs[cluster_key].isin(lin_cats)].copy()
        
        if lin_adata.n_obs < 10:
            print(f"Skipping {lin_name} lineage: not enough cells.")
            continue
            
        try:
            sc.tl.pca(lin_adata)
            sc.pp.neighbors(lin_adata, use_rep='X_pca')
            
            # Find root cell: 优先选择正常的神经元
            root_cat = next((c for c in lin_cats if 'Normal' in c or 'Excitatory Neurons' == c or 'Inhibitory Neurons' == c), lin_cats[0])
            root_cells = lin_adata.obs.index[lin_adata.obs[cluster_key] == root_cat]
            
            if len(root_cells) > 0:
                root_cell_idx = np.where(lin_adata.obs.index == root_cells[0])[0][0]
                lin_adata.uns['iroot'] = root_cell_idx
                
                sc.tl.diffmap(lin_adata)
                sc.tl.dpt(lin_adata)
                
                # 映射回主 adata
                dpt_col = f'dpt_pseudotime_{lin_name}'
                original_adata.obs[dpt_col] = np.nan
                original_adata.obs.loc[lin_adata.obs.index, dpt_col] = lin_adata.obs['dpt_pseudotime']
                
                # Plot DPT for this lineage
                sc.pl.umap(lin_adata, color=['dpt_pseudotime', cluster_key], show=False)
                save_plot(os.path.join(output_dir, f"TSC_DPT_pseudotime_{lin_name}.pdf"))
                plt.close()
                
            else:
                print(f"Could not find root cells for {lin_name}. Skipping DPT.")
        except Exception as e:
            print(f"Error calculating DPT for {lin_name}: {e}")

    # 4. Draw PAGA initialization on UMAP
    sc.pl.paga_compare(adata, show=False)
    save_plot(os.path.join(output_dir, "TSC_PAGA_compare.pdf"))
    plt.close()

def analyze_cnv_dummy(adata, output_dir):
    """
    Placeholder for CNV analysis when infercnvpy is missing.
    """
    print("\n" + "="*50)
    print("CNV Analysis Warning")
    print("="*50)
    print("Could not run advanced CNV analysis (InferCNV) because:")
    print("1. 'infercnvpy' package is not installed.")
    print("2. Gene genomic position data (GTF) is not available in the environment.")
    print("3. 'mygene' or 'pybiomart' are missing for automatic fetching.")
    print("\nRecommendation: To perform CNV analysis, please install 'infercnvpy' and provide a GTF file.")
    print("Creating a placeholder text file in results.")
    
    with open(os.path.join(output_dir, "CNV_Analysis_README.txt"), "w") as f:
        f.write("CNV Analysis Skipped.\n")
        f.write("Missing dependencies: infercnvpy, genomic annotations.\n")
        f.write("Please install 'infercnvpy' and configure genomic positions to run this step.\n")

def main():
    print("Step 19: Trajectory & CNV Analysis Started...")
    
    # 标准化输出目录：轨迹结果与 CNV 结果分别进入各自脚本目录
    traj_dir = PIPELINE_FILES["trajectory_dir"]
    cnv_dir = PIPELINE_FILES["cnv_dir"]
        
    os.makedirs(traj_dir, exist_ok=True)
    os.makedirs(cnv_dir, exist_ok=True)
    
    # Load TSC Data (Single Cell resolution is best for Trajectory)
    tsc_path = PIPELINE_FILES["annotated_h5ad"]
    if not os.path.exists(tsc_path):
        print(f"Error: TSC data not found at {tsc_path}")
        return

    print(f"Loading TSC data from {tsc_path}...")
    adata_tsc = sc.read_h5ad(tsc_path)
    adata_tsc.obs_names_make_unique()
    
    # Run Trajectory
    try:
        # analyze_trajectory will compute dpt_pseudotime and modify adata_tsc in-place
        analyze_trajectory(adata_tsc, traj_dir)
    except Exception as e:
        print(f"Error in Trajectory Analysis: {e}")
        import traceback
        traceback.print_exc()

    # Run CNV (now adata_tsc contains dpt_pseudotime)
    analyze_cnv(adata_tsc, cnv_dir)
    
    print("Step 19 Completed.")

if __name__ == "__main__":
    main()
