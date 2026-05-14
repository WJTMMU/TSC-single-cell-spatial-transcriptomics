# ==============================================================================
# 文件名称: 01_Data_Preprocessing.py
# 功能描述: 数据预处理、聚类与标记基因提取
# 联动说明: 作为整个分析流程的起点，读取原始空转数据，进行质控、降维聚类，并输出聚类结果和每个 cluster 的 top 50 markers 供下一步手动注释使用。
# ==============================================================================

import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import scanpy.external as sce

# 忽略警告
warnings.filterwarnings("ignore")

# 绘图设置
sc.settings.verbosity = 3
# sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white', vector_friendly=True)
plt.rcParams['figure.figsize'] = (8, 8)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# 引入配置文件中的常量
from analysis_config import *
from analysis_contract import stamp_data_contract

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


def load_seekspace_sample(sample_name, base_dir):
    sample_path = os.path.join(base_dir, sample_name)
    matrix_dir = None
    for root, dirs, files in os.walk(sample_path):
        if 'matrix.mtx.gz' in files:
            matrix_dir = root
            break
    
    if not matrix_dir:
        raise FileNotFoundError(f"Could not find matrix.mtx.gz in {sample_path}")
    
    print(f"Loading matrix from: {matrix_dir}")
    adata = sc.read_10x_mtx(matrix_dir, cache=True)
    
    loc_file = os.path.join(matrix_dir, "cell_locations.tsv.gz")
    if os.path.exists(loc_file):
        print(f"Loading locations from: {loc_file}")
        # Fix: header=0 because file has a header
        spatial_df = pd.read_csv(loc_file, sep='\t', header=0, index_col=0)
        # Standardize column names
        spatial_df.columns = ['x', 'y']
        
        common_cells = adata.obs_names.intersection(spatial_df.index)
        print(f"  Matched {len(common_cells)} cells out of {adata.n_obs} (Matrix) and {len(spatial_df)} (Locations)")
        
        if len(common_cells) == 0:
            print("  WARNING: No cells matched! Check barcode formats (e.g. -1 suffix).")
            # Try adding -1 to spatial_df index if adata has it
            if adata.obs_names[0].endswith('-1') and not spatial_df.index[0].endswith('-1'):
                print("  Attempting to add '-1' suffix to location barcodes...")
                spatial_df.index = spatial_df.index + '-1'
                common_cells = adata.obs_names.intersection(spatial_df.index)
                print(f"  Re-matched: {len(common_cells)} cells")
            # Try removing -1 from spatial_df index if adata doesn't have it but spatial_df does
            elif not adata.obs_names[0].endswith('-1') and spatial_df.index[0].endswith('-1'):
                print("  Attempting to remove '-1' suffix from location barcodes...")
                spatial_df.index = spatial_df.index.str.replace(r'-1$', '', regex=True)
                common_cells = adata.obs_names.intersection(spatial_df.index)
                print(f"  Re-matched: {len(common_cells)} cells")

        adata = adata[common_cells].copy()
        spatial_df = spatial_df.loc[common_cells]
        
        adata.obsm['spatial'] = spatial_df[['x', 'y']].values.astype(float)
        adata.obs['spatial_x'] = spatial_df['x'].astype(float)
        adata.obs['spatial_y'] = spatial_df['y'].astype(float)
    else:
        print(f"Warning: Location file not found at {loc_file}")
    
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    return adata

def auto_annotate_clusters(adata, marker_dict, groupby='leiden'):
    sc.tl.dendrogram(adata, groupby=groupby)
    
    cluster_map = {}
    clusters = adata.obs[groupby].unique()
    for cluster in clusters:
        scores = {}
        for cell_type, genes in marker_dict.items():
            valid_genes = [g for g in genes if g in adata.var_names]
            if not valid_genes: 
                continue
            sc.tl.score_genes(adata, gene_list=valid_genes, score_name=f'score_{cell_type}')
            scores[cell_type] = adata.obs[adata.obs[groupby] == cluster][f'score_{cell_type}'].mean()
        
        if scores:
            best_match = max(scores, key=scores.get)
            cluster_map[cluster] = best_match
        else:
            cluster_map[cluster] = "Unknown"
            
    return cluster_map

def main():
    step_dir = get_step_dir("01")
    # 1. 加载数据
    adatas = []
    for sample in TSC_SAMPLES:
        try:
            print(f"Processing {sample}...")
            ad = load_seekspace_sample(sample, TSC_DATA_DIR)
            adatas.append(ad)
        except Exception as e:
            print(f"Error loading {sample}: {e}")

    if not adatas:
        print("No data loaded!")
        return

    adata_tsc = sc.concat(adatas, join='outer', label='sample', keys=TSC_SAMPLES)
    print(f"Combined Data: {adata_tsc.shape}")

    # 2. 质控
    adata_tsc.var['mt'] = adata_tsc.var_names.str.startswith(('MT-', 'mt-'))
    sc.pp.calculate_qc_metrics(adata_tsc, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata_tsc, min_genes=200)
    sc.pp.filter_genes(adata_tsc, min_cells=3)
    adata_tsc = adata_tsc[adata_tsc.obs.pct_counts_mt < 20, :].copy()

    # 3. 预处理与高变基因
    # 确保 X 是整数类型 (concat 可能引入浮点数)
    from scipy import sparse
    if sparse.issparse(adata_tsc.X):
        adata_tsc.X = adata_tsc.X.astype(int)
    else:
        adata_tsc.X = adata_tsc.X.astype(int)

    # 保存原始计数到 layers 和 raw
    adata_tsc.layers['counts'] = adata_tsc.X.copy()
    adata_tsc.raw = adata_tsc # raw 现在保存的是全基因的原始计数
    print(f"Raw counts saved. Max value in raw: {adata_tsc.raw.X.max()}")
    
    sc.pp.normalize_total(adata_tsc, target_sum=1e4)
    sc.pp.log1p(adata_tsc)
    
    print("Selecting HVGs using Seurat method...")
    sc.pp.highly_variable_genes(
        adata_tsc, 
        n_top_genes=2000, 
        batch_key='sample',
        flavor='seurat'
    )
    # adata_tsc.raw = adata_tsc # 移除这行，避免覆盖为归一化数据
    adata_tsc = adata_tsc[:, adata_tsc.var['highly_variable']].copy()
    
    # 备份 Log-Normalized 数据用于差异分析 (避免使用 Scaled data 或 Raw counts)
    adata_tsc.layers['log1p'] = adata_tsc.X.copy()
    stamp_data_contract(adata_tsc, producer_script="01_Data_Preprocessing.py")

    # 4. 降维与整合
    sc.pp.scale(adata_tsc, max_value=10)
    
    if np.isnan(adata_tsc.X).any():
        print("Warning: NaNs found in X. Filling with 0...")
        adata_tsc.X = np.nan_to_num(adata_tsc.X, nan=0.0)

    print("Running PCA...")
    sc.tl.pca(adata_tsc, svd_solver='auto')
    
    if np.isnan(adata_tsc.obsm['X_pca']).any():
        print("Warning: NaNs found in X_pca. Filling with 0...")
        adata_tsc.obsm['X_pca'] = np.nan_to_num(adata_tsc.obsm['X_pca'], nan=0.0)

    print("Running Harmony integration...")
    try:
        sce.pp.harmony_integrate(adata_tsc, 'sample', basis='X_pca', adjusted_basis='X_pca_harmony')
        print("Harmony integration done.")
    except Exception as e:
        print(f"Harmony failed: {e}")
        print("Fallback: Using standard PCA.")
        adata_tsc.obsm['X_pca_harmony'] = adata_tsc.obsm['X_pca'].copy()

    # 5. 聚类
    sc.pp.neighbors(adata_tsc, use_rep='X_pca_harmony', n_neighbors=30)
    sc.tl.leiden(adata_tsc, resolution=0.2, key_added='leiden')
    sc.tl.umap(adata_tsc)

    # 6. 导出 Marker Genes 供手动注释
    print("Calculating marker genes for manual annotation...")
    # 使用 log1p 层进行差异分析，并显式禁用 use_raw
    sc.tl.rank_genes_groups(adata_tsc, 'leiden', method='wilcoxon', layer='log1p', use_raw=False)
    
    # 导出前50个基因
    result = adata_tsc.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    # 构建 DataFrame
    markers_df = pd.DataFrame(
        {group + '_' + key: result[key][group]
        for group in groups for key in ['names', 'scores', 'pvals_adj']}
    ).head(50)
    
    markers_csv_path = PIPELINE_FILES["cluster_markers_top50"]
    markers_df.to_csv(markers_csv_path)
    print(f"Exported top 50 marker genes to {markers_csv_path}")
    
    # 自动注释已移除，转为手动
    # print("Auto-annotating clusters...")
    # ... (原有自动注释代码已移除)

    # 7. 保存到 DIRS['basis']
    # 注意：此时还没有 'celltype' 列，将在 01.5 步骤中添加
    output_file = PIPELINE_FILES["processed_h5ad"]
    adata_tsc.write(output_file)
    print(f"Saved clustering result to {output_file}")
    
    sc.pl.umap(adata_tsc, color=['sample', 'leiden'], show=False)
    plt.savefig(os.path.join(step_dir, "umap_tsc_clustering.png"), bbox_inches='tight')
    print(f"Saved UMAP plot to {step_dir}")
    
    print("\n" + "="*50)
    print("Step 01 Finished!")
    print(f"Please check {markers_csv_path} to identify cell types.")
    print("Then open and edit '01_5_Manual_Annotation.py' to fill in the cell types.")
    print("Finally run '01_5_Manual_Annotation.py' to generate the annotated data for Step 03.")
    print("="*50 + "\n")

if __name__ == "__main__":
    main()
