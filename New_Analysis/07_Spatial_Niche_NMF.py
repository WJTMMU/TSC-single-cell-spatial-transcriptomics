# ==============================================================================
# 文件名称: 07_Spatial_Niche_NMF.py
# 功能描述: 基于 NMF 的空间生态位 (Niche) 分析
# 联动说明: 基于 02 步带有注释的数据，利用非负矩阵分解 (NMF) 提取空间微环境特征，识别特定的空间共定位模式和生态位。
# ==============================================================================

import sys
import os
sys.path.append(os.getcwd())
from analysis_config import *
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import NMF
from analysis_contract import (
    DEFAULT_EXPRESSION_LAYER,
    compute_samplewise_spatial_neighbors,
    get_expression_matrix,
    require_data_contract,
    save_nmf_artifacts,
)

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


# 设置绘图参数 (符合出版要求)
sc.set_figure_params(dpi=300, facecolor='white', vector_friendly=True)
set_nature_style()

def save_plot(filename, fig=None):
    """同时保存 PDF 和 PNG"""
    if fig is None:
        fig = plt.gcf()
    
    base, ext = os.path.splitext(filename)
    # 强制保存 PDF 和 PNG
    fig.savefig(f"{base}.pdf", bbox_inches='tight', dpi=300)
    fig.savefig(f"{base}.png", bbox_inches='tight', dpi=300)
    print(f"Saved plots to {base}.pdf/png")

def get_nmf_input_matrix(adata_view):
    """
    为 NMF 提供非负输入矩阵。
    统一使用数据契约定义的 log1p 表达层，禁止回退到 `.X`。
    """
    return get_expression_matrix(adata_view, layer=DEFAULT_EXPRESSION_LAYER, nonnegative=True)

def analyze_tsc_niche(adata, output_dir):
    """
    分析 TSC (SeekSpace, 单细胞分辨率) 的微环境
    """
    print("\n=== Analyzing TSC (Single-Cell Resolution) Niche ===", flush=True)
    require_data_contract(adata, "07_Spatial_Niche_NMF.py")
    
    # 1. 确定聚类 Key
    cluster_key = "celltype"
    if cluster_key not in adata.obs.columns:
        print(f"Warning: '{cluster_key}' not found in adata.obs. Falling back to 'leiden'.", flush=True)
        cluster_key = "leiden"
    
    print(f"Using cluster key: {cluster_key}", flush=True)

    # 按样本循环计算空间特征
    samples = adata.obs['sample'].unique()
    print(f"Found {len(samples)} samples. Processing independently...")

    for sample_id in samples:
        print(f"\n--- Processing Sample: {sample_id} ---", flush=True)
        # 获取单个样本的数据子集
        adata_sub = adata[adata.obs['sample'] == sample_id].copy()

        # 2. 构建空间邻居图 (按样本独立构建)
        print(f"[{sample_id}] Computing Spatial Neighbors (Generic Coordinates)...", flush=True)
        compute_samplewise_spatial_neighbors(adata_sub, spatial_key="spatial", sample_key="sample")

        # 3. 邻域富集分析 (Neighborhood Enrichment)
        print(f"[{sample_id}] Computing Neighborhood Enrichment...", flush=True)
        try:
            # 增加置换次数以获得更稳定的 P-value (Nature级别通常要求 n_perms >= 1000)
            sq.gr.nhood_enrichment(adata_sub, cluster_key=cluster_key, n_perms=1000)
            
            # 绘图 - 优化配色和标签，使用更高级的 cmap
            plt.figure(figsize=(10, 8))
            # cmap="vlag" 是一种对色弱友好的红蓝发散色板，比 coolwarm 更适合高水平期刊
            sq.pl.nhood_enrichment(adata_sub, cluster_key=cluster_key, method="ward", show=False, cmap="vlag", vmin=-1, vmax=1, annot=True, fmt=".2f")
            plt.title(f"Spatial Neighborhood Enrichment ({sample_id})", fontsize=16, pad=20)
            save_plot(os.path.join(output_dir, f"TSC_{sample_id}_nhood_enrichment.pdf"))
            plt.close()
        except Exception as e:
            print(f"[{sample_id}] Error in Neighborhood Enrichment: {e}", flush=True)

        # 4. 细胞类型交互共现 (Co-occurrence)
        print(f"[{sample_id}] Computing Co-occurrence probability...", flush=True)
        try:
            sq.gr.co_occurrence(adata_sub, cluster_key=cluster_key)
            
            # 绘制所有细胞类型的共现
            all_clusters = adata_sub.obs[cluster_key].unique().tolist()
            
            # 根据细胞类型的数量动态调整画布大小
            n_clusters = len(all_clusters)
            fig_width = min(20, max(12, n_clusters * 2))
            
            plt.figure(figsize=(fig_width, 6))
            sq.pl.co_occurrence(adata_sub, cluster_key=cluster_key, clusters=all_clusters)
            save_plot(os.path.join(output_dir, f"TSC_{sample_id}_co_occurrence.pdf"))
            plt.close()
        except Exception as e:
            print(f"[{sample_id}] Error in Co-occurrence: {e}", flush=True)

        # 5. 空间模块识别 (基于非负表达层的 NMF)
        print(f"[{sample_id}] Running Spatial NMF on Gene Expression...")
        try:
            X = get_nmf_input_matrix(adata_sub)
            n_components = 8
            model = NMF(n_components=n_components, init='nndsvda', random_state=42, max_iter=500)
            W = model.fit_transform(X) # (n_cells, n_components)
            H = model.components_      # (n_components, n_genes)
            save_nmf_artifacts(
                output_dir,
                sample_id,
                adata_sub.obs_names,
                adata_sub.var_names,
                W,
                H,
                source_layer=DEFAULT_EXPRESSION_LAYER,
            )

            # 将因子权重存入 obs
            for i in range(n_components):
                adata_sub.obs[f'NMF_Factor_{i}'] = W[:, i]

            # 绘制 NMF 因子在空间上的分布
            print(f"[{sample_id}] Plotting NMF Factors...")
            
            # 使用更小的 size 和没有边框的散点图，防止高分辨率下单细胞点重叠严重
            # cmap 使用 'magma' 或者 'rocket' 以突显高表达区域
            sc.pl.embedding(adata_sub, basis="spatial", color=[f'NMF_Factor_{i}' for i in range(n_components)], 
                            show=False, cmap="magma", ncols=4, size=5, frameon=False, wspace=0.1)
            save_plot(os.path.join(output_dir, f"TSC_{sample_id}_NMF_Factors_Spatial.pdf"))
            plt.close()
            
            # 提取每个因子的 Top 基因
            feature_names = adata_sub.var_names
            top_genes_dict = {}
            for topic_idx, topic in enumerate(H):
                top_features_ind = topic.argsort()[:-11:-1]
                top_features = feature_names[top_features_ind]
                top_genes_dict[f'NMF_Factor_{topic_idx}'] = top_features.tolist()
                print(f"[{sample_id}] NMF Factor {topic_idx}: {', '.join(top_features)}")

            pd.DataFrame(dict((k, pd.Series(v)) for k, v in top_genes_dict.items())).to_csv(
                os.path.join(output_dir, f"TSC_{sample_id}_NMF_TopGenes.csv"),
                index=False
            )
        except Exception as e:
            print(f"[{sample_id}] Error in NMF: {e}")

        # 释放内存
        del adata_sub
        import gc
        gc.collect()
        
    print("\n--- Processing Combined Samples ---", flush=True)
    try:
        # 为了避免不同样本的空间坐标重叠导致错误计算，我们需要为不同样本的坐标添加偏移量
        # 使得它们在同一个坐标系下但是彼此远离，或者使用 squidpy 的 library_key 参数
        # squidpy.gr.spatial_neighbors 支持 library_key 参数，它会在不同的 library (这里即 sample) 内部独立构建邻居图
        # 这样就可以合并分析而不会产生跨样本的错误边
        print("Computing Spatial Neighbors across all samples (separated by library_key)...", flush=True)
        compute_samplewise_spatial_neighbors(adata, spatial_key="spatial", sample_key="sample")
        
        print("Computing Neighborhood Enrichment for all combined samples...", flush=True)
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, n_perms=1000)
        
        plt.figure(figsize=(10, 8))
        sq.pl.nhood_enrichment(adata, cluster_key=cluster_key, method="ward", show=False, cmap="vlag", vmin=-1, vmax=1, annot=True, fmt=".2f")
        plt.title(f"Spatial Neighborhood Enrichment (Combined All Samples)", fontsize=16, pad=20)
        save_plot(os.path.join(output_dir, f"TSC_Combined_AllSamples_nhood_enrichment.pdf"))
        plt.close()
        print("Successfully generated combined neighborhood enrichment plot.", flush=True)
    except Exception as e:
        print(f"Error in Combined Neighborhood Enrichment: {e}", flush=True)

    print("\n--- Computing Combined Co-occurrence probability ---", flush=True)
    try:
        # 在合并的所有样本上计算共现性
        sq.gr.co_occurrence(adata, cluster_key=cluster_key)
        
        # 绘制所有细胞类型的共现
        all_clusters_combined = adata.obs[cluster_key].unique().tolist()
        
        n_clusters_comb = len(all_clusters_combined)
        fig_width_comb = min(20, max(12, n_clusters_comb * 2))
        
        plt.figure(figsize=(fig_width_comb, 6))
        sq.pl.co_occurrence(adata, cluster_key=cluster_key, clusters=all_clusters_combined)
        save_plot(os.path.join(output_dir, f"TSC_Combined_AllSamples_co_occurrence.pdf"))
        plt.close()
        print("Successfully generated combined co-occurrence plot.", flush=True)
    except Exception as e:
        print(f"Error in Combined Co-occurrence: {e}", flush=True)

    print("\n--- Running Combined NMF using the same contract-defined matrix ---", flush=True)
    try:
        X = get_nmf_input_matrix(adata)
        n_components = 8
        model = NMF(n_components=n_components, init='nndsvda', random_state=42, max_iter=500)
        W = model.fit_transform(X)
        H = model.components_
        save_nmf_artifacts(
            output_dir,
            "Combined_AllSamples",
            adata.obs_names,
            adata.var_names,
            W,
            H,
            source_layer=DEFAULT_EXPRESSION_LAYER,
        )
        print("Saved combined NMF artifacts for downstream reuse.", flush=True)
    except Exception as e:
        print(f"Error in Combined NMF: {e}", flush=True)

    print("\nNiche Analysis for all samples completed.")

def main():
    print("Step 07: Spatial Niche Analysis...")
    output_dir = PIPELINE_FILES["nmf_dir"]

    # --- 1. TSC 分析 ---
    tsc_path = PIPELINE_FILES["annotated_h5ad"]
    if os.path.exists(tsc_path):
        print(f"Loading TSC data from {tsc_path}...")
        adata_tsc = sc.read_h5ad(tsc_path)
        analyze_tsc_niche(adata_tsc, output_dir)
        del adata_tsc
    else:
        print("TSC data not found. Skipping.")

    print("\nStep 07 Done.")

if __name__ == "__main__":
    main()


