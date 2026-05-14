# 15_Spatial_GRN_Squidpy.py
# 06 Module: Spatial Gene Regulatory Network (GRN) Analysis
# Description: Uses Squidpy to build spatial gene co-expression networks
# specifically within pathological microenvironments (e.g., GC/DN niches).

import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

from analysis_config import DIRS, TSC_SAMPLES, PIPELINE_FILES
from analysis_contract import (
    DEFAULT_EXPRESSION_LAYER,
    compute_samplewise_spatial_neighbors,
    get_expression_matrix,
    query_radius_neighbors_by_sample,
    require_data_contract,
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


def main():
    print("==================================================")
    print("06 Module: Spatial Gene Regulatory Network (Squidpy)")
    print("==================================================")
    
    # 1. Load Data
    input_h5ad = PIPELINE_FILES["annotated_h5ad"]
    print(f"Loading preprocessed data from: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    require_data_contract(adata, "15_Spatial_GRN_Squidpy.py")
    
    output_dir = PIPELINE_FILES["grn_dir"]
    print(f"Output directory: {output_dir}")
    
    # 3. Define the Pathological Microenvironment (Niche)
    # We use the GC/DN neighborhood approach (Hypothesis-driven based on scripts 12/14)
    print("Defining Pathological Niche (GC/DN + Immediate Neighbors)...")
    
    # Identify spots that are GC or DN
    # Based on the actual data, the names contain spaces, not underscores
    target_cells = ['Giant Cells', 'Dysmorphic Neurons(Excit)', 'Dysmorphic Neurons(Inhib)']
    if 'cell_type' in adata.obs.columns:
        cell_col = 'cell_type'
    elif 'celltype' in adata.obs.columns:
        cell_col = 'celltype'
    else:
        print(f"Error: Could not find cell type column in {adata.obs.columns}")
        return
        
    target_spots = adata.obs[cell_col].isin(target_cells)
    
    print("Finding immediate physical neighbors of GC/DN...")
    # Query radius: 50 pixels (~50 um depending on resolution)
    # This defines the "immediate neighborhood"
    radius = 100.0 
    
    if int(target_spots.sum()) == 0:
        print("Warning: No GC or DN cells found in the dataset.")
        return

    all_niche_indices = query_radius_neighbors_by_sample(
        adata,
        target_spots.to_numpy(),
        radius=radius,
        spatial_key="spatial",
        sample_key="sample",
    )
    
    # Create mask
    niche_mask = np.zeros(adata.n_obs, dtype=bool)
    niche_mask[all_niche_indices] = True
    adata.obs['is_pathological_niche'] = niche_mask.astype(str)
    
    # 4. Subset to Niche
    adata_niche = adata[adata.obs['is_pathological_niche'] == 'True'].copy()
    print(f"Subsampled {adata_niche.n_obs} spots in the pathological niche out of {adata.n_obs} total.")
    
    if adata_niche.n_obs < 50:
        print("Warning: Too few spots in the defined niche. Exiting.")
        return

    # 5. Spatial Autocorrelation (Moran's I) within the Niche
    print("Computing Spatial Autocorrelation (Moran's I) within the niche...")
    # Recompute spatial neighbors for the subset
    compute_samplewise_spatial_neighbors(adata_niche, spatial_key="spatial", sample_key="sample")
    
    # Use top 2000 highly variable genes for speed, or all if preferred
    if 'highly_variable' in adata_niche.var:
        genes_to_test = adata_niche.var_names[adata_niche.var['highly_variable']][:2000]
    else:
        genes_to_test = adata_niche.var_names[:2000]
        
    sq.gr.spatial_autocorr(
        adata_niche,
        mode="moran",
        genes=genes_to_test,
        n_perms=100,
        n_jobs=-1
    )
    
    # Extract top spatially variable genes
    moran_df = adata_niche.uns['moranI']
    top_spatial_genes = moran_df[moran_df['pval_norm'] < 0.05].sort_values('I', ascending=False).head(200).index.tolist()
    
    print(f"Identified {len(top_spatial_genes)} significant spatially variable genes in the niche.")
    
    # 5. Build Co-expression Network (GRN)
    print("Building Gene Co-expression Network...")
    # Get expression matrix for top genes
    expr_mat = get_expression_matrix(adata_niche, layer=DEFAULT_EXPRESSION_LAYER, nonnegative=False)
    expr_mat = expr_mat[:, [adata_niche.var_names.get_loc(g) for g in top_spatial_genes]]

    expr_df = pd.DataFrame(expr_mat, columns=top_spatial_genes)
    
    # Compute correlation matrix (Spearman for robustness)
    corr_matrix = expr_df.corr(method='spearman')
    
    # Threshold to create edges
    threshold = 0.4
    edges = np.where(np.abs(corr_matrix) > threshold)
    
    G = nx.Graph()
    for i in range(len(edges[0])):
        u = corr_matrix.index[edges[0][i]]
        v = corr_matrix.columns[edges[1][i]]
        weight = corr_matrix.iloc[edges[0][i], edges[1][i]]
        if u != v: # no self-loops
            G.add_edge(u, v, weight=weight)
            
    # 6. Network Analysis (Hub Genes)
    print("Analyzing Network Centrality...")
    degree_dict = dict(G.degree(weight='weight'))
    centrality = nx.eigenvector_centrality(G, max_iter=1000, weight='weight', tol=1e-04)
    
    hub_df = pd.DataFrame({
        'Gene': list(centrality.keys()),
        'Eigenvector_Centrality': list(centrality.values()),
        'Weighted_Degree': [degree_dict[g] for g in centrality.keys()]
    }).sort_values('Eigenvector_Centrality', ascending=False)
    
    hub_df.to_csv(os.path.join(output_dir, "Pathological_Niche_Hub_Genes.csv"), index=False)
    print("Saved Hub Genes to:", os.path.join(output_dir, "Pathological_Niche_Hub_Genes.csv"))
    
    # 7. Visualization
    # 7.1 Plot Spatial Distribution of Top Hub Genes
    top_hubs = hub_df['Gene'].head(4).tolist()
    for sample in adata_niche.obs['sample'].unique():
        print(f"Plotting spatial hub genes for sample {sample}...")
        sample_adata = adata[adata.obs['sample'] == sample].copy()
        
        # Ensure we have spatial coordinates for plotting
        if 'spatial' not in sample_adata.obsm:
            print(f"Warning: No spatial coordinates found for sample {sample}")
            continue
            
        # Plot top 4 hub genes
        try:
            # We must use basis='spatial' for scatter plot if spatial plot fails
            fig, ax = plt.subplots(figsize=(10, 10))
            sc.pl.embedding(sample_adata, basis="spatial", color=top_hubs, cmap='Reds', size=20, show=False)
            plt.suptitle(f"Top Niche Hub Genes: {sample}")
            # Ensure equal aspect ratio per memory
            for a in plt.gcf().axes:
                a.set_aspect('equal')
            plt.savefig(os.path.join(output_dir, f"Spatial_Hub_Genes_{sample}.png"), dpi=300, bbox_inches='tight')
            plt.savefig(os.path.join(output_dir, f"Spatial_Hub_Genes_{sample}.pdf"), bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"Error plotting spatial hub genes for {sample}: {e}")
            
    # 7.2 Plot Network Graph (Top 50 nodes for clarity)
    top_nodes = hub_df['Gene'].head(50).tolist()
    top_5_hubs = hub_df['Gene'].head(5).tolist() # 前5大Hub基因
    sub_G = G.subgraph(top_nodes).copy()
    
    # 稀疏化网络连边：只保留权重最高的连边 (比如 absolute weight > 0.6)
    edges_to_keep = [(u, v) for u, v, d in sub_G.edges(data=True) if abs(d.get('weight', 0)) > 0.6]
    sub_G_sparse = sub_G.edge_subgraph(edges_to_keep).copy()
    
    # 移除孤立节点 (度数为0的节点)
    isolated_nodes = list(nx.isolates(sub_G_sparse))
    sub_G_sparse.remove_nodes_from(isolated_nodes)
    
    # 更新绘图的节点列表，防止原先的前50包含了被删除的节点
    valid_nodes = list(sub_G_sparse.nodes())
    
    plt.figure(figsize=(10, 10), dpi=600)
    # 使用 spring_layout，基于无孤立节点的稀疏网络
    pos = nx.spring_layout(sub_G_sparse, seed=42, k=0.8)
    
    # 节点大小和颜色 (仅对有效节点)
    node_sizes = [centrality[n] * 5000 for n in valid_nodes]
    node_colors = ['#E64B35' if n in top_5_hubs else '#4DBBD5' for n in valid_nodes]
    
    # 绘制节点 (带边框)
    nx.draw_networkx_nodes(sub_G_sparse, pos, node_size=node_sizes, node_color=node_colors, 
                           edgecolors='black', linewidths=1.0, alpha=0.9)
                           
    # 绘制连边 (根据权重调整粗细)
    edges = sub_G_sparse.edges()
    if edges:
        weights = [abs(sub_G_sparse[u][v]['weight']) * 2 for u, v in edges]
        nx.draw_networkx_edges(sub_G_sparse, pos, edgelist=edges, width=weights, edge_color='gray', alpha=0.4)
        
    # 绘制标签
    nx.draw_networkx_labels(sub_G_sparse, pos, font_size=9, font_weight='bold', font_family='Arial', font_color='black')
    
    plt.title("Spatial Gene Regulatory Network\n(Top 50 Hubs, Highlighted Top 5)", fontweight='bold', fontsize=14)
    plt.axis('off')
    plt.savefig(os.path.join(output_dir, "Spatial_GRN_Network.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "Spatial_GRN_Network.pdf"), bbox_inches='tight')
    plt.close()
    
    print("Spatial GRN analysis completed successfully!")

if __name__ == "__main__":
    main()
