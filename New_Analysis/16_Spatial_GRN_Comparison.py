# ==============================================================================
# 文件名称: 16_Spatial_GRN_Comparison.py
# 功能描述: 病理与正常微环境的空间 GRN 对比
# 联动说明: 承接 15 步，专门对比 GC/DN 所在病理生态位与正常生态位在基因调控网络上的差异，找出致痫枢纽基因 (如 TNR)。
# ==============================================================================

import os
import sys
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.spatial import cKDTree

# Add current directory to path to import config
sys.path.append(os.getcwd())
from analysis_config import DIRS, PIPELINE_FILES
from analysis_contract import (
    DEFAULT_EXPRESSION_LAYER,
    compute_samplewise_spatial_neighbors,
    get_expression_matrix,
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


# 绘图设置，保证 PDF 中的文字可以被编辑
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

try:
    from analysis_config import set_nature_style, NPG_COLORS
    set_nature_style()
except ImportError:
    pass

def get_hub_genes(adata_subset, niche_name, output_dir, cell_col='celltype'):
    """
    Compute spatial autocorrelation and build GRN for a specific niche subset.
    Includes GC/DN class balancing correction to ensure rare cells (like GC) 
    are not overshadowed by abundant cells (like DN).
    """
    print(f"\n--- Processing {niche_name} ---")
    
    if adata_subset.n_obs < 50:
        print(f"Warning: Too few spots ({adata_subset.n_obs}) in {niche_name}. Skipping.")
        return None, None
        
    print(f"Subsampled {adata_subset.n_obs} spots.")

    # 1. Spatial Autocorrelation
    print("Computing Spatial Autocorrelation (Moran's I)...")
    compute_samplewise_spatial_neighbors(adata_subset, spatial_key="spatial", sample_key="sample")
    
    # Use top 2000 highly variable genes for speed
    if 'highly_variable' in adata_subset.var:
        base_genes = adata_subset.var_names[adata_subset.var['highly_variable']][:2000].tolist()
    else:
        base_genes = adata_subset.var_names[:2000].tolist()
        
    genes_to_test = list(base_genes)
    
    # --- GC/DN CORRECTION: Feature Selection ---
    # Force inclusion of specific markers to prevent rare signals from being overshadowed
    if niche_name in ['DN_Niche', 'GC_Niche'] and cell_col in adata_subset.obs:
        print(f"Applying correction: Forcing {niche_name} markers into spatial test...")
        gc_mask = adata_subset.obs[cell_col].str.contains('Giant', na=False)
        dn_mask = adata_subset.obs[cell_col].str.contains('Dysmorphic', na=False)
        
        tmp_group = pd.Series('Other', index=adata_subset.obs_names)
        tmp_group[gc_mask] = 'GC'
        tmp_group[dn_mask] = 'DN'
        adata_subset.obs['tmp_group'] = tmp_group.astype('category')
        
        try:
            # Use layer='log1p' if available to avoid raw counts warning
            layer = 'log1p' if 'log1p' in adata_subset.layers else None
            
            target_group = 'GC' if niche_name == 'GC_Niche' else 'DN'
            if target_group in adata_subset.obs['tmp_group'].values:
                sc.tl.rank_genes_groups(adata_subset, groupby='tmp_group', groups=[target_group], reference='Other', method='wilcoxon', layer=layer)
                if target_group in adata_subset.uns['rank_genes_groups']['names'].dtype.names:
                    markers = pd.DataFrame(adata_subset.uns['rank_genes_groups']['names'])[target_group].head(200).tolist()
                    markers = [g for g in markers if g in adata_subset.var_names]
                    genes_to_test.extend(markers)
                    genes_to_test = list(set(genes_to_test))
                    print(f"Added specific markers for {target_group}. Total genes to test: {len(genes_to_test)}")
        except Exception as e:
            print(f"Warning: Could not compute markers for correction: {e}")

    # Using n_jobs=1 to avoid PicklingError and memory limit (WinError 1450) on Windows
    # Multiprocessing with large AnnData objects often fails when serializing data to child processes
    sq.gr.spatial_autocorr(
        adata_subset,
        mode="moran",
        genes=genes_to_test,
        n_perms=100,
        n_jobs=1
    )
    
    moran_df = adata_subset.uns['moranI']
    top_spatial_genes = moran_df[moran_df['pval_norm'] < 0.05].sort_values('I', ascending=False).head(200).index.tolist()
    print(f"Identified {len(top_spatial_genes)} significant spatially variable genes.")

    if not top_spatial_genes:
        print("No significant spatial genes found.")
        return None, None

    # 2. Build Co-expression Network
    print("Building Gene Co-expression Network...")
    expr_mat = get_expression_matrix(adata_subset, layer=DEFAULT_EXPRESSION_LAYER, nonnegative=False)
    expr_mat = expr_mat[:, [adata_subset.var_names.get_loc(g) for g in top_spatial_genes]]

    expr_df = pd.DataFrame(expr_mat, columns=top_spatial_genes)
    
    # --- GC/DN CORRECTION: Co-expression Balancing ---
    # To prevent abundant cells (DN) from washing out rare cells (GC), 
    # we take the max absolute correlation across target, and the global niche.
    if niche_name in ['DN_Niche', 'GC_Niche'] and cell_col in adata_subset.obs:
        print(f"Applying correction: Calculating balanced co-expression network for {niche_name}...")
        gc_mask = adata_subset.obs[cell_col].str.contains('Giant', na=False)
        dn_mask = adata_subset.obs[cell_col].str.contains('Dysmorphic', na=False)
        other_mask = ~(gc_mask | dn_mask)
        
        target_mask = gc_mask if niche_name == 'GC_Niche' else dn_mask
        
        global_corr = expr_df.corr(method='spearman').fillna(0)
        max_corr = global_corr.copy()
        max_abs_corr = global_corr.abs()
        
        for mask, name in zip([target_mask, other_mask], ['Target', 'Other']):
            if mask.sum() > 10: # Only compute if enough cells
                sub_corr = expr_df[mask.values].corr(method='spearman').fillna(0)
                is_larger = sub_corr.abs() > max_abs_corr
                max_abs_corr = max_abs_corr.where(~is_larger, sub_corr.abs())
                max_corr = max_corr.where(~is_larger, sub_corr)
        
        corr_matrix = max_corr
    else:
        corr_matrix = expr_df.corr(method='spearman').fillna(0)
    
    # Threshold
    threshold = 0.4
    edges = np.where(np.abs(corr_matrix) > threshold)
    
    G = nx.Graph()
    for i in range(len(edges[0])):
        u = corr_matrix.index[edges[0][i]]
        v = corr_matrix.columns[edges[1][i]]
        weight = corr_matrix.iloc[edges[0][i], edges[1][i]]
        if u != v:
            G.add_edge(u, v, weight=weight)
            
    if G.number_of_nodes() == 0:
        print("Network is empty after thresholding.")
        return None, None
        
    # 3. Network Centrality
    print("Analyzing Network Centrality...")
    degree_dict = dict(G.degree(weight='weight'))
    centrality = nx.eigenvector_centrality(G, max_iter=1000, weight='weight', tol=1e-04)
    
    hub_df = pd.DataFrame({
        'Gene': list(centrality.keys()),
        'Eigenvector_Centrality': list(centrality.values()),
        'Weighted_Degree': [degree_dict[g] for g in centrality.keys()]
    }).sort_values('Eigenvector_Centrality', ascending=False)
    
    hub_df.to_csv(os.path.join(output_dir, f"GRN_Hub_Genes_{niche_name.replace(' ', '_')}.csv"), index=False)
    
    # Save the correlation matrix so downstream scripts can reconstruct the network
    corr_matrix.to_csv(os.path.join(output_dir, f"{niche_name.replace(' ', '_')}_Corr_Matrix.csv"))
    
    # Plot Network
    top_nodes = hub_df['Gene'].head(50).tolist()
    top_5_hubs = hub_df['Gene'].head(5).tolist() # 前5大Hub基因
    sub_G = G.subgraph(top_nodes).copy()
    
    # 稀疏化网络连边：只保留绝对相关性较高的边 (例如 > 0.6)
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
    
    # 根据网络类别和 Hub 级别设置颜色 (仅对有效节点)
    node_sizes = [centrality[n] * 5000 for n in valid_nodes]
    base_color = '#E64B35' if 'Pathological' in niche_name else '#4DBBD5'
    node_colors = [base_color if n in top_5_hubs else 'lightgray' for n in valid_nodes]
    
    nx.draw_networkx_nodes(sub_G_sparse, pos, node_size=node_sizes, node_color=node_colors, 
                           edgecolors='black', linewidths=1.0, alpha=0.9)
                           
    edges = sub_G_sparse.edges()
    if edges:
        weights = [abs(sub_G_sparse[u][v]['weight']) * 2 for u, v in edges]
        nx.draw_networkx_edges(sub_G_sparse, pos, edgelist=edges, width=weights, edge_color='gray', alpha=0.4)
        
    nx.draw_networkx_labels(sub_G_sparse, pos, font_size=9, font_weight='bold', font_family='Arial', font_color='black')
    
    plt.title(f"Spatial GRN (Top 50 Hubs) - {niche_name}\nHighlighted Top 5 Hubs", fontweight='bold', fontsize=14)
    plt.axis('off')
    plt.savefig(os.path.join(output_dir, f"Network_{niche_name.replace(' ', '_')}.pdf"), bbox_inches='tight')
    plt.close()
    
    return hub_df, top_spatial_genes

def main():
    print("==================================================")
    print("06 Module Add-on: Spatial GRN Pathological vs Normal Comparison")
    print("==================================================")
    
    input_h5ad = PIPELINE_FILES["annotated_h5ad"]
    output_dir = PIPELINE_FILES["grn_comparison_dir"]
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading data from: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    require_data_contract(adata, "16_Spatial_GRN_Comparison.py")
    
    cell_col = 'celltype' if 'celltype' in adata.obs.columns else 'cell_type'
    
    # Define cell categories
    dn_types = [c for c in adata.obs[cell_col].unique() if 'Dysmorphic' in str(c)]
    gc_types = [c for c in adata.obs[cell_col].unique() if 'Giant' in str(c)]
    patho_types = dn_types + gc_types
    neuron_types = [c for c in adata.obs[cell_col].unique() if 'Neuron' in str(c) and 'Dysmorphic' not in str(c)]
    
    dn_spots = adata.obs[cell_col].isin(dn_types)
    gc_spots = adata.obs[cell_col].isin(gc_types)
    neuron_spots = adata.obs[cell_col].isin(neuron_types)
    
    print("Computing physical neighborhoods...")
    coords = adata.obsm['spatial']
    
    # Radii configuration based on script 12
    dn_radius = 100.0
    gc_radius = 100.0  # Changed to 100um as requested
    neuron_radius = 100.0 
    
    all_dn_niche_indices = []
    all_gc_niche_indices = []
    all_normal_niche_indices = []
    
    samples = adata.obs['sample'].unique()
    for sample in samples:
        sample_mask = adata.obs['sample'] == sample
        sample_indices = np.where(sample_mask)[0]
        sample_coords = coords[sample_indices]
        tree = cKDTree(sample_coords)
        
        # 1. Pathological Niches
        sample_dn_spots = dn_spots.iloc[sample_indices]
        sample_gc_spots = gc_spots.iloc[sample_indices]
        
        target_coords_dn = sample_coords[sample_dn_spots]
        target_coords_gc = sample_coords[sample_gc_spots]
        
        if len(target_coords_dn) > 0:
            neighbor_lists_dn = tree.query_ball_point(target_coords_dn, r=dn_radius)
            local_dn_indices = np.unique(np.concatenate(neighbor_lists_dn)).astype(int)
            global_dn_indices = sample_indices[local_dn_indices]
            all_dn_niche_indices.extend(global_dn_indices)
            
        if len(target_coords_gc) > 0:
            neighbor_lists_gc = tree.query_ball_point(target_coords_gc, r=gc_radius)
            local_gc_indices = np.unique(np.concatenate(neighbor_lists_gc)).astype(int)
            global_gc_indices = sample_indices[local_gc_indices]
            all_gc_niche_indices.extend(global_gc_indices)
            
        # 2. Normal Niche (Neurons, but exclude microenvironments with > 5% DN/GC cells)
        sample_neuron_spots = neuron_spots.iloc[sample_indices]
        target_coords_neuron = sample_coords[sample_neuron_spots]
        
        if len(target_coords_neuron) > 0:
            neighbor_lists_n = tree.query_ball_point(target_coords_neuron, r=neuron_radius)
            
            valid_normal_niche_local = []
            valid_neuron_anchors_count = 0
            
            for i, neighbors in enumerate(neighbor_lists_n):
                n_total = len(neighbors)
                if n_total == 0: continue
                
                # Check proportion of DN/GC in this neighborhood
                is_dn = sample_dn_spots.iloc[neighbors].values
                is_gc = sample_gc_spots.iloc[neighbors].values
                n_patho = np.sum(is_dn | is_gc)
                
                if n_patho / n_total <= 0.05:
                    valid_normal_niche_local.extend(neighbors)
                    valid_neuron_anchors_count += 1
                    
            if valid_normal_niche_local:
                local_n_indices = np.unique(valid_normal_niche_local).astype(int)
                global_n_indices = sample_indices[local_n_indices]
                all_normal_niche_indices.extend(global_n_indices)
            
            print(f"  Sample {sample}: Retained {valid_neuron_anchors_count} out of {len(target_coords_neuron)} normal neuron anchors (<= 5% DN/GC).")

    # Convert to sets for logic
    set_dn_niche = set(all_dn_niche_indices)
    set_gc_niche = set(all_gc_niche_indices)
    set_normal_niche = set(all_normal_niche_indices)
    
    dn_mask = np.zeros(adata.n_obs, dtype=bool)
    dn_mask[list(set_dn_niche)] = True
    
    gc_mask = np.zeros(adata.n_obs, dtype=bool)
    gc_mask[list(set_gc_niche)] = True
    
    normal_mask = np.zeros(adata.n_obs, dtype=bool)
    normal_mask[list(set_normal_niche)] = True
    
    # We will create three separate adatas
    adata_dn = adata[dn_mask].copy()
    adata_dn.obs['Niche_Type'] = 'DN_Niche'
    
    adata_gc = adata[gc_mask].copy()
    adata_gc.obs['Niche_Type'] = 'GC_Niche'
    
    adata_normal = adata[normal_mask].copy()
    adata_normal.obs['Niche_Type'] = 'Normal_Neuron_Niche'
    
    print(f"\n--- Final Niche Sizes ---")
    print(f"DN Niche size: {adata_dn.n_obs} spots")
    print(f"GC Niche size: {adata_gc.n_obs} spots")
    print(f"Normal Neuron Niche size (before imputation): {adata_normal.n_obs} spots")
    
    if adata_normal.n_obs < 50:
        print("WARNING: The Normal Neuron Niche has fewer than 50 spots. This is generally insufficient for robust GRN construction.")
    else:
        print("Normal Neuron Niche size is sufficient for GRN construction.")
    
    # --- Spatial KNN Imputation for Normal Niche ---
    print("\n--- Imputing DN/GC cells in the Normal Niche ---")
    import scipy.sparse as sp
    
    is_patho_in_normal = adata_normal.obs[cell_col].isin(patho_types)
    is_neuron_in_normal = adata_normal.obs[cell_col].isin(neuron_types)
    
    if is_patho_in_normal.sum() > 0:
        print(f"Found {is_patho_in_normal.sum()} DN/GC cells in the Normal Niche. Replacing them using Spatial KNN (k=5)...")
        
        X_is_sparse = sp.issparse(adata_normal.X)
        X_dense = adata_normal.X.toarray() if X_is_sparse else adata_normal.X.copy()
        
        has_raw = adata_normal.raw is not None
        if has_raw:
            raw_X_is_sparse = sp.issparse(adata_normal.raw.X)
            raw_X_dense = adata_normal.raw.X.toarray() if raw_X_is_sparse else adata_normal.raw.X.copy()
            
        for sample in samples:
            sample_mask_norm = adata_normal.obs['sample'] == sample
            sample_patho_idx = np.where(is_patho_in_normal & sample_mask_norm)[0]
            sample_neuron_idx = np.where(is_neuron_in_normal & sample_mask_norm)[0]
            
            if len(sample_patho_idx) == 0:
                continue
                
            if len(sample_neuron_idx) == 0:
                print(f"Warning: Sample {sample} has DN/GC cells but no normal neurons for imputation.")
                continue
                
            norm_coords = adata_normal.obsm['spatial']
            neuron_coords = norm_coords[sample_neuron_idx]
            patho_coords = norm_coords[sample_patho_idx]
            
            tree_imp = cKDTree(neuron_coords)
            actual_k = min(5, len(sample_neuron_idx))
            distances, indices = tree_imp.query(patho_coords, k=actual_k)
            
            if actual_k == 1:
                indices = indices.reshape(-1, 1)
                
            for i, p_idx in enumerate(sample_patho_idx):
                n_idx_global = sample_neuron_idx[indices[i]]
                X_dense[p_idx, :] = np.mean(X_dense[n_idx_global, :], axis=0)
                if has_raw:
                    raw_X_dense[p_idx, :] = np.mean(raw_X_dense[n_idx_global, :], axis=0)
                    
        adata_normal.X = sp.csr_matrix(X_dense) if X_is_sparse else X_dense
        if has_raw:
            new_raw_adata = sc.AnnData(X=sp.csr_matrix(raw_X_dense) if raw_X_is_sparse else raw_X_dense,
                                       obs=adata_normal.obs, var=adata_normal.raw.var)
            adata_normal.raw = new_raw_adata
            
        if adata_normal.obs[cell_col].dtype.name == 'category':
            if 'Imputed_Neuron' not in adata_normal.obs[cell_col].cat.categories:
                adata_normal.obs[cell_col] = adata_normal.obs[cell_col].cat.add_categories(['Imputed_Neuron'])
        adata_normal.obs.loc[is_patho_in_normal, cell_col] = 'Imputed_Neuron'
        print("Imputation completed.")
    else:
        print("No DN/GC cells found in the Normal Niche. No imputation needed.")
        
    # Run GRN
    hub_dn, _ = get_hub_genes(adata_dn, 'DN_Niche', output_dir, cell_col=cell_col)
    hub_gc, _ = get_hub_genes(adata_gc, 'GC_Niche', output_dir, cell_col=cell_col)
    hub_n, _ = get_hub_genes(adata_normal, 'Normal_Neuron_Niche', output_dir, cell_col=cell_col)
    
    # Comparison
    if hub_dn is not None and hub_gc is not None and hub_n is not None:
        print("\n--- Comparing Hub Genes ---")
        top_dn = set(hub_dn['Gene'].head(50))
        top_gc = set(hub_gc['Gene'].head(50))
        top_n = set(hub_n['Gene'].head(50))
        
        shared_all = top_dn.intersection(top_gc).intersection(top_n)
        shared_patho = top_dn.intersection(top_gc) - top_n
        unique_dn = top_dn - top_gc - top_n
        unique_gc = top_gc - top_dn - top_n
        unique_n = top_n - top_dn - top_gc
        
        # Venn Diagram (3-way)
        try:
            from matplotlib_venn import venn3
            plt.figure(figsize=(8, 8))
            venn3([top_dn, top_gc, top_n], set_labels=('DN Niche Hubs', 'GC Niche Hubs', 'Normal Neuron Niche Hubs'))
            plt.title("Overlap of Top 50 GRN Hub Genes")
            plt.savefig(os.path.join(output_dir, "Hub_Genes_Venn3.pdf"), bbox_inches='tight')
            plt.close()
        except ImportError:
            print("matplotlib_venn not installed. Skipping Venn diagram.")
            
        with open(os.path.join(output_dir, "Hub_Genes_Comparison.txt"), "w") as f:
            f.write(f"Shared Across All Niches ({len(shared_all)}):\n{', '.join(shared_all)}\n\n")
            f.write(f"Shared between GC and DN only ({len(shared_patho)}):\n{', '.join(shared_patho)}\n\n")
            f.write(f"Unique to DN Niche ({len(unique_dn)}):\n{', '.join(unique_dn)}\n\n")
            f.write(f"Unique to GC Niche ({len(unique_gc)}):\n{', '.join(unique_gc)}\n\n")
            f.write(f"Unique to Normal Niche ({len(unique_n)}):\n{', '.join(unique_n)}\n")
            
        print("Comparison saved to Hub_Genes_Comparison.txt")

    # 4. Map Top Hub Genes back to Spatial Coordinates for each sample
    print("\n--- Mapping Top Hub Genes to Spatial Space ---")
    
    def plot_spatial_hubs(adata_full, hub_df, niche_name, top_n=12):
        if hub_df is None or hub_df.empty:
            return
            
        top_genes = hub_df['Gene'].head(top_n).tolist()
        
        # We need to make sure the genes are actually in the adata.var
        valid_genes = [g for g in top_genes if g in adata_full.var_names]
        if not valid_genes:
            print(f"None of the top genes for {niche_name} are in the main adata.")
            return
            
        samples = adata_full.obs['sample'].unique()
        for sample in samples:
            print(f"Plotting spatial hubs for {niche_name} - Sample {sample}...")
            sample_adata = adata_full[adata_full.obs['sample'] == sample].copy()
            
            if 'spatial' not in sample_adata.obsm:
                print(f"Warning: No spatial coordinates for sample {sample}")
                continue
                
            try:
                fig, ax = plt.subplots(figsize=(15, 15))
                # Plot using sc.pl.embedding with basis='spatial'
                sc.pl.embedding(
                    sample_adata, 
                    basis="spatial", 
                    color=valid_genes, 
                    cmap='Reds', 
                    size=20, 
                    show=False,
                    ncols=4 # Arrange in a grid
                )
                plt.suptitle(f"Top {top_n} Hub Genes: {niche_name} (Sample {sample})", y=1.02, fontsize=16)
                
                # Ensure equal aspect ratio
                for a in plt.gcf().axes:
                    if a.get_label() != '<colorbar>':
                        a.set_aspect('equal')
                        
                safe_niche = niche_name.replace(' ', '_')
                plt.savefig(os.path.join(output_dir, f"Spatial_Hubs_{safe_niche}_{sample}.png"), dpi=300, bbox_inches='tight')
                plt.savefig(os.path.join(output_dir, f"Spatial_Hubs_{safe_niche}_{sample}.pdf"), bbox_inches='tight')
                plt.close('all')
            except Exception as e:
                print(f"Error plotting spatial hubs for {sample}: {e}")

    # Plot for DN Niche
    if hub_dn is not None:
        plot_spatial_hubs(adata, hub_dn, 'DN_Niche', top_n=12)
        
    # Plot for GC Niche
    if hub_gc is not None:
        plot_spatial_hubs(adata, hub_gc, 'GC_Niche', top_n=12)
        
    # Plot for Normal Niche
    if hub_n is not None:
        plot_spatial_hubs(adata, hub_n, 'Normal_Neuron_Niche', top_n=12)

    print("Spatial GRN Comparison completed successfully!")

if __name__ == "__main__":
    main()
