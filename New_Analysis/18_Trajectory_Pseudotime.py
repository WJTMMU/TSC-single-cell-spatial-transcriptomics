# ==============================================================================
# 文件名称: 18_Trajectory_Pseudotime.py
# 功能描述: 拟时序与空间轨迹推断
# 联动说明: 利用 PAGA 和 Diffusion Map 等算法，以正常神经元/OPC为起点，推断 GC 和 DN 细胞的发育/去分化演进轨迹。
# ==============================================================================

import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from analysis_config import DIRS, TSC_SAMPLES, PIPELINE_FILES

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
    print("08 Module: Spatial Trajectory Inference (PAGA & DPT)")
    print("==================================================")
    
    # 1. Load Data
    input_h5ad = PIPELINE_FILES["annotated_h5ad"]
    print(f"Loading preprocessed data from: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    
    output_dir = PIPELINE_FILES["trajectory_dir"]
    print(f"Output directory: {output_dir}")
    
    # Ensure standard cell type annotations are present
    if 'cell_type' in adata.obs.columns:
        cell_col = 'cell_type'
    elif 'celltype' in adata.obs.columns:
        cell_col = 'celltype'
    else:
        print(f"Error: Cell type column not found in {adata.obs.columns}.")
        return
    
    # Define root cells for trajectory (e.g., Progenitors, OPCs, or Normal Neurons)
    # We want to trace the path to Giant Cells and Dysmorphic Neurons
    print("Setting root cell for trajectory...")
    
    # 2. PAGA Trajectory Graph
    # We need to compute neighborhood graph and clustering first if not present
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    
    print("Computing PAGA graph...")
    sc.tl.paga(adata, groups=cell_col)
    
    # Plot PAGA graph
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.paga(adata, plot=False, ax=ax, threshold=0.05, 
               node_size_scale=1.5, edge_width_scale=1.0)
    ax.set_title("PAGA Trajectory Graph across Cell Types")
    plt.savefig(os.path.join(output_dir, "PAGA_Trajectory_Graph.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "PAGA_Trajectory_Graph.pdf"), bbox_inches='tight')
    plt.close()
    
    # 3. Diffusion Pseudotime (DPT)
    print("Computing Diffusion Maps and Pseudotime...")
    sc.tl.diffmap(adata)
    
    # Try to find a sensible root: e.g., 'OPC' or 'Excitatory_Neurons' or generic progenitor
    root_candidates = ['Excitatory Neurons', 'Excitatory_Neurons', 'OPC', 'Astrocytes']
    root_idx = None
    for cand in root_candidates:
        if cand in adata.obs[cell_col].cat.categories:
            root_idx = np.flatnonzero(adata.obs[cell_col] == cand)[0]
            print(f"Selected '{cand}' as the root cell type for pseudotime.")
            break
    
    if root_idx is None:
        print("Warning: Root cell type not found. Using random root.")
        root_idx = 0
        
    adata.uns['iroot'] = root_idx
    
    # Calculate DPT
    sc.tl.dpt(adata)
    
    # Optional: Highlight specific transition path (Excitatory Neurons -> Giant Cells)
    print("Highlighting Excitatory Neurons to Giant Cells transition...")
    
    # Create a focused subset for clearer visualization if both types exist
    target_types = ['Excitatory Neurons', 'Giant Cells', 'Dysmorphic Neurons(Excit)']
    available_targets = [t for t in target_types if t in adata.obs[cell_col].cat.categories]
    
    if len(available_targets) >= 2:
        fig, ax = plt.subplots(figsize=(8, 6))
        # Plot UMAP highlighting only the target transition cells
        sc.pl.umap(adata, color=cell_col, groups=available_targets, 
                   palette='Set1', show=False, ax=ax, title="Path: Excitatory Neurons -> Giant Cells/DN")
        plt.savefig(os.path.join(output_dir, "UMAP_Transition_Path_Highlight.pdf"), bbox_inches='tight')
        plt.close()
    
    # Plot Pseudotime on UMAP
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='dpt_pseudotime', cmap='viridis', show=False, ax=ax)
    ax.set_title("Developmental Pseudotime")
    plt.savefig(os.path.join(output_dir, "UMAP_Pseudotime.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "UMAP_Pseudotime.pdf"), bbox_inches='tight')
    plt.close()
    
    # Plot Pseudotime by Cell Type (Violin Plot)
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # order by pseudotime mean to show the progression from root to end
    type_dpt_mean = adata.obs.groupby(cell_col)['dpt_pseudotime'].mean().sort_values().index
    sns.violinplot(x=cell_col, y='dpt_pseudotime', data=adata.obs, inner='quartile', ax=ax, palette='Set2', order=type_dpt_mean)
    
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_title("Pseudotime Distribution by Cell Type (Ordered by Progression)")
    ax.set_ylabel("Diffusion Pseudotime")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Pseudotime_by_CellType_Violin.png"), dpi=300)
    plt.savefig(os.path.join(output_dir, "Pseudotime_by_CellType_Violin.pdf"))
    plt.close()
    
    # 4. Spatial Pseudotime Visualization (For each sample)
    print("Generating spatial pseudotime maps...")
    for sample in adata.obs['sample'].unique():
        print(f"Plotting spatial pseudotime for sample {sample}...")
        sample_adata = adata[adata.obs['sample'] == sample].copy()
        
        if 'spatial' not in sample_adata.obsm:
            print(f"Warning: No spatial coordinates found for sample {sample}")
            continue
            
        try:
            fig, ax = plt.subplots(figsize=(8, 8))
            # Use embedding with basis='spatial' instead of spatial to avoid Visium-specific requirements
            sc.pl.embedding(sample_adata, basis="spatial", color='dpt_pseudotime', cmap='viridis', 
                          size=20, alpha=0.8, show=False, ax=ax)
            ax.set_title(f"Spatial Pseudotime: {sample}")
            
            # Plot preferences from core memories: equal aspect ratio, PDF+PNG
            ax.set_aspect('equal')
            
            plt.savefig(os.path.join(output_dir, f"Spatial_Pseudotime_{sample}.png"), dpi=300, bbox_inches='tight')
            plt.savefig(os.path.join(output_dir, f"Spatial_Pseudotime_{sample}.pdf"), bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"Error plotting spatial pseudotime for {sample}: {e}")
        
    # 5. Save Trajectory Data
    print("Saving trajectory results...")
    trajectory_df = adata.obs[['sample', cell_col, 'dpt_pseudotime']]
    trajectory_df.to_csv(os.path.join(output_dir, "Trajectory_Pseudotime_Results.csv"))
    
    print(f"Trajectory analysis completed successfully. Results saved to: {output_dir}")

if __name__ == "__main__":
    main()
