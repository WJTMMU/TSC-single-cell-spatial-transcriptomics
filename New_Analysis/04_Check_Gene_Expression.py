# ==============================================================================
# 文件名称: 04_Check_Gene_Expression.py
# 功能描述: 核心基因空间表达可视化
# 联动说明: 承接 02 步带有注释的 h5ad 文件，在空间坐标和组织切片上验证特定标记基因或细胞类型的分布是否符合预期。
# ==============================================================================

import sys
import os
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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


# Add current directory and New_Analysis to path
sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "New_Analysis"))

try:
    from New_Analysis.analysis_config import *
except ImportError:
    try:
        from analysis_config import *
    except ImportError:
        # Fallback if analysis_config is not found
        print("Warning: analysis_config.py not found. Using default paths.")
        DIRS = {'basis': 'Results/Basis'}
        os.makedirs(DIRS['basis'], exist_ok=True)

# Ignore warnings
warnings.filterwarnings("ignore")

# Plotting settings
sc.settings.verbosity = 3
# sc.settings.set_figure_params(dpi=120, facecolor='white', vector_friendly=True)
plt.rcParams['figure.figsize'] = (8, 8)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

def main():
    print("=== Step 04: Check Gene Expression and Spatial Distribution ===")
    step_dir = get_step_dir("04")
    
    # 1. Load Data
    adata_path = PIPELINE_FILES["annotated_h5ad"]
    if not os.path.exists(adata_path):
        print(f"Error: Data file not found at {adata_path}")
        return

    print(f"Loading data from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    print(f"Data loaded: {adata.shape}")
    
    # Ensure output directory for plots exists
    plot_dir = step_dir
    os.makedirs(plot_dir, exist_ok=True)
    print(f"Plots will be saved to: {plot_dir}")

    # 2. Define Genes of Interest
    # Default list + typical markers
    genes_to_check = [
        'GFAP', 'RELN', 'SST', 'VIP', 'PVALB', 'GAD1', 'GAD2',  # Interneurons / Glia
        'SLC17A7', 'NEUROD2', 'SATB2', 'CUX2', 'RORB',          # Excitatory Neurons
        'MKI67', 'TOP2A',                                       # Proliferation
        'CD3D', 'CD3E', 'CD19', 'CD79A', 'CD14', 'CD68',        # Immune cells
        'CLDN5', 'PECAM1', 'ACTA2'                              # Endothelial / Mural
    ]
    
    # Check if we should use raw
    use_raw = False
    if adata.raw is not None:
        print("Using adata.raw for gene expression plotting (includes all genes).")
        use_raw = True
        available_genes = adata.raw.var_names
    else:
        print("adata.raw is None. Using adata.var_names (only HVGs).")
        available_genes = adata.var_names

    # Filter genes
    valid_genes = [g for g in genes_to_check if g in available_genes]
    missing_genes = [g for g in genes_to_check if g not in available_genes]
    
    print(f"\nValid genes found: {len(valid_genes)}")
    if missing_genes:
        print(f"Missing genes: {missing_genes}")

    # 3. Plot Gene Expression (DotPlot)
    if valid_genes and False:  # SKIP for now to avoid permission error and speed up
        print("\nGenerating DotPlot...")
        sc.pl.dotplot(adata, valid_genes, groupby='leiden', standard_scale='var', use_raw=use_raw, show=False)
        plt.savefig(os.path.join(plot_dir, "gene_expression_dotplot.pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(plot_dir, "gene_expression_dotplot.png"), bbox_inches='tight', dpi=300)
        plt.close()
        print("Saved dotplot (PDF & PNG).")
        
        # Plot Violin for top genes
        print("Generating Violin Plots...")
        sc.pl.violin(adata, valid_genes[:5], groupby='leiden', rotation=90, use_raw=use_raw, show=False)
        plt.savefig(os.path.join(plot_dir, "gene_expression_violin_top5.pdf"), bbox_inches='tight')
        plt.savefig(os.path.join(plot_dir, "gene_expression_violin_top5.png"), bbox_inches='tight', dpi=300)
        plt.close()
        print("Saved violin plots (PDF & PNG).")

    # 4. Spatial Coordinate Diagnosis
    print("\n=== Spatial Coordinates Diagnosis ===")
    has_spatial = False
    if 'spatial_x' in adata.obs and 'spatial_y' in adata.obs:
        has_spatial = True
        # Ensure obsm['spatial'] matches obs
        adata.obsm['spatial'] = adata.obs[['spatial_x', 'spatial_y']].values
        
        print("Preview of spatial coordinates (first 5 rows):")
        print(adata.obs[['spatial_x', 'spatial_y']].head())
        
        # Check correlation
        corr = adata.obs['spatial_x'].corr(adata.obs['spatial_y'])
        print(f"Correlation between X and Y: {corr:.4f}")
        
        if abs(corr) > 0.9:
            print("WARNING: High correlation detected! This indicates a diagonal line issue.")
            print("Possible causes: Data reading error, or columns are identical.")
    else:
        print("Error: 'spatial_x' or 'spatial_y' not found in obs.")

    # 5. Plot Spatial Clusters
    if has_spatial:
        samples = adata.obs['sample'].unique()
        
        if True: # Enable spatial clusters plotting
            print("\n=== Plotting Spatial Distribution ===")
            
            # Check if celltype exists, otherwise use leiden
            color_key = 'celltype' if 'celltype' in adata.obs else 'leiden'
            print(f"Using '{color_key}' for coloring.")

            for sample in samples:
                print(f"Plotting Sample: {sample}")
                subset = adata[adata.obs['sample'] == sample]
                
                # Create plot
                fig, ax = plt.subplots(figsize=(6, 6))
                
                # Plot embedding
                sc.pl.embedding(
                    subset, 
                    basis='spatial', 
                    color=color_key, 
                    title=f"Spatial Clusters ({color_key}) - {sample}",
                    frameon=False,
                    s=4,  # Smaller point size as requested
                    ax=ax,
                    show=False
                )
                
                # Force aspect ratio to be equal
                ax.set_aspect('equal', 'box')
                
                # Save
                out_file_png = os.path.join(plot_dir, f"spatial_clusters_{sample}.png")
                out_file_pdf = os.path.join(plot_dir, f"spatial_clusters_{sample}.pdf")
                
                plt.savefig(out_file_pdf, bbox_inches='tight')
                plt.savefig(out_file_png, bbox_inches='tight', dpi=300)
                plt.close()
                print(f"Saved plots to {out_file_pdf} and {out_file_png}")

        # --- New Section: Plot Individual Cell Types (Highlight One, Gray Others) ---
        if 'celltype' in adata.obs:
            print("\n=== Plotting Individual Cell Type Distributions ===")
            dist_plot_dir = os.path.join(plot_dir, "cell_type_distributions")
            os.makedirs(dist_plot_dir, exist_ok=True)
            
            # Define a consistent color palette for cell types
            cell_types = adata.obs['celltype'].cat.categories
            # Use scanpy's default palette if available, otherwise tab20
            if 'celltype_colors' in adata.uns and len(adata.uns['celltype_colors']) == len(cell_types):
                 palette = dict(zip(cell_types, adata.uns['celltype_colors']))
            else:
                 import matplotlib.colors as mcolors
                 cmap = plt.get_cmap('tab20')
                 palette = {ct: cmap(i/len(cell_types)) for i, ct in enumerate(cell_types)}

            for sample in samples:
                print(f"Processing cell type plots for {sample}...")
                subset = adata[adata.obs['sample'] == sample]
                
                # Get coordinates
                x_all = subset.obs['spatial_x'].values
                y_all = subset.obs['spatial_y'].values
                
                for ct in cell_types:
                    # Create figure
                    fig, ax = plt.subplots(figsize=(6, 6))
                    
                    # 1. Plot all cells as background (gray)
                    # Use lighter gray for better contrast
                    ax.scatter(x_all, y_all, c='#e0e0e0', s=4, label='Others', edgecolors='none')
                    
                    # 2. Plot specific cell type
                    mask = subset.obs['celltype'] == ct
                    if np.sum(mask) > 0:
                        ax.scatter(
                            x_all[mask], 
                            y_all[mask], 
                            c=[palette[ct]], 
                            s=4, # Small size as requested
                            label=ct, 
                            edgecolors='none'
                        )
                    
                    ax.set_title(f"{sample} - {ct}")
                    ax.set_aspect('equal', 'box')
                    ax.axis('off')
                    
                    # Save
                    safe_ct_name = str(ct).replace(" ", "_").replace("/", "-")
                    out_name = f"{sample}_{safe_ct_name}"
                    plt.savefig(os.path.join(dist_plot_dir, f"{out_name}.png"), bbox_inches='tight', dpi=300)
                    plt.savefig(os.path.join(dist_plot_dir, f"{out_name}.pdf"), bbox_inches='tight')
                    plt.close()
            print(f"Cell type distribution plots saved to {dist_plot_dir}")
        # ---------------------------------------------------------------------------

    # 6. Overlay on HE Image
        print("\n=== Overlaying on HE Images ===")
        import matplotlib.image as mpimg
        import glob

        for sample in samples:
            print(f"Processing overlay for {sample}...")
            
            # Find the image file
            # Construct path: TSC_DATA_DIR / sample / *aligned_HE_TIMG.png
            search_path = os.path.join(TSC_DATA_DIR, sample, "*_aligned_HE_TIMG.png")
            found_images = glob.glob(search_path)
            
            if not found_images:
                print(f"  Warning: No HE image found for {sample} at {search_path}")
                continue
                
            img_path = found_images[0]
            print(f"  Found image: {img_path}")
            
            # Load image
            try:
                img = mpimg.imread(img_path)
            except Exception as e:
                print(f"  Error reading image: {e}")
                continue
                
            # Get subset of data
            subset = adata[adata.obs['sample'] == sample]
            x = subset.obs['spatial_x'].values
            y = subset.obs['spatial_y'].values
            
            # Plot
            # Create a larger figure for high res overlay
            fig, ax = plt.subplots(figsize=(12, 12))
            
            # Show image
            # Use origin='lower' to flip the image vertically (row 0 at bottom)
            ax.imshow(img, origin='lower')
            
            # Scatter points
            # Map leiden clusters to colors
            groups = subset.obs['leiden'].unique()
            # Sort groups for consistent legend
            try:
                 groups = sorted(groups, key=int)
            except:
                 groups = sorted(groups)

            # Use tab20 colormap
            cmap = plt.get_cmap('tab20')
            colors = [cmap(i) for i in np.linspace(0, 1, len(groups))]
            
            for i, group in enumerate(groups):
                mask = subset.obs['leiden'] == group
                ax.scatter(
                    x[mask], 
                    y[mask], 
                    label=group,
                    color=colors[i % len(colors)],
                    s=3,  # Very small points for overlay
                    alpha=0.7,
                    edgecolors='none'
                )
                
            ax.legend(markerscale=5, loc='center left', bbox_to_anchor=(1, 0.5), title="Cluster")
            ax.set_title(f"Overlay: {sample}")
            ax.axis('off') 
            
            # Save
            out_file_png = os.path.join(plot_dir, f"overlay_HE_{sample}.png")
            out_file_pdf = os.path.join(plot_dir, f"overlay_HE_{sample}.pdf")
            
            plt.savefig(out_file_pdf, bbox_inches='tight')
            plt.savefig(out_file_png, bbox_inches='tight', dpi=300)
            plt.close()
            print(f"  Saved overlay to {out_file_png}")

    print("\nDone. Please check the output folder.")

if __name__ == "__main__":
    main()

