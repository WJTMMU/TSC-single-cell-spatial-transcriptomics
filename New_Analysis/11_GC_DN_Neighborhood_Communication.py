# ==============================================================================
# 文件名称: 11_GC_DN_Neighborhood_Communication.py
# 功能描述: GC/DN 邻域细胞通讯特异性分析
# 联动说明: 在细胞注释和空间生态位分析基础上，专门聚焦于巨细胞 (GC) 和发育不良神经元 (DN) 作为病理锚点，与其空间邻居之间的局部通讯网络。
# ==============================================================================

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.spatial import cKDTree

import scanpy as sc


# Add current directory to path to import config
sys.path.append(os.getcwd())
try:
    from analysis_config import *
except ImportError:
    sys.path.append(r"d:\Data_Analysis\空间转录组分析\ST_analysis\New_Analysis")
    from analysis_config import *

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


# 绘图设置，保证 PDF 中的文字可以被 Illustrator 或 WPS 等软件编辑
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# 忽略警告
warnings.filterwarnings("ignore")

sc.set_figure_params(dpi=300, facecolor='white', vector_friendly=True)

DN_RADIUS_UM = 100.0
GC_RADIUS_UM = 100.0
NORMAL_RADIUS_UM = 100.0
STEP11_DIR = get_step_dir("11")
LOG_FILE = os.path.join(STEP11_DIR, "neighborhood.log")

def log_msg(msg):
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(msg + "\n")
    print(msg, flush=True)

def analyze_neighborhood_communication_dynamic(adata, pixels_per_micron=1.0):
    log_msg(
        f"Finding cells within dynamic radii "
        f"(DN: {DN_RADIUS_UM:.0f} um, GC: {GC_RADIUS_UM:.0f} um, Normal: {NORMAL_RADIUS_UM:.0f} um)..."
    )

    results = {}
    
    # Store compositions across samples for averaging
    all_dn_comps = []
    all_gc_comps = []
    all_normal_comps = []
    
    # Define uniform colors for all cell types across all pie charts
    unique_celltypes = sorted([str(c) for c in adata.obs['celltype'].unique()])
    
    # 尝试从 adata.uns 获取统一的 Nature 风格配色，如果没有则回退
    global_color_dict = {}
    if 'celltype_colors' in adata.uns and 'celltype' in adata.obs:
        cats = adata.obs['celltype'].cat.categories
        colors = adata.uns['celltype_colors']
        for c, color in zip(cats, colors):
            global_color_dict[str(c)] = color
    else:
        colors_palette = sns.color_palette("tab20", n_colors=max(20, len(unique_celltypes)))
        if len(unique_celltypes) > 20:
            # Extend palette if more than 20 cell types
            colors_palette = sns.color_palette("husl", n_colors=len(unique_celltypes))
        global_color_dict = {ctype: mcolors.to_hex(color) for ctype, color in zip(unique_celltypes, colors_palette)}
        
    global_color_dict['Other (<2%)'] = '#d3d3d3' # Gray for others
    
    for sample_id in adata.obs['sample'].unique():
        log_msg(f"\nProcessing Sample: {sample_id}")
        
        # Avoid loading the whole matrix, just get indices
        log_msg("  Getting sample mask...")
        sample_mask = (adata.obs['sample'] == sample_id).values
        sample_indices = np.where(sample_mask)[0]
        
        log_msg("  Getting obs_sub...")
        # Get obs and obsm for this sample only
        obs_sub = adata.obs.iloc[sample_indices].copy()
        
        log_msg("  Checking spatial...")
        if 'spatial' not in adata.obsm:
            log_msg(f"Error: 'spatial' coordinates not found.")
            continue
            
        log_msg("  Getting coords...")
        coords = adata.obsm['spatial'][sample_indices]
        
        # 1. Process DN niche
        dn_types = [c for c in obs_sub['celltype'].unique() if 'Dysmorphic' in str(c)]
        dn_mask = obs_sub['celltype'].isin(dn_types).values
        dn_indices = np.where(dn_mask)[0]
        
        # 2. Process GC niche
        gc_types = [c for c in obs_sub['celltype'].unique() if 'Giant' in str(c)]
        gc_mask = obs_sub['celltype'].isin(gc_types).values
        gc_indices = np.where(gc_mask)[0]
        
        if len(dn_indices) == 0 and len(gc_indices) == 0:
            log_msg(f"No target cells (DN/GC) found in sample {sample_id}. Skipping.")
            continue
            
        log_msg(f"Found {len(dn_indices)} DN cells and {len(gc_indices)} GC cells.")
        
        tree = cKDTree(coords)
        all_neighbor_indices = []
        
        if len(dn_indices) > 0:
            dn_coords = coords[dn_indices]
            dn_radius_pixels = DN_RADIUS_UM * pixels_per_micron
            dn_neighbors_lists = tree.query_ball_point(dn_coords, r=dn_radius_pixels)
            dn_niche_indices = np.unique(np.concatenate(dn_neighbors_lists)).astype(int)
            all_neighbor_indices.extend(dn_niche_indices)
        else:
            dn_niche_indices = np.array([])
            
        if len(gc_indices) > 0:
            gc_coords = coords[gc_indices]
            gc_radius_pixels = GC_RADIUS_UM * pixels_per_micron
            gc_neighbors_lists = tree.query_ball_point(gc_coords, r=gc_radius_pixels)
            gc_niche_indices = np.unique(np.concatenate(gc_neighbors_lists)).astype(int)
            all_neighbor_indices.extend(gc_niche_indices)
        else:
            gc_niche_indices = np.array([])
            
        # 3. Process Normal Neuron Niche
        neuron_types = [c for c in obs_sub['celltype'].unique() if 'Neuron' in str(c) and 'Dysmorphic' not in str(c)]
        neuron_mask = obs_sub['celltype'].isin(neuron_types).values
        neuron_indices = np.where(neuron_mask)[0]
        
        if len(neuron_indices) > 0:
            neuron_coords = coords[neuron_indices]
            neuron_radius_pixels = NORMAL_RADIUS_UM * pixels_per_micron
            neuron_neighbors_lists = tree.query_ball_point(neuron_coords, r=neuron_radius_pixels)
            
            # The user wants the ACTUAL cell proportions in the Normal Neuron Niche, 
            # meaning we DO NOT filter out neighborhoods that contain DN/GC cells.
            # We simply take all cells within 100um of ANY normal neuron.
            normal_niche_indices = np.unique(np.concatenate(neuron_neighbors_lists)).astype(int)
        else:
            normal_niche_indices = np.array([])
        
        unique_neighbor_indices_in_sub = np.unique(all_neighbor_indices).astype(int)
        log_msg(f"Found {len(unique_neighbor_indices_in_sub)} total cells (including targets) within the combined neighborhoods.")
        
        log_msg("  Assigning neighborhood mask...")
        # Create a categorical column to distinguish roles for better visualization
        obs_sub = obs_sub.copy()
        obs_sub['neighborhood_role'] = 'Background'
        obs_sub.iloc[unique_neighbor_indices_in_sub, obs_sub.columns.get_loc('neighborhood_role')] = 'Neighborhood Cells'
        obs_sub.iloc[dn_indices, obs_sub.columns.get_loc('neighborhood_role')] = 'DN (Target)'
        obs_sub.iloc[gc_indices, obs_sub.columns.get_loc('neighborhood_role')] = 'GC (Target)'
        
        # Keep boolean mask for filtering
        obs_sub['in_target_neighborhood'] = False
        obs_sub.iloc[unique_neighbor_indices_in_sub, obs_sub.columns.get_loc('in_target_neighborhood')] = True
        
        # New: Specific masks for DN, GC, and Normal niches
        obs_sub['in_dn_niche'] = False
        if len(dn_niche_indices) > 0:
            obs_sub.iloc[dn_niche_indices, obs_sub.columns.get_loc('in_dn_niche')] = True
            
        obs_sub['in_gc_niche'] = False
        if len(gc_niche_indices) > 0:
            obs_sub.iloc[gc_niche_indices, obs_sub.columns.get_loc('in_gc_niche')] = True
            
        obs_sub['in_normal_niche'] = False
        if len(normal_niche_indices) > 0:
            obs_sub.iloc[normal_niche_indices, obs_sub.columns.get_loc('in_normal_niche')] = True
        
        log_msg("  Making output directory...")
        out_dir = STEP11_DIR
        os.makedirs(out_dir, exist_ok=True)
        
        # Skip spatial plotting to save memory/prevent crashes
        log_msg("  Skipping spatial visualization for now.")
        
        log_msg("  Calculating composition...")
        
        # Calculate combined composition
        neighbor_cells_obs = obs_sub[obs_sub['in_target_neighborhood']]
        composition = neighbor_cells_obs['celltype'].value_counts()
        composition_df = pd.DataFrame({'Count': composition, 'Percentage': composition / len(neighbor_cells_obs) * 100})
        composition_df.to_csv(os.path.join(out_dir, f"{sample_id}_Combined_Neighborhood_Composition.csv"))
        
        # Calculate DN Niche composition
        if len(dn_niche_indices) > 0:
            dn_niche_obs = obs_sub[obs_sub['in_dn_niche']]
            dn_comp = dn_niche_obs['celltype'].value_counts(normalize=True) * 100
            dn_comp.name = sample_id
            all_dn_comps.append(dn_comp)
            
            # Local saving
            dn_comp_counts = dn_niche_obs['celltype'].value_counts()
            dn_comp_df = pd.DataFrame({'Count': dn_comp_counts, 'Percentage': dn_comp})
            dn_comp_df.to_csv(os.path.join(out_dir, f"{sample_id}_DN_Niche_Composition.csv"))
            
        # Calculate GC Niche composition
        if len(gc_niche_indices) > 0:
            gc_niche_obs = obs_sub[obs_sub['in_gc_niche']]
            gc_comp = gc_niche_obs['celltype'].value_counts(normalize=True) * 100
            gc_comp.name = sample_id
            all_gc_comps.append(gc_comp)
            
            # Local saving
            gc_comp_counts = gc_niche_obs['celltype'].value_counts()
            gc_comp_df = pd.DataFrame({'Count': gc_comp_counts, 'Percentage': gc_comp})
            gc_comp_df.to_csv(os.path.join(out_dir, f"{sample_id}_GC_Niche_Composition.csv"))
            
        # Calculate Normal Niche composition
        if len(normal_niche_indices) > 0:
            normal_niche_obs = obs_sub[obs_sub['in_normal_niche']]
            normal_comp = normal_niche_obs['celltype'].value_counts(normalize=True) * 100
            normal_comp.name = sample_id
            all_normal_comps.append(normal_comp)
            
            # Local saving
            normal_comp_counts = normal_niche_obs['celltype'].value_counts()
            normal_comp_df = pd.DataFrame({'Count': normal_comp_counts, 'Percentage': normal_comp})
            normal_comp_df.to_csv(os.path.join(out_dir, f"{sample_id}_Normal_Niche_Composition.csv"))
            
        # Visualize DN vs GC Niche Proportions (Per Sample)
        if len(dn_niche_indices) > 0 and len(gc_niche_indices) > 0:
            try:
                plot_df = pd.DataFrame({
                    'DN_Niche': dn_comp,
                    'GC_Niche': gc_comp
                }).fillna(0)
                
                # Plot side-by-side bar chart
                fig, ax = plt.subplots(figsize=(12, 6))
                plot_df.plot(kind='bar', ax=ax, color=['#1f77b4', '#d62728'], width=0.8)
                
                plt.title(f"Cell Type Proportions in Pathological Niches ({sample_id})", fontsize=14)
                plt.xlabel("Cell Type", fontsize=12)
                plt.ylabel("Percentage (%)", fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.legend(title="Niche Type")
                plt.tight_layout()
                
                fig.savefig(os.path.join(out_dir, f"{sample_id}_Niche_Proportions_Comparison.pdf"), bbox_inches='tight', dpi=300)
                fig.savefig(os.path.join(out_dir, f"{sample_id}_Niche_Proportions_Comparison.png"), bbox_inches='tight', dpi=300)
                plt.close(fig)
                log_msg(f"Saved niche proportion comparison plot for {sample_id}")
            except Exception as e:
                log_msg(f"Error plotting niche proportions for {sample_id}: {e}")
        
        # Map back to global indices
        global_neighbor_indices = sample_indices[unique_neighbor_indices_in_sub]
        results[sample_id] = global_neighbor_indices

    # ==========================================
    # Global Average Visualizations across samples
    # ==========================================
    log_msg("\n--- Generating Global Averaged Visualizations ---")
    if all_dn_comps and all_gc_comps:
        # Calculate mean across samples
        df_dn_all = pd.DataFrame(all_dn_comps).fillna(0)
        df_gc_all = pd.DataFrame(all_gc_comps).fillna(0)
        
        mean_dn_comp = df_dn_all.mean(axis=0)
        mean_gc_comp = df_gc_all.mean(axis=0)
        
        # Calculate Normal mean if available
        if all_normal_comps:
            df_normal_all = pd.DataFrame(all_normal_comps).fillna(0)
            mean_normal_comp = df_normal_all.mean(axis=0)
            avg_df = pd.DataFrame({'DN_Niche_Avg': mean_dn_comp, 'GC_Niche_Avg': mean_gc_comp, 'Normal_Niche_Avg': mean_normal_comp}).fillna(0)
        else:
            avg_df = pd.DataFrame({'DN_Niche_Avg': mean_dn_comp, 'GC_Niche_Avg': mean_gc_comp}).fillna(0)
        
        # Save average data
        avg_df.to_csv(os.path.join(out_dir, "AllSamples_Avg_Niche_Proportions.csv"))
        
        # 1. Combined Bar Chart
        fig, ax = plt.subplots(figsize=(12, 6))
        
        bar_colors = ['#1f77b4', '#d62728'] if not all_normal_comps else ['#1f77b4', '#d62728', '#2ca02c']
        avg_df.plot(kind='bar', ax=ax, color=bar_colors, width=0.8)
        
        plt.title("Average Cell Type Proportions in Niches (All Samples)", fontsize=14)
        plt.xlabel("Cell Type", fontsize=12)
        plt.ylabel("Average Percentage (%)", fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.legend(title="Niche Type")
        plt.tight_layout()
        fig.savefig(os.path.join(out_dir, "AllSamples_Avg_Niche_Proportions_Bar.pdf"), bbox_inches='tight', dpi=300)
        fig.savefig(os.path.join(out_dir, "AllSamples_Avg_Niche_Proportions_Bar.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)
        
        # 2. Pie Charts
        def plot_pie(data_series, title, filename):
            # Filter out very small categories for cleaner pie chart, group into 'Other'
            threshold = 2.0
            mask = data_series > threshold
            tail = data_series[~mask].sum()
            data_to_plot = data_series[mask].copy()
            if tail > 0:
                data_to_plot['Other (<2%)'] = tail
                
            fig, ax = plt.subplots(figsize=(8, 8))
            
            # Map global colors to current slice
            pie_colors = [global_color_dict.get(str(idx), '#d3d3d3') for idx in data_to_plot.index]
            
            ax.pie(data_to_plot, labels=data_to_plot.index, autopct='%1.1f%%', 
                   startangle=140, colors=pie_colors, pctdistance=0.85,
                   textprops={'fontsize': 10})
            
            # Draw circle in center to make it a donut chart (looks more modern)
            centre_circle = plt.Circle((0,0), 0.50, fc='white')
            fig.gca().add_artist(centre_circle)
            
            plt.title(title, fontsize=16)
            plt.tight_layout()
            fig.savefig(os.path.join(out_dir, filename + ".pdf"), bbox_inches='tight', dpi=300)
            fig.savefig(os.path.join(out_dir, filename + ".png"), bbox_inches='tight', dpi=300)
            plt.close(fig)

        plot_pie(mean_dn_comp, "Average Cell Composition: DN Niche (100μm)", "AllSamples_Avg_DN_Niche_Pie")
        plot_pie(mean_gc_comp, "Average Cell Composition: GC Niche (200μm)", "AllSamples_Avg_GC_Niche_Pie")
        
        if all_normal_comps:
            plot_pie(mean_normal_comp, "Average Cell Composition: Normal Neuron Niche (100μm)", "AllSamples_Avg_Normal_Niche_Pie")
            
        log_msg("Saved global average bar and pie charts.")

        
    return results

def main():
    open(LOG_FILE, "w", encoding="utf-8").close()
    log_msg("Step 11: Targeted Neighborhood Cell-Cell Communication Analysis")
    
    tsc_path = PIPELINE_FILES["annotated_h5ad"]
    if not os.path.exists(tsc_path):
        log_msg(f"Error: {tsc_path} not found.")
        return
        
    log_msg("Loading data...")
    adata = sc.read_h5ad(tsc_path, backed='r')
    
    target_types = [c for c in adata.obs['celltype'].unique() if 'Giant' in str(c) or 'Dysmorphic' in str(c)]
    log_msg(f"Identified target cell types: {target_types}")
    
    neighbor_indices_dict = analyze_neighborhood_communication_dynamic(adata, pixels_per_micron=1.0)
    
    log_msg("\nNext step: You can run CellPhoneDB, COMMOT, or Squidpy ligand-receptor analysis")
    log_msg("specifically on these extracted neighborhood AnnData objects to find what interacts with GC/DN.")
    
    all_global_indices = []
    for sample_id, indices in neighbor_indices_dict.items():
        all_global_indices.extend(indices)
        
    if all_global_indices:
        log_msg("Loading subset of data for export...")
        
        # 释放一些不需要的变量以节省内存
        import gc
        del neighbor_indices_dict
        gc.collect()
        
        # Subset the backed AnnData using indices
        adata_all_neighbors = adata[all_global_indices].to_memory()
        
        # Explicitly modify obs_names and avoid relying on views
        # Make them globally unique by combining sample and original barcode
        new_obs_names = [f"{sample}_{barcode}" for sample, barcode in zip(adata_all_neighbors.obs['sample'], adata_all_neighbors.obs_names)]
        
        # It's possible the original barcodes already had sample prefixes, so we use a robust counter to guarantee uniqueness
        unique_obs_names = []
        seen = set()
        for name in new_obs_names:
            if name in seen:
                i = 1
                while f"{name}_{i}" in seen:
                    i += 1
                unique_name = f"{name}_{i}"
            else:
                unique_name = name
            seen.add(unique_name)
            unique_obs_names.append(unique_name)
            
        # Force re-assignment to ensure uniqueness
        import anndata as ad
        # 创建新的 AnnData 时，手动把 raw 的 X 提取出来，防止 Subset 后引用错位导致崩溃
        # scipy mmwrite 不能直接写入某些特殊格式的视图，必须拷贝为独立的稀疏矩阵
        new_raw = None
        if adata_all_neighbors.raw is not None:
            raw_X = adata_all_neighbors.raw.X.copy()
            new_raw = ad.AnnData(X=raw_X, var=adata_all_neighbors.raw.var)
            
        adata_all_neighbors = ad.AnnData(
            X=adata_all_neighbors.X.copy(),
            obs=adata_all_neighbors.obs,
            var=adata_all_neighbors.var,
            layers={k: v.copy() for k, v in adata_all_neighbors.layers.items()} if adata_all_neighbors.layers else None
        )
        if new_raw is not None:
            adata_all_neighbors.raw = new_raw
            
        adata_all_neighbors.obs_names = unique_obs_names
        
        out_dir = PIPELINE_FILES["communication_input_dir"]
        os.makedirs(out_dir, exist_ok=True)

        params_df = pd.DataFrame([
            {"Parameter": "pixels_per_micron", "Value": 1.0},
            {"Parameter": "dn_radius_um", "Value": DN_RADIUS_UM},
            {"Parameter": "gc_radius_um", "Value": GC_RADIUS_UM},
            {"Parameter": "normal_radius_um", "Value": NORMAL_RADIUS_UM},
        ])
        params_df.to_csv(os.path.join(out_dir, "neighborhood_parameters.csv"), index=False)
        
        meta = pd.DataFrame({
            'Cell': adata_all_neighbors.obs_names,
            'cell_type': adata_all_neighbors.obs['celltype']
        })
        meta.to_csv(os.path.join(out_dir, "meta.txt"), sep='\t', index=False)
        
        layer_to_use = 'counts' if 'counts' in adata_all_neighbors.layers else None
        
        # Instead of large dense CSV, use Matrix Market format for R memory efficiency
        import scipy.io as io
        import scipy.sparse as sp
        
        # Free memory before matrix extraction
        del adata
        gc.collect()
        
        # 彻底转换为标准的 scipy.sparse 格式以防止任何 AnnData 或 DataFrame 的视图问题
        # 使用 raw 或者 X 作为输出，保证 CellChat 能读取到原始 counts 或者归一化后的数据
        # 14 脚本中如果直接跑 CellChat，通常推荐使用 raw counts 然后在其内部 normalize
        # 确保只获取当前子集的 raw 数据，而不是全局的 raw
        if layer_to_use is not None:
            log_msg(f"  Using '{layer_to_use}' layer for MTX export...")
            matrix_to_save = adata_all_neighbors.layers[layer_to_use].copy()
            var_names_to_save = adata_all_neighbors.var_names.tolist()
        elif adata_all_neighbors.raw is not None:
            log_msg("  Using raw counts for MTX export...")
            matrix_to_save = adata_all_neighbors.raw.X.copy()
            var_names_to_save = adata_all_neighbors.raw.var_names.tolist()
        else:
            log_msg("  Using .X for MTX export...")
            matrix_to_save = adata_all_neighbors.X.copy()
            var_names_to_save = adata_all_neighbors.var_names.tolist()
        
        # 强制将矩阵转为 csr_matrix 并剔除任何 AnnData 包装
        if not sp.issparse(matrix_to_save):
            matrix_to_save = sp.csr_matrix(matrix_to_save)
        else:
            # 即使是稀疏矩阵，也强制转一次 csr_matrix 确保它是标准格式
            matrix_to_save = sp.csr_matrix(matrix_to_save)
            
        mtx_path = os.path.join(out_dir, "counts.mtx")
        log_msg(f"  Writing MTX file to: {mtx_path}")
        try:
            # 必须写成 matrix_to_save.T，因为 Seurat/CellChat 需要 Gene x Cell 的矩阵
            # mmwrite 需要标准的 scipy sparse matrix，不能是 AnnData views
            log_msg(f"  Matrix shape (Gene x Cell) before writing: {matrix_to_save.T.shape}")
            with open(mtx_path, 'wb') as f:
                io.mmwrite(f, matrix_to_save.T)
            log_msg(f"  Successfully wrote {mtx_path}")
            
            # Verify file exists and its size
            import time
            if os.path.exists(mtx_path):
                time.sleep(1) # wait for file system sync
                size = os.path.getsize(mtx_path)
                log_msg(f"  VERIFICATION: File exists on disk, size: {size} bytes")
                if size == 0:
                    log_msg("  WARNING: File size is 0 bytes! Writing might have failed silently.")
            else:
                log_msg("  ERROR: File does not exist on disk after mmwrite!")
                
        except Exception as e:
            log_msg(f"  ERROR writing MTX: {str(e)}")
            import traceback
            log_msg(traceback.format_exc())
        
        with open(os.path.join(out_dir, "genes.tsv"), "w") as f:
            for g in var_names_to_save:
                f.write(f"{g}\n")
                
        with open(os.path.join(out_dir, "barcodes.tsv"), "w") as f:
            for b in adata_all_neighbors.obs_names:
                f.write(f"{b}\n")
                
        log_msg(f"\nExported neighborhood data to Matrix Market format at: {out_dir}")
        log_msg("You can now run CellChat or NicheNet on these files to find localized interactions.")
        log_msg("Done!")

if __name__ == "__main__":
    main()

