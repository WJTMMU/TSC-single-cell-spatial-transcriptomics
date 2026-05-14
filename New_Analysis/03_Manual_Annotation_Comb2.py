# ==============================================================================
# 文件名称: 03_Manual_Annotation_Comb2.py
# 功能描述: 手动细胞类型注释
# 联动说明: 承接 01 步输出的聚类和 marker 数据，由用户手动指定细胞类型字典，将聚类编号映射为具体的细胞类型（如 DN, GC 等），并保存为带有细胞类型注释的 h5ad 文件。
# ==============================================================================

import sys
import os
sys.path.append(os.getcwd())
from analysis_config import NPG_COLORS, set_nature_style, DIRS, PIPELINE_FILES, get_step_dir
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse
import numpy as np

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


# ==============================================================================
# 04 手动注释增强版脚本
# 使用方法：
# 1. 运行 01_Data_Preprocessing.py 生成 cluster_markers_top50.csv
# 2. 打开 csv 文件，分析每个 Cluster 的 Marker 基因
# 3. 在下方 manual_annotation 字典中填写 Cluster ID 对应的细胞类型
# 4. 运行本脚本：python 03_Manual_Annotation_Comb2.py
# ==============================================================================

# Plotting parameters
set_nature_style()
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600

def main():
    print("Step 03: Manual Annotation (Enhanced Visualization)...")
    step_dir = get_step_dir("03")
    
    # 1. 加载 01 步生成的聚类结果
    tsc_path = PIPELINE_FILES["processed_h5ad"]
    if not os.path.exists(tsc_path):
        print(f"Error: Data not found at {tsc_path}. Run 01_Data_Preprocessing.py first.")
        return
    
    print(f"Loading {tsc_path}...")
    adata = sc.read_h5ad(tsc_path)
    print(f"Loaded data with {len(adata.obs['leiden'].unique())} clusters.")
    print(f"Current clusters: {sorted(adata.obs['leiden'].unique())}")

    # ==========================================================================
    # Merge Cluster 0 and 9
    # ==========================================================================
    print("Merging cluster 9 into cluster 0 as requested...")
    if '9' in adata.obs['leiden'].cat.categories:
        # Create a mapping dictionary
        cluster_map = {c: c for c in adata.obs['leiden'].cat.categories}
        cluster_map['9'] = '0'
        cluster_map['12'] = '0'
        
        
        # Apply mapping
        adata.obs['leiden'] = adata.obs['leiden'].map(cluster_map).astype('category')
        print("Merged cluster 9 into 0 successfully.")
        print(f"New clusters: {sorted(adata.obs['leiden'].unique())}")
    else:
        print("Warning: Cluster 9 not found in data.")

    # ==========================================================================
    # TODO: 请在此处修改注释字典！
    # 格式：'Cluster_ID': 'Cell_Type_Name'
    # ==========================================================================
    manual_annotation = {
        '0': 'Dysmorphic Neurons(Excit)',
        '1': 'Oligodendrocytes',
        '2': 'Astrocytes',
        '3': 'Excitatory Neurons',
        '4': 'Microglia',
        '5': 'Endothelial Cells',
        '6': 'Dysmorphic Neurons(Inhib)',
        '7': 'Inhibitory Neurons',
        '8': 'OPCs',
        # '9': 'Dendritic_cells', # Merged into 0
        '10': 'Giant Cells',
        '11': 'T Cells',  
    }
    # ==========================================================================

    print("Applying manual annotation...")
    # 检查是否有未注释的 Cluster
    all_clusters = set(adata.obs['leiden'].unique())
    annotated_clusters = set(manual_annotation.keys())
    
    missing = all_clusters - annotated_clusters
    if missing:
        print(f"Warning: The following clusters are missing from your annotation map: {missing}")
        print("They will be labeled as 'Unknown'.")
    
    # 应用映射
    # Convert to string first to avoid "Cannot setitem on a Categorical with a new category" error
    adata.obs['celltype'] = adata.obs['leiden'].astype(str).map(manual_annotation).fillna('Unknown')
    
    # 4. 设置细胞类型顺序 (Reorder categories)
    target_order = [
        'Excitatory Neurons',
        'Dysmorphic Neurons(Excit)',
        'Inhibitory Neurons',
        'Dysmorphic Neurons(Inhib)',
        'Giant Cells',
        'Astrocytes',
        'Endothelial Cells',
        'Microglia',
        'OPCs',
        'Oligodendrocytes',
        'T Cells'
    ]
    
    # Filter out types that might not exist in the data to avoid errors
    existing_types = adata.obs['celltype'].unique()
    valid_order = [t for t in target_order if t in existing_types]
    
    # Append any remaining types not in target_order
    remaining = [t for t in existing_types if t not in valid_order]
    final_order = valid_order + remaining
    
    print(f"Reordering cell types to: {final_order}")
    adata.obs['celltype'] = adata.obs['celltype'].astype('category')
    adata.obs['celltype'] = adata.obs['celltype'].cat.reorder_categories(final_order)
    
    # --- 设定 Nature 风格的全局细胞颜色 (保存在 adata.uns 中) ---
    adata.uns['celltype_colors'] = [NPG_COLORS[i % len(NPG_COLORS)] for i in range(len(final_order))]
    
    # 打印统计
    print("Cell type counts:")
    print(adata.obs['celltype'].value_counts())

    # ==========================================================================
    # 处理样本匿名化映射 (Sample 1, Sample 2, ...)
    # ==========================================================================
    sample_col_original = 'sample' if 'sample' in adata.obs.columns else 'batch'
    sample_col = 'plot_sample'
    
    if sample_col_original in adata.obs.columns:
        unique_samples = adata.obs[sample_col_original].unique()
        sample_mapping = {samp: f"Sample {i+1}" for i, samp in enumerate(unique_samples)}
        adata.obs[sample_col] = adata.obs[sample_col_original].map(sample_mapping)
        print(f"Sample name mapping applied: {sample_mapping}")
    else:
        adata.obs[sample_col] = 'Sample 1'

    # ==========================================================================
    # 6. 新增需求：合并所有神经元类为 'Neurons' 大类 (celltype_major)
    # ==========================================================================
    print("\nGenerating major cell type annotation (Merging all Neurons)...")
    major_map = {
        'Excitatory Neurons': 'Neurons',
        'Dysmorphic Neurons(Excit)': 'Neurons',
        'Inhibitory Neurons': 'Neurons',
        'Dysmorphic Neurons(Inhib)': 'Neurons',
        'Giant Cells': 'Neurons',
    }
    
    # 复制原始细胞类型，然后应用映射（如果不在字典中则保持原样）
    adata.obs['celltype_major'] = adata.obs['celltype'].astype(str)
    adata.obs['celltype_major'] = adata.obs['celltype_major'].apply(lambda x: major_map.get(x, x))
    
    # 设置大类的顺序
    major_order = ['Neurons', 'Astrocytes', 'Endothelial Cells', 'Microglia', 'OPCs', 'Oligodendrocytes', 'T Cells']
    valid_major_order = [t for t in major_order if t in adata.obs['celltype_major'].unique()]
    adata.obs['celltype_major'] = adata.obs['celltype_major'].astype('category').cat.reorder_categories(valid_major_order)
    
    # 设定大类的 Nature 风格配色
    adata.uns['celltype_major_colors'] = [NPG_COLORS[i % len(NPG_COLORS)] for i in range(len(valid_major_order))]
    
    print("Major Cell type counts:")
    print(adata.obs['celltype_major'].value_counts())

    # ==========================================================================
    # 绘制大类 (Major Cell Types) 比例的 Stacked Barplot (合并 Per-sample 和 Overall)
    # ==========================================================================
    print("Generating stacked barplot for major celltypes (combined)...")
    
    # 1. Per-sample proportion
    counts_major = adata.obs.groupby([sample_col, 'celltype_major'], observed=False).size().unstack(fill_value=0)
    proportions_major = counts_major.div(counts_major.sum(axis=1), axis=0)
    
    # 2. Overall proportion (总平均)
    counts_major_total = adata.obs['celltype_major'].value_counts()
    prop_major_total = counts_major_total / counts_major_total.sum()
    prop_major_total = prop_major_total.to_frame(name="Overall").T
    
    # 合并 DataFrame
    combined_proportions_major = pd.concat([proportions_major, prop_major_total])
    
    # 强制按照 valid_major_order 排序，并反转
    ordered_major_cols = [c for c in valid_major_order if c in combined_proportions_major.columns][::-1]
    combined_proportions_major = combined_proportions_major[ordered_major_cols]
    
    # 颜色列表也要相应反转，保持颜色和细胞类型的映射一致
    original_major_colors = adata.uns['celltype_major_colors'][:len(ordered_major_cols)]
    major_colors = original_major_colors[::-1]

    fig_major_bar, ax_major = plt.subplots(figsize=(1.5 * len(combined_proportions_major) + 2, 6))
    combined_proportions_major.plot(
        kind='bar', stacked=True, width=0.8, ax=ax_major,
        color=major_colors,
        edgecolor='black', linewidth=0.5
    )
    ax_major.set_ylabel("Fraction of Cells", fontsize=14, fontweight='bold')
    ax_major.set_xlabel("", fontsize=14)
    ax_major.set_title("Major Cell Type Proportions", fontsize=16, fontweight='bold')
    ax_major.spines['top'].set_visible(False)
    ax_major.spines['right'].set_visible(False)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    
    # 将图例的顺序再次反转回来
    handles, labels = ax_major.get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], title='Major Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=10)
    
    plt.tight_layout()
    barplot_major_path = os.path.join(step_dir, "celltype_major_proportions_combined.pdf")
    plt.savefig(barplot_major_path, bbox_inches='tight', dpi=600)
    plt.close(fig_major_bar)
    print(f"Saved major celltype combined stacked barplot to {barplot_major_path}")

    # 保存带有大类注释的数据
    adata.write(PIPELINE_FILES["annotated_h5ad"])
    print(f"Saved canonical annotated data to {PIPELINE_FILES['annotated_h5ad']}")

    # 绘制大类的 UMAP (Seurat Style: 图在左，图例在右，去掉坐标轴边框，并在聚类中心显示名称)
    plt.figure()
    sc.pl.umap(adata, color=['celltype_major'], legend_loc='right margin', frameon=False, show=False)
    # 获取当前的 axes，然后添加 on data 的标签
    ax = plt.gca()
    sc.pl.umap(adata, color=['celltype_major'], legend_loc='on data', frameon=False, show=False, ax=ax, title="")
    major_plot_path = os.path.join(step_dir, "umap_tsc_annotated_major.png")
    plt.savefig(major_plot_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Saved major cell type UMAP to {major_plot_path}")

    # 绘制大类的 DotPlot
    print("Generating DotPlot for major cell types...")
    major_marker_dict = {
        'Neurons': ['RBFOX3', 'SLC17A7', 'GAD1', 'GAD2'], # 通用神经元标志物
        'Astrocytes': ['GFAP', 'AQP4', 'SLC1A2'],
        'Endothelial Cells': ['PECAM1', 'FLT1', 'VWF'],
        'Microglia': ['CSF1R', 'C3', 'P2RY12'],
        'OPCs': ['PDGFRA', 'CSPG4', 'VCAN'],
        'Oligodendrocytes': ['MOBP', 'ST18'],
        'T Cells': ['CD3E', 'CD247']
    }
    
    flat_major_markers = {}
    for group, genes in major_marker_dict.items():
        valid_genes = [g for g in genes if g in adata.var_names or (adata.raw is not None and g in adata.raw.var_names)]
        if valid_genes:
            flat_major_markers[group] = valid_genes
            
    if flat_major_markers:
        dp_major = sc.pl.dotplot(
            adata, 
            flat_major_markers, 
            groupby='celltype_major', 
            standard_scale='var', 
            cmap='Blues',           # 换用蓝色渐变区分
            colorbar_title='Scaled Mean Expression',
            show=False,
            return_fig=True
        )
        dp_major.style(cmap='Blues', dot_edge_color='black', dot_edge_lw=0.5)
        
        fig_major = dp_major.get_axes()['mainplot_ax'].figure
        main_ax_major = dp_major.get_axes()['mainplot_ax']
        
        labels_major = [item.get_text() for item in main_ax_major.get_xticklabels()]
        main_ax_major.set_xticks(range(len(labels_major)))
        main_ax_major.set_xticklabels(labels_major, rotation=45, ha='right', rotation_mode='anchor')
        
        if 'group_ax' in dp_major.get_axes():
            group_ax_major = dp_major.get_axes()['group_ax']
            for text_obj in group_ax_major.texts:
                text_obj.set_rotation(-45)
                text_obj.set_ha('left')
                text_obj.set_va('bottom')
                
        fig_major.tight_layout()
        dotplot_major_path = os.path.join(step_dir, "annotation_dotplot_major.pdf")
        fig_major.savefig(dotplot_major_path, bbox_inches='tight')
        plt.close(fig_major)
        print(f"Saved major cell type dotplot to {dotplot_major_path}")

    # 3. 绘图验证 (UMAP)
    # 画第一个图（带右侧图例）
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    sc.pl.umap(adata, color=['leiden'], legend_loc='right margin', frameon=False, show=False, ax=axes[0])
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data', frameon=False, show=False, ax=axes[0], title="leiden")
    
    sc.pl.umap(adata, color=['celltype'], legend_loc='right margin', frameon=False, show=False, ax=axes[1])
    sc.pl.umap(adata, color=['celltype'], legend_loc='on data', frameon=False, show=False, ax=axes[1], title="celltype")
    
    plt.tight_layout()
    plot_path = os.path.join(step_dir, "umap_tsc_annotated.png")
    plt.savefig(plot_path, bbox_inches='tight')
    print(f"Saved annotated UMAP plot to {plot_path}")
    
    # 5. 绘制气泡图 (DotPlot) 验证 Marker
    print("Generating DotPlot for marker verification...")
    
    # Define markers for each type (adjust as needed)
    # 优化点 1: 为共享特征和区分特征添加分组标签 (Marker Groups)
    # 优化点 2: 补充正常神经元与DN的差异Marker，以及GC/星胶/内皮的区分Marker
    marker_dict = {
        'Excitatory (Pan-Markers)': ['SLC17A7', 'SATB2', 'NPTX1', 'SLC25A22'], # 正常的兴奋性特征
        'DN (Pan-Markers)': ['NEFL', 'NEFM', 'NEFH'], # DN的通用异常形态标志物
        'Inhibitory (Pan-Markers)': ['GAD1', 'GAD2', 'DLX6-AS1'], # 抑制性特征
        'Giant Cells': ['FN1', 'COL1A1'],  # GC特有的强基质特征 (添加COL1A1区分星胶)
        'GC/Astro/Endo (Shared)': ['VIM', 'SLC7A11'], # GC、星胶、内皮共享的MAPK/间质/应激特征
        'Astrocytes': ['GFAP', 'AQP4', 'SLC1A2'],   # 星胶经典特征
        'Endothelial Cells': ['PECAM1', 'FLT1', 'VWF'], # 内皮经典特征
        'Microglia': ['CSF1R', 'C3','P2RY12'],
        'OPCs': ['PDGFRA', 'CSPG4','VCAN'],
        'Oligodendrocytes': ['MOBP', 'ST18'],
        'T Cells': ['CD3E', 'CD247']
    }
    
    # Flatten list and filter available genes (保持字典的层级结构用于分组显示)
    flat_markers = {}
    for group, genes in marker_dict.items():
        valid_genes = [g for g in genes if g in adata.var_names or (adata.raw is not None and g in adata.raw.var_names)]
        if valid_genes:
            flat_markers[group] = valid_genes
    
    if flat_markers:
        # 优化点 3: 使用分组字典传入 dotplot，这样横轴上方会有漂亮的分组标签框
        # 优化点 4: 更换颜色图 (cmap)，使用更具区分度的颜色（如 'Reds' 或 'viridis'，去掉默认的底层颜色）
        # 优化点 5: 添加 dendrogram=False 确保按照我们设定的顺序排列
        # 优化点 6: 调整字体旋转角度
        dp = sc.pl.dotplot(
            adata, 
            flat_markers, 
            groupby='celltype', 
            standard_scale='var', 
            cmap='Reds',           # 换用红色渐变，视觉上更干净
            colorbar_title='Scaled Mean Expression',
            show=False,
            return_fig=True
        )
        
        # Scanpy's dotplot object has built-in methods for styling
        dp.style(cmap='Reds', dot_edge_color='black', dot_edge_lw=0.5)
        
        # Force the rotation by applying it directly to the figure axes
        fig = dp.get_axes()['mainplot_ax'].figure
        
        axes_dict = dp.get_axes()
        main_ax = axes_dict['mainplot_ax']
        
        # Get labels and set them with rotation
        # 底部的基因名恢复为逆时针45度 (rotation=45, ha='right')
        labels = [item.get_text() for item in main_ax.get_xticklabels()]
        main_ax.set_xticks(range(len(labels)))
        main_ax.set_xticklabels(labels, rotation=45, ha='right', rotation_mode='anchor')
        
        if 'group_ax' in axes_dict:
            group_ax = axes_dict['group_ax']
            
            # Scanpy's group labels are actually stored as Text objects in the axes,
            # not as standard xticklabels. We need to iterate through the text objects.
            for text_obj in group_ax.texts:
                text_obj.set_rotation(-45)
                # Adjust position slightly down and to the right to avoid overlapping the brackets
                text_obj.set_ha('left')
                text_obj.set_va('bottom')
                
        fig.tight_layout()

        dotplot_path = os.path.join(step_dir, "annotation_dotplot_optimized.pdf")
        fig.savefig(dotplot_path, bbox_inches='tight')
        print(f"Saved optimized annotation dotplot to {dotplot_path}")
    else:
        print("Warning: No standard marker genes found in dataset.")
        
    # ==========================================================================
    # Neuron Subtypes 专属可视化：组合 2 (空间与密度强调)
    # ==========================================================================
    # 1. 2D 密度图 (展示过渡态)
    # 1.5 纯色底版 Neuron UMAP (配合密度图使用)
    # 2. Matrixplot / Complex Heatmap (全局差异基因)
    # 3. 坡度连线图 (样本间比例变化)
    # ==========================================================================
    print("\nGenerating Neuron subtypes specific visualizations (Combination 2)...")
    adata_neurons = adata[adata.obs['celltype_major'] == 'Neurons'].copy()
    adata_neurons.obs['celltype'] = adata_neurons.obs['celltype'].cat.remove_unused_categories()
    
    print("Recomputing PCA, Neighbors, and UMAP for Neuron subtypes...")
    sc.tl.pca(adata_neurons, svd_solver='arpack')
    sc.pp.neighbors(adata_neurons, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_neurons)
    
    # 1.5 纯色底版 Neuron UMAP
    print("Generating plain UMAP for Neurons (no labels, single color)...")
    # 获取大类中 'Neurons' 对应的颜色
    neuron_color_idx = valid_major_order.index('Neurons')
    neuron_color = adata.uns['celltype_major_colors'][neuron_color_idx]
    
    # 创建一个临时列用于统一上色
    adata_neurons.obs['neuron_base'] = 'Neurons'
    adata_neurons.uns['neuron_base_colors'] = [neuron_color]
    
    plt.figure()
    sc.pl.umap(adata_neurons, color='neuron_base', legend_loc=None, frameon=False, show=False, title="")
    neuron_umap_path = os.path.join(step_dir, "umap_tsc_neurons_base.pdf")
    plt.savefig(neuron_umap_path, bbox_inches='tight', dpi=600)
    plt.close()
    print(f"Saved plain Neuron UMAP to {neuron_umap_path}")

    # 1. 2D 密度图
    print("Generating Embedding Density for Neuron subtypes...")
    sc.tl.embedding_density(adata_neurons, basis='umap', groupby='celltype')
    plt.figure()
    sc.pl.embedding_density(adata_neurons, basis='umap', key='umap_density_celltype', show=False)
    density_path = os.path.join(step_dir, "density_tsc_neurons_comb2.pdf")
    plt.savefig(density_path, bbox_inches='tight', dpi=600)
    plt.close()
    print(f"Saved Neuron density UMAP to {density_path}")

    # 2. Matrixplot (矩阵热图)
    print("Generating Matrixplot for Neuron subtypes...")
    neuron_marker_dict = {
        'Excitatory': ['SLC17A7', 'SATB2', 'NPTX1', 'SLC25A22'],
        'Inhibitory': ['GAD1', 'GAD2', 'DLX6-AS1'],
        'DN Markers': ['NEFL', 'NEFM', 'NEFH'],
        'GC Markers': ['FN1', 'COL1A1', 'VIM']
    }
    
    flat_neuron_markers = {}
    for group, genes in neuron_marker_dict.items():
        valid_genes = [g for g in genes if g in adata_neurons.var_names or (adata_neurons.raw is not None and g in adata_neurons.raw.var_names)]
        if valid_genes:
            flat_neuron_markers[group] = valid_genes
            
    if flat_neuron_markers:
        fig_matrix = sc.pl.matrixplot(
            adata_neurons, 
            flat_neuron_markers, 
            groupby='celltype', 
            standard_scale='var', 
            cmap='viridis',
            colorbar_title='Scaled Mean Expression',
            show=False,
            return_fig=True
        )
        matrix_path = os.path.join(step_dir, "annotation_matrixplot_neurons_comb2.pdf")
        fig_matrix.savefig(matrix_path, bbox_inches='tight')
        print(f"Saved Neuron subtypes matrixplot to {matrix_path}")
    
    # 3. 坡度图 / 连线折线图 (样本间比例变化 + Overall)
    print("Generating Line Plot (Slopegraph) for Neuron subtypes proportions...")
    valid_neuron_order = [c for c in final_order if c in adata_neurons.obs['celltype'].unique()]
    counts_neu = adata_neurons.obs.groupby([sample_col, 'celltype'], observed=False).size().unstack(fill_value=0)
    proportions_neu = counts_neu.div(counts_neu.sum(axis=1), axis=0)
    
    # 增加 Overall (总平均)
    counts_neu_total = adata_neurons.obs['celltype'].value_counts()
    prop_neu_total = counts_neu_total / counts_neu_total.sum()
    prop_neu_total = prop_neu_total.to_frame(name="Overall").T
    
    # 合并样本和总平均
    proportions_neu_combined = pd.concat([proportions_neu, prop_neu_total])
    proportions_neu_combined = proportions_neu_combined[valid_neuron_order]
    
    color_dict = dict(zip(final_order, adata.uns['celltype_colors']))
    neu_colors = [color_dict[c] for c in valid_neuron_order]
    
    fig_line, ax_line = plt.subplots(figsize=(7, 4))
    
    # 画折线
    proportions_neu_combined.plot(kind='line', marker='o', markersize=8, linewidth=2.5, ax=ax_line, color=neu_colors)
    
    ax_line.set_ylabel("Fraction within Neurons", fontsize=14, fontweight='bold')
    ax_line.set_xlabel("", fontsize=14, fontweight='bold') # 去掉 Sample 标签，因为包含了 Overall
    ax_line.set_title("Neuron Subtypes Dynamics & Overall", fontsize=16, fontweight='bold')
    ax_line.spines['top'].set_visible(False)
    ax_line.spines['right'].set_visible(False)
    
    # 增加一条垂直虚线分隔样本和Overall，增加层次感
    ax_line.axvline(x=len(proportions_neu.index) - 0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    plt.xticks(range(len(proportions_neu_combined.index)), proportions_neu_combined.index, rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    
    plt.legend(title='Neuron Subtype', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=10)
    plt.tight_layout()
    
    lineplot_path = os.path.join(step_dir, "celltype_neurons_lineplot_comb2.pdf")
    plt.savefig(lineplot_path, bbox_inches='tight', dpi=600)
    plt.close(fig_line)
    print(f"Saved neuron subtypes line plot to {lineplot_path}")
    
    # ==========================================================================
    # 4. 方案 A：细胞质量/转录活性的横向小提琴图 (Transcriptional Activity)
    # ==========================================================================
    print("Generating horizontal violin plot for transcriptional activity (Plan A)...")
    # 计算每个细胞检测到的基因数 (n_genes)
    if 'n_genes_by_counts' in adata_neurons.obs.columns:
        n_genes_col = 'n_genes_by_counts'
    else:
        # 手动计算检测到的基因数
        if scipy.sparse.issparse(adata_neurons.X):
            adata_neurons.obs['n_genes'] = (adata_neurons.X > 0).sum(axis=1).A1
        else:
            adata_neurons.obs['n_genes'] = (adata_neurons.X > 0).sum(axis=1)
        n_genes_col = 'n_genes'

    # 按照 valid_neuron_order 排序 (为了和下方折线图的图例从上到下严格一致，这里不需要反转)
    # 因为 seaborn 的横向小提琴图默认也是从上往下画第一项
    df_qc = adata_neurons.obs[[n_genes_col, 'celltype']].copy()
    df_qc['celltype'] = pd.Categorical(df_qc['celltype'], categories=valid_neuron_order, ordered=True)

    fig_qc, ax_qc = plt.subplots(figsize=(8, 4))
    sns.violinplot(
        data=df_qc, x=n_genes_col, y='celltype',
        palette=neu_colors, ax=ax_qc, orient='h',
        inner='box', linewidth=1, scale='width', saturation=1.0 # 强制饱和度为 1.0，防止颜色变暗
    )
    
    # 隐藏 y 轴的 label，因为和下面的图例共用了
    ax_qc.set_ylabel("")
    ax_qc.set_xlabel("Number of Detected Genes per Cell", fontsize=14, fontweight='bold')
    ax_qc.set_title("Transcriptional Activity across Neuron Subtypes", fontsize=16, fontweight='bold')
    ax_qc.spines['top'].set_visible(False)
    ax_qc.spines['right'].set_visible(False)
    
    # 既然要和下面共享图例，我们可以把 Y 轴的文字也去掉，只留刻度，或者连刻度也去掉
    # 这里我们只保留刻度线，去掉文字，让它看起来完全是一个图例的延伸
    ax_qc.set_yticklabels([])
    
    plt.xticks(fontsize=12)
    plt.tight_layout()
    
    qc_path = os.path.join(step_dir, "celltype_neurons_activity_violin_PlanA.pdf")
    plt.savefig(qc_path, bbox_inches='tight', dpi=600)
    plt.close(fig_qc)
    print(f"Saved transcriptional activity violin plot to {qc_path}")

    # ==========================================================================
    # 5. 方案 B：细胞周期比例的横向分组柱状图 (Cell Cycle Phase - Grouped Barplot)
    # ==========================================================================
    print("Generating grouped horizontal barplot for cell cycle phases (Plan B)...")
    
    # 简化版人类细胞周期基因集 (Tirosh et al, 2015)
    s_genes = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
    g2m_genes = ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2', 'CBX5', 'CENPA']

    # 过滤出真实存在于你数据中的基因
    s_genes_valid = [g for g in s_genes if g in adata_neurons.var_names]
    g2m_genes_valid = [g for g in g2m_genes if g in adata_neurons.var_names]

    if len(s_genes_valid) > 5 and len(g2m_genes_valid) > 5:
        # 使用 scanpy 内置函数进行细胞周期打分
        sc.tl.score_genes_cell_cycle(adata_neurons, s_genes=s_genes_valid, g2m_genes=g2m_genes_valid)

        # 统计每个亚群处于 G1, S, G2M 期的比例
        cc_counts = adata_neurons.obs.groupby(['celltype', 'phase'], observed=False).size().unstack(fill_value=0)
        cc_props = cc_counts.div(cc_counts.sum(axis=1), axis=0)
        cc_props = cc_props.reindex(valid_neuron_order)

        # 绘制并排横向柱状图 (Grouped Horizontal Barplot, NOT stacked)
        # plot(kind='barh', stacked=False) 会自动将各个 phase 并排展示
        fig_cc, ax_cc = plt.subplots(figsize=(8, 4.5))
        
        # 自定义 G1, S, G2M 的颜色 (使用 NPG 里的蓝、橙、绿)
        phase_colors = [NPG_COLORS[1], NPG_COLORS[0], NPG_COLORS[2]]
        
        cc_props.plot(
            kind='barh', stacked=False, ax=ax_cc,
            color=phase_colors,
            edgecolor='black', linewidth=0.5, width=0.7
        )

        ax_cc.set_xlabel("Fraction of Cells", fontsize=14, fontweight='bold')
        ax_cc.set_ylabel("", fontsize=14)
        ax_cc.set_title("Cell Cycle Phase Distribution", fontsize=16, fontweight='bold')
        ax_cc.spines['top'].set_visible(False)
        ax_cc.spines['right'].set_visible(False)
        
        # 反转 y 轴，让 Excitatory 在最上面
        ax_cc.invert_yaxis()
        
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12, fontweight='bold')
        plt.legend(title='Phase', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=10)
        plt.tight_layout()

        cc_path = os.path.join(step_dir, "celltype_neurons_cellcycle_grouped_PlanB.pdf")
        plt.savefig(cc_path, bbox_inches='tight', dpi=600)
        plt.close(fig_cc)
        print(f"Saved grouped cell cycle barplot to {cc_path}")
    else:
        print("Not enough cell cycle genes found in dataset to run Plan B.")
    
    # 结尾的大类 Nested Donut Chart 不在 comb2 中强制输出，或者保留也可以
    # 在 comb2 中我们把它也保留下来，因为大类的可视化不受影响
    # ==========================================================================
    # 高级可视化：双层嵌套环形图 (Nested Donut Chart) 展示细胞层级比例
    # ==========================================================================
    print("Generating nested donut chart for hierarchical cell proportions...")
    
    # 统计各个细分亚群的数量
    cell_counts = adata.obs['celltype'].value_counts()
    
    # 准备外环数据 (Sub-types) 和 内环数据 (Major-types)
    inner_labels = []
    inner_sizes = []
    inner_colors = []
    
    outer_labels = []
    outer_sizes = []
    outer_colors = []
    
    # 获取颜色映射字典
    major_color_dict = dict(zip(valid_major_order, adata.uns['celltype_major_colors']))
    sub_color_dict = dict(zip(final_order, adata.uns['celltype_colors']))
    
    for major_type in valid_major_order:
        # 获取属于该大类的所有亚型细胞数量
        subtypes_in_major = adata[adata.obs['celltype_major'] == major_type].obs['celltype'].value_counts()
        subtypes_in_major = subtypes_in_major[subtypes_in_major > 0] # 过滤掉 0 数量的
        
        # 按照 final_order 排序亚型
        ordered_subs = [t for t in final_order if t in subtypes_in_major.index]
        
        major_total = 0
        for sub in ordered_subs:
            count = subtypes_in_major[sub]
            major_total += count
            outer_labels.append(sub)
            outer_sizes.append(count)
            outer_colors.append(sub_color_dict[sub])
            
        if major_total > 0:
            inner_labels.append(major_type)
            inner_sizes.append(major_total)
            inner_colors.append(major_color_dict[major_type])

    # 画双层环形图
    fig_donut, ax_donut = plt.subplots(figsize=(10, 10), subplot_kw=dict(aspect="equal"))
    
    # 外环 (细分亚型)
    wedges_out, texts_out = ax_donut.pie(
        outer_sizes, radius=1.0, colors=outer_colors,
        wedgeprops=dict(width=0.3, edgecolor='w', linewidth=1.5),
        startangle=90, counterclock=False
    )
    
    # 内环 (大类)
    wedges_in, texts_in = ax_donut.pie(
        inner_sizes, radius=0.7, colors=inner_colors,
        wedgeprops=dict(width=0.3, edgecolor='w', linewidth=2),
        startangle=90, counterclock=False
    )
    
    # 添加引线和标签 (外环)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="none", alpha=0.8)
    kw = dict(arrowprops=dict(arrowstyle="-", color="0.2", lw=1),
              bbox=bbox_props, zorder=0, va="center")
    
    for i, p in enumerate(wedges_out):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = f"angle,angleA=0,angleB={ang}"
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        
        # 计算百分比
        percentage = f"{(outer_sizes[i]/sum(outer_sizes))*100:.1f}%"
        label_text = f"{outer_labels[i]}\n({percentage})"
        
        ax_donut.annotate(label_text, xy=(x, y), xytext=(1.2*np.sign(x), 1.2*y),
                          horizontalalignment=horizontalalignment, **kw, fontsize=10)
                          
    # 在中心添加总细胞数
    total_cells = sum(inner_sizes)
    ax_donut.text(0, 0, f"Total Cells\n{total_cells}", ha='center', va='center', fontsize=14, fontweight='bold')
    
    plt.title("Hierarchical Cell Type Proportions", fontsize=18, fontweight='bold', pad=30)
    
    donut_path = os.path.join(step_dir, "celltype_hierarchical_donut.pdf")
    plt.savefig(donut_path, bbox_inches='tight', dpi=600)
    plt.close(fig_donut)
    print(f"Saved hierarchical donut chart to {donut_path}")

    print("\nDone! You can now proceed to the next analysis step.")

if __name__ == "__main__":
    main()

