# ==============================================================================
# 文件名称: 02_Manual_Annotation.py
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
# 01.5 手动注释脚本
# 使用方法：
# 1. 运行 01_TSC_SeekSpace_Workflow.py 生成 cluster_markers_top50.csv
# 2. 打开 csv 文件，分析每个 Cluster 的 Marker 基因
# 3. 在下方 manual_annotation 字典中填写 Cluster ID 对应的细胞类型
# 4. 运行本脚本：python 01_5_Manual_Annotation.py
# ==============================================================================

# Plotting parameters
set_nature_style()
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600

def main():
    print("Step 01.5: Manual Annotation...")
    step_dir = get_step_dir("02")
    
    # 1. 加载 01 步生成的聚类结果
    tsc_path = PIPELINE_FILES["processed_h5ad"]
    if not os.path.exists(tsc_path):
        print(f"Error: Data not found at {tsc_path}. Run 01_TSC_SeekSpace_Workflow.py first.")
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
    # 绘制大类 (Major Cell Types) 比例的 Stacked Barplot
    # ==========================================================================
    print("Generating stacked barplot for major celltypes...")
    
    # 1. Per-sample proportion
    counts_major = adata.obs.groupby([sample_col, 'celltype_major'], observed=False).size().unstack(fill_value=0)
    proportions_major = counts_major.div(counts_major.sum(axis=1), axis=0)
    
    # 强制按照 valid_major_order 排序
    # Pandas 的 plot(stacked=True) 是从 DataFrame 的第一列开始画（画在最下面）
    # 通常图例是从上往下排列，所以我们想要的第一项（比如 Neurons）应该画在最上面
    # 所以我们需要把列的顺序反转 (::-1)
    ordered_major_cols = [c for c in valid_major_order if c in proportions_major.columns][::-1]
    proportions_major = proportions_major[ordered_major_cols]
    
    # 颜色列表也要相应反转，保持颜色和细胞类型的映射一致
    original_major_colors = adata.uns['celltype_major_colors'][:len(ordered_major_cols)]
    major_colors = original_major_colors[::-1]

    fig_major_bar, ax_major = plt.subplots(figsize=(1.5 * len(proportions_major) + 2, 6))
    proportions_major.plot(
        kind='bar', stacked=True, width=0.8, ax=ax_major,
        color=major_colors,
        edgecolor='black', linewidth=0.5
    )
    ax_major.set_ylabel("Fraction of Cells", fontsize=14, fontweight='bold')
    ax_major.set_xlabel("", fontsize=14)
    ax_major.set_title("Major Cell Type (Per Sample)", fontsize=16, fontweight='bold')
    ax_major.spines['top'].set_visible(False)
    ax_major.spines['right'].set_visible(False)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    
    # 因为我们在画图时反转了顺序，所以 pandas 自动生成的图例也会反转（最底下的画在最下面）
    # 我们需要将图例的顺序再次反转回来，让它从上到下的顺序和 valid_major_order 一致
    handles, labels = ax_major.get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], title='Major Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=10)
    
    plt.tight_layout()
    barplot_major_path = os.path.join(step_dir, "celltype_major_proportions_per_sample.pdf")
    plt.savefig(barplot_major_path, bbox_inches='tight', dpi=600)
    plt.close(fig_major_bar)

    # 2. Overall proportion (总平均)
    counts_major_total = adata.obs['celltype_major'].value_counts()
    prop_major_total = counts_major_total / counts_major_total.sum()
    prop_major_total = prop_major_total[ordered_major_cols].to_frame(name="Overall")
    
    fig_major_tot, ax_major_tot = plt.subplots(figsize=(2.5, 6))
    prop_major_total.T.plot(
        kind='bar', stacked=True, width=0.6, ax=ax_major_tot,
        color=major_colors,
        edgecolor='black', linewidth=0.5
    )
    ax_major_tot.set_ylabel("Fraction of Cells", fontsize=14, fontweight='bold')
    ax_major_tot.set_xlabel("", fontsize=14)
    ax_major_tot.set_title("Major Cell Type (Overall)", fontsize=16, fontweight='bold')
    ax_major_tot.spines['top'].set_visible(False)
    ax_major_tot.spines['right'].set_visible(False)
    plt.xticks(rotation=0, fontsize=12)
    plt.yticks(fontsize=12)
    
    # 同样反转总平均图的图例
    handles_tot, labels_tot = ax_major_tot.get_legend_handles_labels()
    plt.legend(handles_tot[::-1], labels_tot[::-1], title='Major Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=10)
    
    plt.tight_layout()
    barplot_major_tot_path = os.path.join(step_dir, "celltype_major_proportions_overall.pdf")
    plt.savefig(barplot_major_tot_path, bbox_inches='tight', dpi=600)
    plt.close(fig_major_tot)
    print(f"Saved major celltype stacked barplots to {step_dir}")

    # 保存带有大类注释的数据
    preview_h5ad = os.path.join(step_dir, "adata_tsc_annotated_basic.h5ad")
    adata.write(preview_h5ad)
    print(f"Saved basic annotated preview to {preview_h5ad}")

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
    # 7. 新增需求：提取 Neurons 大类进行单独的亚型展示 (UMAP 和 DotPlot)
    # ==========================================================================
    print("\nGenerating Neuron subtypes specific visualizations...")
    adata_neurons = adata[adata.obs['celltype_major'] == 'Neurons'].copy()
    
    # 过滤掉 adata_neurons 中 'celltype' 未使用的 categories
    adata_neurons.obs['celltype'] = adata_neurons.obs['celltype'].cat.remove_unused_categories()
    
    # 重新计算 Neuron 亚群的 PCA, Neighbors 和 UMAP，使降维图更紧凑，避免其他细胞类型的空白留白
    print("Recomputing PCA, Neighbors, and UMAP for Neuron subtypes...")
    sc.tl.pca(adata_neurons, svd_solver='arpack')
    sc.pp.neighbors(adata_neurons, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_neurons)
    
    # 绘制 Neuron 亚类的 UMAP
    plt.figure()
    sc.pl.umap(adata_neurons, color=['celltype'], legend_loc='right margin', frameon=False, show=False)
    ax_neu = plt.gca()
    sc.pl.umap(adata_neurons, color=['celltype'], legend_loc='on data', frameon=False, show=False, ax=ax_neu, title="")
    neuron_umap_path = os.path.join(step_dir, "umap_tsc_neurons_only.png")
    plt.savefig(neuron_umap_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Saved Neuron subtypes UMAP to {neuron_umap_path}")

    # 绘制 Neuron 亚类的 DotPlot
    print("Generating DotPlot for Neuron subtypes...")
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
        dp_neuron = sc.pl.dotplot(
            adata_neurons, 
            flat_neuron_markers, 
            groupby='celltype', 
            standard_scale='var', 
            cmap='Oranges',           # 换用橙色渐变区分
            colorbar_title='Scaled Mean Expression',
            show=False,
            return_fig=True
        )
        dp_neuron.style(cmap='Oranges', dot_edge_color='black', dot_edge_lw=0.5)
        
        fig_neuron = dp_neuron.get_axes()['mainplot_ax'].figure
        main_ax_neuron = dp_neuron.get_axes()['mainplot_ax']
        
        labels_neuron = [item.get_text() for item in main_ax_neuron.get_xticklabels()]
        main_ax_neuron.set_xticks(range(len(labels_neuron)))
        main_ax_neuron.set_xticklabels(labels_neuron, rotation=45, ha='right', rotation_mode='anchor')
        
        if 'group_ax' in dp_neuron.get_axes():
            group_ax_neuron = dp_neuron.get_axes()['group_ax']
            for text_obj in group_ax_neuron.texts:
                text_obj.set_rotation(-45)
                text_obj.set_ha('left')
                text_obj.set_va('bottom')
                
        fig_neuron.tight_layout()
        dotplot_neuron_path = os.path.join(step_dir, "annotation_dotplot_neurons_only.pdf")
        fig_neuron.savefig(dotplot_neuron_path, bbox_inches='tight')
        plt.close(fig_neuron)
        print(f"Saved Neuron subtypes dotplot to {dotplot_neuron_path}")
    
    # ==========================================================================
    # 绘制神经元亚型 (Neuron Subtypes) 比例的 Stacked Barplot
    # ==========================================================================
    print("Generating stacked barplot for Neuron subtypes...")
    
    # 提取 Neuron 特有的颜色映射 (对应 final_order)
    color_dict = dict(zip(final_order, adata.uns['celltype_colors']))
    
    # 确定亚型的顺序 (只取存在于 adata_neurons 中的)
    valid_neuron_order = [c for c in final_order if c in adata_neurons.obs['celltype'].unique()]
    
    # 1. Per-sample proportion
    counts_neu = adata_neurons.obs.groupby([sample_col, 'celltype'], observed=False).size().unstack(fill_value=0)
    proportions_neu = counts_neu.div(counts_neu.sum(axis=1), axis=0)
    
    # 反转列顺序，使第一项在最上方
    ordered_neu_cols = [c for c in valid_neuron_order if c in proportions_neu.columns][::-1]
    proportions_neu = proportions_neu[ordered_neu_cols]
    
    # 颜色顺序也跟着反转
    neu_colors = [color_dict[c] for c in ordered_neu_cols]

    fig_neu_bar, ax_neu = plt.subplots(figsize=(1.5 * len(proportions_neu) + 2, 6))
    proportions_neu.plot(
        kind='bar', stacked=True, width=0.8, ax=ax_neu,
        color=neu_colors,
        edgecolor='black', linewidth=0.5
    )
    ax_neu.set_ylabel("Fraction within Neurons", fontsize=14, fontweight='bold')
    ax_neu.set_xlabel("", fontsize=14)
    ax_neu.set_title("Neuron Subtypes (Per Sample)", fontsize=16, fontweight='bold')
    ax_neu.spines['top'].set_visible(False)
    ax_neu.spines['right'].set_visible(False)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    
    # 反转图例顺序
    handles_neu, labels_neu = ax_neu.get_legend_handles_labels()
    plt.legend(handles_neu[::-1], labels_neu[::-1], title='Neuron Subtype', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=10)
    
    plt.tight_layout()
    barplot_neu_path = os.path.join(step_dir, "celltype_neurons_proportions_per_sample.pdf")
    plt.savefig(barplot_neu_path, bbox_inches='tight', dpi=600)
    plt.close(fig_neu_bar)

    # 2. Overall proportion (总平均)
    counts_neu_total = adata_neurons.obs['celltype'].value_counts()
    prop_neu_total = counts_neu_total / counts_neu_total.sum()
    prop_neu_total = prop_neu_total[ordered_neu_cols].to_frame(name="Overall")
    
    fig_neu_tot, ax_neu_tot = plt.subplots(figsize=(2.5, 6))
    prop_neu_total.T.plot(
        kind='bar', stacked=True, width=0.6, ax=ax_neu_tot,
        color=neu_colors,
        edgecolor='black', linewidth=0.5
    )
    ax_neu_tot.set_ylabel("Fraction within Neurons", fontsize=14, fontweight='bold')
    ax_neu_tot.set_xlabel("", fontsize=14)
    ax_neu_tot.set_title("Neuron Subtypes (Overall)", fontsize=16, fontweight='bold')
    ax_neu_tot.spines['top'].set_visible(False)
    ax_neu_tot.spines['right'].set_visible(False)
    plt.xticks(rotation=0, fontsize=12)
    plt.yticks(fontsize=12)
    
    # 同样反转图例
    handles_neu_tot, labels_neu_tot = ax_neu_tot.get_legend_handles_labels()
    plt.legend(handles_neu_tot[::-1], labels_neu_tot[::-1], title='Neuron Subtype', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=10)
    
    plt.tight_layout()
    barplot_neu_tot_path = os.path.join(step_dir, "celltype_neurons_proportions_overall.pdf")
    plt.savefig(barplot_neu_tot_path, bbox_inches='tight', dpi=600)
    plt.close(fig_neu_tot)
    print(f"Saved neuron subtypes stacked barplots to {step_dir}")

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
