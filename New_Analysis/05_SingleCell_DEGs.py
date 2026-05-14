# ==============================================================================
# 文件名称: 04_SingleCell_DEGs.py
# 功能描述: 单细胞水平差异基因分析 (DEGs)
# 联动说明: 承接 02 步的注释结果，针对病理细胞（GC、DN）与正常细胞（Astrocytes、正常神经元）进行差异表达分析，输出差异基因列表。
# ==============================================================================

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import gseapy as gp

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


# Add current directory to path to import config
sys.path.append(os.getcwd())
try:
    from analysis_config import *
except ImportError:
    sys.path.append(r"d:\Data_Analysis\空间转录组分析\ST_analysis\New_Analysis")
    from analysis_config import *

# 绘图设置
sc.set_figure_params(dpi=300, facecolor='white', vector_friendly=True)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

def save_plot(filename, fig=None):
    """同时保存 PDF 和 PNG"""
    if fig is None:
        fig = plt.gcf()
    base, ext = os.path.splitext(filename)
    fig.savefig(f"{base}.pdf", bbox_inches='tight', dpi=300)
    fig.savefig(f"{base}.png", bbox_inches='tight', dpi=300)
    print(f"Saved plots to {base}.pdf/png")

def run_and_plot_enrichment(genes_list, title_prefix, out_dir, comparison_name):
    """
    运行 GO/KEGG 富集分析并绘制条形图
    """
    if not genes_list or len(genes_list) < 3:
        print(f"Not enough genes for enrichment analysis in {comparison_name}")
        return

    print(f"Running enrichment analysis for {comparison_name} ({len(genes_list)} genes)...")
    
    try:
        enr = gp.enrichr(gene_list=genes_list,
                         gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                         organism='human',
                         outdir=None)
        
        res = enr.results
        if res is None or res.empty:
            print(f"No enrichment results for {comparison_name}")
            return
            
        # Filter significant results
        sig_res = res[res['Adjusted P-value'] < 0.05].copy()
        
        if sig_res.empty:
            print(f"No significant enrichment terms (Adj P < 0.05) for {comparison_name}")
            # Fallback to nominal p-value just to show something if needed, but better to keep strict
            return
            
        # Save full results
        sig_res.to_csv(os.path.join(out_dir, f"Enrichment_{comparison_name}.csv"), index=False)
        
        # Plot top 10 terms combined
        top_terms = sig_res.sort_values('Adjusted P-value').head(15)
        
        # Format Term names for better display
        top_terms['Short_Term'] = top_terms['Term'].apply(lambda x: (x[:45] + '...') if len(x) > 48 else x)
        top_terms['Log_P'] = -np.log10(top_terms['Adjusted P-value'])
        
        plt.figure(figsize=(8, 6))
        # Color by database
        palette = {'GO_Biological_Process_2023': '#1f77b4', 'KEGG_2021_Human': '#ff7f0e'}
        
        sns.barplot(data=top_terms, x='Log_P', y='Short_Term', hue='Gene_set', dodge=False, palette=palette)
        
        plt.title(f'Top Enriched Terms: {title_prefix}', fontsize=14, pad=15)
        plt.xlabel('-Log10(Adjusted P-value)', fontsize=12)
        plt.ylabel('')
        
        # Move legend
        plt.legend(title='Database', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        save_plot(os.path.join(out_dir, f"Enrichment_Barplot_{comparison_name}.pdf"))
        plt.close()
        
    except Exception as e:
        print(f"Enrichment analysis failed for {comparison_name}: {e}")

def plot_custom_volcano(df, title, out_path, top_n=10, fc_thresh=1.0, p_thresh=0.05, highlight_up=None):
    """
    绘制 Nature 级别的自定义火山图。
    使用 adjustText 解决标签重叠问题，并优化配色和点大小。
    可以传入 highlight_up 列表，强制高亮这些上调基因。
    """
    df['log_pval'] = -np.log10(df['pvals_adj'] + 1e-300) # 防止 log(0)
    
    # 分类点
    conditions = [
        (df['logfoldchanges'] >= fc_thresh) & (df['pvals_adj'] < p_thresh),
        (df['logfoldchanges'] <= -fc_thresh) & (df['pvals_adj'] < p_thresh)
    ]
    choices = ['Up', 'Down']
    df['Significance'] = np.select(conditions, choices, default='Not Sig')
    
    fig, ax = plt.subplots(figsize=(7, 6))
    
    # 绘制背景点
    sns.scatterplot(
        data=df[df['Significance'] == 'Not Sig'], 
        x='logfoldchanges', y='log_pval', 
        color='grey', alpha=0.3, s=15, edgecolor='none', ax=ax
    )
    
    # 绘制上调点
    sns.scatterplot(
        data=df[df['Significance'] == 'Up'], 
        x='logfoldchanges', y='log_pval', 
        color='#D62728', alpha=0.8, s=30, edgecolor='none', ax=ax # 经典红
    )
    
    # 绘制下调点
    sns.scatterplot(
        data=df[df['Significance'] == 'Down'], 
        x='logfoldchanges', y='log_pval', 
        color='#1F77B4', alpha=0.8, s=30, edgecolor='none', ax=ax # 经典蓝
    )
    
    # 添加阈值线
    ax.axvline(fc_thresh, ls='--', color='black', alpha=0.5, lw=1)
    ax.axvline(-fc_thresh, ls='--', color='black', alpha=0.5, lw=1)
    ax.axhline(-np.log10(p_thresh), ls='--', color='black', alpha=0.5, lw=1)
    
    # 添加标签
    texts = []
    
    if highlight_up is not None:
        # 如果指定了高亮上调基因，就只画这些
        up_genes = df[(df['Significance'] == 'Up') & (df['gene'].isin(highlight_up))]
    else:
        # 否则选取得分最高的前 N 个上调基因
        up_genes = df[df['Significance'] == 'Up'].sort_values('logfoldchanges', ascending=False).head(top_n)
        
    # 下调基因正常画前 N 个
    down_genes = df[df['Significance'] == 'Down'].sort_values('logfoldchanges', ascending=True).head(top_n)
    
    for _, row in pd.concat([up_genes, down_genes]).iterrows():
        texts.append(ax.text(row['logfoldchanges'], row['log_pval'], row['gene'], 
                             fontsize=9, fontstyle='italic'))
    
    # 使用 adjust_text 自动避让标签
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
        
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10(Adj. P-value)', fontsize=12)
    
    # 移除顶部和右侧边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    save_plot(out_path, fig=fig)
    plt.close(fig)

def main():
    print("Step 11: Single-Cell Resolution Differential Expression Analysis...")
    
    # 建立输出目录
    out_dir = get_step_dir("05")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # 1. 加载已经处理好、注释好的单细胞级别空间数据
    tsc_path = PIPELINE_FILES["annotated_h5ad"]
    if not os.path.exists(tsc_path):
        print(f"Error: Annotated data not found at {tsc_path}. Please run previous steps.")
        return
        
    print(f"Loading data from {tsc_path}...")
    adata = sc.read_h5ad(tsc_path)
    
    # 确保 'celltype' 列存在
    if 'celltype' not in adata.obs.columns:
        print("Error: 'celltype' column not found. Have you run 01_5_Manual_Annotation.py?")
        return
        
    print(f"Cell types available: {adata.obs['celltype'].unique().tolist()}")

    # 确保使用归一化后的数据进行差异分析
    if 'log1p' in adata.layers:
        print("Using 'log1p' layer for DE analysis.")
        layer_to_use = 'log1p'
    else:
        print("Warning: 'log1p' layer not found. Using .X.")
        layer_to_use = None

    # ==========================================
    # 分析 0: 细胞类型比例可视化 (单柱堆叠图)
    # ==========================================
    print("\n--- Plotting Cell Type Proportions ---")
    
    # 强制应用指定的细胞类型顺序
    target_order = [
        'Excitatory Neurons',
        'Dysmorphic Neurons(Excit)',
        'Inhibitory Neurons',
        'Dysmorphic Neurons(Inhib)',
        'Giant Cells',
        'Astrocytes',
        'Microglia',
        'OPCs',
        'Oligodendrocytes',
        'Endothelial Cells',
        'T Cells'
    ]
    
    # 获取数据中实际存在的细胞类型并按照 target_order 排序
    existing_types = adata.obs['celltype'].unique()
    ordered_types = [t for t in target_order if t in existing_types]
    # 添加可能存在但不在 target_order 中的类型
    remaining_types = [t for t in existing_types if t not in ordered_types]
    final_order = ordered_types + remaining_types
    
    # 将细胞类型转为按顺序排列的类别
    adata.obs['celltype'] = adata.obs['celltype'].astype('category')
    adata.obs['celltype'] = adata.obs['celltype'].cat.reorder_categories(final_order)
    
    cell_counts = adata.obs['celltype'].value_counts(sort=False) # 保持类别顺序，不按数量排序
    cell_props = adata.obs['celltype'].value_counts(normalize=True, sort=False) * 100
    
    # 保存细胞比例表格
    df_props = pd.DataFrame({
        'CellType': cell_counts.index,
        'Count': cell_counts.values,
        'Percentage': cell_props.values
    })
    df_props.to_csv(os.path.join(out_dir, "CellType_Proportions.csv"), index=False)
    print("Saved CellType_Proportions.csv")

    # 绘制单柱 100% 堆叠图
    plt.figure(figsize=(4, 8))
    
    # 累积的底部位置
    bottom = 0
    colors = sns.color_palette('tab20', n_colors=len(cell_props))
    if len(cell_props) > 20:
        colors = sns.color_palette('husl', n_colors=len(cell_props))
        
    for i, (ctype, prop) in enumerate(cell_props.items()):
        plt.bar('TSC Overall', prop, bottom=bottom, color=colors[i], edgecolor='white', label=f"{ctype} ({prop:.1f}%)")
        bottom += prop
        
    plt.title("Cell Type Proportions in TSC", fontsize=14, pad=15)
    plt.ylabel("Percentage (%)", fontsize=12)
    plt.ylim(0, 100)
    
    # 将图例放在图的右侧
    plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    # 移除上下左右边框以保持美观
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    plt.tight_layout()
    save_plot(os.path.join(out_dir, "CellType_Proportions_Stacked.pdf"))
    plt.close()
    
    # ------------------------------------------
    # 补充分析: 正常神经元与DN中兴奋性/抑制性比例对比
    # ------------------------------------------
    print("\n--- Calculating Excitatory/Inhibitory Ratio (Normal vs DN) ---")
    # 获取所有的正常神经元和DN
    neuron_types = [c for c in adata.obs['celltype'].unique() if 'Neuron' in str(c) and 'Dysmorphic' not in str(c)]
    dn_types = [c for c in adata.obs['celltype'].unique() if 'Dysmorphic' in str(c)]
    
    def get_ex_in_counts(cell_type_list):
        ex_count = 0
        in_count = 0
        for ct in cell_type_list:
            count = adata.obs[adata.obs['celltype'] == ct].shape[0]
            if 'Excit' in ct:
                ex_count += count
            elif 'Inhib' in ct:
                in_count += count
        return ex_count, in_count

    normal_ex, normal_in = get_ex_in_counts(neuron_types)
    dn_ex, dn_in = get_ex_in_counts(dn_types)
    
    if normal_in > 0 and dn_in > 0:
        normal_total = normal_ex + normal_in
        dn_total = dn_ex + dn_in
        
        normal_ratio = normal_ex / normal_in
        dn_ratio = dn_ex / dn_in
        
        normal_ex_pct = normal_ex / normal_total * 100
        normal_in_pct = normal_in / normal_total * 100
        dn_ex_pct = dn_ex / dn_total * 100
        dn_in_pct = dn_in / dn_total * 100
        
        # 统计学差异检验: Fisher Exact Test 或 卡方检验
        import scipy.stats as stats
        contingency_table = [[normal_ex, normal_in], [dn_ex, dn_in]]
        # 样本量大时使用卡方检验更合适
        chi2, p_val, dof, expected = stats.chi2_contingency(contingency_table)
        
        ratio_df = pd.DataFrame({
            'Group': ['Normal Neurons', 'Dysmorphic Neurons (DN)'],
            'Excitatory_Count': [normal_ex, dn_ex],
            'Inhibitory_Count': [normal_in, dn_in],
            'Total_Count': [normal_total, dn_total],
            'Excitatory_Pct': [normal_ex_pct, dn_ex_pct],
            'Inhibitory_Pct': [normal_in_pct, dn_in_pct],
            'E_I_Ratio': [normal_ratio, dn_ratio]
        })
        ratio_df.to_csv(os.path.join(out_dir, "Excitatory_Inhibitory_Ratio_Comparison.csv"), index=False)
        print(f"Saved Excitatory_Inhibitory_Ratio_Comparison.csv (Chi-square p-value: {p_val:.2e})")
        
        # 绘制 100% 堆叠柱状图
        plt.figure(figsize=(6, 7))
        bar_width = 0.6
        
        # 绘制抑制性神经元 (底部)
        p1 = plt.bar(ratio_df['Group'], ratio_df['Inhibitory_Pct'], color='#1f77b4', edgecolor='white', width=bar_width, label='Inhibitory')
        # 绘制兴奋性神经元 (顶部)
        p2 = plt.bar(ratio_df['Group'], ratio_df['Excitatory_Pct'], bottom=ratio_df['Inhibitory_Pct'], color='#d62728', edgecolor='white', width=bar_width, label='Excitatory')
        
        # 添加具体数量文本
        for i, row in ratio_df.iterrows():
            # 抑制性文本 (居中于下半部分)
            plt.text(i, row['Inhibitory_Pct'] / 2, f"n={row['Inhibitory_Count']}", ha='center', va='center', color='white', fontweight='bold')
            # 兴奋性文本 (居中于上半部分)
            plt.text(i, row['Inhibitory_Pct'] + row['Excitatory_Pct'] / 2, f"n={row['Excitatory_Count']}", ha='center', va='center', color='white', fontweight='bold')
            # 顶部显示 E/I Ratio
            plt.text(i, 102, f"E/I Ratio: {row['E_I_Ratio']:.2f}", ha='center', va='bottom', fontweight='bold')
            
        # 添加 P-value 显著性标记
        sig_text = "ns"
        if p_val < 0.001: sig_text = "***"
        elif p_val < 0.01: sig_text = "**"
        elif p_val < 0.05: sig_text = "*"
        
        # 画显著性连线
        y_max = 108
        plt.plot([0, 0, 1, 1], [y_max, y_max+2, y_max+2, y_max], lw=1.2, c='black')
        plt.text(0.5, y_max+3, sig_text + (f" (p={p_val:.1e})" if p_val < 0.05 else f" (p={p_val:.2f})"), ha='center', va='bottom', color='black', fontsize=10)
        
        plt.title("Excitatory / Inhibitory Ratio Imbalance", fontsize=14, pad=25)
        plt.ylabel("Percentage (%)", fontsize=12)
        plt.ylim(0, 120) # 留出空间给标签和连线
        plt.legend(loc='upper right', bbox_to_anchor=(1.35, 1))
        
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        save_plot(os.path.join(out_dir, "Barplot_EI_Ratio_Stacked.pdf"))
        plt.close()
    else:
        print("Warning: Could not calculate E/I ratio because Inhibitory count is 0 in one of the groups.")

    # 用于保存 Analysis 1 的结果，方便后续对比
    dn_up_genes = set()
    dn_down_genes = set()

    # ==========================================
    # 分析 1a: 兴奋性异型神经元 (DN Excitatory) vs 兴奋性正常神经元 (Normal Excitatory)
    # ==========================================
    print("\n--- DE Analysis 1a: Dysmorphic Neurons (Excitatory) vs Normal Excitatory Neurons ---")
    
    dn_excit_types = [c for c in adata.obs['celltype'].unique() if 'Dysmorphic' in str(c) and 'Excit' in str(c)]
    normal_excit_types = [c for c in adata.obs['celltype'].unique() if 'Neuron' in str(c) and 'Excitatory' in str(c) and 'Dysmorphic' not in str(c)]
    
    print(f"Identified DN Excitatory types: {dn_excit_types}")
    print(f"Identified Normal Excitatory types: {normal_excit_types}")
    
    if dn_excit_types and normal_excit_types:
        adata.obs['DE_Group_Excit'] = 'Other'
        adata.obs.loc[adata.obs['celltype'].isin(dn_excit_types), 'DE_Group_Excit'] = 'DN_Excit'
        adata.obs.loc[adata.obs['celltype'].isin(normal_excit_types), 'DE_Group_Excit'] = 'Normal_Excit'
        
        adata_excit = adata[adata.obs['DE_Group_Excit'].isin(['DN_Excit', 'Normal_Excit'])].copy()
        
        sc.tl.rank_genes_groups(
            adata_excit, 
            groupby='DE_Group_Excit', 
            groups=['DN_Excit'], 
            reference='Normal_Excit', 
            method='wilcoxon',
            layer=layer_to_use,
            use_raw=False
        )
        
        result = adata_excit.uns['rank_genes_groups']
        df_excit_de = pd.DataFrame({
            'gene': result['names']['DN_Excit'],
            'logfoldchanges': result['logfoldchanges']['DN_Excit'],
            'pvals': result['pvals']['DN_Excit'],
            'pvals_adj': result['pvals_adj']['DN_Excit'],
            'scores': result['scores']['DN_Excit']
        })
        df_excit_de.to_csv(os.path.join(out_dir, "DE_DN_Excit_vs_Normal_Excit.csv"), index=False)
        print("Saved DE_DN_Excit_vs_Normal_Excit.csv")
        
        # 提取选出的核心基因并在上调中高亮
        core_genes = [
            'EPHA3', 'EPHA5', 'EPHA6', 'EFNA5', 'PLXNA4', 'LAMB1', 'TNR',
            'GABRA5', 'GABRG3', 'SLC35F3', 
            'CNTN3', 'CNTN4', 'CNTN5', 'CNTNAP2', 'CNTNAP5', 
            'CDH9', 'CDH12', 'CDH18', 'CDH22',
            'ROBO1', 'ROBO2', 'LRP1B', 'LRP8', 'NCAM2'
        ]
        
        # 筛选在上调基因中存在于 core_genes 的基因，并按 FC 排序
        up_df = df_excit_de[(df_excit_de['logfoldchanges'] >= 0.5) & (df_excit_de['pvals_adj'] < 0.05)]
        core_up = up_df[up_df['gene'].isin(core_genes)].sort_values('logfoldchanges', ascending=False)
        
        highlight_list = core_up['gene'].tolist()
        if 'TNR' in highlight_list:
            highlight_list.remove('TNR')
        
        # 提取前 10 个，并确保包含 TNR
        top_10_up = ['TNR'] + highlight_list[:9] if 'TNR' in core_up['gene'].values else highlight_list[:10]
        
        plot_custom_volcano(
            df_excit_de, 
            title="DN (Excitatory) vs Normal Excitatory Neurons", 
            out_path=os.path.join(out_dir, "Volcano_Custom_DN_Excit_vs_Normal_Excit.pdf"),
            top_n=10, fc_thresh=0.5,
            highlight_up=top_10_up
        )
        
        sig_genes_excit = df_excit_de[df_excit_de['pvals_adj'] < 0.05]
        excit_up_list = sig_genes_excit[sig_genes_excit['logfoldchanges'] > 0.5]['gene'].tolist()
        excit_down_list = sig_genes_excit[sig_genes_excit['logfoldchanges'] < -0.5]['gene'].tolist()
        
        dn_up_genes.update(excit_up_list)
        dn_down_genes.update(excit_down_list)
        
        run_and_plot_enrichment(excit_up_list, "Upregulated in DN Excit (vs Normal Excit)", out_dir, "DN_Excit_vs_Normal_Up")
        run_and_plot_enrichment(excit_down_list, "Downregulated in DN Excit (vs Normal Excit)", out_dir, "DN_Excit_vs_Normal_Down")
        
        top_up = sig_genes_excit.sort_values('logfoldchanges', ascending=False).head(10)['gene'].tolist()
        top_down = sig_genes_excit.sort_values('logfoldchanges', ascending=True).head(10)['gene'].tolist()
        plot_genes = top_up + top_down
        
        if plot_genes:
            sc.pl.dotplot(adata_excit, plot_genes, groupby='DE_Group_Excit', 
                          standard_scale='var', cmap='Reds', show=False)
            save_plot(os.path.join(out_dir, "Dotplot_DN_Excit_vs_Normal_Excit.pdf"))
            plt.close()
        del adata_excit

    # ==========================================
    # 分析 1b: 抑制性异型神经元 (DN Inhibitory) vs 抑制性正常神经元 (Normal Inhibitory)
    # ==========================================
    print("\n--- DE Analysis 1b: Dysmorphic Neurons (Inhibitory) vs Normal Inhibitory Neurons ---")
    
    dn_inhib_types = [c for c in adata.obs['celltype'].unique() if 'Dysmorphic' in str(c) and 'Inhib' in str(c)]
    normal_inhib_types = [c for c in adata.obs['celltype'].unique() if 'Neuron' in str(c) and 'Inhibitory' in str(c) and 'Dysmorphic' not in str(c)]
    
    print(f"Identified DN Inhibitory types: {dn_inhib_types}")
    print(f"Identified Normal Inhibitory types: {normal_inhib_types}")
    
    if dn_inhib_types and normal_inhib_types:
        adata.obs['DE_Group_Inhib'] = 'Other'
        adata.obs.loc[adata.obs['celltype'].isin(dn_inhib_types), 'DE_Group_Inhib'] = 'DN_Inhib'
        adata.obs.loc[adata.obs['celltype'].isin(normal_inhib_types), 'DE_Group_Inhib'] = 'Normal_Inhib'
        
        adata_inhib = adata[adata.obs['DE_Group_Inhib'].isin(['DN_Inhib', 'Normal_Inhib'])].copy()
        
        sc.tl.rank_genes_groups(
            adata_inhib, 
            groupby='DE_Group_Inhib', 
            groups=['DN_Inhib'], 
            reference='Normal_Inhib', 
            method='wilcoxon',
            layer=layer_to_use,
            use_raw=False
        )
        
        result = adata_inhib.uns['rank_genes_groups']
        df_inhib_de = pd.DataFrame({
            'gene': result['names']['DN_Inhib'],
            'logfoldchanges': result['logfoldchanges']['DN_Inhib'],
            'pvals': result['pvals']['DN_Inhib'],
            'pvals_adj': result['pvals_adj']['DN_Inhib'],
            'scores': result['scores']['DN_Inhib']
        })
        df_inhib_de.to_csv(os.path.join(out_dir, "DE_DN_Inhib_vs_Normal_Inhib.csv"), index=False)
        print("Saved DE_DN_Inhib_vs_Normal_Inhib.csv")
        
        plot_custom_volcano(
            df_inhib_de, 
            title="DN (Inhibitory) vs Normal Inhibitory Neurons", 
            out_path=os.path.join(out_dir, "Volcano_Custom_DN_Inhib_vs_Normal_Inhib.pdf"),
            top_n=12, fc_thresh=0.5
        )
        
        sig_genes_inhib = df_inhib_de[df_inhib_de['pvals_adj'] < 0.05]
        inhib_up_list = sig_genes_inhib[sig_genes_inhib['logfoldchanges'] > 0.5]['gene'].tolist()
        inhib_down_list = sig_genes_inhib[sig_genes_inhib['logfoldchanges'] < -0.5]['gene'].tolist()
        
        dn_up_genes.update(inhib_up_list)
        dn_down_genes.update(inhib_down_list)
        
        run_and_plot_enrichment(inhib_up_list, "Upregulated in DN Inhib (vs Normal Inhib)", out_dir, "DN_Inhib_vs_Normal_Up")
        run_and_plot_enrichment(inhib_down_list, "Downregulated in DN Inhib (vs Normal Inhib)", out_dir, "DN_Inhib_vs_Normal_Down")
        
        top_up = sig_genes_inhib.sort_values('logfoldchanges', ascending=False).head(10)['gene'].tolist()
        top_down = sig_genes_inhib.sort_values('logfoldchanges', ascending=True).head(10)['gene'].tolist()
        plot_genes = top_up + top_down
        
        if plot_genes:
            sc.pl.dotplot(adata_inhib, plot_genes, groupby='DE_Group_Inhib', 
                          standard_scale='var', cmap='Reds', show=False)
            save_plot(os.path.join(out_dir, "Dotplot_DN_Inhib_vs_Normal_Inhib.pdf"))
            plt.close()
        del adata_inhib

    # 获取所有 normal neurons，用于后续的 GC 比较
    normal_neuron_types = [c for c in adata.obs['celltype'].unique() if 'Neuron' in str(c) and 'Dysmorphic' not in str(c)]

    # ==========================================
    # 分析 2: 巨细胞 (Giant Cells) vs 胶质细胞 (Astrocytes)
    # ==========================================
    print("\n--- DE Analysis 2: Giant Cells vs Astrocytes (Glial context) ---")
    
    gc_types = [c for c in adata.obs['celltype'].unique() if 'Giant' in str(c)]
    glial_types = [c for c in adata.obs['celltype'].unique() if 'Astrocyte' in str(c)] # 也可以加上 OPCs, Oligodendrocytes 等
    
    print(f"Identified GC types: {gc_types}")
    print(f"Identified Reference Glial types: {glial_types}")
    
    if gc_types and glial_types:
        adata.obs['DE_Group_GC_Astro'] = 'Other'
        adata.obs.loc[adata.obs['celltype'].isin(gc_types), 'DE_Group_GC_Astro'] = 'Giant_Cells'
        adata.obs.loc[adata.obs['celltype'].isin(glial_types), 'DE_Group_GC_Astro'] = 'Astrocytes'
        
        adata_gc_astro = adata[adata.obs['DE_Group_GC_Astro'].isin(['Giant_Cells', 'Astrocytes'])].copy()
        
        sc.tl.rank_genes_groups(
            adata_gc_astro, 
            groupby='DE_Group_GC_Astro', 
            groups=['Giant_Cells'], 
            reference='Astrocytes', 
            method='wilcoxon',
            layer=layer_to_use,
            use_raw=False
        )
        
        result = adata_gc_astro.uns['rank_genes_groups']
        df_gc_astro = pd.DataFrame({
            'gene': result['names']['Giant_Cells'],
            'logfoldchanges': result['logfoldchanges']['Giant_Cells'],
            'pvals_adj': result['pvals_adj']['Giant_Cells']
        })
        df_gc_astro.to_csv(os.path.join(out_dir, "DE_GiantCells_vs_Astrocytes.csv"), index=False)
        
        plot_custom_volcano(
            df_gc_astro, 
            title="Giant Cells vs Astrocytes", 
            out_path=os.path.join(out_dir, "Volcano_Custom_GiantCells_vs_Astrocytes.pdf"),
            top_n=12, fc_thresh=0.5
        )
        
        sig_genes = df_gc_astro[df_gc_astro['pvals_adj'] < 0.05]
        
        # 分离上调和下调基因进行富集分析
        up_genes_list = sig_genes[sig_genes['logfoldchanges'] > 0.5]['gene'].tolist()
        down_genes_list = sig_genes[sig_genes['logfoldchanges'] < -0.5]['gene'].tolist()
        
        run_and_plot_enrichment(up_genes_list, "Upregulated in GC (vs Astro)", out_dir, "GC_vs_Astro_Up")
        run_and_plot_enrichment(down_genes_list, "Downregulated in GC (vs Astro)", out_dir, "GC_vs_Astro_Down")
        
        top_up = sig_genes.sort_values('logfoldchanges', ascending=False).head(10)['gene'].tolist()
        top_down = sig_genes.sort_values('logfoldchanges', ascending=True).head(10)['gene'].tolist()
        plot_genes = top_up + top_down
        
        if plot_genes:
            sc.pl.dotplot(adata_gc_astro, plot_genes, groupby='DE_Group_GC_Astro', 
                          standard_scale='var', cmap='Reds', show=False)
            save_plot(os.path.join(out_dir, "Dotplot_GiantCells_vs_Astrocytes.pdf"))
            plt.close()
        
        del adata_gc_astro

    # ==========================================
    # 分析 3: 巨细胞 (Giant Cells) vs 正常神经元
    # ==========================================
    print("\n--- DE Analysis 3: Giant Cells vs Normal Neurons ---")
    
    if gc_types and normal_neuron_types:
        adata.obs['DE_Group_GC_Neuron'] = 'Other'
        adata.obs.loc[adata.obs['celltype'].isin(gc_types), 'DE_Group_GC_Neuron'] = 'Giant_Cells'
        adata.obs.loc[adata.obs['celltype'].isin(normal_neuron_types), 'DE_Group_GC_Neuron'] = 'Normal_Neuron'
        
        adata_gc_neuron = adata[adata.obs['DE_Group_GC_Neuron'].isin(['Giant_Cells', 'Normal_Neuron'])].copy()
        
        sc.tl.rank_genes_groups(
            adata_gc_neuron, 
            groupby='DE_Group_GC_Neuron', 
            groups=['Giant_Cells'], 
            reference='Normal_Neuron', 
            method='wilcoxon',
            layer=layer_to_use,
            use_raw=False
        )
        
        result = adata_gc_neuron.uns['rank_genes_groups']
        df_gc_neuron = pd.DataFrame({
            'gene': result['names']['Giant_Cells'],
            'logfoldchanges': result['logfoldchanges']['Giant_Cells'],
            'pvals_adj': result['pvals_adj']['Giant_Cells']
        })
        df_gc_neuron.to_csv(os.path.join(out_dir, "DE_GiantCells_vs_NormalNeuron.csv"), index=False)
        
        plot_custom_volcano(
            df_gc_neuron, 
            title="Giant Cells vs Normal Neurons", 
            out_path=os.path.join(out_dir, "Volcano_Custom_GiantCells_vs_NormalNeuron.pdf"),
            top_n=12, fc_thresh=0.5
        )
        
        sig_genes_gc = df_gc_neuron[df_gc_neuron['pvals_adj'] < 0.05]
        gc_up_genes = set(sig_genes_gc[sig_genes_gc['logfoldchanges'] > 0.5]['gene'])
        gc_down_genes = set(sig_genes_gc[sig_genes_gc['logfoldchanges'] < -0.5]['gene'])
        
        run_and_plot_enrichment(list(gc_up_genes), "Upregulated in GC (vs Normal Neuron)", out_dir, "GC_vs_Normal_Up")
        run_and_plot_enrichment(list(gc_down_genes), "Downregulated in GC (vs Normal Neuron)", out_dir, "GC_vs_Normal_Down")
        
        # --- Venn Diagram: GC vs Neuron AND DN vs Neuron ---
        print("Generating Venn diagram for overlapping DEGs (GC/DN vs Normal Neurons)...")
        try:
            from matplotlib_venn import venn2
            fig, axes = plt.subplots(1, 2, figsize=(10, 5))
            
            venn2([gc_up_genes, dn_up_genes], set_labels=('GC Up (vs Neuron)', 'DN Up (vs Neuron)'), ax=axes[0])
            axes[0].set_title('Upregulated Genes Overlap')
            
            venn2([gc_down_genes, dn_down_genes], set_labels=('GC Down', 'DN Down'), ax=axes[1])
            axes[1].set_title('Downregulated Genes Overlap')
            
            plt.tight_layout()
            save_plot(os.path.join(out_dir, "Venn_GC_DN_vs_Neuron_Overlap.pdf"))
            plt.close()
        except ImportError:
            print("matplotlib_venn not installed. Skipping Venn diagram plotting.")
            
        # Save intersecting genes
        intersect_up = gc_up_genes.intersection(dn_up_genes)
        intersect_down = gc_down_genes.intersection(dn_down_genes)
        with open(os.path.join(out_dir, "Overlap_GC_DN_vs_Neuron_Genes.txt"), "w") as f:
            f.write("Common Upregulated Genes (GC and DN vs Normal Neuron):\n")
            f.write(", ".join(intersect_up) + "\n\n")
            f.write("Common Downregulated Genes (GC and DN vs Normal Neuron):\n")
            f.write(", ".join(intersect_down) + "\n")
            
        # --- Enrichment Analysis for Overlapping DEGs ---
        print("Running enrichment analysis for overlapping DEGs...")
        run_and_plot_enrichment(list(intersect_up), "Common Upregulated in GC & DN", out_dir, "Overlap_GC_DN_Up")
        run_and_plot_enrichment(list(intersect_down), "Common Downregulated in GC & DN", out_dir, "Overlap_GC_DN_Down")
            
        del adata_gc_neuron

    # ==========================================
    # 分析 4: 巨细胞 (Giant Cells) 与 星形胶质细胞 (Astrocytes) 的相似性分析
    # ==========================================
    print("\n--- Analysis 4: Expression Similarity between Giant Cells and Astrocytes ---")
    if gc_types and 'glial_types' in locals() and glial_types:
        gc_cells = adata[adata.obs['celltype'].isin(gc_types)]
        astro_cells = adata[adata.obs['celltype'].isin(glial_types)]
        
        if layer_to_use and layer_to_use in gc_cells.layers:
            gc_mean = np.asarray(gc_cells.layers[layer_to_use].mean(axis=0)).flatten()
            astro_mean = np.asarray(astro_cells.layers[layer_to_use].mean(axis=0)).flatten()
        else:
            gc_mean = np.asarray(gc_cells.X.mean(axis=0)).flatten()
            astro_mean = np.asarray(astro_cells.X.mean(axis=0)).flatten()
            
        df_sim = pd.DataFrame({'GC_Mean': gc_mean, 'Astro_Mean': astro_mean}, index=adata.var_names)
        
        # Filter out genes that are rarely expressed to reduce noise
        df_sim = df_sim[(df_sim['GC_Mean'] > 0.1) | (df_sim['Astro_Mean'] > 0.1)]
        
        # Calculate Pearson Correlation
        corr = df_sim['GC_Mean'].corr(df_sim['Astro_Mean'])
        
        plt.figure(figsize=(7, 6))
        sns.scatterplot(data=df_sim, x='Astro_Mean', y='GC_Mean', alpha=0.3, s=10, color='grey', edgecolor='none')
        
        # ---------------------------------------------------------
        # NEW LOGIC: Find genes close to y=x and highly expressed
        # ---------------------------------------------------------
        # 1. Expression level threshold (both must be high)
        expr_thresh = 0.5
        
        # 2. Distance to y=x line (absolute log2 fold change should be close to 0)
        # We use log2(GC/Astro) ratio. A ratio close to 0 means they are similar.
        # Add small epsilon to prevent div by zero
        df_sim['Log2FC_Ratio'] = np.log2((df_sim['GC_Mean'] + 1e-5) / (df_sim['Astro_Mean'] + 1e-5))
        
        # Filter candidates based on expression and fold-change
        mask_candidates = (df_sim['GC_Mean'] > expr_thresh) & \
                          (df_sim['Astro_Mean'] > expr_thresh) & \
                          (df_sim['Log2FC_Ratio'].abs() < 0.5) # Fold change within ~1.4x
                          
        top_shared = df_sim[mask_candidates].copy()
        
        # Sort by average expression to get the most prominent shared markers
        top_shared['Avg_Expr'] = (top_shared['GC_Mean'] + top_shared['Astro_Mean']) / 2
        top_shared = top_shared.sort_values('Avg_Expr', ascending=False)
        
        # Highlight top shared highly expressed genes (markers)
        top_plot = top_shared.head(25)
        
        sns.scatterplot(data=top_plot, x='Astro_Mean', y='GC_Mean', color='#D62728', s=30, edgecolor='none')
        
        # Add y=x reference line
        max_val = max(df_sim['GC_Mean'].max(), df_sim['Astro_Mean'].max())
        plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1)
        
        texts = []
        for gene, row in top_plot.iterrows():
            texts.append(plt.text(row['Astro_Mean'], row['GC_Mean'], gene, fontsize=9, fontstyle='italic'))
        
        if texts:
            adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
            
        plt.title(f"Expression Similarity: GC vs Astrocytes\nPearson R = {corr:.3f} | {len(top_shared)} Similar Genes", fontsize=14)
        plt.xlabel("Astrocytes Average Expression", fontsize=12)
        plt.ylabel("Giant Cells Average Expression", fontsize=12)
        
        # Remove top/right spines
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        save_plot(os.path.join(out_dir, "Scatter_Similarity_GC_vs_Astro.pdf"))
        plt.close()
        
        # Save shared markers list
        top_shared.to_csv(os.path.join(out_dir, "Shared_Markers_GC_and_Astrocytes.csv"))
        
        # --- Enrichment Analysis for Similarity Genes ---
        print(f"Running enrichment analysis for {len(top_shared)} shared similarity genes...")
        shared_genes_list = top_shared.index.tolist()
        
        if len(shared_genes_list) >= 3:
            run_and_plot_enrichment(shared_genes_list, "Similar Genes (GC & Astro)", out_dir, "Similar_GC_Astro")

    # ==========================================
    # 分析 5: 巨细胞 (Giant Cells) 与 内皮细胞 (Endothelial Cells) 的相似性分析
    # ==========================================
    print("\n--- Analysis 5: Expression Similarity between Giant Cells and Endothelial Cells ---")
    endo_types = [c for c in adata.obs['celltype'].unique() if 'Endothelial' in str(c) or 'Endo' in str(c)]
    
    print(f"Identified Endothelial types: {endo_types}")
    
    if gc_types and endo_types:
        gc_cells = adata[adata.obs['celltype'].isin(gc_types)]
        endo_cells = adata[adata.obs['celltype'].isin(endo_types)]
        
        if layer_to_use and layer_to_use in gc_cells.layers:
            gc_mean = np.asarray(gc_cells.layers[layer_to_use].mean(axis=0)).flatten()
            endo_mean = np.asarray(endo_cells.layers[layer_to_use].mean(axis=0)).flatten()
        else:
            gc_mean = np.asarray(gc_cells.X.mean(axis=0)).flatten()
            endo_mean = np.asarray(endo_cells.X.mean(axis=0)).flatten()
            
        df_sim_endo = pd.DataFrame({'GC_Mean': gc_mean, 'Endo_Mean': endo_mean}, index=adata.var_names)
        
        df_sim_endo = df_sim_endo[(df_sim_endo['GC_Mean'] > 0.1) | (df_sim_endo['Endo_Mean'] > 0.1)]
        corr_endo = df_sim_endo['GC_Mean'].corr(df_sim_endo['Endo_Mean'])
        
        plt.figure(figsize=(7, 6))
        sns.scatterplot(data=df_sim_endo, x='Endo_Mean', y='GC_Mean', alpha=0.3, s=10, color='grey', edgecolor='none')
        
        expr_thresh = 0.5
        df_sim_endo['Log2FC_Ratio'] = np.log2((df_sim_endo['GC_Mean'] + 1e-5) / (df_sim_endo['Endo_Mean'] + 1e-5))
        
        mask_candidates_endo = (df_sim_endo['GC_Mean'] > expr_thresh) & \
                               (df_sim_endo['Endo_Mean'] > expr_thresh) & \
                               (df_sim_endo['Log2FC_Ratio'].abs() < 0.5)
                               
        top_shared_endo = df_sim_endo[mask_candidates_endo].copy()
        
        top_shared_endo['Avg_Expr'] = (top_shared_endo['GC_Mean'] + top_shared_endo['Endo_Mean']) / 2
        top_shared_endo = top_shared_endo.sort_values('Avg_Expr', ascending=False)
        
        top_plot_endo = top_shared_endo.head(25)
        sns.scatterplot(data=top_plot_endo, x='Endo_Mean', y='GC_Mean', color='#1F77B4', s=30, edgecolor='none')
        
        max_val_endo = max(df_sim_endo['GC_Mean'].max(), df_sim_endo['Endo_Mean'].max())
        plt.plot([0, max_val_endo], [0, max_val_endo], 'k--', alpha=0.5, linewidth=1)
        
        texts_endo = []
        for gene, row in top_plot_endo.iterrows():
            texts_endo.append(plt.text(row['Endo_Mean'], row['GC_Mean'], gene, fontsize=9, fontstyle='italic'))
        
        if texts_endo:
            adjust_text(texts_endo, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
            
        plt.title(f"Expression Similarity: GC vs Endothelial Cells\nPearson R = {corr_endo:.3f} | {len(top_shared_endo)} Similar Genes", fontsize=14)
        plt.xlabel("Endothelial Cells Average Expression", fontsize=12)
        plt.ylabel("Giant Cells Average Expression", fontsize=12)
        
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        save_plot(os.path.join(out_dir, "Scatter_Similarity_GC_vs_Endo.pdf"))
        plt.close()
        
        top_shared_endo.to_csv(os.path.join(out_dir, "Shared_Markers_GC_and_Endothelial.csv"))
        
        print(f"Running enrichment analysis for {len(top_shared_endo)} shared similarity genes (GC vs Endo)...")
        shared_genes_list_endo = top_shared_endo.index.tolist()
        
        if len(shared_genes_list_endo) >= 3:
            run_and_plot_enrichment(shared_genes_list_endo, "Similar Genes (GC & Endo)", out_dir, "Similar_GC_Endo")

    print(f"\nAll DE analyses completed successfully! Results are in {out_dir}")

if __name__ == "__main__":
    main()
