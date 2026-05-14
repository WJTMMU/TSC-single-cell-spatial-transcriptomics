# ==============================================================================
# 文件名称: 08_NMF_Downstream.py
# 功能描述: NMF 空间生态位下游分析
# 联动说明: 承接 07 步 NMF 空间生态位分析，进一步解析各个生态位 Factor 的基因特征并进行富集分析，确定病理生态位的生物学属性。
# ==============================================================================

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp # Added gseapy for enrichment analysis
from analysis_contract import (
    DEFAULT_EXPRESSION_LAYER,
    get_expression_matrix,
    load_nmf_artifacts,
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


# Plotting parameters
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# Add current directory to path to import config
sys.path.append(os.getcwd())
try:
    from analysis_config import *
except ImportError:
    sys.path.append(r"d:\Data_Analysis\空间转录组分析\ST_analysis\New_Analysis")
    from analysis_config import *


def align_saved_nmf_to_adata(adata_subset, output_dir, sample_id):
    cell_weights, gene_weights, metadata, _ = load_nmf_artifacts(output_dir, sample_id)
    cell_weights = cell_weights.reindex(adata_subset.obs_names)
    gene_weights = gene_weights.reindex(adata_subset.var_names)

    if cell_weights.isna().any().any():
        raise RuntimeError(f"NMF cell weights for '{sample_id}' do not align with AnnData obs_names.")
    if gene_weights.isna().any().any():
        raise RuntimeError(f"NMF gene weights for '{sample_id}' do not align with AnnData var_names.")

    factor_cols = cell_weights.columns.tolist()
    return cell_weights, gene_weights, metadata, factor_cols

def analyze_nmf_downstream():
    print("Step 08: NMF Downstream Analysis (Interpreting Niches)")
    
    tsc_path = PIPELINE_FILES["annotated_h5ad"]
    if not os.path.exists(tsc_path):
        print(f"Error: {tsc_path} not found.")
        return
        
    print("Loading data...", flush=True)
    adata = sc.read_h5ad(tsc_path)
    require_data_contract(adata, "08_NMF_Downstream.py")
    
    out_dir = PIPELINE_FILES["nmf_downstream_dir"]
    os.makedirs(out_dir, exist_ok=True)

    import scipy.stats as stats

    def process_one_dataset(dataset_id, adata_subset, plot_spatial):
        print(f"\nProcessing Sample: {dataset_id}")
        cell_weights, gene_weights, metadata, factor_cols = align_saved_nmf_to_adata(adata_subset, PIPELINE_FILES["nmf_dir"], dataset_id)
        print(
            f"  Loaded NMF artifacts from Step 07 "
            f"(layer={metadata.get('source_layer')}, n_factors={len(factor_cols)})."
        )

        expr_mat = get_expression_matrix(adata_subset, layer=DEFAULT_EXPRESSION_LAYER, nonnegative=False)
        top_genes_dict = {}

        for factor in factor_cols:
            top_features_raw = gene_weights[factor].sort_values(ascending=False).head(200).index.tolist()
            factor_weights = cell_weights[factor].to_numpy()
            sig_genes = []

            for gene in top_features_raw:
                g_idx = adata_subset.var_names.get_loc(gene)
                g_expr = expr_mat[:, g_idx]
                corr, pval = stats.pearsonr(g_expr, factor_weights)
                if pval < 0.05 and corr > 0:
                    sig_genes.append(gene)

            top_genes_dict[factor] = sig_genes

        max_len = max((len(v) for v in top_genes_dict.values()), default=0)
        save_dict = {k: v + [np.nan] * (max_len - len(v)) for k, v in top_genes_dict.items()}
        pd.DataFrame(save_dict).to_csv(os.path.join(out_dir, f"{dataset_id}_NMF_Top200_Sig_Genes.csv"), index=False)
        print(f"  Saved Top Significant Genes per NMF Factor to {dataset_id}_NMF_Top200_Sig_Genes.csv")

        print(f"  Running GO and KEGG Enrichment Analysis for {dataset_id}...")
        enrichment_results = []
        for factor, genes in top_genes_dict.items():
            valid_genes = [g for g in genes if pd.notna(g)]
            if len(valid_genes) < 3:
                print(f"    Skipping {factor}: Not enough significant genes ({len(valid_genes)}).")
                continue

            try:
                enr = gp.enrichr(
                    gene_list=valid_genes,
                    gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                    organism='human',
                    outdir=None,
                )

                res = enr.results
                sig_res = res[res['Adjusted P-value'] < 0.05].copy()
                for _, row in sig_res.iterrows():
                    enrichment_results.append({
                        'Factor': factor,
                        'Database': row['Gene_set'],
                        'Term': row['Term'],
                        'P-value': row['P-value'],
                        'Adjusted P-value': row['Adjusted P-value'],
                        'Overlap': row['Overlap'],
                        'Genes': row['Genes']
                    })
            except Exception as e:
                print(f"    Failed to run enrichment for {factor}: {e}")

        if enrichment_results:
            df_enrichment = pd.DataFrame(enrichment_results)
            df_enrichment.to_csv(os.path.join(out_dir, f"{dataset_id}_NMF_Factors_Enrichment.csv"), index=False)
            print(f"  Saved GO/KEGG Enrichment to {dataset_id}_NMF_Factors_Enrichment.csv")

            print(f"  Plotting Enrichment Results for {dataset_id}...")
            for factor in df_enrichment['Factor'].unique():
                factor_df = df_enrichment[df_enrichment['Factor'] == factor].copy()
                if factor_df.empty:
                    continue

                factor_df['Clean_Term'] = factor_df['Term'].apply(lambda x: x.split(' (GO:')[0] if '(GO:' in x else x)
                factor_df = factor_df.sort_values('P-value', ascending=True).head(10)

                plt.figure(figsize=(10, max(4, len(factor_df) * 0.4)))
                sns.barplot(
                    x=-np.log10(factor_df['Adjusted P-value']),
                    y='Clean_Term',
                    hue='Database',
                    data=factor_df,
                    dodge=False
                )
                plt.title(f"{dataset_id} - {factor} Enrichment (Top 10)")
                plt.xlabel("-log10(Adjusted P-value)")
                plt.ylabel("")
                plt.tight_layout()
                plt.savefig(os.path.join(out_dir, f"{dataset_id}_{factor}_Enrichment.pdf"))
                plt.savefig(os.path.join(out_dir, f"{dataset_id}_{factor}_Enrichment.png"), dpi=300)
                plt.close()

        for factor in factor_cols:
            adata_subset.obs[factor] = cell_weights[factor].to_numpy()

        df_factors = adata_subset.obs[['celltype'] + factor_cols]
        mean_weights = df_factors.groupby('celltype')[factor_cols].mean()
        mean_weights_norm = mean_weights.div(mean_weights.sum(axis=1), axis=0)
        mean_weights.to_csv(os.path.join(out_dir, f"{dataset_id}_NMF_CellType_Weights.csv"))

        plt.figure(figsize=(10, 8))
        sns.heatmap(mean_weights_norm, cmap='Reds', annot=False)
        plt.title(f"{dataset_id} - NMF Factor Weights by Cell Type")
        plt.ylabel("Cell Type")
        plt.xlabel("NMF Factor")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{dataset_id}_NMF_CellType_Heatmap.pdf"))
        plt.savefig(os.path.join(out_dir, f"{dataset_id}_NMF_CellType_Heatmap.png"))
        plt.close()
        print(f"  Saved CellType-Factor Heatmap to {dataset_id}_NMF_CellType_Heatmap.png")

        if plot_spatial:
            print(f"  Plotting spatial overlay for factors...")
            for factor in factor_cols:
                fig, ax = plt.subplots(figsize=(6, 6))
                sc.pl.embedding(
                    adata_subset,
                    basis='spatial',
                    color=factor,
                    cmap='magma',
                    show=False,
                    ax=ax,
                    size=5,
                    frameon=False
                )
                ax.set_title(f"{dataset_id} - {factor.replace('_', ' ')} Spatial Distribution")
                ax.set_aspect('equal')
                plt.savefig(os.path.join(out_dir, f"{dataset_id}_{factor}_Spatial.pdf"), bbox_inches='tight')
                plt.savefig(os.path.join(out_dir, f"{dataset_id}_{factor}_Spatial.png"), bbox_inches='tight', dpi=300)
                plt.close()

    for sample_id in adata.obs['sample'].unique():
        adata_sub = adata[adata.obs['sample'] == sample_id].copy()
        process_one_dataset(sample_id, adata_sub, plot_spatial=True)

    print("\n--- Processing Combined All Samples ---")
    process_one_dataset("Combined_AllSamples", adata.copy(), plot_spatial=False)

if __name__ == "__main__":
    analyze_nmf_downstream()

