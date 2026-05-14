import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

sys.path.append(os.getcwd())
from analysis_config import DIRS, PIPELINE_FILES

# ---------------------------------------------------------
# 配置
# ---------------------------------------------------------
NPG_COLORS = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"]

def get_grn_network(hub_file, corr_file, top_n=50, threshold=0.6):
    """
    复现 17 文件的逻辑：读取 Hub 基因和相关性矩阵，返回过滤孤立节点后的 DiGraph 和有效节点列表
    """
    if not os.path.exists(hub_file) or not os.path.exists(corr_file):
        print(f"Error: Missing {hub_file} or {corr_file}")
        sys.exit(1)
        
    df_hub = pd.read_csv(hub_file)
    top_genes = df_hub['Gene'].head(top_n).tolist()
    
    corr_matrix = pd.read_csv(corr_file, index_col=0)
    sub_corr = corr_matrix.loc[top_genes, top_genes]
    
    edges = np.where(np.abs(sub_corr) > threshold)
    
    G = nx.Graph()
    G.add_nodes_from(top_genes)
    
    for i in range(len(edges[0])):
        u = sub_corr.index[edges[0][i]]
        v = sub_corr.columns[edges[1][i]]
        if u != v:
            weight = sub_corr.iloc[edges[0][i], edges[1][i]]
            G.add_edge(u, v, weight=weight)
            
    # 剔除度数为 0 的孤立节点
    isolated_nodes = list(nx.isolates(G))
    G.remove_nodes_from(isolated_nodes)
    
    valid_genes = list(G.nodes())
    
    # 转换为有向图用于后续统一的绘制逻辑 (虽然共表达是无向的，但这里为了兼容后续绘图代码，使用 DiGraph)
    # 因为共表达是对称的，我们在绘图时只取一边，或者把原来的逻辑改写
    return G, valid_genes

def main():
    print("==================================================")
    print("21_Targeted_Network_Visualization (Using 17_GRN Data directly)")
    print("==================================================")
    
    out_dir = PIPELINE_FILES["targeted_scenic_dir"]
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. 28 个双重验证核心靶点 (DN Niche)
    dn_core_genes = [
        'EPHA3', 'EPHA5', 'EPHA6', 'EFNA5', 'PLXNA4', 'LAMB1', 'TNR',
        'GABRA5', 'GABRG3', 'SLC35F3', 
        'CNTN3', 'CNTN4', 'CNTN5', 'CNTNAP2', 'CNTNAP5', 
        'CDH9', 'CDH12', 'CDH18', 'CDH22',
        'ROBO1', 'ROBO2', 'LRP1B', 'LRP8', 'NCAM2'
    ]
    
    # 2. 从 17 文件读取并复原 GC 的网络和节点
    gc_hub_file = os.path.join(PIPELINE_FILES["grn_comparison_dir"], "GRN_Hub_Genes_GC_Niche.csv")
    gc_corr_file = os.path.join(PIPELINE_FILES["grn_comparison_dir"], "GC_Niche_Corr_Matrix.csv")
    
    # 也把 DN 的网络读取出来，用于验证那 28 个基因的连接
    dn_hub_file = os.path.join(PIPELINE_FILES["grn_comparison_dir"], "GRN_Hub_Genes_DN_Niche.csv")
    dn_corr_file = os.path.join(PIPELINE_FILES["grn_comparison_dir"], "DN_Niche_Corr_Matrix.csv")

    print("Loading GRN for GC Niche...")
    G_gc_full, gc_valid_genes = get_grn_network(gc_hub_file, gc_corr_file)
    print(f"GC Niche: Retained {len(gc_valid_genes)} valid connected genes from Top 50.")
    
    print("Loading GRN for DN Niche...")
    G_dn_full, dn_valid_genes = get_grn_network(dn_hub_file, dn_corr_file)
    
    # 筛选出 DN 中属于这 28 个 core genes 且有连接的节点
    valid_dn_targets = [g for g in dn_core_genes if g in dn_valid_genes]
    print(f"DN Niche: {len(valid_dn_targets)}/28 core genes are present in the top 50 connected GRN.")
    
    analysis_targets = {
        "Excitatory_DN_Core_Targets": (valid_dn_targets, G_dn_full),
        "Giant_Cells_GRN_Targets": (gc_valid_genes, G_gc_full)
    }
    
    # ---------------------------------------------------------
    # 循环绘制靶向网络图
    # ---------------------------------------------------------
    for label, (key_molecules, G_full) in analysis_targets.items():
        print(f"\nStarting analysis for {label}")
        
        if not key_molecules:
            print(f"No valid targets found for {label}. Skipping.")
            continue
            
        # 从全网络中提取只包含 target 及其直接邻居的子图
        # 如果只想画 targets 之间的网络，就直接 subgraph
        sub_G = G_full.subgraph(key_molecules).copy()
        
        # 剔除在 target 内部子图中没有连边的孤立节点
        isolated = list(nx.isolates(sub_G))
        sub_G.remove_nodes_from(isolated)
        
        if len(sub_G.nodes) == 0:
            print(f"No connections found among the targets in {label}. Skipping.")
            continue
            
        plt.figure(figsize=(10, 10), dpi=600)
        
        pos = nx.spring_layout(sub_G, seed=42, k=0.8)
        
        # 节点颜色分配与大小调整
        node_colors = []
        node_sizes = []
        
        for node in sub_G.nodes():
            node_colors.append(NPG_COLORS[0]) # 红色
            # 根据度数设置大小
            degree = sub_G.degree(node)
            node_sizes.append(1000 + degree * 200) 
                
        # 绘制节点
        nx.draw_networkx_nodes(sub_G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.9, edgecolors='black', linewidths=1.5)
        
        # 绘制标签
        nx.draw_networkx_labels(sub_G, pos, font_size=10, font_family="Arial", font_weight='bold', font_color='black')
        
        # 绘制连边
        edges = sub_G.edges()
        weights = [abs(sub_G[u][v]['weight']) * 3 for u, v in edges]
        nx.draw_networkx_edges(sub_G, pos, edgelist=edges, width=weights, edge_color="#A9A9A9", alpha=0.6)
        
        plt.title(f"Targeted Co-expression Network: {label}", fontsize=14, fontweight='bold', pad=20)
        plt.axis("off")
        
        output_pdf = os.path.join(out_dir, f"Targeted_Network_{label}.pdf")
        output_png = os.path.join(out_dir, f"Targeted_Network_{label}.png")
        plt.savefig(output_pdf, format="pdf", bbox_inches="tight")
        plt.savefig(output_png, format="png", bbox_inches="tight")
        plt.close()
        
        print(f"Saved {output_pdf}")

if __name__ == "__main__":
    main()
