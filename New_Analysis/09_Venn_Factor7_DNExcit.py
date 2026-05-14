import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from analysis_config import PIPELINE_FILES, get_step_dir

# --- 路径配置 ---
nmf_file = os.path.join(PIPELINE_FILES["nmf_downstream_dir"], "Combined_AllSamples_NMF_Top200_Sig_Genes.csv")
deg_file = os.path.join(PIPELINE_FILES["deg_dir"], "DE_DN_Excit_vs_Normal_Excit.csv")
out_dir = get_step_dir("09")

# 确保输出目录存在
os.makedirs(out_dir, exist_ok=True)

print("Loading NMF Factor 7 genes...")
df_nmf = pd.read_csv(nmf_file)
if 'Factor_7' in df_nmf.columns:
    # 提取 Factor_7 的基因并去除 NaN
    factor7_genes = set(df_nmf['Factor_7'].dropna().astype(str))
    print(f"Number of genes in Factor 7: {len(factor7_genes)}")
else:
    print("Error: 'Factor_7' not found in NMF file.")
    factor7_genes = set()

print("Loading Excitatory DN DEGs...")
df_deg = pd.read_csv(deg_file)
# 过滤条件：pvals_adj < 0.05
sig_degs = df_deg[df_deg['pvals_adj'] < 0.05].copy()
# 提取上调和下调的基因 (DN Excit 相对于 Normal Excit)
dn_up_genes = set(sig_degs[sig_degs['logfoldchanges'] > 0.5]['gene'].astype(str))
dn_down_genes = set(sig_degs[sig_degs['logfoldchanges'] < -0.5]['gene'].astype(str))
dn_all_sig_genes = dn_up_genes.union(dn_down_genes)

print(f"Number of Excitatory DN Upregulated genes: {len(dn_up_genes)}")
print(f"Number of Excitatory DN Downregulated genes: {len(dn_down_genes)}")
print(f"Number of All Significant DEGs: {len(dn_all_sig_genes)}")

# --- 求交集 ---
intersect_up = factor7_genes.intersection(dn_up_genes)
intersect_down = factor7_genes.intersection(dn_down_genes)
intersect_all = factor7_genes.intersection(dn_all_sig_genes)

print(f"Intersection with Up-regulated: {len(intersect_up)}")
print(f"Intersection with Down-regulated: {len(intersect_down)}")

# 保存交集基因到文本文件
out_txt = os.path.join(out_dir, "Overlap_Factor7_and_DNExcit_DEGs.txt")
with open(out_txt, "w") as f:
    f.write("Intersection between NMF Factor 7 and Excitatory DN Upregulated DEGs:\n")
    f.write(", ".join(sorted(list(intersect_up))) + "\n\n")
    f.write("Intersection between NMF Factor 7 and Excitatory DN Downregulated DEGs:\n")
    f.write(", ".join(sorted(list(intersect_down))) + "\n")
print(f"Saved intersection gene list to {out_txt}")

# --- 绘制 Venn 图 ---
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'sans-serif']
plt.figure(figsize=(10, 5))

# 1. Factor 7 vs DN Excit Up
plt.subplot(1, 2, 1)
v1 = venn2([factor7_genes, dn_up_genes], set_labels=('NMF Factor 7', 'DN(Excit) Up-DEGs'), set_colors=('#4DBBD5', '#E64B35'))
plt.title("Factor 7 vs DN Upregulated")

# 2. Factor 7 vs DN Excit Down
plt.subplot(1, 2, 2)
v2 = venn2([factor7_genes, dn_down_genes], set_labels=('NMF Factor 7', 'DN(Excit) Down-DEGs'), set_colors=('#4DBBD5', '#00A087'))
plt.title("Factor 7 vs DN Downregulated")

plt.tight_layout()
out_pdf = os.path.join(out_dir, "Venn_Factor7_vs_DNExcit.pdf")
plt.savefig(out_pdf, bbox_inches='tight', dpi=600)
plt.close()
print(f"Saved Venn diagram to {out_pdf}")
