import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial import cKDTree
import warnings
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

# ==============================================================================
# 全局随机种子设置 (Global Seed Setting)
# ==============================================================================
import random
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
sc.settings.seed = SEED
# ==============================================================================

# 配置路径
sys.path.append(os.getcwd())
try:
    from analysis_config import *
except ImportError:
    sys.path.append(r"d:\Data_Analysis\空间转录组分析\ST_analysis\New_Analysis")
    from analysis_config import *

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
warnings.filterwarnings("ignore")

def main():
    print("Step 13: Spatial Visualization of GC/DN Cell Communication (6 Modes)")
    
    # 1. 路径设置
    adata_path = PIPELINE_FILES["annotated_h5ad"]
    cc_path = os.path.join(PIPELINE_FILES["cellchat_dir"], "cellchat_inferred_network.csv")
    out_dir_base = PIPELINE_FILES["spatial_comm_viz_dir"]
    os.makedirs(out_dir_base, exist_ok=True)
    
    if not os.path.exists(adata_path) or not os.path.exists(cc_path):
        print(f"Error: Missing input files. Need:\n{adata_path}\n{cc_path}")
        return
        
    print("Loading AnnData...")
    adata = sc.read_h5ad(adata_path)
    
    print("Loading CellChat Network...")
    df_net = pd.read_csv(cc_path)
    
    # 获取原始表达量用于查询 Ligand 和 Receptor
    if adata.raw is not None:
        expr_adata = adata.raw.to_adata()
    else:
        expr_adata = adata
        
    samples = adata.obs['sample'].unique()
    
    # 预加载高清 H&E 图片以提高性能
    print("Pre-loading HD HE images...")
    import glob
    import json
    import matplotlib.image as mpimg
    from matplotlib.widgets import Slider, Button, TextBox
    he_images = {}
    for sample in samples:
        found_he = glob.glob(os.path.join(TSC_DATA_DIR, sample, "*HE*.png"))
        if found_he:
            print(f"  Loading {os.path.basename(found_he[0])} for {sample}...")
            try:
                he_images[sample] = mpimg.imread(found_he[0])
            except Exception as e:
                print(f"  Error loading {found_he[0]}: {e}")
        else:
            print(f"  No HD HE image found for {sample}")
            
    base_extent = [0, 55050, 0, 19906]
    
    # === 提取全局颜色映射 ===
    global_celltypes = adata.obs['celltype'].unique().tolist()
    global_color_map = {}
    
    if 'celltype_colors' in adata.uns:
        cats = adata.obs['celltype'].cat.categories
        colors = adata.uns['celltype_colors']
        for c, color in zip(cats, colors):
            global_color_map[c] = color
    else:
        # Fallback to NPG colors
        npg_colors = [
            '#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', 
            '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85',
            '#2E2A2B', '#E5A5C2', '#586E6E', '#C0C0C0', '#D4E2E0'
        ]
        for i, c in enumerate(global_celltypes):
            global_color_map[c] = npg_colors[i % len(npg_colors)]
            
    # 定义 6 种分析模式
    modes = [
        {'name': 'GC_as_Source', 'type': 'Giant Cells', 'role': 'source'},
        {'name': 'GC_as_Receiver', 'type': 'Giant Cells', 'role': 'target'},
        {'name': 'DN_Ex_as_Source', 'type': 'Dysmorphic Neurons(Excit)', 'role': 'source'},
        {'name': 'DN_Ex_as_Receiver', 'type': 'Dysmorphic Neurons(Excit)', 'role': 'target'},
        {'name': 'DN_In_as_Source', 'type': 'Dysmorphic Neurons(Inhib)', 'role': 'source'},
        {'name': 'DN_In_as_Receiver', 'type': 'Dysmorphic Neurons(Inhib)', 'role': 'target'}
    ]
    
    for mode in modes:
        mode_name = mode['name']
        focus_type = mode['type']
        role = mode['role']
        
        print(f"\n========================================================")
        print(f"Processing Mode: {mode_name} (Focus: {focus_type} as {role})")
        print(f"========================================================")
        
        mode_dir = os.path.join(out_dir_base, mode_name)
        os.makedirs(mode_dir, exist_ok=True)
        
        # 筛选通讯对
        if role == 'source':
            mask = (df_net['source'] == focus_type) & (df_net['target'] != focus_type)
            groupby_cols = ['source', 'target']
        else:
            mask = (df_net['target'] == focus_type) & (df_net['source'] != focus_type)
            groupby_cols = ['target', 'source']
            
        df_filtered = df_net[mask].copy()
        
        if df_filtered.empty:
            print(f"No significant interactions found for {mode_name}.")
            continue
            
        # 提取通讯概率最显著 (Top 1) 的一条通讯，并取总体 Prob 前 15
        idx_to_keep = df_filtered.groupby(groupby_cols)['prob'].idxmax()
        df_top = df_filtered.loc[idx_to_keep].sort_values('prob', ascending=False).head(15)
        
        for idx, row in df_top.iterrows():
            source = row['source']
            target = row['target']
            ligand = row['ligand']
            receptor = row['receptor']
            prob = row['prob']
            
            print(f"\n  Interaction: {source} ({ligand}) -> {target} ({receptor}) [Prob: {prob:.4f}]")
            
            if ligand not in expr_adata.var_names or receptor not in expr_adata.var_names:
                print(f"    Warning: {ligand} or {receptor} not in data. Skipping.")
                continue
                
            for sample in samples:
                sample_mask = adata.obs['sample'] == sample
                sample_indices = np.where(sample_mask)[0]
                
                adata_sub = adata[sample_indices].copy()
                expr_sub = expr_adata[sample_indices].copy()
                coords = adata_sub.obsm['spatial']
                
                l_expr = expr_sub[:, ligand].X
                r_expr = expr_sub[:, receptor].X
                
                if hasattr(l_expr, "todense"): l_expr = l_expr.todense()
                if hasattr(r_expr, "todense"): r_expr = r_expr.todense()
                l_expr = np.asarray(l_expr).flatten()
                r_expr = np.asarray(r_expr).flatten()
                
                source_mask = (adata_sub.obs['celltype'] == source).values & (l_expr > 0)
                source_idx = np.where(source_mask)[0]
                
                target_mask = (adata_sub.obs['celltype'] == target).values & (r_expr > 0)
                target_idx = np.where(target_mask)[0]
                
                if len(source_idx) == 0 or len(target_idx) == 0:
                    print(f"    [{sample}] Missing expressing source or target cells. Skipping.")
                    continue
                    
                distance_threshold = 100.0
                tree_target = cKDTree(coords[target_idx])
                neighbors = tree_target.query_ball_point(coords[source_idx], r=distance_threshold)
                
                pairs = []
                for i, n_list in enumerate(neighbors):
                    for j in n_list:
                        pairs.append((source_idx[i], target_idx[j]))
                        
                if len(pairs) == 0:
                    print(f"    [{sample}] No spatial pairs within {distance_threshold} um. Skipping.")
                    continue
                    
                # --- 寻找互作热点微环境进行局部放大展示 ---
                candidate_indices = list(set([s for s, t in pairs] + [t for s, t in pairs]))
                if len(candidate_indices) == 0:
                    continue
                    
                fov_radius = 150.0  # 300x300um square -> half-width is 150um
                
                # 计算每个 candidate cell 覆盖的 local pairs 数量
                cell_pair_counts = []
                for c_idx in candidate_indices:
                    cx, cy = coords[c_idx]
                    x_min, x_max = cx - fov_radius, cx + fov_radius
                    y_min, y_max = cy - fov_radius, cy + fov_radius
                    
                    local_count = 0
                    for s_i, t_i in pairs:
                        sx, sy = coords[s_i]
                        tx, ty = coords[t_i]
                        if (x_min <= sx <= x_max and y_min <= sy <= y_max) or \
                           (x_min <= tx <= x_max and y_min <= ty <= y_max):
                            local_count += 1
                            
                    cell_pair_counts.append((c_idx, local_count))
                    
                # 按照包含的通讯对数量降序排序
                cell_pair_counts.sort(key=lambda x: x[1], reverse=True)
                
                # 贪心算法选取不重叠（或少重叠）的 top FOVs
                selected_centers = []
                covered_cells = set()
                
                for c_idx, count in cell_pair_counts:
                    if c_idx in covered_cells:
                        continue
                        
                    cx, cy = coords[c_idx]
                    selected_centers.append(coords[c_idx])
                    
                    x_min, x_max = cx - fov_radius, cx + fov_radius
                    y_min, y_max = cy - fov_radius, cy + fov_radius
                    
                    for other_c in candidate_indices:
                        ox, oy = coords[other_c]
                        if x_min <= ox <= x_max and y_min <= oy <= y_max:
                            covered_cells.add(other_c)
                            
                # 如果是 GC 相关的模式，我们画出所有找到的有效 FOVs
                # 如果是 DN 相关的模式，FOVs 数量可能高达几百个，为了防止崩溃，最多只画前 5 个
                is_gc_mode = ('GC' in mode_name)
                max_fovs_to_draw = len(selected_centers) if is_gc_mode else min(5, len(selected_centers))
                
                selected_centers = selected_centers[:max_fovs_to_draw]
                print(f"    [{sample}] Found {len(cell_pair_counts)} candidate FOVs. Drawing top {len(selected_centers)} FOVs.")
                
                # 我们只在绘制第一个 FOV 时进行交互式对齐 (因为所有 FOV 共享同一张 HE 图像和坐标偏移)
                # 记录一次对齐结果，后续 FOV 直接复用
                
                for fov_idx, center_coord in enumerate(selected_centers):
                    print(f"    --- Drawing FOV {fov_idx+1}/{len(selected_centers)} ---")
                    
                    x_min, x_max = center_coord[0] - fov_radius, center_coord[0] + fov_radius
                    y_min, y_max = center_coord[1] - fov_radius, center_coord[1] + fov_radius
                    
                    local_pairs = []
                    for s_i, t_i in pairs:
                        sx, sy = coords[s_i]
                        tx, ty = coords[t_i]
                        if (x_min <= sx <= x_max and y_min <= sy <= y_max) or \
                           (x_min <= tx <= x_max and y_min <= ty <= y_max):
                            local_pairs.append((s_i, t_i))
                            
                    if len(local_pairs) == 0:
                        continue
                    
                    # ==============================================================================
                    # 局部微环境交互式微调 (Interactive Local Alignment)
                    # ==============================================================================
                    has_hd_image = sample in he_images
                    if has_hd_image:
                        img_hd = he_images[sample]
                    else:
                        img_hd = None
                
                    offset_file = os.path.join(out_dir_base, f"alignment_offset_{sample}_{ligand}_{receptor}_FOV{fov_idx+1}.json")
                    
                    local_offset = {'dx': 0.0, 'dy': 0.0}
                    if has_hd_image:
                        if os.path.exists(offset_file):
                            with open(offset_file, 'r') as f:
                                local_offset = json.load(f)
                            print(f"    Loaded saved alignment offset for {sample} ({ligand}-{receptor}) FOV{fov_idx+1}: {local_offset}")
                        else:
                            print(f"\n    [Interactive Local Alignment] Launching window for {sample} ({ligand}-{receptor}) FOV{fov_idx+1}...")
                            print("      --> Use 'Zoom (x)' slider to switch between 1x (Global), 10x (overview) and 40x (cell-level).")
                            print("      --> First align the macro tissue shape at 1x, then zoom in to align dots precisely.")
                            print("      --> Click 'Save & Continue' when done.")
                        
                            # 提升交互式窗口的 DPI 和画布大小，以最大化呈现高清图细节
                            fig_align, ax_align = plt.subplots(figsize=(14, 14), dpi=150)
                        # 给底部留出更多空间放三个滑块
                            plt.subplots_adjust(bottom=0.30)
                    
                            # 使用 antialiased 保证缩放时的边缘平滑度
                            img_obj = ax_align.imshow(img_hd, extent=base_extent, aspect='equal', interpolation='antialiased')
                    
                            # 在交互界面中，为了支持 1x 全局宏观缩放，我们需要画出所有细胞
                        
                            # 绘制用于对齐的局部散点 (Source 和 Target) - 只画当前 FOV 范围内的红色/蓝色目标点，保持焦点清晰
                            local_source_indices = list(set([s for s, t in local_pairs]))
                            local_target_indices = list(set([t for s, t in local_pairs]))
                            c_source = coords[local_source_indices]
                            c_target = coords[local_target_indices]
                            special_indices_align = set(local_source_indices + local_target_indices)
                    
                            # 绘制切片内的所有其他背景细胞
                            other_indices_align = [idx for idx in range(len(coords)) if idx not in special_indices_align]
                            if len(other_indices_align) > 0:
                                other_types_align = adata_sub.obs['celltype'].iloc[other_indices_align].values
                                other_coords_align = coords[other_indices_align]
                                unique_other_types_align = np.unique(other_types_align)
                                for ct in unique_other_types_align:
                                    ct_mask = other_types_align == ct
                                    # 背景点画小一点 (s=15)
                                    ax_align.scatter(other_coords_align[ct_mask, 0], other_coords_align[ct_mask, 1], 
                                                     c=global_color_map.get(ct, 'grey'), s=15, alpha=0.35, zorder=1, edgecolors='none', label=f'Env: {ct}')
                    
                            # 将当前 FOV 的 Source 和 Target 画在最上层，并且更大更醒目
                            ax_align.scatter(c_source[:, 0], c_source[:, 1], c='#DC143C', s=150, edgecolor='white', linewidth=1.5, zorder=4, label=f'Source: {source}')
                            ax_align.scatter(c_target[:, 0], c_target[:, 1], c='#00008B', s=150, edgecolor='white', linewidth=1.5, zorder=4, label=f'Target: {target}')
                    
                            # 在全局画布上画一个红框指示当前的 FOV
                            import matplotlib.patches as patches
                            fov_rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, 
                                                         linewidth=2, edgecolor='red', facecolor='none', linestyle='--', zorder=5)
                            ax_align.add_patch(fov_rect)
                    
                            # 初始视图设为全局宏观 (1X) 或者较低倍数 (例如 10x，即展示原本 FOV 的 4 倍范围)
                            init_zoom = 1.0  # 1X for global overview
                            def set_zoom_limits(zoom_level):
                                # zoom_level: 1.0 (Global), 10.0 to 40.0. 
                                # 40x means showing the exact fov_radius (150).
                                # 10x means showing 4 times the fov_radius (600).
                                # 1x (Global) means showing the whole tissue section
                                if zoom_level <= 1.0:
                                    ax_align.set_xlim(base_extent[0], base_extent[1])
                                    ax_align.set_ylim(base_extent[3], base_extent[2]) # y_max, y_min
                                else:
                                    current_radius = fov_radius * (40.0 / zoom_level)
                                    ax_align.set_xlim(center_coord[0] - current_radius, center_coord[0] + current_radius)
                                    ax_align.set_ylim(center_coord[1] + current_radius, center_coord[1] - current_radius)
                    
                            set_zoom_limits(init_zoom)
                    
                            ax_align.set_title(f"Align HE to Dots: {sample} ({ligand}-{receptor}) FOV{fov_idx+1}\nRed: {source} | Blue: {target}", fontweight='bold', fontsize=12)
                    
                            # 显示图例
                            ax_align.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), framealpha=0.8, fontsize=9)
                    
                            axcolor = 'lightgoldenrodyellow'
                            ax_zoom = plt.axes([0.15, 0.20, 0.65, 0.03], facecolor=axcolor)
                            ax_dx = plt.axes([0.15, 0.15, 0.65, 0.03], facecolor=axcolor)
                            ax_dy = plt.axes([0.15, 0.10, 0.65, 0.03], facecolor=axcolor)
                    
                            # 添加类似显微镜倍数的缩放滑块 (1x 到 40x)
                            s_zoom = Slider(ax_zoom, 'Magnification (x)', 1.0, 40.0, valinit=init_zoom, valstep=1.0)

                    
                            # 局部微调步长更精细，范围扩大至 -1000 到 1000
                            s_dx = Slider(ax_dx, 'HE X Offset', -1000.0, 1000.0, valinit=0.0, valstep=1.0)
                            s_dy = Slider(ax_dy, 'HE Y Offset', -1000.0, 1000.0, valinit=0.0, valstep=1.0)
                    
                            # ==========================================
                            # 添加自定义步长和方向按钮控制
                            # ==========================================
                            # 步长输入框
                            ax_step = plt.axes([0.15, 0.03, 0.1, 0.04])
                            text_step = TextBox(ax_step, 'Step Size: ', initial='10')
                    
                            # 方向按钮 (左、右、上、下)
                            ax_left = plt.axes([0.35, 0.03, 0.05, 0.04])
                            ax_right = plt.axes([0.45, 0.03, 0.05, 0.04])
                            ax_up = plt.axes([0.40, 0.05, 0.05, 0.04])
                            ax_down = plt.axes([0.40, 0.01, 0.05, 0.04])
                    
                            btn_left = Button(ax_left, '< Left')
                            btn_right = Button(ax_right, 'Right >')
                            btn_up = Button(ax_up, '^ Up')
                            btn_down = Button(ax_down, 'v Down')
                    
                            def get_step():
                                try:
                                    return float(text_step.text)
                                except ValueError:
                                    return 10.0
                            
                            def move_left(event):
                                # 如果想要图片往左走，我们需要减小它的 x 坐标
                                s_dx.set_val(s_dx.val - get_step())
                            def move_right(event):
                                # 如果想要图片往右走，我们需要增加它的 x 坐标
                                s_dx.set_val(s_dx.val + get_step())
                            def move_up(event):
                                # 注意：我们之前把 imshow 的 origin 设为了默认的 upper，即 0 在最上面，19906 在最下面
                                # 所以 Y 轴的数值是往下递增的。
                                # 如果想要图片在视觉上“往上走”，实际上是要减小它的 Y 坐标值。
                                s_dy.set_val(s_dy.val - get_step())
                            def move_down(event):
                                # 同理，图片“往下走”，Y 坐标值要增加。
                                s_dy.set_val(s_dy.val + get_step())
                        
                            btn_left.on_clicked(move_left)
                            btn_right.on_clicked(move_right)
                            btn_up.on_clicked(move_up)
                            btn_down.on_clicked(move_down)
                            # ==========================================
                    
                            offset_dict = {'dx': 0.0, 'dy': 0.0, 'saved': False}
                    
                            def update(val):
                                dx = s_dx.val
                                dy = s_dy.val
                                zoom = s_zoom.val
                        
                                # 更新底图位置
                                new_extent = [base_extent[0] + dx, base_extent[1] + dx,
                                              base_extent[2] + dy, base_extent[3] + dy]
                                img_obj.set_extent(new_extent)
                        
                                # 更新视野缩放 (Zoom)
                                set_zoom_limits(zoom)
                        
                                fig_align.canvas.draw_idle()
                        
                            s_dx.on_changed(update)
                            s_dy.on_changed(update)
                            s_zoom.on_changed(update)
                    
                            saveax = plt.axes([0.8, 0.05, 0.15, 0.04])
                            button = Button(saveax, 'Save & Continue', color='lightblue', hovercolor='0.975')
                    
                            def save_and_close(event):
                                offset_dict['dx'] = s_dx.val
                                offset_dict['dy'] = s_dy.val
                                offset_dict['saved'] = True
                                plt.close(fig_align)
                        
                            button.on_clicked(save_and_close)
                            plt.show()
                    
                            if offset_dict['saved']:
                                local_offset = {'dx': offset_dict['dx'], 'dy': offset_dict['dy']}
                                with open(offset_file, 'w') as f:
                                    json.dump(local_offset, f)
                                print(f"    Saved alignment offset: {local_offset}")
                            else:
                                print(f"    Alignment skipped. Using 0 offset.")
                            

                    current_extent = [base_extent[0] + local_offset['dx'], base_extent[1] + local_offset['dx'],
                                      base_extent[2] + local_offset['dy'], base_extent[3] + local_offset['dy']]

                    # --- 绘图 (全局切片 + 10x大视野纯HE + 局部放大带点) ---
                    # 将最终保存的主画布大小改为 3 列，分别展示 宏观全切片 -> 中观10X组织形态 -> 微观40X细胞互作
                    fig, axes = plt.subplots(1, 3, figsize=(36, 12), gridspec_kw={'width_ratios': [1, 1, 1]}, dpi=300)
                    ax_global = axes[0]
                    ax_mid = axes[1]
                    ax_local = axes[2]
                
                    # ----------------------------------------------------
                    # 图 1：全局空间切片 (Global)
                    # ----------------------------------------------------
                    if has_hd_image:
                        # 统一使用 'none' 或 'nearest' 确保极高分辨率下不被涂抹模糊
                        ax_global.imshow(img_hd, extent=current_extent, aspect='equal', alpha=1.0, interpolation='none')
                    else:
                        ax_global.set_facecolor('#1a1a1a')
                    
                    all_coords_img = coords
                
                    for ct in global_celltypes:
                        ct_mask = adata_sub.obs['celltype'] == ct
                        if np.sum(ct_mask) > 0:
                            ax_global.scatter(all_coords_img[ct_mask, 0], all_coords_img[ct_mask, 1], 
                                              c=global_color_map.get(ct, 'grey'), s=5, alpha=0.3, edgecolors='none')
                
                    circle_radius = fov_radius
                    center_x = center_coord[0]
                    center_y = center_coord[1]
                
                    circle = Circle((center_x, center_y), circle_radius, fill=False, edgecolor='red', linewidth=3, linestyle='--')
                    ax_global.add_patch(circle)
                
                    ax_global.set_title(f"Global Tissue Section\n({sample})", fontsize=16, fontweight='bold', fontname='Arial')
                    ax_global.axis('off')
                    ax_global.set_aspect('equal')
                
                    # ----------------------------------------------------
                    # 图 2：中观 10X 纯组织大视野 (Mid-view: 10X HE)
                    # ----------------------------------------------------
                    if has_hd_image:
                        ax_mid.imshow(img_hd, extent=current_extent, aspect='equal', alpha=1.0, interpolation='none')
                    else:
                        ax_mid.set_facecolor('#1a1a1a')
                    
                    # 10X 视野，相当于 fov_radius 的 4 倍
                    mid_radius = fov_radius * 4.0
                    mid_x_min = center_x - mid_radius
                    mid_x_max = center_x + mid_radius
                    mid_y_min = center_y - mid_radius
                    mid_y_max = center_y + mid_radius
                
                    ax_mid.set_xlim(mid_x_min, mid_x_max)
                    ax_mid.set_ylim(mid_y_max, mid_y_min)  # y_max 在前，维持自上而下
                
                    # 画一个红框指示真正的 40X 放大区域在哪里
                    import matplotlib.patches as patches
                    rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, 
                                             linewidth=3, edgecolor='red', facecolor='none', linestyle='--')
                    ax_mid.add_patch(rect)
                
                    # 添加 10X 视图的比例尺 (假设 100 像素 = 27微米，画个长一点的，比如 500 像素 = 135 微米)
                    mid_sb_len = 500.0
                    mid_sb_x2 = mid_x_max - 0.05 * (mid_x_max - mid_x_min)
                    mid_sb_x1 = mid_sb_x2 - mid_sb_len
                    mid_sb_y = mid_y_max - 0.05 * (mid_y_max - mid_y_min)
                
                    import matplotlib.patheffects as pe
                    ax_mid.plot([mid_sb_x1, mid_sb_x2], [mid_sb_y, mid_sb_y], color='white', linewidth=5, zorder=10,
                                path_effects=[pe.Stroke(linewidth=7, foreground='black'), pe.Normal()])
                
                    ax_mid.text((mid_sb_x1 + mid_sb_x2)/2, mid_sb_y - 0.02 * (mid_y_max - mid_y_min), r'135 $\mu$m', color='white', 
                                fontsize=15, fontweight='bold', ha='center', va='bottom', zorder=10,
                                path_effects=[pe.withStroke(linewidth=3, foreground="black")])
                
                    ax_mid.set_title(f"10X Overview (No Cells)\nNiche Locator", fontsize=16, fontweight='bold', fontname='Arial')
                    ax_mid.axis('off')
                    ax_mid.set_aspect('equal')

                    # ----------------------------------------------------
                    # 图 3：微观 40X 局部放大微环境带通讯散点 (Local Niche)
                    # ----------------------------------------------------
                    if has_hd_image:
                        # 移除半透明 alpha=0.5，改为 1.0；插值改为 'none' 保证极致清晰
                        ax_local.imshow(img_hd, extent=current_extent, aspect='equal', alpha=1.0, interpolation='none')
                    else:
                        ax_local.set_facecolor('#1a1a1a')
                
                    img_x_min, img_x_max = x_min, x_max
                    img_y_min, img_y_max = y_min, y_max
                
                    in_fov_mask = (all_coords_img[:, 0] >= img_x_min) & (all_coords_img[:, 0] <= img_x_max) & \
                                  (all_coords_img[:, 1] >= img_y_min) & (all_coords_img[:, 1] <= img_y_max)
                    in_fov_indices = np.where(in_fov_mask)[0]
                
                    local_source_indices = list(set([s for s, t in local_pairs]))
                    local_target_indices = list(set([t for s, t in local_pairs]))
                    special_indices = set(local_source_indices + local_target_indices)
                
                    source_color = '#DC143C'
                    target_color = '#00008B'
                    link_color = '#FFD700'
                
                    other_indices = [idx for idx in in_fov_indices if idx not in special_indices]
                    if len(other_indices) > 0:
                        other_types = adata_sub.obs['celltype'].iloc[other_indices].values
                        other_coords = all_coords_img[other_indices]
                        unique_other_types = np.unique(other_types)
                        for ct in unique_other_types:
                            ct_mask = other_types == ct
                            ax_local.scatter(other_coords[ct_mask, 0], other_coords[ct_mask, 1], 
                                             c=global_color_map.get(ct, 'grey'), s=80, alpha=0.35, zorder=1, label=f'Env: {ct}', edgecolors='none')
                    
                    c_source = coords[local_source_indices]
                    # 点画大一些，原来是 250，放大到 600
                    ax_local.scatter(c_source[:, 0], c_source[:, 1], c=source_color, s=600, edgecolor='white', linewidth=3.0, zorder=4)
                
                    c_target = coords[local_target_indices]
                    # 点画大一些，原来是 250，放大到 600
                    ax_local.scatter(c_target[:, 0], c_target[:, 1], c=target_color, s=600, edgecolor='white', linewidth=3.0, zorder=4)
                
                    lines = []
                    for s_i, t_i in local_pairs:
                        p1 = coords[s_i]
                        p2 = coords[t_i]
                        lines.append([p1, p2])
                    
                    # 线画粗一些，原来是 3.5，放大到 8.0
                    lc = LineCollection(lines, colors=link_color, linewidths=8.0, alpha=0.9, zorder=2)
                    ax_local.add_collection(lc)
                
                    # 绘制 Scale Bar (由于 100 像素实际上约等于 18~20 微米，修改标签以反映真实物理尺度)
                    scale_bar_len_plot = 100.0
                    # 放置在右下角：距离右边缘 5%，距离下边缘 5%
                    sb_x2 = img_x_max - 0.05 * (img_x_max - img_x_min)
                    sb_x1 = sb_x2 - scale_bar_len_plot
                    sb_y = img_y_min + 0.05 * (img_y_max - img_y_min)
                
                    # 绘制白色带黑边的线条作为 Scale Bar 增加可见度
                    import matplotlib.patheffects as pe
                    ax_local.plot([sb_x1, sb_x2], [sb_y, sb_y], color='white', linewidth=4, zorder=10,
                                  path_effects=[pe.Stroke(linewidth=6, foreground='black'), pe.Normal()])
                
                    text_y = sb_y + 0.02 * (img_y_max - img_y_min)
                    # 将原来的 100 \mu m 改为约等于其实际物理长度的 20 \mu m
                    ax_local.text((sb_x1 + sb_x2)/2, text_y, r'20 $\mu$m', color='white', 
                                  fontsize=14, fontweight='bold', ha='center', va='bottom', zorder=10,
                                  path_effects=[pe.withStroke(linewidth=3, foreground="black")])
                
                    ax_local.set_xlim(img_x_min, img_x_max)
                    # 使用正常的 y_min 到 y_max，不颠倒 Y 轴
                    ax_local.set_ylim(img_y_max, img_y_min)
                    ax_local.set_aspect('equal')
                
                    ax_local.set_title(f"Niche: {ligand}-{receptor}", fontsize=16, fontweight='bold', fontname='Arial')
                    ax_local.axis('off')
                
                    # 整理 Legend
                    custom_lines = []
                    custom_lines.append(Line2D([0], [0], color='none', label='Communication Pair'))
                    custom_lines.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=source_color, markeredgecolor='white', markersize=14, label=f'Source: {source}'))
                    custom_lines.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=target_color, markeredgecolor='white', markersize=14, label=f'Target: {target}'))
                    custom_lines.append(Line2D([0], [0], color=link_color, lw=3.0, label=f'{ligand} \u2192 {receptor}'))
                
                    custom_lines.append(Line2D([0], [0], color='none', label=' '))
                    custom_lines.append(Line2D([0], [0], color='none', label='Microenvironment Cells'))
                
                    if len(other_indices) > 0:
                        for ct in unique_other_types:
                            custom_lines.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=global_color_map.get(ct, 'grey'), markersize=12, alpha=0.7, label=ct))
                        
                    leg = ax_local.legend(handles=custom_lines, loc='center left', bbox_to_anchor=(1.02, 0.5), 
                                          frameon=False, fontsize=13, prop={'family': 'Arial', 'size': 13})
                
                    for text in leg.get_texts():
                        if text.get_text() in ['Communication Pair', 'Microenvironment Cells']:
                            text.set_fontweight('bold')
                
                    plt.tight_layout()
                
                    # 保存为样本名加通路名 (强制极高分辨率 600 DPI)
                    filename = f"{sample}_{ligand}_{receptor}"
                    if len(selected_centers) > 1:
                        filename += f"_FOV{fov_idx+1}"
                        
                    plt.savefig(os.path.join(mode_dir, f"{filename}.pdf"), bbox_inches='tight', dpi=600)
                    plt.savefig(os.path.join(mode_dir, f"{filename}.png"), bbox_inches='tight', dpi=600)
                    plt.close()
                
                    print(f"    [{sample}] Saved {filename} to {mode_name}")
                
                    # ==================================================================
                    # 绘制纯高清 H&E 切片 (不画散点和连线，单独保存)
                    # ==================================================================
                    if sample in he_images:
                        try:
                            he_dir = os.path.join(mode_dir, "HE_Views")
                            os.makedirs(he_dir, exist_ok=True)
                        
                            img_he = he_images[sample]
                        
                            # 增大画布尺寸并提高清晰度
                            fig_he, axes_he = plt.subplots(1, 2, figsize=(24, 12), gridspec_kw={'width_ratios': [1, 1]})
                            ax_g_he = axes_he[0]
                            ax_l_he = axes_he[1]
                        
                            # 1. 全景图
                            # 根据对齐验证测试 (蓝点对得准)：
                            # 使用默认的 origin (upper)，并且直接将 Y 坐标与图像对齐，无需计算任何翻转 (flipped_y)
                            img_extent = current_extent
                        
                            # 在全景图上画红框/红圈
                            ax_g_he.imshow(img_he, extent=img_extent, aspect='equal', interpolation='bicubic')
                            # 直接使用原始的 center_coord[1]，因为坐标系天然对齐
                            circle_he = Circle((center_coord[0], center_coord[1]), fov_radius, fill=False, edgecolor='red', linewidth=3, linestyle='--')
                            ax_g_he.add_patch(circle_he)
                            ax_g_he.set_title(f"Global HE Section\n({sample})", fontsize=16, fontweight='bold', fontname='Arial')
                            ax_g_he.axis('off')
                            ax_g_he.set_aspect('equal')
                        
                            # 2. 局部放大图
                            ax_l_he.imshow(img_he, extent=img_extent, aspect='equal', interpolation='bicubic')
                        
                            # 局部放大范围：直接使用原始的 x_min, x_max, y_min, y_max
                            ax_l_he.set_xlim(x_min, x_max)
                            # 因为图像默认 origin='upper' (0 在上面)，所以设置 ylim 时应该按照 (y_max, y_min) 的顺序
                            # 这样能让 0 保持在上面，与全景图和空转坐标系的方向严格一致
                            ax_l_he.set_ylim(y_max, y_min)
                        
                            # 添加 Scale Bar (反映约 20um 的真实物理尺度)
                            sb_x2_he = x_max - 0.05 * (x_max - x_min)
                            sb_x1_he = sb_x2_he - 100.0
                            # 因为 Y 轴大值在下面，所以我们要把比例尺放在 y_max (下方) 附近
                            sb_y_he = y_max - 0.05 * (y_max - y_min)
                        
                            import matplotlib.patheffects as pe
                            ax_l_he.plot([sb_x1_he, sb_x2_he], [sb_y_he, sb_y_he], color='white', linewidth=4, zorder=10,
                                          path_effects=[pe.Stroke(linewidth=6, foreground='black'), pe.Normal()])
                        
                            # 文字放在比例尺稍微靠上一点的位置 (即 Y 坐标更小一点)
                            text_y_he = sb_y_he - 0.02 * (y_max - y_min)
                            ax_l_he.text((sb_x1_he + sb_x2_he)/2, text_y_he, r'20 $\mu$m', color='white', 
                                          fontsize=14, fontweight='bold', ha='center', va='bottom', zorder=10,
                                          path_effects=[pe.withStroke(linewidth=3, foreground="black")])
                        
                            ax_l_he.set_title(f"HD HE Niche: {ligand}-{receptor}", fontsize=16, fontweight='bold', fontname='Arial')
                            ax_l_he.axis('off')
                            ax_l_he.set_aspect('equal')
                        
                            plt.tight_layout()
                        
                            # 仅保存为 PNG，使用 600 DPI 来保证极高清晰度
                            he_filename = f"{sample}_{ligand}_{receptor}_HE"
                            if len(selected_centers) > 1:
                                he_filename += f"_FOV{fov_idx+1}"
                                
                            plt.savefig(os.path.join(he_dir, f"{he_filename}.png"), bbox_inches='tight', dpi=600)
                            plt.close(fig_he)
                        
                            print(f"    [{sample}] Saved HD HE View to {mode_name}/HE_Views")
                        except Exception as e:
                            print(f"    [{sample}] Error generating HE view: {e}")
                
    print("\nLocal spatial visualization of communications finished!")

if __name__ == "__main__":
    main()

