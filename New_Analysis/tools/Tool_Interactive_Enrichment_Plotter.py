# ==============================================================================
# 文件名称: Tool_Interactive_Enrichment_Plotter.py
# 功能描述: 工具脚本: 交互式富集结果绘图 (Python)
# 联动说明: 独立工具，提供图形化界面 (GUI)，方便用户读取任意富集分析 CSV 文件并交互式生成可编辑的高清 PDF 图表。
# ==============================================================================

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.font import Font

# --- 导入全局 Nature 风格配置 ---
try:
    sys.path.append(os.getcwd())
    from analysis_config import NPG_COLORS, set_nature_style
    set_nature_style()
except:
    # Fallback 如果没找到 config 文件
    NPG_COLORS = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"]
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'sans-serif']

plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600

def clean_term_name(term):
    """清理通路名称，去除 GO ID 等后缀，并限制长度"""
    term = str(term)
    if '(GO:' in term:
        term = term.split(' (GO:')[0]
    if len(term) > 50:
        term = term[:47] + '...'
    return term

def plot_enrichment(df, term_col, pval_col, db_col, count_col, custom_title, out_path_prefix, plot_style):
    """核心绘图函数，支持多种 Nature 风格可视化"""
    if df.empty:
        messagebox.showinfo("Info", "没有可绘制的数据！")
        return

    df = df.copy()
    df['Clean_Term'] = df[term_col].apply(clean_term_name)
    df['Log_P'] = -np.log10(pd.to_numeric(df[pval_col], errors='coerce'))
    df = df.dropna(subset=['Log_P'])
    
    # 按照 P-value 排序 (越小排在越上面，即 Log_P 越大排在上面)
    # matplotlib 画图是从下往上画，所以这里按 Log_P 升序排列，这样最大的在最上面
    df = df.sort_values('Log_P', ascending=True)

    fig_height = max(4.0, len(df) * 0.35)
    fig, ax = plt.subplots(figsize=(8, fig_height))
    
    # =========================================================
    # 风格 1: 经典多色柱状图 (Classic Barplot by Database)
    # =========================================================
    if plot_style == "1. 经典多色柱状图 (按数据库分类)":
        if db_col and db_col in df.columns and db_col != "None":
            unique_dbs = df[db_col].unique()
            plot_palette = {db: NPG_COLORS[i % len(NPG_COLORS)] for i, db in enumerate(unique_dbs)}
            hue_col = db_col
        else:
            df['Database'] = 'Enrichment'
            hue_col = 'Database'
            plot_palette = {'Enrichment': NPG_COLORS[0]}
            
        sns.barplot(
            data=df, x='Log_P', y='Clean_Term', hue=hue_col, 
            dodge=False, palette=plot_palette, ax=ax, 
            edgecolor='black', linewidth=0.8
        )
        plt.legend(title='Database', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        
    # =========================================================
    # 风格 2: 渐变色柱状图 (Gradient Barplot)
    # =========================================================
    elif plot_style == "2. 渐变色柱状图 (按 P 值渐变)":
        norm = plt.Normalize(df['Log_P'].min(), df['Log_P'].max())
        sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=norm)
        sm.set_array([])
        
        ax.barh(df['Clean_Term'], df['Log_P'], color=sm.to_rgba(df['Log_P']), edgecolor='black', linewidth=0.8, height=0.7)
        cbar = fig.colorbar(sm, ax=ax, shrink=0.5, aspect=15)
        cbar.set_label(f'-Log10({pval_col})', rotation=270, labelpad=15)
        
    # =========================================================
    # 风格 3: 气泡图 (Dotplot / Bubble Plot)
    # =========================================================
    elif plot_style == "3. 气泡图 (Dotplot)":
        # 气泡大小映射
        if count_col and count_col in df.columns and count_col != "None":
            sizes = pd.to_numeric(df[count_col], errors='coerce').fillna(1) * 20
            size_label = count_col
        else:
            sizes = df['Log_P'] * 20
            size_label = "-Log10(P-value)"
            
        scatter = ax.scatter(
            df['Log_P'], df['Clean_Term'], 
            s=sizes, c=df['Log_P'], cmap="viridis", 
            alpha=0.9, edgecolors='black', linewidth=0.5, zorder=3
        )
        
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=15)
        cbar.set_label(f'-Log10({pval_col})', rotation=270, labelpad=15)
        
        # 添加网格线辅助视觉
        ax.grid(True, linestyle='--', alpha=0.4, axis='y', zorder=0)
        
        # 图例
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, num=4, func=lambda s: s/20)
        ax.legend(handles, labels, title=size_label, bbox_to_anchor=(1.25, 1), loc='upper left', frameon=False)

    # =========================================================
    # 风格 4: 棒棒糖图 (Lollipop Plot)
    # =========================================================
    elif plot_style == "4. 棒棒糖图 (Lollipop Plot)":
        # 画横线
        ax.hlines(y=df['Clean_Term'], xmin=0, xmax=df['Log_P'], color='gray', alpha=0.7, linewidth=2, zorder=1)
        # 画圆点
        scatter = ax.scatter(
            df['Log_P'], df['Clean_Term'], 
            s=120, c=df['Log_P'], cmap="coolwarm", 
            alpha=1.0, edgecolors='black', linewidth=0.8, zorder=2
        )
        
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=15)
        cbar.set_label(f'-Log10({pval_col})', rotation=270, labelpad=15)

    # =========================================================
    # 风格 5: 极坐标柱状图 (Circular Barplot)
    # =========================================================
    elif plot_style == "5. 极坐标柱状图 (Circular Barplot)":
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
        
        # 排序：对于极坐标，将值降序排列（大的在外，小的在内圈比较好看，或者顺时针排列）
        df = df.sort_values('Log_P', ascending=False).reset_index(drop=True)
        N = len(df)
        
        # 计算角度
        theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
        width = 2 * np.pi / N * 0.8
        radii = df['Log_P'].values
        
        # 颜色映射
        norm = plt.Normalize(radii.min(), radii.max())
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        colors = sm.to_rgba(radii)
        
        bars = ax.bar(theta, radii, width=width, bottom=radii.max()*0.2, color=colors, edgecolor='white', alpha=0.9)
        
        # 标签处理
        ax.set_xticks(theta)
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.spines['polar'].set_visible(False)
        
        # 添加带有旋转角度的文本标签
        for t, r, label in zip(theta, radii, df['Clean_Term']):
            rotation = np.degrees(t)
            alignment = 'left'
            # 保证文字不颠倒
            if rotation >= 90 and rotation < 270:
                rotation += 180
                alignment = 'right'
            
            ax.text(t, r + radii.max()*0.25, label, 
                    ha=alignment, va='center', rotation=rotation, rotation_mode='anchor', 
                    fontsize=9, color='#333333')
                    
        cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.1)
        cbar.set_label(f'-Log10({pval_col})')

    # =========================================================
    # 风格 6: 南丁格尔玫瑰图 (Rose Chart)
    # =========================================================
    elif plot_style == "6. 南丁格尔玫瑰图 (Rose Chart)":
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
        
        df = df.sort_values('Log_P', ascending=False).reset_index(drop=True)
        N = len(df)
        
        theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
        width = 2 * np.pi / N # 玫瑰图宽度填满
        radii = df['Log_P'].values
        
        norm = plt.Normalize(radii.min(), radii.max())
        sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
        sm.set_array([])
        colors = sm.to_rgba(radii)
        
        # 玫瑰图通常从圆心开始，宽度占满
        bars = ax.bar(theta, radii, width=width, bottom=0.0, color=colors, edgecolor='white', linewidth=1, alpha=0.85)
        
        ax.set_xticks(theta)
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.spines['polar'].set_visible(False)
        
        for t, r, label in zip(theta, radii, df['Clean_Term']):
            rotation = np.degrees(t)
            alignment = 'left'
            if rotation >= 90 and rotation < 270:
                rotation += 180
                alignment = 'right'
            
            # 文字稍往外偏移
            ax.text(t, r + radii.max()*0.05, label, 
                    ha=alignment, va='center', rotation=rotation, rotation_mode='anchor', 
                    fontsize=9, color='#333333')

        cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.1)
        cbar.set_label(f'-Log10({pval_col})')

    # =========================================================
    # 风格 7: 矩形树图 (Treemap)
    # =========================================================
    elif plot_style == "7. 矩形树图 (Treemap)":
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # 对于矩形树图，通常用面积代表基因数，颜色代表 P-value
        if count_col and count_col in df.columns and count_col != "None":
            sizes = pd.to_numeric(df[count_col], errors='coerce').fillna(1).values
        else:
            # 备用方案：用 P-value 代表大小
            sizes = df['Log_P'].values
            
        # 颜色映射
        norm = plt.Normalize(df['Log_P'].min(), df['Log_P'].max())
        sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)
        sm.set_array([])
        colors = sm.to_rgba(df['Log_P'])
        
        # 使用 squarify 绘制
        try:
            import squarify
            labels = [f"{t}\n({s})" for t, s in zip(df['Clean_Term'], np.round(df['Log_P'], 2))]
            squarify.plot(sizes=sizes, label=labels, color=colors, alpha=0.8, 
                          ax=ax, edgecolor="white", linewidth=2, text_kwargs={'fontsize':8, 'weight':'bold', 'color':'white'})
        except ImportError:
            messagebox.showerror("Error", "需要安装 squarify 库 (pip install squarify) 才能绘制矩形树图！")
            return
            
        ax.axis('off') # 去除坐标轴
        cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
        cbar.set_label(f'-Log10({pval_col})')

    # =========================================================
    # 风格 8: 放射状棒棒糖图 (Radial Lollipop)
    # =========================================================
    elif plot_style == "8. 放射状棒棒糖图 (Radial Lollipop)":
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
        
        df = df.sort_values('Log_P', ascending=False).reset_index(drop=True)
        N = len(df)
        
        theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
        radii = df['Log_P'].values
        
        # 颜色和大小映射
        norm = plt.Normalize(radii.min(), radii.max())
        sm = plt.cm.ScalarMappable(cmap="YlGnBu", norm=norm)
        sm.set_array([])
        colors = sm.to_rgba(radii)
        
        if count_col and count_col in df.columns and count_col != "None":
            sizes = pd.to_numeric(df[count_col], errors='coerce').fillna(1) * 30
        else:
            sizes = radii * 30
            
        # 绘制放射线段 (从底座 0.2 开始到 r)
        bottom = radii.max()*0.1
        for t, r in zip(theta, radii):
            ax.plot([t, t], [bottom, r + bottom], color='gray', linestyle='--', linewidth=1, alpha=0.6)
            
        # 绘制末端气泡
        ax.scatter(theta, radii + bottom, s=sizes, color=colors, alpha=0.9, edgecolor='black', zorder=10)
        
        # 画个中心的圈作为原点
        ax.fill_between(np.linspace(0, 2*np.pi, 100), 0, bottom, color='lightgray', alpha=0.3)
        
        ax.set_xticks(theta)
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.spines['polar'].set_visible(False)
        
        for t, r, label in zip(theta, radii, df['Clean_Term']):
            rotation = np.degrees(t)
            alignment = 'left'
            if rotation >= 90 and rotation < 270:
                rotation += 180
                alignment = 'right'
            
            ax.text(t, r + bottom + radii.max()*0.1, label, 
                    ha=alignment, va='center', rotation=rotation, rotation_mode='anchor', 
                    fontsize=9, color='#333333')

        cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.1)
        cbar.set_label(f'-Log10({pval_col})')

    # --- 全局修饰 ---
    ax.set_title(custom_title, fontsize=15, fontweight='bold', pad=20)
    ax.set_xlabel(f'-Log10({pval_col})', fontsize=13, fontweight='bold')
    ax.set_ylabel('')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    out_pdf = f"{out_path_prefix}.pdf"
    plt.savefig(out_pdf, bbox_inches='tight', dpi=600)
    plt.close()
    
    messagebox.showinfo("Success", f"绘图成功！\n已保存至:\n{out_pdf}")


class EnrichmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Universal Interactive Enrichment Plotter (Nature Style)")
        self.root.geometry("1000x950")
        self.root.minsize(900, 800)
        self.root.configure(bg="white")
        
        self.style = ttk.Style()
        self.style.theme_use('clam')
        
        self.base_font = Font(family="Microsoft YaHei", size=11)
        self.title_font = Font(family="Microsoft YaHei", size=12, weight="bold")
        self.btn_font = Font(family="Microsoft YaHei", size=12, weight="bold")
        
        self.style.configure('.', font=self.base_font, background="white", foreground="black")
        self.style.configure('TFrame', background="white")
        self.style.configure('TLabel', background="white")
        self.style.configure('TLabelframe', background="white")
        self.style.configure('TLabelframe.Label', font=self.title_font, foreground="#333333", background="white")
        self.style.configure('TRadiobutton', background="white")
        self.style.configure('TButton', font=self.base_font, padding=6)
        self.style.configure('TEntry', font=self.base_font)
        self.style.configure('TCombobox', font=self.base_font)
        self.style.configure('TNotebook', background="white")
        self.style.configure('TNotebook.Tab', background="white")
        
        # 让下拉列表本身的字体也能动态跟随
        self.root.option_add('*TCombobox*Listbox.font', self.base_font)
        
        self.default_dir = r"d:\Data_Analysis\空间转录组分析\ST_analysis\New_Analysis"
        self.csv_path = tk.StringVar()
        self.selected_sheet = None
        self.df = None
        
        self.main_container = ttk.Frame(self.root, padding="15")
        self.main_container.pack(fill="both", expand=True)
        
        self.plot_styles = [
            "1. 经典多色柱状图 (按数据库分类)",
            "2. 渐变色柱状图 (按 P 值渐变)",
            "3. 气泡图 (Dotplot)",
            "4. 棒棒糖图 (Lollipop Plot)",
            "5. 极坐标柱状图 (Circular Barplot)",
            "6. 南丁格尔玫瑰图 (Rose Chart)",
            "7. 矩形树图 (Treemap)",
            "8. 放射状棒棒糖图 (Radial Lollipop)"
        ]
        self.selected_style = tk.StringVar(value=self.plot_styles[0])
        
        self.create_widgets()
        self.draw_preview() # 初始化示意图
        
        # 绑定窗口缩放事件，实现动态字体
        self.last_width = 1000
        self.last_height = 950
        self.root.bind('<Configure>', self.on_resize)

    def on_resize(self, event):
        if event.widget == self.root:
            w, h = event.width, event.height
            if w == self.last_width and h == self.last_height:
                return
            self.last_width, self.last_height = w, h
            
            # 计算缩放比例 (以初始 1000x950 为基准)
            scale = min(w / 1000, h / 950)
            
            # 动态更新字体大小 (保持下限)
            self.base_font.config(size=max(9, int(11 * scale)))
            self.title_font.config(size=max(10, int(12 * scale)))
            self.btn_font.config(size=max(10, int(12 * scale)))

    def create_widgets(self):
        # 1. 文件选择
        frame_file = ttk.LabelFrame(self.main_container, text="1. 选择富集分析数据文件 (CSV / Excel)", padding="10")
        frame_file.pack(fill="x", pady=(0, 10))
        frame_file.columnconfigure(1, weight=1)
        
        ttk.Label(frame_file, text="文件路径:").grid(row=0, column=0, padx=(0, 10))
        ttk.Entry(frame_file, textvariable=self.csv_path, state='readonly', font=self.base_font).grid(row=0, column=1, sticky="ew", padx=(0, 10))
        ttk.Button(frame_file, text="📂 浏览...", command=self.browse_file).grid(row=0, column=2)

        # 2. 列名映射
        frame_cols = ttk.LabelFrame(self.main_container, text="2. 列名映射配置 (智能识别)", padding="10")
        frame_cols.pack(fill="x", pady=(0, 10))
        frame_cols.columnconfigure(1, weight=1)
        frame_cols.columnconfigure(3, weight=1)

        ttk.Label(frame_cols, text="Term (通路):").grid(row=0, column=0, sticky="e", pady=5)
        self.cb_term = ttk.Combobox(frame_cols, state="readonly", font=self.base_font)
        self.cb_term.grid(row=0, column=1, sticky="ew", padx=5)

        ttk.Label(frame_cols, text="P-value (显著性):").grid(row=0, column=2, sticky="e", pady=5)
        self.cb_pval = ttk.Combobox(frame_cols, state="readonly", font=self.base_font)
        self.cb_pval.grid(row=0, column=3, sticky="ew", padx=5)

        ttk.Label(frame_cols, text="Database (分组,可选):").grid(row=1, column=0, sticky="e", pady=5)
        self.cb_db = ttk.Combobox(frame_cols, state="readonly", font=self.base_font)
        self.cb_db.grid(row=1, column=1, sticky="ew", padx=5)
        
        ttk.Label(frame_cols, text="Count (气泡大小,可选):").grid(row=1, column=2, sticky="e", pady=5)
        self.cb_count = ttk.Combobox(frame_cols, state="readonly", font=self.base_font)
        self.cb_count.grid(row=1, column=3, sticky="ew", padx=5)
        
        tk.Button(frame_cols, text="🔄 加载数据", bg="#2196F3", fg="white", font=self.title_font, command=self.load_data).grid(row=2, column=0, columnspan=4, pady=10)

        # 3. 绘图风格选择与预览
        frame_style = ttk.LabelFrame(self.main_container, text="3. 选择 Nature 可视化风格", padding="10")
        frame_style.pack(fill="x", pady=(0, 10))
        
        style_left = ttk.Frame(frame_style)
        style_left.pack(side="left", fill="y", padx=(0, 20))
        
        for style in self.plot_styles:
            ttk.Radiobutton(style_left, text=style, value=style, variable=self.selected_style, command=self.draw_preview).pack(anchor="w", pady=5)
            
        style_right = ttk.Frame(frame_style)
        style_right.pack(side="right", fill="both", expand=True)
        ttk.Label(style_right, text="风格示意图 Preview:", font=self.title_font, foreground="#E64B35").pack(anchor="n")
        self.preview_canvas = tk.Canvas(style_right, width=300, height=150, bg="white", highlightthickness=1, highlightbackground="#cccccc")
        self.preview_canvas.pack(pady=5)

        # 4. 绘图模式选择 (Notebook)
        self.notebook = ttk.Notebook(self.main_container)
        self.notebook.pack(fill="both", expand=True, pady=(0, 10))

        # --- Tab A: Top N ---
        self.tab_top = ttk.Frame(self.notebook, padding="20")
        self.notebook.add(self.tab_top, text=" 📊 模式 A: 按 P-value 取 Top N ")
        
        ttk.Label(self.tab_top, text="显示前 N 个通路 (Top N):").grid(row=0, column=0, sticky="e", pady=10)
        self.entry_top_n = ttk.Entry(self.tab_top, font=self.base_font)
        self.entry_top_n.insert(0, "15")
        self.entry_top_n.grid(row=0, column=1, sticky="w", pady=10)

        ttk.Label(self.tab_top, text="P-value 过滤阈值:").grid(row=1, column=0, sticky="e", pady=10)
        self.entry_pval = ttk.Entry(self.tab_top, font=self.base_font)
        self.entry_pval.insert(0, "0.05")
        self.entry_pval.grid(row=1, column=1, sticky="w", pady=10)
        
        tk.Button(self.tab_top, text="🎨 生成图表 (Top N)", bg="#4CAF50", fg="white", font=self.btn_font, command=self.run_plot_top).grid(row=2, column=0, columnspan=2, pady=20, ipadx=20)

        # --- Tab B: Custom ---
        self.tab_custom = ttk.Frame(self.notebook, padding="15")
        self.notebook.add(self.tab_custom, text=" 🎯 模式 B: 自定义勾选通路 ")
        
        frame_search = ttk.Frame(self.tab_custom)
        frame_search.pack(fill="x", pady=(0, 5))
        ttk.Label(frame_search, text="🔍 搜索: ").pack(side="left")
        self.search_var = tk.StringVar()
        self.search_var.trace("w", self.filter_listbox)
        ttk.Entry(frame_search, textvariable=self.search_var, font=self.base_font).pack(side="left", fill="x", expand=True)
        
        frame_list = ttk.Frame(self.tab_custom)
        frame_list.pack(fill="both", expand=True)
        scrollbar = ttk.Scrollbar(frame_list)
        scrollbar.pack(side="right", fill="y")
        self.listbox = tk.Listbox(frame_list, selectmode=tk.MULTIPLE, yscrollcommand=scrollbar.set, font=self.base_font)
        self.listbox.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=self.listbox.yview)
        
        tk.Button(self.tab_custom, text="🎨 生成图表 (已勾选)", bg="#FF9800", fg="white", font=self.btn_font, command=self.run_plot_custom).pack(pady=10, ipadx=20)

        # 5. 输出设置
        frame_out = ttk.LabelFrame(self.main_container, text="4. 输出设置", padding="10")
        frame_out.pack(fill="x")
        ttk.Label(frame_out, text="自定义图表标题:").pack(side="left")
        self.entry_title = ttk.Entry(frame_out, width=50, font=self.base_font)
        self.entry_title.pack(side="left", padx=10)

    def draw_preview(self):
        """在 Canvas 上绘制对应风格的简单示意图"""
        c = self.preview_canvas
        c.delete("all")
        style = self.selected_style.get()
        
        # 绘制模拟的坐标轴 (非极坐标时使用)
        if "5" not in style and "6" not in style and "7" not in style and "8" not in style:
            c.create_line(40, 20, 40, 130, width=2, fill="#333") # Y轴
            c.create_line(40, 130, 280, 130, width=2, fill="#333") # X轴
            c.create_text(160, 145, text="-Log10(P-value)", font=("Arial", 9, "bold"), fill="#333")
            
        y_positions = [110, 80, 50, 20]
        lengths = [200, 150, 100, 60]
        
        if "1" in style: # Classic Barplot
            colors = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488"]
            for y, l, col in zip(y_positions, lengths, colors):
                c.create_rectangle(42, y-8, 42+l, y+8, fill=col, outline="black")
                c.create_text(35, y, text="Term", anchor="e", font=("Arial", 8))
                
        elif "2" in style: # Gradient Barplot
            colors = ["#800026", "#bd0026", "#fd8d3c", "#fed976"]
            for y, l, col in zip(y_positions, lengths, colors):
                c.create_rectangle(42, y-8, 42+l, y+8, fill=col, outline="black")
                c.create_text(35, y, text="Term", anchor="e", font=("Arial", 8))
                
        elif "3" in style: # Dotplot
            sizes = [12, 9, 6, 4]
            colors = ["#fde725", "#41b6c4", "#225ea8", "#081d58"]
            for y, l, s, col in zip(y_positions, lengths, sizes, colors):
                c.create_line(40, y, 280, y, dash=(2,2), fill="#ccc")
                c.create_oval(42+l-s, y-s, 42+l+s, y+s, fill=col, outline="black")
                c.create_text(35, y, text="Term", anchor="e", font=("Arial", 8))
                
        elif "4" in style: # Lollipop
            colors = ["#d73027", "#f46d43", "#74add1", "#4575b4"]
            for y, l, col in zip(y_positions, lengths, colors):
                c.create_line(42, y, 42+l, y, width=2, fill="gray")
                c.create_oval(42+l-6, y-6, 42+l+6, y+6, fill=col, outline="black")
                c.create_text(35, y, text="Term", anchor="e", font=("Arial", 8))
                
        elif "5" in style: # Circular Barplot
            c.create_oval(110, 35, 190, 115, outline="#ccc", dash=(2,2)) # 内圈
            c.create_arc(110, 35, 190, 115, start=0, extent=60, fill="#21908C", outline="white")
            c.create_arc(100, 25, 200, 125, start=60, extent=60, fill="#3B528B", outline="white")
            c.create_arc(90, 15, 210, 135, start=120, extent=60, fill="#5DC863", outline="white")
            c.create_arc(80, 5, 220, 145, start=180, extent=60, fill="#FDE725", outline="white")
            c.create_text(150, 75, text="Circular\nBarplot", font=("Arial", 8), fill="#333", justify="center")
            
        elif "6" in style: # Rose Chart
            c.create_arc(70, -5, 230, 155, start=0, extent=45, fill="#F0F921", outline="white")
            c.create_arc(90, 15, 210, 135, start=45, extent=45, fill="#FCA636", outline="white")
            c.create_arc(110, 35, 190, 115, start=90, extent=45, fill="#E16462", outline="white")
            c.create_arc(120, 45, 180, 105, start=135, extent=45, fill="#B12A90", outline="white")
            c.create_arc(130, 55, 170, 95, start=180, extent=45, fill="#6A00A8", outline="white")
            c.create_arc(135, 60, 165, 90, start=225, extent=45, fill="#0D0887", outline="white")
            c.create_text(150, 75, text="Rose", font=("Arial", 8, "bold"), fill="white")
            
        elif "7" in style: # Treemap
            c.create_rectangle(50, 20, 150, 130, fill="#9E0142", outline="white", width=2)
            c.create_rectangle(150, 20, 250, 80, fill="#D53E4F", outline="white", width=2)
            c.create_rectangle(150, 80, 210, 130, fill="#F46D43", outline="white", width=2)
            c.create_rectangle(210, 80, 250, 110, fill="#FDAE61", outline="white", width=2)
            c.create_rectangle(210, 110, 250, 130, fill="#FEE08B", outline="white", width=2)
            c.create_text(100, 75, text="Term A", fill="white", font=("Arial", 9, "bold"))
            c.create_text(200, 50, text="Term B", fill="white", font=("Arial", 8))
            
        elif "8" in style: # Radial Lollipop
            cx, cy = 150, 75
            c.create_oval(cx-10, cy-10, cx+10, cy+10, fill="#eee", outline="#ccc")
            angles = [0, 45, 90, 135, 180, 225, 270, 315]
            lengths = [60, 50, 40, 30, 20, 30, 40, 50]
            sizes = [10, 8, 6, 5, 4, 5, 6, 8]
            colors = ["#081d58", "#225ea8", "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1", "#ffffd9", "#41b6c4"]
            for a, l, s, col in zip(angles, lengths, sizes, colors):
                rad = np.radians(a)
                end_x = cx + l * np.cos(rad)
                end_y = cy - l * np.sin(rad)
                start_x = cx + 10 * np.cos(rad)
                start_y = cy - 10 * np.sin(rad)
                c.create_line(start_x, start_y, end_x, end_y, dash=(2,2), fill="gray")
                c.create_oval(end_x-s, end_y-s, end_x+s, end_y+s, fill=col, outline="black")

    def prompt_sheet_selection(self, file_path, sheet_names):
        top = tk.Toplevel(self.root)
        top.title("选择 Excel Sheet")
        top.geometry("350x150")
        top.transient(self.root)
        top.grab_set()
        top.configure(bg="white")
        
        ttk.Label(top, text="检测到多个 Sheet，请选择要读取的表：", font=self.base_font).pack(pady=15)
        cb = ttk.Combobox(top, values=sheet_names, state="readonly", font=self.base_font)
        cb.set(sheet_names[0])
        cb.pack(pady=5, fill="x", padx=30)
        
        def on_confirm():
            self.selected_sheet = cb.get()
            top.destroy()
            self.update_columns(file_path)
            
        ttk.Button(top, text="确定", command=on_confirm).pack(pady=15)

    def browse_file(self):
        filetypes = (("数据文件 (CSV/Excel)", "*.csv *.xlsx *.xls"), ("CSV 文件", "*.csv"), ("Excel 文件", "*.xlsx *.xls"), ("All", "*.*"))
        file_path = filedialog.askopenfilename(initialdir=self.default_dir, filetypes=filetypes)
        if file_path:
            self.csv_path.set(file_path)
            self.selected_sheet = None
            
            if file_path.lower().endswith(('.xlsx', '.xls')):
                try:
                    xl = pd.ExcelFile(file_path)
                    sheet_names = xl.sheet_names
                    if len(sheet_names) > 1:
                        self.prompt_sheet_selection(file_path, sheet_names)
                    else:
                        self.selected_sheet = sheet_names[0]
                        self.update_columns(file_path)
                except Exception as e:
                    messagebox.showerror("Error", f"读取 Excel 失败:\n{e}")
            else:
                self.update_columns(file_path)

    def update_columns(self, file_path):
        try:
            if file_path.lower().endswith(('.xlsx', '.xls')):
                cols = pd.read_excel(file_path, sheet_name=self.selected_sheet, nrows=0).columns.tolist()
            else:
                cols = pd.read_csv(file_path, nrows=0).columns.tolist()
                
            for cb in [self.cb_term, self.cb_pval, self.cb_db, self.cb_count]:
                cb.config(values=["None"] + cols)
            
            self.cb_term.set(next((c for c in cols if c.lower() in ['term', 'description', 'pathway']), cols[0] if cols else "None"))
            self.cb_pval.set(next((c for c in cols if 'adj' in c.lower() or 'padj' in c.lower() or 'p-value' in c.lower()), cols[1] if len(cols)>1 else "None"))
            self.cb_db.set(next((c for c in cols if c.lower() in ['database', 'source']), "None"))
            self.cb_count.set(next((c for c in cols if c.lower() in ['count', 'genecount']), "None"))
            
            self.df = None
            self.listbox.delete(0, tk.END)
        except Exception as e:
            messagebox.showerror("Error", f"读取列名失败:\n{e}")

    def load_data(self):
        if not self.csv_path.get(): return
        try:
            file_path = self.csv_path.get()
            if file_path.lower().endswith(('.xlsx', '.xls')):
                self.df = pd.read_excel(file_path, sheet_name=self.selected_sheet)
            else:
                self.df = pd.read_csv(file_path)
                
            self.all_terms = self.df[self.cb_term.get()].astype(str).unique().tolist()
            self.filter_listbox()
            sheet_info = f" (Sheet: {self.selected_sheet})" if self.selected_sheet else ""
            messagebox.showinfo("Success", f"加载成功，共 {len(self.df)} 条数据。{sheet_info}")
        except Exception as e:
            messagebox.showerror("Error", f"加载失败:\n{e}")

    def filter_listbox(self, *args):
        search = self.search_var.get().lower()
        self.listbox.delete(0, tk.END)
        for t in self.all_terms:
            if search in t.lower(): self.listbox.insert(tk.END, t)

    def _get_params(self):
        if self.df is None:
            messagebox.showwarning("Warning", "请先加载数据！")
            return None
        term, pval, db, count = self.cb_term.get(), self.cb_pval.get(), self.cb_db.get(), self.cb_count.get()
        title = self.entry_title.get() or "Enrichment Analysis"
        prefix = os.path.splitext(self.csv_path.get())[0] + "_Plot"
        return term, pval, db, count, title, prefix, self.selected_style.get()

    def run_plot_top(self):
        params = self._get_params()
        if not params: return
        term, pval, db, count, title, prefix, style = params
        
        try:
            top_n = int(self.entry_top_n.get())
            thresh = float(self.entry_pval.get())
            plot_df = self.df[pd.to_numeric(self.df[pval], errors='coerce') < thresh].copy()
            plot_df = plot_df.sort_values(pval).head(top_n)
            plot_enrichment(plot_df, term, pval, db, count, f"{title} (Top {top_n})", f"{prefix}_Top{top_n}", style)
        except Exception as e:
            messagebox.showerror("Error", f"绘图失败:\n{e}")

    def run_plot_custom(self):
        params = self._get_params()
        if not params: return
        term, pval, db, count, title, prefix, style = params
        
        sel = [self.listbox.get(i) for i in self.listbox.curselection()]
        if not sel: return
        plot_df = self.df[self.df[term].isin(sel)].copy().drop_duplicates(subset=[term])
        plot_enrichment(plot_df, term, pval, db, count, f"{title} (Custom)", f"{prefix}_Custom", style)

if __name__ == "__main__":
    root = tk.Tk()
    app = EnrichmentApp(root)
    try:
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    except: pass
    root.mainloop()
