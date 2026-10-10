# 比较同一冻结状态在结直肠与回肠的 marker 区分度和供者覆盖；缺失 Glial 保留为 NA。
# 在本机 Python 环境运行，读取描述性比较表；生成可编辑图件和面板对齐记录。
# 队列和抽样不同，这张图不是准确率比较，也不做成对患者统计。
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Step1. 读取全状态表，顺序沿用模型谱系与状态编号
project = Path(__file__).resolve().parents[3]
root = project / 'IBD_ecotype_analysis/06_validation/GSE214695_scRNA_recovery_20261009'
data = pd.read_csv(root / '53状态_结肠回肠描述性比较.tsv', sep='\t')
types = ['B', 'CD4_T', 'CD8_T', 'Endothelial', 'Enteroendocrine', 'Fibroblast', 'Glial', 'Goblet', 'ILC', 'Mast', 'Monocyte_Macrophage', 'Pericyte']
data['TypeOrder'] = data.CellType.map(dict(zip(types, range(len(types)))))
data = data.sort_values(['TypeOrder', 'State'])
short = {'CD4_T': 'CD4 T', 'CD8_T': 'CD8 T', 'Enteroendocrine': 'Endocrine', 'Monocyte_Macrophage': 'Mono/Mac'}
labels = [short.get(t, t) + ' ' + s for t, s in zip(data.CellType, data.State)]
sys.path.insert(0, str(Path.home() / '.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 6, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.dpi': 600})

# Step2. 两个互补面板：自身 marker 优势与具有分配细胞的供者比例
width_mm = 183
fig, axes = plt.subplots(1, 2, figsize=(width_mm / 25.4, 220 / 25.4), sharey=True)
fig.subplots_adjust(left=0.20, right=0.985, bottom=0.10, top=0.87, wspace=0.26)
y = np.arange(len(data))
for ax in axes:
    ax.set_ylim(len(data) - 0.5, -0.5)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(axis='y', length=0)
    for boundary in np.flatnonzero(data.CellType.ne(data.CellType.shift()).to_numpy())[1:]:
        ax.axhline(boundary - 0.5, color='#DDDDDD', lw=0.5)
axes[0].set_yticks(y, labels)
for suffix, color, offset, display in [('Colon', '#0072B2', -0.16, 'Colon/rectum: GSE214695'), ('Ileum', '#E69F00', 0.16, 'Terminal ileum: IBDverse')]:
    contrast = data['OwnMinusBestOther_InternalOnly_' + suffix].to_numpy()
    donor_fraction = data[suffix + 'DonorSupportFraction'].to_numpy()
    axes[0].scatter(contrast, y + offset, s=9, color=color, label=display, edgecolor='none')
    axes[1].scatter(donor_fraction, y + offset, s=9, color=color, edgecolor='none')
axes[0].axvline(0, color='#777777', lw=0.7, linestyle='--')
axes[0].set_xlabel('Own marker mean Z minus best other state')
axes[1].set_xlim(-0.025, 1.025)
axes[1].set_xticks([0, 0.5, 1])
axes[1].set_xlabel('Supporting donors / donors with lineage')
axes[0].set_title('a  Internal marker separation', loc='left', pad=9)
axes[1].set_title('b  Donor coverage', loc='left', pad=9)
fig.text(0.20, 0.975, 'Same frozen states, independent intestinal cohorts', fontsize=9, fontweight='bold')
fig.text(0.20, 0.947, 'All 53 states shown; Glial unavailable in IBDverse. Missing values are not filled with zero.', fontsize=7)
fig.legend(*axes[0].get_legend_handles_labels(), loc='upper left', bbox_to_anchor=(0.20, 0.924), ncol=2, frameon=False, fontsize=6.5, borderaxespad=0)
fig.text(0.20, 0.02, 'Different cohorts and cell caps; marker separation is descriptive and internally derived, not accuracy.', fontsize=6.5)
fig.canvas.draw()
require_matplotlib_panel_alignment(fig, json_out=root / 'figure_qa/Fig6_alignment.json', panel_ids=['a', 'b'])
fig.savefig(root / 'figures/Fig6_colon_ileum_comparison.pdf')
fig.savefig(root / 'figures/Fig6_colon_ileum_comparison.svg')
fig.savefig(root / 'figures/Fig6_colon_ileum_comparison.png', dpi=600)
data.to_csv(root / 'figures/Fig6_source_data.tsv', sep='\t', index=False)
plt.close(fig)
print('Saved Fig6')
