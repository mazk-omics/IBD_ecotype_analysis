# 参照Cell2021按细胞类型×数据集展示验证表现，同时补充绝对误差和RNA分母的影响。
# Python/matplotlib运行；只画已有官方CIBERSORTx结果，保留全部20行、不可计算和未验证状态。
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Step1. 图件契约：趋势恢复不等于准确校准，三个面板回答不同问题
out = Path(__file__).resolve().parents[1]
root = out.parent
types = (root / '00_audit/frozen_types.txt').read_text().splitlines()
table = pd.read_csv(out / 'EXTERNAL_RECOMPUTED_METRICS.tsv', sep='\t')
table = table.loc[table.method.eq('CIBERSORTx') & table.cell_type.isin(types)]
figdir = out / 'figures'
sys.path.insert(0, str(Path.home() / '.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'legend.frameon': False})
contract = {'question': 'Which frozen Joint20 CIBERSORTx types track known external composition, and where does calibration or coverage fail?', 'claim': 'Recovery is type- and estimand-dependent; whole-matrix or correlation-only claims are insufficient.', 'archetype': 'quantitative grid', 'backend': 'Python only', 'panel_a': 'Cell-count Pearson across mixtures: original-style type/dataset recovery unit, 34 eligible pairs, 33 finite correlations and 1 all-zero failure', 'panel_b': 'Cell-count MAE: calibration evidence including high-r Plasma bias', 'panel_c': 'Full-gene RNA MAE: explicit estimand comparison rather than exchanging truths silently', 'points': '100 designed mixtures per cohort, 1000 draws each, cells and donors reused; no biological CI or hypothesis test', 'coverage': 'All 20 rows kept; NE means no positive truth; IBDverse ILC Pearson NA because all predictions zero; no failure removed', 'exports': '183 x 185 mm; vector PDF/SVG plus 350 dpi PNG preview; final glyphs at least 6 pt', 'original_comparison': 'Cell2021 STAR Methods and Figure S1B assessment unit verified from cached full text; original supplemental page layout not re-rendered or claimed as an exact facsimile'}
(figdir / 'FIGURE_CONTRACT.json').write_text(json.dumps(contract, indent=2), encoding='utf-8')

# Step2. 全20类共用坐标行，用颜色和形状区分两个固定外部数据集
fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 185 / 25.4), sharey=True)
fig.subplots_adjust(left=0.29, right=0.98, top=0.86, bottom=0.20, wspace=0.40)
positions = np.arange(len(types))
display = [x.replace('_', ' ') for x in types]
panels = [('Cell_count_SECONDARY', 'Pearson_r', 'Composition tracking', 'Cell-count Pearson r'), ('Cell_count_SECONDARY', 'MAE_pp', 'Cell-count calibration', 'Cell-count MAE (pp)'), ('Full_gene_RNA_SECONDARY', 'MAE_pp', 'RNA-contribution calibration', 'Full-gene RNA MAE (pp)')]
for ax, (scale, metric, title, xlabel) in zip(axes, panels):
    for name, color, marker, offset in [('GSE214695', '#397D9A', 'o', -0.13), ('IBDverse', '#A47758', 's', 0.13)]:
        block = table.loc[table.dataset.eq(name) & table.truth_scale.eq(scale)].set_index('cell_type').loc[types]
        valid = block.evaluable & block[metric].notna()
        ax.scatter(block.loc[valid, metric], positions[valid] + offset, color=color, marker=marker, s=12, linewidths=0, label=name, zorder=3)
        if metric == 'Pearson_r':
            missing = block.evaluable & block[metric].isna()
            for y in positions[missing]:
                ax.text(0.03, y + offset, 'NA', color=color, fontsize=6, ha='left', va='center')
    ax.set_title(title, fontsize=7, pad=16)
    ax.set_xlabel(xlabel)
    ax.set_ylim(len(types) - 0.4, -0.6)
    ax.grid(axis='x', color='#E8E8E8', linewidth=0.5, zorder=0)
axes[0].set_yticks(positions, display)
axes[0].set_xlim(-0.06, 1.06)
axes[0].set_xticks([0, 0.5, 1])
axes[1].set_xlim(-0.7, 30)
axes[1].set_xticks([0, 10, 20, 30])
axes[2].set_xlim(-0.15, 6)
axes[2].set_xticks([0, 2, 4, 6])
paneth = types.index('Paneth')
for letter, ax in zip('abc', axes):
    ax.annotate(letter, xy=(0, 1), xycoords='axes fraction', xytext=(-8, 8), textcoords='offset points', fontsize=8, weight='bold', annotation_clip=False)
    ax.text(0.08, paneth, 'NE', transform=ax.get_yaxis_transform(), fontsize=6, color='#777777', va='center')
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.62, 0.965), ncol=2, fontsize=7)
fig.text(0.29, 0.105, '100 mixtures per cohort; 1000 cell draws per mixture. No prediction renormalization.', fontsize=6)
fig.text(0.29, 0.078, 'NE: no positive external truth. NA: IBDverse ILC predictions are all zero.', fontsize=6)
fig.text(0.29, 0.051, '34 eligible type/cohort pairs; all 20 parent types retained. Descriptive statistics; no biological CI.', fontsize=6)
require_matplotlib_panel_alignment(fig, json_out=figdir / 'SIGNATURE_VALIDATION.alignment.json', strict=True)
fig.savefig(figdir / 'SIGNATURE_VALIDATION.pdf')
fig.savefig(figdir / 'SIGNATURE_VALIDATION.svg')
fig.savefig(figdir / 'SIGNATURE_VALIDATION.png', dpi=350)
plt.close(fig)
print('One compact, three-panel external signature validation figure exported.')
