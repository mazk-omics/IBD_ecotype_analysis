# 绘制实际两版n4的top10重合与逐状态留出结果，完整保留53状态，缺失不当成失败。
# 本机cibersortx_review环境运行；只读取核查通过的marker表，导出PDF/PNG/SVG和源表。
# 两张单面板图分别回答“排名改动多少”和“共同marker能否对应冻结状态”，不是细胞身份准确率。
from pathlib import Path
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path.home() / '.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment

# Step1. 使用固定53状态顺序，所有不可评价值保持NA
out = Path(__file__).resolve().parents[1]
figures = out / '08_figures'
states = pd.read_csv(out / '00_audit/frozen_states.tsv', sep='\t').sort_values(['Ecotype', 'CellType', 'State'])
order = states.CellState.to_list()
comparison = pd.read_csv(out / '02_external_markers/n4_disease_vs_region_comparison.tsv', sep='\t')
top10 = comparison[comparison.HeldOut.eq('All')].set_index('CellState').reindex(order)[['Top10Jaccard']]
top10.columns = ['Disease vs Region3']
loocv = pd.read_csv(out / '02_external_markers/marker_LOOCV_per_dataset.tsv', sep='\t')
loocv['DisplayConcordance'] = loocv.Concordant.where(loocv.Status.eq('EVALUABLE'))
concordance = loocv.pivot(index='CellState', columns=['HeldOut', 'Variant'], values='DisplayConcordance').reindex(order)
concordance = concordance.reindex(columns=pd.MultiIndex.from_product([['Colon_GSE214695', 'Ileum_IBDverse'], ['DiseaseGroup', 'Region3']]))
concordance.columns = ['Colon / Disease', 'Colon / Region3', 'Ileum / Disease', 'Ileum / Region3']
comparison.to_csv(figures / 'Figure_B_top10_source.tsv', sep='\t', index=False)
loocv.to_csv(figures / 'Figure_B_LOOCV_source.tsv', sep='\t', index=False)
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 7, 'legend.frameon': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False})

# Step2. 各图单独绘制，固定色标且图例明确“灰色不可评价”；不添加显著性星号
for name, matrix, title, label in [('Figure_B_top10_overlap', top10, 'n4 definition: top10 marker overlap', 'Top10 Jaccard'), ('Figure_B_LOOCV', concordance, 'Marker LOOCV across two external datasets', 'LOOCV concordance (0 or 1)')]:
    fig, ax = plt.subplots(figsize=(7.35, 9.0))
    fig.subplots_adjust(left=0.47, right=0.86, bottom=0.14, top=0.92)
    cmap = plt.get_cmap('Blues').copy()
    cmap.set_bad('#dddddd')
    image = ax.imshow(matrix.to_numpy(), cmap=cmap, vmin=0, vmax=1, aspect='auto', interpolation='none')
    ax.set_yticks(np.arange(53), labels=order, fontsize=6.4)
    ax.set_xticks(np.arange(matrix.shape[1]), labels=matrix.columns, rotation=45 if matrix.shape[1]>1 else 0, ha='right' if matrix.shape[1]>1 else 'center', rotation_mode='anchor', fontsize=6.5)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=8, pad=10)
    color_ax = fig.add_axes([0.90, 0.32, 0.016, 0.30])
    colorbar = fig.colorbar(image, cax=color_ax)
    colorbar.set_label(label, fontsize=7)
    if name.endswith('LOOCV'):
        colorbar.set_ticks([0, 1])
        colorbar.set_ticklabels(['0', '1'])
    fig.text(0.47, 0.025, 'Gray = not evaluable; LOOCV: 0 mismatch / 1 concordant\nTwo cohorts; conditional marker correspondence, not cell-identity accuracy', fontsize=6)
    require_matplotlib_panel_alignment(fig, axes=[ax], json_out=figures / (name+'.alignment.json'), tolerance_pt=1.5, strict=True)
    fig.savefig(figures / (name+'.pdf'))
    fig.savefig(figures / (name+'.svg'))
    fig.savefig(figures / (name+'.png'), dpi=350)
    plt.close(fig)
print('N4_COMPARISON_FIGURES_COMPLETE')
