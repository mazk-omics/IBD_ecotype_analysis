# 每个外部数据集仅一页核心验证图，区分RNA贡献和cell-count误差，保留全部可评价类型。
# 本机Python/matplotlib运行；读实际输出表，导出PDF/PNG并保留渲染对齐证据。
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Step1. 读取两套外部验证的真实数值
project = Path(__file__).resolve().parents[2]
out = Path((project / '06_validation/fraction_sensitivity_current_path.txt').read_text(encoding='utf-8'))
sys.path.insert(0, str(Path.home() / '.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
qa_notes = ['# 外部验证图件契约', '', '每个数据集一页：左图各类型MAE揭示RNA分母与细胞数估计量的差异；右图Pearson/Spearman回答跨混合物排序恢复程度。', '100份模拟，每份1000次原始细胞抽样。参考固定Full20；图显示全部作者可评价类型，未对预测比例重归一化。', '描述统计，无误差区间或独立生物学重复推断；100个混合物复用来源细胞/供者。Unsupported-type质量另见数值表。', '采用183mm定宽、白背景、最小6pt字体、编辑性PDF/SVG；不制造模拟图示数据。']
for dataset in ['GSE214695', 'IBDverse']:
    folder = out / '08_external_pseudobulk' / dataset
    table = pd.read_csv(folder / '04_accuracy/accuracy_metrics_all_scales.tsv', sep='\t')
    table = table.loc[table.method.eq('BayesPrism_final') & table.evaluable & ~table.cell_type.str.startswith('POOLED')]
    RNA = table.loc[table.truth_scale.eq('Retained_gene_RNA_PRIMARY')].set_index('cell_type')
    cell = table.loc[table.truth_scale.eq('Cell_count_SECONDARY')].set_index('cell_type').loc[RNA.index]
    RNA = RNA.sort_values('MAE_pp', ascending=False)
    cell = cell.loc[RNA.index]
    positions = np.arange(len(RNA))
    fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 160 / 25.4), sharey=True)
    fig.subplots_adjust(left=0.30, right=0.97, top=0.90, bottom=0.16, wspace=0.38)
    axes[0].plot(RNA.MAE_pp, positions - 0.10, 'o', color='#397D9A', markersize=3, label='Retained-gene RNA truth')
    axes[0].plot(cell.MAE_pp, positions + 0.10, 'o', color='#A47758', markersize=3, label='Cell-count truth')
    axes[0].set_yticks(positions, [t.replace('_', ' ') for t in RNA.index])
    axes[0].invert_yaxis()
    axes[0].set_xlim(left=0)
    axes[0].set_xlabel('Mean absolute error (pp)')
    axes[0].set_title(dataset + ' | final theta', fontsize=8, pad=15)
    axes[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), fontsize=6, frameon=False)
    axes[1].plot(RNA.Pearson_r, positions - 0.10, 'o', color='#397D9A', markersize=3, label='Pearson r')
    axes[1].plot(RNA.Spearman_rho, positions + 0.10, 'o', color='#A47758', markersize=3, label='Spearman rho')
    axes[1].set_xlim(-1.05, 1.05)
    axes[1].set_xticks([-1, -0.5, 0, 0.5, 1])
    axes[1].set_xlabel('Correlation against RNA truth')
    axes[1].set_title('100 mixtures × 1000 draws', fontsize=8, pad=15)
    axes[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), fontsize=6, frameon=False)
    for letter, ax in zip(['a', 'b'], axes):
        ax.annotate(letter, xy=(0, 1), xycoords='axes fraction', xytext=(-12, 8), textcoords='offset points', fontsize=8, weight='bold', annotation_clip=False)
    filename = dataset + '_external_accuracy'
    require_matplotlib_panel_alignment(fig, json_out=out / '06_figures' / (filename + '.alignment.json'), strict=True)
    fig.savefig(out / '06_figures' / (filename + '.pdf'))
    fig.savefig(out / '06_figures' / (filename + '.svg'))
    fig.savefig(out / '06_figures' / (filename + '.png'), dpi=350)
    plt.close(fig)
(out / '06_figures/external_figure_contract.md').write_text('\n'.join(qa_notes) + '\n', encoding='utf-8')
print('External figures complete: two pages')
