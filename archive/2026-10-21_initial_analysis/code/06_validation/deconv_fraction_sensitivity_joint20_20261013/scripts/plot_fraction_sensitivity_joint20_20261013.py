# 用完整20类数据绘制真实bulk一致性和自然模拟准确性，报告图分页面保留全部观察值。
# 本机Python/matplotlib运行；读本轮数值表，输出PDF/SVG/PNG及最终绘图区对齐证据。
from pathlib import Path
import sys
import string
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Step1. 读取图件数据，固定单位、颜色与导出规范
project = Path(__file__).resolve().parents[2]
out = Path((project / '06_validation/fraction_sensitivity_current_path.txt').read_text(encoding='utf-8'))
figdir = out / '06_figures'
sys.path.insert(0, str(Path.home() / '.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
types = (out / '00_audit/frozen_types.txt').read_text().splitlines()
display = {t: t.replace('_', ' ') for t in types}
display.update({'Monocyte_Macrophage': 'Monocyte / macrophage', 'Stem_TA_progenitor': 'Stem / TA progenitor'})
metrics = pd.read_csv(out / '04_MSCCR120_concordance/type_concordance_metrics.tsv', sep='\t')
main = metrics.loc[metrics.group.eq('All120')].set_index('cell_type').loc[types]
paired = pd.read_csv(out / '04_MSCCR120_concordance/paired_fractions_long.tsv', sep='\t')
contracts = ['# Figure contract', '', 'Question: where do matched-reference fraction estimates agree, and where are material paired discrepancies?', 'Real-bulk figures measure concordance only. Natural-mixture figures use exact paired mixtures and frozen full-gene UMI truth.', 'Archetype: quantitative grids. Backend: Python exclusively. All 20 types and all requested observations included.', 'Exports: editable PDF/SVG plus 350-dpi PNG; 183-mm report pages; minimum text 6 pt; alignment gate before every export.', 'Scatter: sample-level pairing. Bland–Altman: signed differences, bias and descriptive 95% limits of agreement. Difference forest: patient-bootstrap 95% interval of bias and median absolute differences.', 'Simulation: development benchmark; no donor-independent confidence interval is claimed because cells/donors recur across mixtures.', 'No significance tests or p-value thresholds. No cross-method renormalization. Near-zero threshold is 0.01 percentage points.']
(figdir / 'figure_contract_and_legends.md').write_text('\n'.join(contracts) + '\n', encoding='utf-8')
export_figures = []

# Step2. 20类Pearson/Spearman与三个肠段分层热图
for filename, data, labels, title in [
    ('real_concordance_heatmap', main[['Pearson_r', 'Spearman_rho']].to_numpy(), ['Pearson r', 'Spearman rho'], 'Real-bulk fraction concordance | n = 120 patients'),
    ('regional_concordance_heatmap', np.column_stack([metrics.loc[metrics.group.eq('Region3:' + r)].set_index('cell_type').loc[types, 'Spearman_rho'] for r in ['Colon', 'Rectum', 'Ileum']]), ['Colon (n=40)', 'Rectum (n=40)', 'Ileum (n=40)'], 'Region-specific concordance | Spearman rho')]:
    fig, ax = plt.subplots(figsize=(183 / 25.4, 165 / 25.4))
    fig.subplots_adjust(left=0.34, right=0.78, bottom=0.10, top=0.90)
    image = ax.imshow(data, vmin=-1, vmax=1, cmap='RdBu_r', aspect='auto', interpolation='none')
    ax.set_yticks(range(20), [display[t] for t in types])
    ax.set_xticks(range(len(labels)), labels)
    ax.set_title(title, pad=14)
    cax = fig.add_axes([0.83, 0.24, 0.025, 0.48])
    fig.colorbar(image, cax=cax, label='Correlation')
    require_matplotlib_panel_alignment(fig, json_out=figdir / (filename + '.alignment.json'), exclude_axes=[cax], strict=True)
    export_figures.append((filename, fig))

# Step3. 类型偏差及配对绝对差异，不能用相关性替代数值一致性
fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 175 / 25.4), sharey=True)
fig.subplots_adjust(left=0.31, right=0.97, bottom=0.10, top=0.91, wspace=0.30)
positions = np.arange(20)
axes[0].errorbar(main.mean_difference_pp, positions, xerr=np.vstack([main.mean_difference_pp - main.mean_difference_pp_CI_low, main.mean_difference_pp_CI_high - main.mean_difference_pp]), fmt='o', color='#397D9A', markersize=3, linewidth=0.8)
axes[0].axvline(0, color='#A8A8A8', linestyle='--', linewidth=0.6)
axes[0].set_yticks(positions, [display[t] for t in types])
axes[0].invert_yaxis()
axes[0].set_xlabel('BayesPrism − CIBERSORTx (pp)')
axes[0].set_title('Mean difference | 95% patient CI', pad=15)
axes[1].plot(main.median_absolute_difference_pp, positions, 'o', color='#A47758', markersize=3)
axes[1].set_xlim(left=0)
axes[1].set_xlabel('Median absolute difference (pp)')
axes[1].set_title('Paired discrepancy', pad=15)
for i, ax in enumerate(axes):
    ax.annotate(string.ascii_lowercase[i], xy=(0, 1), xycoords='axes fraction', xytext=(-12, 8), textcoords='offset points', fontsize=8, weight='bold', annotation_clip=False)
require_matplotlib_panel_alignment(fig, json_out=figdir / 'real_paired_difference.alignment.json', strict=True)
export_figures.append(('real_paired_difference', fig))

# Step4. 每个类型散点和Bland–Altman，四页覆盖全部20类
for plot_kind in ['scatter', 'Bland_Altman']:
    for page in range(4):
        page_types = types[page * 5:(page + 1) * 5]
        fig, axes = plt.subplots(3, 2, figsize=(183 / 25.4, 220 / 25.4))
        fig.subplots_adjust(left=0.12, right=0.96, bottom=0.08, top=0.91, wspace=0.34, hspace=0.55)
        fig.delaxes(axes.flat[-1])
        for index, (cell_type, ax) in enumerate(zip(page_types, axes.flat)):
            block = paired.loc[paired.cell_type.eq(cell_type)]
            x = block.CIBERSORTx_Job20.to_numpy() * 100
            y = block.BayesPrism_final.to_numpy() * 100
            if plot_kind == 'scatter':
                bound = max(x.max(), y.max()) * 1.08
                ax.scatter(x, y, s=7, color='#397D9A', alpha=0.65, linewidths=0)
                ax.plot([0, bound], [0, bound], '--', color='#A8A8A8', linewidth=0.6)
                ax.set_xlim(0, bound)
                ax.set_ylim(0, bound)
                ax.set_xlabel('CIBERSORTx Job20 (%)')
                ax.set_ylabel('BayesPrism final (%)')
                ax.set_title(display[cell_type] + '\nr = %.2f; rho = %.2f' % (main.loc[cell_type, 'Pearson_r'], main.loc[cell_type, 'Spearman_rho']), pad=8)
            else:
                average = (x + y) / 2
                delta = y - x
                ax.scatter(average, delta, s=7, color='#397D9A', alpha=0.65, linewidths=0)
                ax.axhline(main.loc[cell_type, 'BA_bias_pp'], color='#A47758', linewidth=0.8)
                for name in ['BA_lower_LoA_pp', 'BA_upper_LoA_pp']:
                    ax.axhline(main.loc[cell_type, name], color='#A8A8A8', linestyle='--', linewidth=0.6)
                ax.set_xlabel('Mean of methods (%)')
                ax.set_ylabel('BayesPrism − CIBERSORTx (pp)')
                ax.set_title(display[cell_type], pad=8)
            ax.annotate(string.ascii_lowercase[index], xy=(0, 1), xycoords='axes fraction', xytext=(-13, 10), textcoords='offset points', fontsize=8, weight='bold', annotation_clip=False)
        filename = 'real_' + plot_kind + '_page' + str(page + 1)
        require_matplotlib_panel_alignment(fig, json_out=figdir / (filename + '.alignment.json'), strict=True)
        export_figures.append((filename, fig))

# Step5. 自然模拟误差，按配对差异回答准确性差别
if (out / '05_existing_pseudobulk_benchmark/paired_type_accuracy_difference.tsv').exists():
    accuracy = pd.read_csv(out / '05_existing_pseudobulk_benchmark/paired_type_accuracy_difference.tsv', sep='\t').set_index('cell_type').loc[types]
    fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 175 / 25.4), sharey=True)
    fig.subplots_adjust(left=0.31, right=0.97, bottom=0.11, top=0.89, wspace=0.32)
    axes[0].plot(accuracy.BP_MAE_pp, positions - 0.10, 'o', markersize=3, color='#397D9A', label='BayesPrism final')
    axes[0].plot(accuracy.CX_MAE_pp, positions + 0.10, 'o', markersize=3, color='#A47758', label='CIBERSORTx Job15')
    axes[0].set_yticks(positions, [display[t] for t in types])
    axes[0].invert_yaxis()
    axes[0].set_xlim(left=0)
    axes[0].set_xlabel('MAE against frozen UMI truth (pp)')
    axes[0].set_title('Natural mixtures | n = 28', pad=16)
    axes[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), fontsize=6)
    axes[1].plot(accuracy.mean_paired_absolute_error_difference_pp, positions, 'o', markersize=3, color='#397D9A')
    axes[1].axvline(0, color='#A8A8A8', linestyle='--', linewidth=0.6)
    axes[1].set_xlabel('Paired absolute error: BP − CX (pp)')
    axes[1].set_title('Negative: lower BayesPrism error', pad=16)
    for i, ax in enumerate(axes):
        ax.annotate(string.ascii_lowercase[i], xy=(0, 1), xycoords='axes fraction', xytext=(-12, 8), textcoords='offset points', fontsize=8, weight='bold', annotation_clip=False)
    require_matplotlib_panel_alignment(fig, json_out=figdir / 'natural28_paired_accuracy.alignment.json', strict=True)
    export_figures.append(('natural28_paired_accuracy', fig))

# Step6. 保留完整图页，图外说明列出统计量和重复单位
for filename, fig in export_figures:
    fig.savefig(figdir / (filename + '.pdf'))
    fig.savefig(figdir / (filename + '.svg'))
    fig.savefig(figdir / (filename + '.png'), dpi=350)
    plt.close(fig)
print('Rendered:', len(export_figures), 'PDF/SVG/PNG figures')
