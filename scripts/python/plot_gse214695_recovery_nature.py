# 按 Cell 2021 的状态 marker 验证思路展示独立 GSE214695 小规模 recovery，采用 nature-figure 的版式审计。
# 在本地 .tools/cibersortx_review 环境运行；读取官方汇总表与既有 GSE72819 审阅，输出 PDF/SVG/PNG 和逐图源数据。
# 图只总结既有分配和官方 Z 分数，不重新拟合模型；缺失值与未分配均明确保留。
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Step1. 读入结果和 nature-figure 面板对齐检查
project = Path(__file__).resolve().parents[3]
root = project / 'IBD_ecotype_analysis/06_validation/GSE214695_scRNA_recovery_20261009'
raw = root / 'raw_server'
audit = raw / 'audit'
figures = root / 'figures'
figures.mkdir(exist_ok=True)
qa = root / 'figure_qa'
qa.mkdir(exist_ok=True)
sys.path.insert(0, str(Path.home() / '.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 7, 'axes.linewidth': 0.6, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.dpi': 600, 'figure.facecolor': 'white'})
width_mm = 183  # 双栏出版宽度；每图高度随状态和基因数量变化
states = pd.read_csv(audit / 'state_recovery_summary.tsv', sep='\t')
cells = pd.read_csv(audit / 'cell_assignments_all.tsv', sep='\t')
donors = pd.read_csv(audit / 'donor_state_assignments.tsv', sep='\t')
markers = pd.read_csv(audit / 'top5_marker_expression.tsv', sep='\t')
bulk = pd.read_csv(project / 'IBD_ecotype_analysis/06_validation/GSE72819_recovery_review_20261009/53状态注释与恢复审阅.tsv', sep='\t')
type_order = ['B', 'CD4_T', 'CD8_T', 'Endothelial', 'Enteroendocrine', 'Fibroblast', 'Glial', 'Goblet', 'ILC', 'Mast', 'Monocyte_Macrophage', 'Pericyte']
short = {'CD4_T': 'CD4 T', 'CD8_T': 'CD8 T', 'Enteroendocrine': 'Endocrine', 'Monocyte_Macrophage': 'Mono/Mac'}
state_colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#56B4E9', '#D55E00', '#8A854B']
stratum_order = ['HC_nonIBD', 'CD_inflamed', 'UC_inflamed']
stratum_colors = ['#0072B2', '#009E73', '#D55E00']
states['TypeOrder'] = states.CellType.map(dict(zip(type_order, range(len(type_order)))))
states = states.sort_values(['TypeOrder', 'State']).reset_index(drop=True)
states = states.merge(bulk[['CellState', 'WinnerSamples', 'InitialProgram']], on='CellState', how='left', validate='one_to_one')
states.to_csv(figures / 'Fig1_source_data.tsv', sep='\t', index=False)
print('Figure inputs:', len(states), 'states;', len(cells), 'cells;', cells.Donor.nunique(), 'donors')

# Step2. 恢复总览：Z 分数、细胞支持和供者支持与 bulk 胜出数并列，分母分别标注
fig, axes = plt.subplots(1, 4, figsize=(width_mm / 25.4, 200 / 25.4), sharey=True)
fig.subplots_adjust(left=0.20, right=0.985, bottom=0.105, top=0.88, wspace=0.24)
y = np.arange(len(states))
labels = [short.get(t, t) + ' ' + s for t, s in zip(states.CellType, states.State)]
for ax in axes:
    ax.set_ylim(len(states) - 0.5, -0.5)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(axis='y', length=0)
    for boundary in np.flatnonzero(states.CellType.ne(states.CellType.shift()).to_numpy())[1:]:
        ax.axhline(boundary - 0.5, color='#DDDDDD', lw=0.5, zorder=0)
axes[0].set_yticks(y, labels)
z_values = states.Z.to_numpy()
finite = np.isfinite(z_values)
colors = np.where(z_values > 1.96, '#0072B2', '#B4B4B4')
axes[0].barh(y[finite], z_values[finite], height=0.72, color=colors[finite], edgecolor='none')
axes[0].axvline(1.96, color='#D55E00', lw=0.7, linestyle='--')
for row in np.flatnonzero(~finite):
    axes[0].text(np.nanmax(z_values) * 0.08, row, 'NA', color='#777777', va='center', fontsize=6)
axes[0].set_xlabel('Recovery Z\n(30 shuffles)')
axes[1].barh(y, states.AssignedCells, height=0.72, color='#56B4E9', edgecolor='none')
axes[1].set_xlabel('Assigned cells')
axes[2].barh(y, states.AssignedDonors, height=0.72, color='#009E73', edgecolor='none')
axes[2].set_xlim(0, 19)
axes[2].set_xticks([0, 6, 12, 18])
axes[2].set_xlabel('Supporting donors')
axes[3].barh(y, states.WinnerSamples, height=0.72, color='#888888', edgecolor='none')
axes[3].set_xlim(0, 75)
axes[3].set_xticks([0, 25, 50, 73])
axes[3].set_xlabel('Bulk winners /73')
for ax, letter, title in zip(axes, 'abcd', ['Official recovery Z', 'Single-cell support', 'Donor coverage', 'GSE72819 context']):
    ax.set_title(letter + '  ' + title, loc='left', pad=9)
fig.text(0.20, 0.972, 'Independent colon/rectum recovery | frozen Region3 CD+UC discovery', fontsize=9, fontweight='bold')
fig.text(0.20, 0.948, f'{cells.Donor.nunique()} donors; {len(cells):,} cells; {cells.CellType.nunique()} lineages. Dashed line: nominal Z = 1.96.', fontsize=7)
fig.text(0.20, 0.017, 'Assignment is not proof of validity. Cell-level Z values are exploratory; no multiple-testing claim.', fontsize=6.5)
fig.canvas.draw()
require_matplotlib_panel_alignment(fig, json_out=qa / 'Fig1_alignment.json', panel_ids=['a', 'b', 'c', 'd'])
fig.savefig(figures / 'Fig1_recovery_overview.pdf')
fig.savefig(figures / 'Fig1_recovery_overview.svg')
fig.savefig(figures / 'Fig1_recovery_overview.png', dpi=600)
plt.close(fig)

# Step3. 全部类型的冻结 marker 对照；每状态固定展示 MaxFC 排序前 3 个，不按 recovery 表现择优
cmap = plt.get_cmap('RdBu_r').copy()
cmap.set_bad('#DDDDDD')
for cell_type in type_order:
    info = markers[markers.CellType.eq(cell_type)].copy()
    state_ids = states.loc[states.CellType.eq(cell_type), 'CellState'].to_list()
    rows = info[['MarkerState', 'Gene']].drop_duplicates().groupby('MarkerState', sort=False).head(3)
    row_index = pd.MultiIndex.from_frame(rows)
    disc = info.pivot(index=['MarkerState', 'Gene'], columns='TargetState', values='DiscoveryMeanZ').reindex(index=row_index, columns=state_ids)
    rec = info.pivot(index=['MarkerState', 'Gene'], columns='TargetState', values='RecoveryMeanZ').reindex(index=row_index, columns=state_ids)
    plot_labels = [s.rsplit('_', 1)[1] + ' | ' + g for s, g in row_index]
    height_mm = max(95, 4.1 * len(rows) + 45)
    fig, axes = plt.subplots(1, 2, figsize=(width_mm / 25.4, height_mm / 25.4), sharey=True)
    fig.subplots_adjust(left=0.27, right=0.87, bottom=0.18, top=0.82, wspace=0.18)
    for ax, data, letter, title in zip(axes, [disc, rec], 'ab', ['Frozen discovery', 'GSE214695 recovery']):
        image = ax.imshow(data.to_numpy(), aspect='auto', cmap=cmap, vmin=-2, vmax=2, interpolation='nearest')
        ax.set_xticks(np.arange(len(state_ids)), [s.rsplit('_', 1)[1] for s in state_ids])
        ax.set_yticks(np.arange(len(rows)), plot_labels)
        ax.set_xlabel('Assigned state')
        ax.set_title(title, loc='left', pad=9)
        ax.text(-0.05, 1.10, letter, transform=ax.transAxes, fontweight='bold', fontsize=9)
        ax.tick_params(length=0)
        ax.spines[:].set_visible(False)
        for boundary in np.flatnonzero(rows.MarkerState.ne(rows.MarkerState.shift()).to_numpy())[1:]:
            ax.axhline(boundary - 0.5, color='white', lw=1)
    cax = fig.add_axes([0.90, 0.32, 0.018, 0.36])
    fig.colorbar(image, cax=cax, label='Mean gene Z', ticks=[-2, 0, 2])
    fig.text(0.27, 0.94, short.get(cell_type, cell_type) + ' | predefined marker programs', fontweight='bold', fontsize=9)
    fig.text(0.27, 0.065, 'Gene scaling is within lineage and dataset; gray columns have no assigned cells.', fontsize=6.5)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=qa / f'Fig2_{cell_type}_alignment.json', axes=axes, panel_ids=['a', 'b'], exclude_axes=[cax])
    fig.savefig(figures / f'Fig2_{cell_type}_marker_recovery.pdf')
    fig.savefig(figures / f'Fig2_{cell_type}_marker_recovery.svg')
    fig.savefig(figures / f'Fig2_{cell_type}_marker_recovery.png', dpi=600)
    info.to_csv(figures / f'Fig2_{cell_type}_source_data.tsv', sep='\t', index=False)
    plt.close(fig)

# Step4. 每位供者的条件状态比例；小规模抽样不用于估计活检整体细胞比例
wide = donors.pivot(index='CellState', columns='Donor', values='AssignedFraction')
donor_info = cells[['Donor', 'Stratum']].drop_duplicates().set_index('Donor')
ordered_donors = sorted(donor_info.index, key=lambda x: (stratum_order.index(donor_info.loc[x, 'Stratum']), x))
wide = wide.reindex(index=states.CellState, columns=ordered_donors)
wide.to_csv(figures / 'Fig3_donor_fraction_source_data.tsv', sep='\t')
fig, ax = plt.subplots(figsize=(width_mm / 25.4, 190 / 25.4))
fig.subplots_adjust(left=0.20, right=0.87, bottom=0.105, top=0.90)
fraction_cmap = plt.get_cmap('viridis').copy()
fraction_cmap.set_bad('#DDDDDD')
image = ax.imshow(wide, aspect='auto', vmin=0, vmax=1, cmap=fraction_cmap, interpolation='nearest')
ax.set_yticks(y, labels)
ax.set_xticks(np.arange(len(ordered_donors)), [f'D{i + 1:02}' for i in range(len(ordered_donors))], rotation=90, rotation_mode='anchor', ha='right', va='center')
ax.tick_params(length=0)
ax.spines[:].set_visible(False)
ax.set_xlabel('Donor (identities recorded in source data)')
for boundary in np.flatnonzero(states.CellType.ne(states.CellType.shift()).to_numpy())[1:]:
    ax.axhline(boundary - 0.5, color='white', lw=0.8)
for i, (name, color) in enumerate(zip(['HC / non-tumor', 'CD / inflamed', 'UC / inflamed'], stratum_colors)):
    ax.plot([i * 6 - 0.5, (i + 1) * 6 - 0.5], [-0.95, -0.95], color=color, lw=3, clip_on=False)
    ax.text((i + 0.5) * 6 - 0.5, -1.8, name, fontsize=6.5, ha='center', va='bottom')
cax = fig.add_axes([0.905, 0.31, 0.018, 0.36])
fig.colorbar(image, cax=cax, label='Fraction within lineage', ticks=[0, 0.5, 1])
fig.text(0.20, 0.965, 'State support across independent donors', fontweight='bold', fontsize=9)
fig.text(0.20, 0.025, 'At most 100 cells per donor and lineage. Gray: lineage absent; unassigned cells remain in denominators.', fontsize=6.5)
fig.canvas.draw()
require_matplotlib_panel_alignment(fig, json_out=qa / 'Fig3_alignment.json', axes=[ax], panel_ids=['a'], exclude_axes=[cax])
fig.savefig(figures / 'Fig3_donor_state_support.pdf')
fig.savefig(figures / 'Fig3_donor_state_support.svg')
fig.savefig(figures / 'Fig3_donor_state_support.png', dpi=600)
plt.close(fig)

# Step5. 预先指定谱系 marker 的表达及检出比例；部分基因也参与投影，不能称作完全独立验证
diagnostic = pd.read_csv(raw / 'input/diagnostic_marker_with_annotation.tsv', sep='\t', index_col=0)
gene_order = ['MS4A1', 'CD22', 'PAX5', 'CD3D', 'CD3E', 'TCF7', 'NKG7', 'TPSAB1', 'TPSB2', 'KIT', 'CPA3', 'CSF3R', 'CXCR2', 'S100A9', 'CSF1R', 'MS4A7', 'CHGA', 'CHGB', 'TPH1', 'MUC2', 'TFF3', 'MKI67', 'TOP2A', 'NOS2', 'DUOX2', 'LCN2', 'COL1A1', 'DCN', 'PECAM1', 'RGS5', 'SOX10', 'S100B', 'MPZ']
gene_order = [g for g in gene_order if g in diagnostic]
means = diagnostic.groupby(['CellType', 'Donor'])[gene_order].mean().groupby('CellType').mean().reindex(type_order)
detected = diagnostic[gene_order].gt(0).join(diagnostic[['CellType', 'Donor']]).groupby(['CellType', 'Donor']).mean().groupby('CellType').mean().reindex(type_order)
log_means = np.log2(means + 1)  # pseudocount=1；保留零表达，不删除未检出基因
relative = log_means.sub(log_means.mean(axis=0), axis=1).div(log_means.std(axis=0, ddof=1).replace(0, np.nan), axis=1).fillna(0)
means.to_csv(figures / 'Fig4_donor_balanced_marker_mean_CPM.tsv', sep='\t')
detected.to_csv(figures / 'Fig4_donor_balanced_marker_detection_fraction.tsv', sep='\t')
fig, ax = plt.subplots(figsize=(width_mm / 25.4, 100 / 25.4))
fig.subplots_adjust(left=0.19, right=0.88, bottom=0.24, top=0.83)
xx, yy = np.meshgrid(np.arange(len(gene_order)), np.arange(len(type_order)))
dots = ax.scatter(xx.ravel(), yy.ravel(), s=detected.to_numpy().ravel() * 28, c=relative.to_numpy().ravel(), cmap='RdBu_r', vmin=-2, vmax=2, edgecolor='none')
ax.set_yticks(np.arange(len(type_order)), [short.get(t, t) for t in type_order])
ax.set_xticks(np.arange(len(gene_order)), gene_order, rotation=90, rotation_mode='anchor', ha='right', va='center')
ax.set_ylim(len(type_order) - 0.5, -0.5)
ax.set_xlim(-0.5, len(gene_order) - 0.5)
ax.tick_params(length=0)
ax.spines[:].set_visible(False)
cax = fig.add_axes([0.91, 0.34, 0.016, 0.35])
fig.colorbar(dots, cax=cax, label='Relative expression Z', ticks=[-2, 0, 2])
fig.text(0.19, 0.94, 'Lineage marker checks | author-labeled cells', fontweight='bold', fontsize=9)
fig.text(0.19, 0.875, 'Color: log2(mean CPM + 1), standardized across lineages for each gene.', fontsize=6.5)
fig.text(0.19, 0.065, 'Dot area: detection fraction (0–100%). Both summaries give each contributing donor equal weight.', fontsize=6.5)
fig.canvas.draw()
require_matplotlib_panel_alignment(fig, json_out=qa / 'Fig4_alignment.json', axes=[ax], panel_ids=['a'], exclude_axes=[cax])
fig.savefig(figures / 'Fig4_lineage_marker_checks.pdf')
fig.savefig(figures / 'Fig4_lineage_marker_checks.svg')
fig.savefig(figures / 'Fig4_lineage_marker_checks.png', dpi=600)
plt.close(fig)

# Step6. 展示官方生态型丰度，不能据此声称完成了论文的生态型共现置换验证
eco = pd.read_csv(audit / 'ecotype_abundance.txt', sep='\t', index_col=0)
eco_assignment = pd.read_csv(audit / 'ecotype_assignment.txt', sep='\t', index_col=0)
sample_info = cells[['Sample', 'Donor', 'Stratum']].drop_duplicates().set_index('Sample')
sample_info['Ecotype'] = sample_info.index.map(eco_assignment.set_index('ID').Ecotype).fillna('Unassigned')
eco_order = [f'E{i}' for i in range(1, 8)]
sample_order = sorted(sample_info.index, key=lambda x: ((eco_order + ['Unassigned']).index(sample_info.loc[x, 'Ecotype']), x))
assert set(sample_info.index) == set(eco.columns), '官方丰度列必须覆盖全部活检样本'
eco = eco.reindex(index=eco_order, columns=sample_order)
eco.to_csv(figures / 'Fig5_ecotype_abundance_source_data.tsv', sep='\t')
sample_info.loc[sample_order].to_csv(figures / 'Fig5_sample_annotation.tsv', sep='\t')
fig, ax = plt.subplots(figsize=(width_mm / 25.4, 100 / 25.4))
fig.subplots_adjust(left=0.12, right=0.86, bottom=0.25, top=0.80)
image = ax.imshow(eco, aspect='auto', cmap=fraction_cmap, vmin=0, vmax=1, interpolation='nearest')
ax.set_yticks(np.arange(7), eco_order)
ax.set_xticks(np.arange(len(sample_order)), [f'D{ordered_donors.index(sample_info.loc[x, "Donor"]) + 1:02}' for x in sample_order], rotation=90, rotation_mode='anchor', ha='right', va='center')
ax.set_xlabel('Donor; ordered by official ecotype assignment')
ax.spines[:].set_visible(False)
ax.tick_params(length=0)
cax = fig.add_axes([0.90, 0.32, 0.018, 0.33])
fig.colorbar(image, cax=cax, label='Official ecotype abundance', ticks=[0, 0.5, 1])
fig.text(0.12, 0.94, 'Ecotype projection | all 12 model lineages', fontsize=9, fontweight='bold')
counts = sample_info.Ecotype.value_counts().reindex(eco_order + ['Unassigned'], fill_value=0)
fig.text(0.12, 0.87, 'Assigned samples: ' + '; '.join(f'{k}: {v}' for k, v in counts.items()), fontsize=6.5)
fig.text(0.12, 0.09, 'Projection is conditional on sampled lineages. This figure does not test multicellular co-occurrence.', fontsize=6.5)
fig.text(0.12, 0.04, 'No clinical association, no spatial validation, and no ecotype significance stars are inferred.', fontsize=6.5)
fig.canvas.draw()
require_matplotlib_panel_alignment(fig, json_out=qa / 'Fig5_alignment.json', axes=[ax], panel_ids=['a'], exclude_axes=[cax])
fig.savefig(figures / 'Fig5_ecotype_projection.pdf')
fig.savefig(figures / 'Fig5_ecotype_projection.svg')
fig.savefig(figures / 'Fig5_ecotype_projection.png', dpi=600)
plt.close(fig)
print('Saved recovery figures:', figures)
