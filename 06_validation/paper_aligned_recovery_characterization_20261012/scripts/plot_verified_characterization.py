# 展示已核验的恢复、共关联、肠段检验及Hallmark结果；缺失统计保持灰色，不填0。
# 本机cibersortx_review环境运行，绘图后另做PDF文字和碰撞审计；读统计表，写PDF/PNG/SVG及源表。
# 所有图均为BayesPrism路线的敏感性分析；Z支持不等于谱系命名正确。
from pathlib import Path
import json
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path.home() / '.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment

# Step1. 固定字体和冻结状态顺序，保留完整53状态；单图只表达一个证据层面
out = Path(__file__).resolve().parents[1]
figures = out / '08_figures'
states = pd.read_csv(out / '00_audit/frozen_states.tsv', sep='\t').sort_values(['Ecotype', 'CellType', 'State'])
order = states.CellState.to_list()
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.7})
datasets = ['Colon_GSE214695', 'Ileum_IBDverse']
recovery = pd.read_csv(out / '01_state_recovery/state_recovery_evidence_by_dataset.tsv', sep='\t')
recovery.to_csv(figures / 'Figure_A_source.tsv', sep='\t', index=False)
matrix_a = recovery.pivot(index='CellState', columns='Dataset', values='OfficialZ').reindex(index=order, columns=datasets)
region = pd.read_csv(out / '04_state_phenotypes/state_region_fisher.tsv', sep='\t')
matrix_c = np.log2(region.pivot(index='CellState', columns='Region3', values='OddsRatio').reindex(index=order, columns=['Ileum', 'Colon', 'Rectum'])).clip(-4, 4)
q_c = region.pivot(index='CellState', columns='Region3', values='Q').reindex(index=order, columns=matrix_c.columns)
region.to_csv(figures / 'Figure_C_source.tsv', sep='\t', index=False)
features = pd.read_csv(out / '06_ecotype_features/ecotype_continuous_features.tsv', sep='\t')
matrix_f = features.pivot(index='Feature', columns='Ecotype', values='MeanDifference').sort_index()
q_f = features.pivot(index='Feature', columns='Ecotype', values='Q').reindex(index=matrix_f.index, columns=matrix_f.columns)
features.to_csv(figures / 'Figure_F_source.tsv', sep='\t', index=False)

# Step2. 恢复Z、Region3、Hallmark各自一张热图；白色星号只来自实际BH q<0.05
feature_limit = np.nanmax(np.abs(matrix_f.to_numpy()))
heatmaps = [('Figure_A_recovery_Z', matrix_a, None, 'Official state recovery Z', 'Official Z; symmetric-log color; 30 shuffles', 'coolwarm', matrix_a.min().min(), matrix_a.max().max()), ('Figure_C_region_Fisher', matrix_c, q_c, 'Discovery Region3 representation', 'log2 odds ratio; display capped at +/-4', 'RdBu_r', -4, 4), ('Figure_F_Hallmark', matrix_f, q_f, 'Ecotype-associated bulk Hallmark scores', 'Mean score difference (ecotype - others)', 'RdBu_r', -feature_limit, feature_limit)]
for name, matrix, q, title, label, palette, vmin, vmax in heatmaps:
    fig, ax = plt.subplots(figsize=(7.35, 9.0))
    fig.subplots_adjust(left=0.47, right=0.87, bottom=0.11, top=0.92)
    cmap = plt.get_cmap(palette).copy()
    cmap.set_bad('#dddddd')
    if name.startswith('Figure_A'):
        image = ax.imshow(matrix.to_numpy(), aspect='auto', interpolation='none', cmap=cmap, norm=mpl.colors.SymLogNorm(linthresh=2, vmin=vmin, vmax=vmax))
    else:
        image = ax.imshow(matrix.to_numpy(), aspect='auto', interpolation='none', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yticks(np.arange(len(matrix)), labels=[x.replace('HALLMARK_', '') for x in matrix.index], fontsize=6.4)
    ax.set_xticks(np.arange(matrix.shape[1]), labels=['Colon (18 donors)', 'Ileum (30 donors)'] if name.startswith('Figure_A') else matrix.columns, rotation=45, ha='right', rotation_mode='anchor', fontsize=7)
    if name.startswith('Figure_A'):
        ax.set_xticks(np.arange(2), labels=['Colon (18 donors)', 'Ileum (30 donors)'], rotation=0, ha='center', fontsize=6.5)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=8, pad=10)
    if q is not None:
        for i, j in np.argwhere(q.to_numpy()<0.05):
            ax.text(j, i, '*', ha='center', va='center', fontsize=6, color='black')
    color_ax = fig.add_axes([0.90, 0.30, 0.016, 0.33])
    colorbar = fig.colorbar(image, cax=color_ax)
    colorbar.set_label(label, fontsize=7)
    if name.startswith('Figure_A'):
        colorbar.set_ticks([0, 2, 10, 30, 100])
        colorbar.set_ticklabels(['0', '2', '10', '30', '100'])
    fig.text(0.47, 0.025, 'BayesPrism sensitivity analysis; gray = not evaluable\n* BH q<0.05 where shown; biopsy/cell dependence remains', ha='left', fontsize=6)
    require_matplotlib_panel_alignment(fig, axes=[ax], json_out=figures / (name+'.alignment.json'), tolerance_pt=1.5, strict=True)
    fig.savefig(figures / (name+'.pdf'))
    fig.savefig(figures / (name+'.svg'))
    fig.savefig(figures / (name+'.png'), dpi=350)
    plt.close(fig)

# Step3. 共关联只展示已核验的四矩阵平均，保留所有undefined边；按冻结生态型排列
for dataset in datasets:
    c = pd.read_csv(out / '05_ecotype_coassociation' / (dataset+'_C_four_matrix_mean.tsv'), sep='\t', index_col=0)
    present = [state for state in order if state in c.index]
    c = c.loc[present, present]
    c.to_csv(figures / ('Figure_D_'+dataset+'_source.tsv'), sep='\t')
    fig, ax = plt.subplots(figsize=(8.0, 8.5))
    fig.subplots_adjust(left=0.21, right=0.87, bottom=0.24, top=0.90)
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#dddddd')
    image = ax.imshow(c.to_numpy(), cmap=cmap, vmin=-1, vmax=1, interpolation='none', aspect='auto')
    ax.set_xticks(np.arange(len(c)), labels=c.columns, rotation=90, rotation_mode='anchor', ha='right', va='center', fontsize=5.7)
    ax.set_yticks(np.arange(len(c)), labels=c.index, fontsize=5.7)
    ax.tick_params(length=0)
    ax.set_title(dataset+': fixed ecotype state order', fontsize=8)
    cumulative = 0
    for ecotype in states.Ecotype.unique():
        size = sum(states.loc[states.Ecotype.eq(ecotype), 'CellState'].isin(present))
        cumulative += size
        if cumulative<len(c):
            ax.axhline(cumulative-0.5, color='#888888', linewidth=0.5)
            ax.axvline(cumulative-0.5, color='#888888', linewidth=0.5)
    color_ax = fig.add_axes([0.895, 0.34, 0.015, 0.30])
    fig.colorbar(image, cax=color_ax).set_label('Mean of C1-C4 Spearman', fontsize=7)
    fig.text(0.21, 0.025, 'Gray = undefined; no zero imputation. Missing lineages are not recovery failures.', fontsize=6)
    name = 'Figure_D_'+dataset
    require_matplotlib_panel_alignment(fig, axes=[ax], json_out=figures / (name+'.alignment.json'), tolerance_pt=1.5, strict=True)
    fig.savefig(figures / (name+'.pdf'))
    fig.savefig(figures / (name+'.svg'))
    fig.savefig(figures / (name+'.png'), dpi=350)
    plt.close(fig)

# Step4. 共关联Z与零分布是生态型层面的证据，不沿用cell-state显著阈值来宣称生态型恢复
coassoc = pd.read_csv(out / '05_ecotype_coassociation/ecotype_coassociation_1000permutations.tsv', sep='\t')
coassoc.to_csv(figures / 'Figure_E_source.tsv', sep='\t', index=False)
fig, ax = plt.subplots(figsize=(5.8, 3.3))
fig.subplots_adjust(left=0.13, right=0.95, bottom=0.24, top=0.88)
for dataset, color, offset, label in zip(datasets, ['#407e9c', '#b37a52'], [-0.10, 0.10], ['Colon / GSE214695', 'Ileum / IBDverse']):
    values = coassoc[coassoc.Dataset.eq(dataset)].set_index('Ecotype').reindex(['E'+str(i) for i in range(1, 8)])
    ax.scatter(values.Z, np.arange(7)+offset, color=color, s=24, label=label)
ax.axvline(0, color='#999999', linewidth=0.6)
ax.set_yticks(np.arange(7), labels=['E'+str(i) for i in range(1, 8)])
ax.invert_yaxis()
ax.set_xlabel('Ecotype co-association Z (1,000 row permutations)')
ax.set_title('Frozen ecotype co-association in two external atlases', fontsize=8)
ax.legend(loc='lower right', fontsize=6.5)
fig.text(0.13, 0.035, 'P-value tail unresolved; no significance stars. Missing-member coverage is in source data.', fontsize=6)
name = 'Figure_E_ecotype_Z'
require_matplotlib_panel_alignment(fig, axes=[ax], json_out=figures / (name+'.alignment.json'), tolerance_pt=1.5, strict=True)
fig.savefig(figures / (name+'.pdf'))
fig.savefig(figures / (name+'.svg'))
fig.savefig(figures / (name+'.png'), dpi=350)
plt.close(fig)
(figures / 'figure_contract.json').write_text(json.dumps({'Backend': 'Python', 'Archetype': 'Quantitative single-panel evidence sequence', 'Claim': 'Frozen recovery assignments, marker validation, co-association and phenotype association are distinct evidence layers', 'FontFloor': 5, 'Gray': 'Missing or undefined, not zero', 'Statistics': 'Actual output tables only; no fabricated LR/global detection panels', 'Export': ['PDF', 'SVG', 'PNG'], 'PublicationStatus': 'Research report; not journal submission clearance'}, indent=2))
print('FIGURES_RENDERED: A, C, D(two atlases), E, F')
