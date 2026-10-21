# 在本轮结果目录的Python环境绘图，输入全部来自实际统计和患者评分。
# 输出PNG/PDF/SVG与源数据；RISK输入溯源保留，不把关联写成预测性能。
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator, NullFormatter
from matplotlib.colors import TwoSlopeNorm

# Step1. 固定显示约定，统计家庭和效应编码不随绘图改变。
out = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(out / 'scripts/figure_qa'))
from audit_panel_alignment import require_matplotlib_panel_alignment
figdir = out / '04_cross_cohort_summary/figures'
figdir.mkdir(exist_ok=True)
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'], 'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .6, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'legend.frameon': False})
data = pd.read_csv(out / '04_cross_cohort_summary/all_associations.tsv', sep='\t')
membership = pd.read_csv(out / 'audit/frozen_state_ecotype_membership.tsv', sep='\t')
endpoints = ['ETR_W10_clinical_remission', 'IFX_W4_6_endoscopic_healing', 'VDZ_W6_endoscopic_healing', 'VDZ_W12_endoscopic_healing', 'VDZ_W52_endoscopic_healing', 'RISK_3yr_B2_B3']
titles = {'ETR_W10_clinical_remission': 'ETR · Week 10 remission', 'IFX_W4_6_endoscopic_healing': 'IFX · Week 4–6 healing', 'VDZ_W6_endoscopic_healing': 'VDZ · Week 6 healing', 'VDZ_W12_endoscopic_healing': 'VDZ · Week 12 healing', 'VDZ_W52_endoscopic_healing': 'VDZ · Week 52 healing', 'RISK_3yr_B2_B3': 'RISK* · 3-year progression'}
blue = '#286A8D'
neutral = '#8B9399'

# Step2. 每个终点显示全部7生态型，使用完整区间而非只显示显著结果。
for stem, visits in [('Fig1_UC_primary_ecotype_associations', endpoints[:3]), ('Fig2_progression_and_later_VDZ', [endpoints[5], endpoints[3], endpoints[4]]), ('Fig_S1_GSE72819_ecotypes', [endpoints[0]]), ('Fig_S2_RISK_ecotypes', [endpoints[5]]), ('Fig_S3_GSE73661_ecotypes', endpoints[1:5])]:
    targetdir = figdir
    if len(visits) == 1:
        fig, axis = plt.subplots(figsize=(3.54, 2.8))
        axes = [axis]
        fig.subplots_adjust(left=.16, right=.96, bottom=.24, top=.76)
        targetdir = out / ('01_GSE72819_treatment_outcome' if '72819' in stem else '02_RISK_disease_progression') / 'figures'
    elif len(visits) == 4:
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
        axes = axes.ravel()
        fig.subplots_adjust(left=.07, right=.97, bottom=.13, top=.86, wspace=.4, hspace=1.0)
        targetdir = out / '03_GSE73661_mucosal_healing/figures'
    else:
        fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8))
        fig.subplots_adjust(left=.065, right=.98, bottom=.24, top=.76, wspace=.45)
    targetdir.mkdir(exist_ok=True)
    source = []
    for panel, (ax, visit) in enumerate(zip(axes, visits)):
        d = data.loc[data.Endpoint.eq(visit) & data.Family.eq('Ecotype')].set_index('Feature').reindex(['E' + str(i) for i in range(1, 8)])
        source.append(d.reset_index())
        both = d.Q_Firth.lt(.05) & d.Q_Wilcoxon.lt(.05)
        for i, (_, row) in enumerate(d.iterrows()):
            rank_only = row.Q_Wilcoxon < .05 and row.Q_Firth >= .05
            color = blue if both.iloc[i] else ('#B6782E' if rank_only else neutral)
            ax.errorbar(row.OR, i, xerr=[[row.OR - row.CI_low], [row.CI_high - row.OR]], fmt='D' if rank_only else 'o', ms=3.2, color=color, elinewidth=.9, capsize=2, capthick=.7, zorder=3)
        ax.axvline(1, color='#C8CDD0', lw=.7, zorder=0)
        ax.set_xscale('log')
        ax.set_yticks(range(7), d.index)
        ax.set_ylim(6.6, -.6)
        lo = d.CI_low.min() / 1.4
        hi = d.CI_high.max() * 1.4
        ax.set_xlim(lo, hi)
        assert (d.CI_low > 0).all() and (d.OR > 0).all()
        ticks = [10.0 ** i for i in range(-5, 4)] if hi / lo > 30 else [.25, .5, 1, 2, 4, 8]
        ticks = [tick for tick in ticks if lo < tick < hi]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel('OR per 1 SD (log scale)', labelpad=5)
        n = int(d.N.iloc[0])
        events = int(d.Events.iloc[0])
        ax.set_title(titles[visit] + f'\n{events} events / {n} patients', pad=11)
        ax.annotate(chr(97 + panel), xy=(0, 1), xycoords='axes fraction', xytext=(-14, 38), textcoords='offset points', fontsize=8, fontweight='bold', annotation_clip=False)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=str(targetdir / (stem + '.alignment.json')), overlay_svg=str(targetdir / (stem + '.alignment.svg')), strict=True, tolerance_pt=1.5)
    fig.savefig(targetdir / (stem + '.png'), dpi=600)
    fig.savefig(targetdir / (stem + '.pdf'))
    fig.savefig(targetdir / (stem + '.svg'))
    pd.concat(source).to_csv(targetdir / (stem + '_source.tsv'), sep='\t', index=False)
    plt.close(fig)

# Step3. 所有53状态按冻结父谱系排序，不按显著性删行或换顺序。
parents = ['B', 'CD4_T', 'CD8_T', 'Endothelial', 'Enteroendocrine', 'Fibroblast', 'Glial', 'Goblet', 'ILC', 'Mast', 'Monocyte_Macrophage', 'Pericyte']
state_order = [state for parent in parents for state in membership.loc[membership.CellType.eq(parent), 'ID']]
states = data.loc[data.Family.eq('State')].copy()
assert states.OR.gt(0).all()
states['Log2OR'] = np.log2(states.OR)
matrix = states.pivot(index='Feature', columns='Endpoint', values='Log2OR').reindex(index=state_order, columns=endpoints)
q1 = states.pivot(index='Feature', columns='Endpoint', values='Q_Firth').reindex_like(matrix)
q2 = states.pivot(index='Feature', columns='Endpoint', values='Q_Wilcoxon').reindex_like(matrix)
fig, ax = plt.subplots(figsize=(7.2, 9.0))
fig.subplots_adjust(left=.33, right=.97, bottom=.105, top=.91)
image = ax.imshow(matrix, cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2), aspect='auto', interpolation='none')
ax.set_xticks(range(6), ['ETR\nW10', 'IFX\nW4–6', 'VDZ\nW6', 'VDZ\nW12', 'VDZ\nW52', 'RISK*\n3-year'])
ax.xaxis.tick_top()
ax.tick_params(axis='x', length=0, pad=6)
ax.set_yticks(range(53), state_order)
ax.tick_params(axis='y', length=0, labelsize=6, pad=5)
ax.set_title('All 53 frozen cell states', pad=30, fontsize=8)
for i in range(53):
    for j in range(6):
        if q1.iloc[i, j] < .05 and q2.iloc[i, j] < .05:
            ax.text(j, i, '*', ha='center', va='center', fontsize=8, color='black')
        elif q2.iloc[i, j] < .05:
            ax.text(j, i, '^', ha='center', va='center', fontsize=8, color='black')
        elif q1.iloc[i, j] < .05:
            ax.text(j, i, '+', ha='center', va='center', fontsize=8, color='black')
boundaries = np.cumsum([int(membership.CellType.eq(p).sum()) for p in parents])[:-1]
for boundary in boundaries:
    ax.axhline(boundary - .5, color='white', lw=.9)
for spine in ax.spines.values():
    spine.set_visible(False)
cax = fig.add_axes([.43, .045, .43, .014])
bar = fig.colorbar(image, cax=cax, orientation='horizontal', extend='both')
bar.set_label('log2 OR per 1 SD (display saturated at ±2)', fontsize=6)
bar.ax.tick_params(labelsize=6, length=2)
fig.canvas.draw()
stem = 'Fig3_all_53_state_associations'
require_matplotlib_panel_alignment(fig, exclude_axes=[cax], json_out=str(figdir / (stem + '.alignment.json')), strict=True)
fig.savefig(figdir / (stem + '.png'), dpi=600)
fig.savefig(figdir / (stem + '.pdf'))
fig.savefig(figdir / (stem + '.svg'))
states.to_csv(figdir / (stem + '_source.tsv'), sep='\t', index=False)
plt.close(fig)

# Step4. 展示3个队列内生态型候选的全部患者，抖动只分离重叠点。
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9))
fig.subplots_adjust(left=.08, right=.98, bottom=.24, top=.76, wspace=.48)
rng = np.random.default_rng(20261015)
sources = []
for panel, (ax, visit, feature, folder) in enumerate(zip(axes, [endpoints[1], endpoints[2], endpoints[5]], ['E3', 'E2', 'E2'], ['03_GSE73661_mucosal_healing', '03_GSE73661_mucosal_healing', '02_RISK_disease_progression'])):
    patients = pd.read_csv(out / folder / (visit + '_analysis_patients.tsv'), sep='\t')
    scores = pd.read_csv(out / folder / 'continuous_ecotype_scores.tsv', sep='\t', index_col=0)
    patients['Score'] = scores.loc[patients.SampleID, feature].to_numpy()
    patients['Feature'] = feature
    patients['Endpoint'] = visit
    sources.append(patients[['SampleID', 'PatientID', 'Outcome', 'Score', 'Feature', 'Endpoint']])
    for group, color in [(0, neutral), (1, blue)]:
        values = patients.loc[patients.Outcome.eq(group), 'Score'].to_numpy()
        ax.scatter(group + rng.uniform(-.12, .12, len(values)), values, s=9, color=color, alpha=.75, edgecolors='none', zorder=3)
        median = np.median(values)
        ax.plot([group - .18, group + .18], [median, median], color='black', lw=1.1, zorder=4)
    labels = ['No healing', 'Healing'] if panel < 2 else ['No progression', 'Progression']
    labels = [label + f'\n(n = {int(patients.Outcome.eq(i).sum())})' for i, label in enumerate(labels)]
    ax.set_xticks([0, 1], labels)
    ax.set_xlim(-.45, 1.45)
    ax.set_ylabel(feature + ' relative abundance')
    ax.set_ylim(-.025, max(.15, patients.Score.max() * 1.12))
    ax.set_title(titles[visit] + '\n' + feature + f' · N = {len(patients)}', pad=11)
    ax.annotate(chr(97 + panel), xy=(0, 1), xycoords='axes fraction', xytext=(-24, 38), textcoords='offset points', fontsize=8, fontweight='bold', annotation_clip=False)
fig.canvas.draw()
stem = 'Fig4_candidate_patient_distributions'
require_matplotlib_panel_alignment(fig, json_out=str(figdir / (stem + '.alignment.json')), overlay_svg=str(figdir / (stem + '.alignment.svg')), strict=True, tolerance_pt=1.5)
fig.savefig(figdir / (stem + '.png'), dpi=600)
fig.savefig(figdir / (stem + '.pdf'))
fig.savefig(figdir / (stem + '.svg'))
pd.concat(sources).to_csv(figdir / (stem + '_source.tsv'), sep='\t', index=False)
plt.close(fig)
print('PLOTS_COMPLETE', flush=True)
