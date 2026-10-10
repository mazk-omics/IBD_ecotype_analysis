# 将冻结53状态、7生态型和两队列恢复证据整理成可逐项审阅的报告与图件。
# 本机cibersortx_review环境运行；只读已有统计，不重算恢复或修改模型。
# 写入展示完善_20261010；源表、分母、判断来源和未完成方法均保留。
from pathlib import Path
import hashlib
import html
import json
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Step1. 读入既有结果，记录本轮使用的文件校验值
project = Path(__file__).resolve().parents[3] / 'IBD_ecotype_analysis'
root = project / '06_validation/ecotyper_original_method_annotation_20261010'
prior = project / '06_validation/recovery_interpretation_20261010'
out = root / '展示完善_20261010'
out.mkdir(exist_ok=True)
source_files = [root / name for name in ['53_states_preliminary_literature_evidence.tsv', 'ecotype_coassociation_both_datasets.tsv', 'state_literature_marker_correspondence.tsv', 'marker_eligibility_both_datasets.tsv', 'inputs/frozen_ecotypes.tsv', 'inputs/literature_marker_panels.tsv']]
source_files += [prior / name for name in ['各状态前三作者亚型_仅描述.tsv', '53状态_数据库基因重叠证据.tsv', '文献方法与数据库来源.tsv']]
source_files += [root / (ds + '_C_four_matrix_mean.tsv') for ds in ['Colon_GSE214695', 'Ileum_IBDverse']]
before = {str(p.relative_to(project)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
states = pd.read_csv(source_files[0], sep='\t')
coassociation = pd.read_csv(source_files[1], sep='\t')
literature = pd.read_csv(source_files[2], sep='\t')
eligibility = pd.read_csv(source_files[3], sep='\t')
ecotypes = pd.read_csv(source_files[4], sep='\t')
authors = pd.read_csv(prior / '各状态前三作者亚型_仅描述.tsv', sep='\t')
hallmark = pd.read_csv(prior / '53状态_数据库基因重叠证据.tsv', sep='\t')
sources = pd.read_csv(prior / '文献方法与数据库来源.tsv', sep='\t')
sources.loc[sources.SourceID.eq('M1'), 'ActuallyUsed'] = '官方reference-guided recovery；人工文献程序对应；两队列F1/F2四套Spearman与CE成员内外共关联'
sources.loc[sources.SourceID.eq('M1'), 'Limit'] = '原状态恢复30次打乱；CE各1,000次行置换；NA合并为明示项目选择；七指标外部marker/整转录组GSEA及全局百万次检验未完成'
authors['CellState'] = authors.CellType + '_' + authors.State
assert len(states) == 53 and states.CellState.is_unique and len(ecotypes) == 7
assert len(coassociation) == 14 and len(literature) == 1060
print('输入：53状态、7生态型、两队列14条共关联结果；20程序×53状态文献对应。')

# Step2. 保留数量、供者覆盖、marker资格和表达方向，不合成“准确率”
for ds, suffix in [('Colon_GSE214695', 'Colon'), ('Ileum_IBDverse', 'Ileum')]:
    qualified = eligibility[eligibility.Dataset.eq(ds)].set_index('CellState')
    states['MarkerEligible_' + suffix] = states.CellState.map(qualified.PaperMarkerEligible)
    states['LineagePresent_' + suffix] = states.CellState.map(qualified.LineagePresent)
states['资格范围'] = np.select([states.MarkerEligible_Colon & states.MarkerEligible_Ileum, states.MarkerEligible_Colon | states.MarkerEligible_Ileum], ['两队列均≥5', '仅一个队列≥5'], default='两队列均不足/不可评估')
states['同源marker方向'] = '至少一队列无可计算读出'
both_delta = states.OwnMinusBestOther_InternalOnly_Colon.notna() & states.OwnMinusBestOther_InternalOnly_Ileum.notna()
states.loc[both_delta, '同源marker方向'] = '至少一队列≤0'
both_positive = both_delta & states.OwnMinusBestOther_InternalOnly_Colon.gt(0) & states.OwnMinusBestOther_InternalOnly_Ileum.gt(0)
states.loc[both_positive, '同源marker方向'] = '两队列均>0'
states['需保留谱系疑点'] = states.CellState.isin(['Glial_S02', 'Glial_S03', 'Mast_S01', 'Mast_S03', 'Mast_S07'])
states['证据身份'] = '候选程序；未完成原文七指标外部marker与整转录组GSEA'
states.to_csv(out / '53状态_完整审阅表.tsv', sep='\t', index=False)
literature = literature.merge(states[['CellState', 'CellType']], on='CellState', validate='many_to_one')
literature['对应范围'] = ['共享程序' if line == 'Any_shared_program' else '父谱系相符' if cell in line.split(';') else '跨父谱系读出' for cell, line in zip(literature.CellType, literature.ExpectedLineage)]
literature.to_csv(out / '文献程序_逐项对应及范围.tsv', sep='\t', index=False)
print('外部marker资格：', states.groupby('资格范围').size().to_dict(), '；两队列同源marker均正差值：', int(both_positive.sum()))
ce_judgments = {
    'E1': '两队列共同变化证据较强；炎症/干扰素候选组合优先解释，Mast_S01与稀有EEC成员仍待核对。',
    'E2': '整体共同变化证据较弱；Fibroblast_S02的FRC样程序不能代表整个E2已恢复。',
    'E3': '两队列共同变化证据较强；基质、结构调控和周期程序混合，不宜统称增殖生态型。',
    'E4': '两队列方向不一致；回肠仅2/3成员可计算，负Z须结合稀有Mast状态和覆盖变化阅读。',
    'E5': '支持主要来自结直肠；回肠缺胶质，不能外推为两个部位共同成立的生态型。',
    'E6': '两队列共同变化证据较弱；成员被分配到不等于组合关系得到恢复。',
    'E7': '两队列共同变化证据较弱；保留成员定义，暂不作确定生态型功能命名。'
}
coassociation['初步判断'] = coassociation.Ecotype.map(ce_judgments)
coassociation = coassociation.merge(ecotypes[['Ecotype', 'InitialInterpretation', 'Members']], on='Ecotype', validate='many_to_one')
coassociation.to_csv(out / '7生态型_共关联与覆盖审阅表.tsv', sep='\t', index=False)
ecotypes.to_csv(out / '7生态型_冻结组成与最高分配描述.tsv', sep='\t', index=False)
sources.to_csv(out / '来源索引.tsv', sep='\t', index=False)
state_order = [s for members in ecotypes.Members for s in members.split(';')]
assert len(state_order) == 53 and set(state_order) == set(states.CellState) and len(set(state_order)) == 53
states.set_index('CellState').loc[state_order, ['Ecotype', 'CellType']].to_csv(out / '图3_矩阵行列顺序.tsv', sep='\t')

# Step3. 固定尺寸绘图；不添加没有计算过的P值或生物学置信区间
sys.path.insert(0, str(Path.home() / '.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .6, 'savefig.facecolor': 'white'})
datasets = ['Colon_GSE214695', 'Ileum_IBDverse']
colors = ['#356C8C', '#C27643']
ce_names = ['E' + str(i) for i in range(1, 8)]
y = np.arange(7)
fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 112/25.4), gridspec_kw={'width_ratios': [1.65, 1]})
fig.subplots_adjust(left=.10, right=.97, bottom=.20, top=.77, wspace=.35)
fig.suptitle('Ecotype coassociation: support and coverage', x=.10, ha='left', y=.98, fontsize=11)
fig.text(.10, .90, 'Colon: 18 donors    |    Ileum: 30 donors    |    Frozen discovery: 7 ecotypes', fontsize=8)
for k, ds in enumerate(datasets):
    data = coassociation[coassociation.Dataset.eq(ds)].set_index('Ecotype').loc[ce_names]
    pos = y + (k-.5)*.42
    axes[0].scatter(data.Z, pos, s=32, color=colors[k], marker=['o', 's'][k], label=['Colon', 'Ileum'][k], zorder=3)
    axes[1].scatter(data.DefinedReferenceMembers / data.FrozenMembers, pos, s=32, color=colors[k], marker=['o', 's'][k], zorder=3)
    for yi, (_, row) in zip(pos, data.iterrows()):
        axes[1].text(1.10, yi, f'{int(row.DefinedReferenceMembers)}/{int(row.FrozenMembers)}', va='center', fontsize=7, color=colors[k])
axes[0].axvline(0, color='#8E98A0', lw=.7, zorder=1)
axes[0].set(xlim=(-5.5, 9), xlabel='Permutation-standardized coassociation Z')
axes[1].set(xlim=(.55, 1.26), xticks=[.6, .8, 1], xlabel='Computable / frozen members')
axes[0].set_title('a  Common variation', loc='left', pad=13)
axes[1].set_title('b  Coverage boundary', loc='left', pad=13)
for ax in axes:
    ax.set(yticks=y, yticklabels=ce_names, ylim=(6.6, -.6))
axes[0].grid(axis='y', color='#E8ECEF', lw=.5, zorder=0)
fig.legend(*axes[0].get_legend_handles_labels(), loc='upper left', bbox_to_anchor=(.10, .87), ncol=2, fontsize=8, frameon=False, borderaxespad=0)
fig.text(.10, .055, '1,000 row permutations per cohort. Z is not the official cell-state recovery Z.\nPoint estimates; no biological confidence intervals. Coverage uses finite member contrasts.', fontsize=7, linespacing=1.6)
require_matplotlib_panel_alignment(fig, json_out=out / '图1.alignment.json', strict=True)
fig.savefig(out / '图1_生态型共关联与覆盖.pdf')
fig.savefig(out / '图1_生态型共关联与覆盖.svg')
fig.savefig(out / '图1_生态型共关联与覆盖.png', dpi=300)
plt.close(fig)

# 图2保留全部53状态；缺失谱系与不足5细胞分别显示。
display = states.sort_values(['CellType', 'State']).reset_index(drop=True)
n_matrix = display[['AssignedCells_Colon', 'AssignedCells_Ileum']].to_numpy(float)
delta_matrix = display[['OwnMinusBestOther_InternalOnly_Colon', 'OwnMinusBestOther_InternalOnly_Ileum']].to_numpy(float)
fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 230/25.4))
fig.subplots_adjust(left=.34, right=.96, bottom=.125, top=.86, wspace=.38)
fig.suptitle('All 53 states: cell coverage and marker contrast', x=.04, ha='left', y=.985, fontsize=10)
fig.text(.04, .945, 'Colors and values answer different questions; no combined validation score.', fontsize=7.5)
counts_cmap = mpl.colormaps['Blues'].copy()
counts_cmap.set_bad('#E5E7E9')
delta_cmap = mpl.colormaps['RdBu_r'].copy()
delta_cmap.set_bad('#E5E7E9')
im0 = axes[0].imshow(np.log10(1 + n_matrix), aspect='auto', cmap=counts_cmap, vmin=0, vmax=2.7, interpolation='none')
im1 = axes[1].imshow(delta_matrix, aspect='auto', cmap=delta_cmap, vmin=-.36, vmax=.36, interpolation='none')
axes[0].set_title('a  Assigned cells', loc='left', pad=13)
axes[1].set_title('b  Own minus other', loc='left', pad=13)
for ax in axes:
    ax.set(xticks=[0, 1], xticklabels=['Colon', 'Ileum'], yticks=np.arange(53))
    ax.tick_params(axis='both', length=0)
axes[0].set_yticklabels(display.CellState, fontsize=7)
axes[1].set_yticklabels([])
for i in range(53):
    for j in range(2):
        n = n_matrix[i, j]
        axes[0].text(j, i, 'NA' if not np.isfinite(n) else str(int(n)) + ('*' if n < 5 else ''), ha='center', va='center', fontsize=7, color='white' if n > 120 else '#202B34')
        d = delta_matrix[i, j]
        axes[1].text(j, i, 'NA' if not np.isfinite(d) else f'{d:+.2f}', ha='center', va='center', fontsize=7, color='white' if abs(d) > .24 else '#202B34')
for i in range(1, 53):
    if display.CellType.iloc[i] != display.CellType.iloc[i-1]:
        for ax in axes:
            ax.axhline(i-.5, color='white', lw=1.1)
cb0 = fig.add_axes([.34, .072, .24, .009])
cb1 = fig.add_axes([.71, .072, .24, .009])
fig.colorbar(im0, cax=cb0, orientation='horizontal', ticks=[0, 1, 2])
fig.colorbar(im1, cax=cb1, orientation='horizontal', ticks=[-.3, 0, .3])
cb0.set_xlabel('Color: log10(1 + cells)', fontsize=7, labelpad=2)
cb1.set_xlabel('Marker-score difference', fontsize=7, labelpad=2)
fig.text(.04, .008, '* <5 cells: ineligible for that cohort’s external marker step. NA: not evaluable.\nShared projection genes: marker contrast is not independent accuracy.', fontsize=7)
require_matplotlib_panel_alignment(fig, json_out=out / '图2.alignment.json', exclude_axes=[cb0, cb1], strict=True)
fig.savefig(out / '图2_53状态证据总览.pdf')
fig.savefig(out / '图2_53状态证据总览.svg')
fig.savefig(out / '图2_53状态证据总览.png', dpi=300)
plt.close(fig)
display.to_csv(out / '图2_源数据及显示顺序.tsv', sep='\t', index=False)

# 图3显示同一冻结顺序；回肠47状态矩阵扩至53行列仅为保留缺失，不改计算。
fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 112/25.4))
fig.subplots_adjust(left=.065, right=.975, bottom=.28, top=.80, wspace=.20)
fig.suptitle('State coassociation within frozen ecotype groups', x=.065, ha='left', y=.98, fontsize=10)
fig.text(.065, .90, 'Same 53-state order and color scale; gray = undefined, including the diagonal.', fontsize=7.5)
corr_cmap = mpl.colormaps['RdBu_r'].copy()
corr_cmap.set_bad('#D9DDE1')
sizes = ecotypes.MemberStates.to_numpy(int)
boundaries = np.cumsum(sizes)
centers = boundaries - sizes/2 - .5
for ax, ds, letter in zip(axes, datasets, ['a', 'b']):
    corr = pd.read_csv(root / (ds + '_C_four_matrix_mean.tsv'), sep='\t', index_col=0).reindex(index=state_order, columns=state_order)
    corr.to_csv(out / (ds + '_图3显示矩阵.tsv'), sep='\t')
    im = ax.imshow(corr, cmap=corr_cmap, vmin=-1, vmax=1, interpolation='none', aspect='equal')
    for edge in boundaries[:-1] - .5:
        ax.axhline(edge, color='#5C646C', lw=.35)
        ax.axvline(edge, color='#5C646C', lw=.35)
    ax.set(xticks=centers, xticklabels=ce_names, yticks=centers, yticklabels=ce_names)
    ax.tick_params(length=0)
    ax.set_title(letter + '  ' + ('Colon: 53 states, 18 donors' if ds.startswith('Colon') else 'Ileum: 47 states, 30 donors'), loc='left', pad=13)
colorbar = fig.add_axes([.30, .16, .40, .024])
fig.colorbar(im, cax=colorbar, orientation='horizontal', ticks=[-1, -.5, 0, .5, 1])
colorbar.set_xlabel('Mean of four finite Spearman correlations', fontsize=7, labelpad=3)
fig.text(.065, .022, 'Borders show discovery membership, not confirmed recovered clusters.\nIleum Glial rows/columns remain NA; row-order table identifies every state.', fontsize=7)
require_matplotlib_panel_alignment(fig, json_out=out / '图3.alignment.json', exclude_axes=[colorbar], strict=True)
fig.savefig(out / '图3_生态型成员关联结构.pdf')
fig.savefig(out / '图3_生态型成员关联结构.svg')
fig.savefig(out / '图3_生态型成员关联结构.png', dpi=300)
plt.close(fig)
print('完成三张图及逐项源数据。')

# Step4. 简明报告给阅读顺序，完整卡片给具体证据
references = {
    'Luca2021': ('Luca et al., Cell 2021', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/'),
    'Steen2021': ('Steen et al., Cancer Cell 2021', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9205168/'),
    'Martin2019': ('Martin et al., Cell 2019', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7060942/'),
    'Kinchen2018': ('Kinchen et al., Cell 2018', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6176871/'),
    'Smillie2019': ('Smillie et al., Cell 2019', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6662628/'),
    'Elmentaite2021': ('Elmentaite et al., Nature 2021', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8426186/')
}
summary = ['# Cell state与生态型：两队列证据审阅', '',
    '2026-10-10展示完善版。对象固定为MSCCR_CORE12_CTSPEC_REGION_IBD：CD+UC discovery、Region3组内归一化、BayesPrism presort；53状态、12父类型、7生态型。此次新增展示、证据串联与阅读说明，没有重新训练模型或改动原始数据。所有名称均为候选程序，不是最终细胞亚型。', '',
    '## 建议阅读顺序', '',
    '1. 图1：哪些生态型成员仍共同变化？E1、E3相对较强；同时看可计算成员覆盖，不能把回肠子集说成完整53状态模型验证。',
    '2. 图2：每个状态有没有足够细胞、marker方向如何？保留全部53状态，不合成未经验证的可信度总分。',
    '3. 图3：汇总统计背后的关联结构。灰色包含缺失谱系、无法计算的边和对角，不是低表达。',
    '4. HTML阅读页：按父类型、生态型或基因筛选53状态，展开具体文献、作者标签和疑点。', '',
    '## 两个恢复队列的评估范围', '',
    '| 项目 | 结直肠GSE214695 | 回肠IBDverse |', '|---|---|---|',
    '| 选中样本/供者 | 18 | 30 |', '| 抽样输入细胞 | 9,752 | 4,818 |',
    '| 项目抽样上限 | 每样本×父类型100 | 每供者×父类型20 |',
    '| 可输入状态 | 53 | 47；无Glial注释 |',
    '| 满足原文≥5细胞外部marker资格 | 48/53 | 44/47 |', '',
    '按父类型设上限的小规模恢复不能代表整块黏膜细胞比例，也不能跨队列比较生态型患病率。结直肠HC来自非肿瘤区域，并非无条件的健康黏膜金标准；来源索引Q3。', '',
    f'按冻结53状态，42个两队列均满足外部marker资格，8个只满足一个，3个均不足或不可评估。共同可输入47状态中，{int(both_positive.sum())}个两队列同源marker差值均>0；这个方向检查不是独立准确率。', '',
    '## 生态型逐项判断', '',
    '| 生态型 | 结直肠Z；可计算成员 | 回肠Z；可计算成员 | 初步判断 |', '|---|---|---|---|']
for ce in ce_names:
    c = coassociation[coassociation.Dataset.eq(datasets[0]) & coassociation.Ecotype.eq(ce)].iloc[0]
    i = coassociation[coassociation.Dataset.eq(datasets[1]) & coassociation.Ecotype.eq(ce)].iloc[0]
    summary.append(f'| {ce} | {c.Z:.2f}；{c.DefinedReferenceMembers}/{c.FrozenMembers} | {i.Z:.2f}；{i.DefinedReferenceMembers}/{i.FrozenMembers} | {ce_judgments[ce]} |')
summary += ['', 'Z为1,000次行置换的生态型共关联Z，不是官方状态恢复Z；未定义尾部P值，也未作7生态型多重检验显著性声明。回肠E4的零分布可计算成员为1–3，观测仅2；负Z不直接等同真实生物学排斥。E5也需注意零分布覆盖变化。', '',
    '![图1](图1_生态型共关联与覆盖.png)', '', '![图3](图3_生态型成员关联结构.png)', '',
    '## 优先解释的状态及边界', '',
    '| 状态 | 可核对的依据 | 候选判断与限制 |', '|---|---|---|',
    '| Fibroblast_S02（E2） | CCL21/FDCSP/CCL19/CLU：Kinchen读出4/7，前20含前三项；两队列marker正差值 | FRC样/淋巴组织基质程序；不等同论文S4，不证明TLS或整个E2恢复 |',
    '| CD4_T_S01（E7） | KLRB1/CD69/PPP1R15A：Martin驻留程序3/31 | 驻留相关候选，不足以确认完整TRM亚型 |',
    '| CD4_T_S03、CD8_T_S02（E2） | LEF1/CCR7/SELL：初始/中央记忆程序3/7，另有部分细胞毒读出 | 低分化/记忆相关候选；CXCL13不能单独定义Tfh或耗竭T |',
    '| Goblet_S01（E3） | 周期程序9/49，命中不在前20；前20另有FOXM1/E2F2等 | 共享周期候选；面板来自T细胞，不转用T亚型名，不代表全部杯状细胞增殖 |',
    '| Glial_S05（E3） | CDH19/MPZ/FOXD3：胶质身份3/7；结直肠marker差值负，回肠不可评估 | 身份marker支持与状态区分不足同时保留 |',
    '| Mast_S01、Mast_S03、Mast_S07 | 粒细胞/淋巴/髓系样读出与谱系解释有疑点 | 保留冻结ID，不能自动改名或认定污染 |', '',
    '对应来源：[Kinchen 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6176871/)、[Martin 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC7060942/)、[Elmentaite 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8426186/)。逐项基因、排序、面板长度和适用谱系见文献对应TSV。', '',
    '## 全部53状态审阅索引', '', '![图2](图2_53状态证据总览.png)', '',
    '| 状态 | 生态型 | 候选程序 | 结直肠细胞/供者 | 回肠细胞/供者 | 外部marker资格 |', '|---|---|---|---|---|---|']
for row in display.itertuples():
    ileum_count = '不可评估' if pd.isna(row.AssignedCells_Ileum) else f'{int(row.AssignedCells_Ileum)}/{int(row.AssignedDonors_Ileum)}'
    scope = states.loc[states.CellState.eq(row.CellState), '资格范围'].iloc[0]
    summary.append(f'| {row.CellState} | {row.Ecotype} | {row.InitialProgram} | {row.AssignedCells_Colon}/{row.AssignedDonors_Colon} | {ileum_count} | {scope} |')
method_steps = [
    '读取冻结状态、成员、两队列既有汇总与四套Spearman平均矩阵，核对53状态、7生态型及成员ID。',
    '按原文≥5分配细胞规则标记外部marker资格，另列供者数；不删除发现状态。',
    '分别展示自身发现marker差值、文献20程序对应、作者标签供者等权构成；它们不是独立分类准确率。',
    '生态型复用F1/F2及零值置NA的四套Spearman：四套均有限才保留平均边。比较成员内外相关，使用保存的1,000次行置换Z，同时记录成员及零分布覆盖。',
    '16个MSigDB Human Hallmark仅作基因交集，无NES/P/FDR；完整集合留在原表，数据库发布版本未标明。文献短marker对应也不称为富集。',
    '此次仅输出表、图和报告，核对原表SHA256、图件布局与字号；没有新的服务器拟合或临床关联。'
]
pending = [
    '同批恢复细胞的整转录组外部DE与原文七指标排序尚未完成；IBD没有“肿瘤类型数”，需要明示适配。',
    '作者亚型整转录组排序与外部状态前50marker的fgsea（10,000置换）尚未完成；短面板交集不能替代。',
    '按数据集留出的前10marker验证、全局CE百万次检验、配体受体与生态型正式Hallmark比较尚未完成。',
    '现有官方状态恢复采用30次打乱，论文为100；本轮CE的1,000次是另一检验，不能视为补足了状态恢复。'
]
summary += ['', '## 实际操作步骤与方法身份', '']
summary += [f'{k}. {step}' for k, step in enumerate(method_steps, 1)]
summary += ['', '按[Luca et al., Cell 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/)程序/marker人工注释思路及单细胞生态型共关联公式整理；[Steen et al., Cancer Cell 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC9205168/)补充亚型比较思路。既有环境保留Seurat3.1.3/fgsea1.14.0，当前展示不重新调用它们。原文未详述的NA合并采用已记录的项目选择，不能称为全部步骤“一字不差”的完整数值复现。', '', '## 仍需完成的证据', '']
summary += ['- ' + item for item in pending]
summary += ['', '当前适合汇报“冻结模型的双队列初步恢复、共同变化与候选程序解释”，尚不适合宣布53个IBD细胞亚型已最终定义。', '', '## 文件入口', '', '- 双击结果审阅.html：筛选、展开全部状态证据。', '- 53状态_完整审阅表.tsv：可筛选的证据总表。', '- 7生态型_共关联与覆盖审阅表.tsv：两队列统计和覆盖。', '- 图1–图3：各自PDF/SVG/PNG；图2源数据、图3行序及矩阵另存。', '- 前序原始计算、方法边界和环境记录仍在上一级目录。', '', '## 主要来源', '']
summary += [f'- [{title}]({url})' for title, url in references.values()]
(out / '先读我_结果审阅.md').write_text('\n'.join(summary) + '\n', encoding='utf8')

# Step5. 生成完整证据卡；跨父类型读出只作为疑点，不自动改变命名
escape = html.escape
cards = []
for row in display.itertuples():
    relevant = literature[literature.CellState.eq(row.CellState) & literature.Hits.gt(0)].sort_values(['Top20Hits', 'Hits'], ascending=False)
    preferred = relevant[relevant['对应范围'].isin(['共享程序', '父谱系相符'])]
    cross = relevant[relevant['对应范围'].eq('跨父谱系读出')]
    marker_lines = [f'<li><b>{escape(e.Program)}</b>：{e.Hits}/{e.PanelSize}，前20命中{e.Top20Hits}；{escape(str(e.HitGenes))}。<a href="{e.URL}">{escape(e.Publication)}</a>；{escape(e.PanelStatus)}</li>' for e in preferred.itertuples()]
    cross_lines = [f'<li>{escape(e.Program)}：{e.Hits}/{e.PanelSize}；{escape(str(e.HitGenes))}。<a href="{e.URL}">{escape(e.Publication)}</a></li>' for e in cross.itertuples()]
    metrics = []
    for suffix, title in [('Colon', '结直肠'), ('Ileum', '回肠')]:
        n = getattr(row, 'AssignedCells_' + suffix)
        d = getattr(row, 'AssignedDonors_' + suffix)
        z = getattr(row, 'Z_' + suffix)
        delta = getattr(row, 'OwnMinusBestOther_InternalOnly_' + suffix)
        n_text = '不可评估：无父谱系输入' if pd.isna(n) else f'{int(n)}细胞 / {int(d)}供者'
        z_text = 'NA' if pd.isna(z) else f'{z:.2f}'
        delta_text = 'NA' if pd.isna(delta) else f'{delta:+.3f}'
        qualified = getattr(row, 'MarkerEligible_' + suffix)
        metrics.append(f'<tr><th>{title}</th><td>{n_text}</td><td>{z_text}</td><td>{delta_text}</td><td>{"满足≥5" if qualified else "不满足/不可评估"}</td></tr>')
    author_lines = []
    for ds, title in zip(datasets, ['结直肠', '回肠']):
        evidence = authors[authors.Dataset.eq(ds) & authors.CellState.eq(row.CellState)].sort_values('DonorBalancedFractionWithinState', ascending=False).head(3)
        text = '；'.join(f'{escape(e.AuthorSubtype)} {e.DonorBalancedFractionWithinState:.1%}（{e.Cells}细胞，{e.LabelDonors}供者）' for e in evidence.itertuples())
        author_lines.append(f'<li>{title}：{text or "没有可对应的分配细胞"}。</li>')
    hall = hallmark[hallmark.CellState.eq(row.CellState) & hallmark.OverlapCount.gt(0)].sort_values('OverlapCount', ascending=False).head(3)
    hall_lines = [f'<li>{escape(e.GeneSet)}：{e.OverlapCount}/{e.PanelGeneSetSize}；{escape(str(e.OverlapGenes))}</li>' for e in hall.itertuples()]
    scope = states.loc[states.CellState.eq(row.CellState), '资格范围'].iloc[0]
    flag = '<span class="warn">保留谱系疑点</span>' if row.需保留谱系疑点 else ''
    search = row.CellState + ' ' + row.Judgment + ' ' + row.Top20Genes + ' ' + str(row.LiteratureCorrespondence)
    cards.append(f'''<details class="state" id="{row.CellState}" data-type="{row.CellType}" data-ce="{row.Ecotype}" data-search="{escape(search.lower(), quote=True)}">
<summary><strong>{row.CellState}</strong><span class="badge">{row.Ecotype}</span>{escape(row.InitialProgram)} {flag}</summary>
<div class="body"><p><b>当前判断：</b>{escape(row.Judgment)}</p><p class="muted">候选程序，未最终定义。原人工支持级别“{escape(row.ProgramConfidence)}”不是统计显著性；外部marker资格：{scope}。</p>
<table><thead><tr><th>队列</th><th>分配覆盖</th><th>官方状态Z</th><th>同源marker差值</th><th>外部marker资格</th></tr></thead><tbody>{''.join(metrics)}</tbody></table>
<p class="muted">三类统计分别阅读；表达差值可能共享恢复投影基因，不能当作独立准确率。</p>
<h4>前20发现关联基因</h4><p class="genes">{escape(row.Top20Genes.replace(',', ', '))}</p>
<h4>父谱系相符或共享程序的文献读出</h4><ul>{''.join(marker_lines) or '<li>当前20程序没有相符读出命中，不等于没有基因表达。</li>'}</ul>
<details><summary>跨父谱系命中：提示共享程序或待审疑点</summary><ul>{''.join(cross_lines) or '<li>没有命中。</li>'}</ul></details>
<h4>定向Hallmark前三交集（非富集）</h4><ul>{''.join(hall_lines) or '<li>16集无交集。</li>'}</ul><p class="muted">MSigDB Human Hallmark，前三按命中数展示；长度与状态marker数不同，不据此比较富集强弱。无NES/P/FDR，发布版本未标明；完整16集在原表。</p>
<h4>作者标签的描述性对应</h4><ul>{''.join(author_lines)}</ul><p class="muted">供者等权的状态内组成，最多前三标签；不是准确率，不重聚类。宽泛父类型标签不能证明状态亚型。</p>
<p class="limit"><b>具体疑点：</b>{escape(row.ReviewReason)}</p><p>原记录来源：{escape(row.SourceIDs)}。页面末索引说明来源与实际用途；文献逐基因来源和适用范围另存TSV。</p></div></details>''')
ce_cards = []
for row in ecotypes.itertuples():
    c = coassociation[coassociation.Dataset.eq(datasets[0]) & coassociation.Ecotype.eq(row.Ecotype)].iloc[0]
    i = coassociation[coassociation.Dataset.eq(datasets[1]) & coassociation.Ecotype.eq(row.Ecotype)].iloc[0]
    links = ', '.join(f'<a href="#{s}">{s}</a>' for s in row.Members.split(';'))
    ce_cards.append(f'''<article class="ecotype"><h3>{row.Ecotype} · {escape(row.InitialInterpretation)}</h3><p>{escape(ce_judgments[row.Ecotype])}</p>
<table><thead><tr><th>队列</th><th>AvgDif</th><th>共关联Z</th><th>可计算/冻结成员</th><th>零分布成员范围</th><th>最高生态型分配样本</th></tr></thead><tbody>
<tr><th>结直肠</th><td>{c.AvgDif:+.3f}</td><td>{c.Z:+.2f}</td><td>{c.DefinedReferenceMembers}/{c.FrozenMembers}</td><td>{c.MinNullDefinedMembers}–{c.MaxNullDefinedMembers}</td><td>{row.ColonOfficialAssignedSamples}/18</td></tr>
<tr><th>回肠</th><td>{i.AvgDif:+.3f}</td><td>{i.Z:+.2f}</td><td>{i.DefinedReferenceMembers}/{i.FrozenMembers}</td><td>{i.MinNullDefinedMembers}–{i.MaxNullDefinedMembers}</td><td>{row.IleumOfficialAssignedSamples}/30</td></tr></tbody></table>
<p class="muted">最高分配只是样本描述，不能代替共关联或跨队列当作组织患病率。AvgDif为成员内相关减成员外相关的成员均值。</p><p><b>冻结成员：</b>{links}</p></article>''')
type_options = ''.join(f'<option>{escape(t)}</option>' for t in sorted(states.CellType.unique()))
ce_options = ''.join(f'<option>{ce}</option>' for ce in ce_names)
reference_html = ''.join(f'<li><a href="{url}">{title}</a></li>' for title, url in references.values())
source_html = sources.to_html(index=False, escape=True, classes='sources')
method_html = ''.join(f'<li>{escape(step)}</li>' for step in method_steps)
pending_html = ''.join(f'<li>{escape(step)}</li>' for step in pending)
page = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Cell state与生态型证据审阅</title>
<style>
*{{box-sizing:border-box}}body{{margin:0;background:#f4f6f7;color:#213341;font:16px/1.75 "Microsoft YaHei",Arial,sans-serif}}main{{max-width:1160px;margin:auto;padding:35px 28px 70px}}h1{{font-size:31px;margin:0 0 10px}}h2{{font-size:25px;margin-top:40px;border-bottom:1px solid #ccd7dd;padding-bottom:8px}}h3{{font-size:19px;margin:0 0 10px}}h4{{margin-bottom:7px}}a{{color:#216789;text-decoration:none}}a:hover{{text-decoration:underline}}p{{margin:10px 0}}.muted{{color:#526776;font-size:14px}}.notice{{background:#e6eff3;border-left:4px solid #356c8c;padding:16px 20px}}nav{{display:flex;gap:20px;flex-wrap:wrap;padding:18px 0}}.tiles{{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:22px 0}}.tile{{background:white;padding:17px;border:1px solid #d6e0e5;border-radius:8px}}.tile strong{{display:block;font-size:28px}}figure{{margin:24px 0;background:white;border:1px solid #d6e0e5;padding:20px}}figure img{{display:block;max-width:100%;height:auto;margin:auto}}figcaption{{font-size:14px;color:#526776;padding-top:15px}}figure.tall img{{max-width:800px;width:100%}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:9px 12px;border-bottom:1px solid #dbe3e7;text-align:left;vertical-align:top}}thead{{background:#eef3f6}}.ecotype,.state{{background:white;border:1px solid #d6e0e5;border-radius:8px;margin:15px 0;padding:20px}}.state{{padding:0}}.state>summary{{padding:16px 20px;cursor:pointer}}.body{{padding:0 22px 22px}}.badge{{display:inline-block;color:#356c8c;background:#e5eef4;border-radius:4px;margin:0 10px;padding:0 8px;font-size:14px}}.warn{{background:#fff0d7;color:#7b5817;padding:3px 8px;border-radius:4px;font-size:13px}}.limit{{background:#fff4e3;padding:10px 13px}}.genes{{font-family:Consolas,monospace;font-size:14px;overflow-wrap:anywhere}}.filters{{position:sticky;top:0;background:#f4f6f7ef;padding:16px 0;display:flex;gap:10px;flex-wrap:wrap;z-index:2}}input,select,button{{font:inherit;padding:7px 10px;border:1px solid #bbcbd3;border-radius:5px}}input{{flex:1;min-width:220px}}button{{background:white;cursor:pointer}}details[hidden]{{display:none}}details:target{{border:2px solid #c27643;scroll-margin-top:90px}}.sources{{font-size:13px}}.sources td{{overflow-wrap:anywhere}}@media(max-width:700px){{main{{padding:20px 14px}}.tiles{{grid-template-columns:1fr}}table{{display:block;overflow-x:auto}}}}@media print{{body{{background:white;font-size:12px}}nav,.filters{{display:none}}main{{padding:0;max-width:none}}details,.ecotype,figure{{break-inside:avoid}}}}
</style><main><p class="muted">冻结Region3 · CD+UC discovery · BayesPrism presort · 2026-10-10</p>
<h1>Cell state与生态型：两队列证据审阅</h1><p>逐项区分恢复分配、marker、作者标签、文献和生态型共同变化。</p>
<nav><a href="#overview">阅读入口</a><a href="#ecotypes">7生态型</a><a href="#states">53状态卡</a><a href="#methods">方法与未完成步骤</a><a href="#sources">来源</a><a href="先读我_结果审阅.md">Markdown报告</a></nav>
<div class="notice">本轮完善展示，不重新拟合或修改原始数据。全部状态保留冻结ID；名称为候选程序。E1、E3共同变化相对较强，但不代表每个成员的谱系和功能命名已被验证。</div>
<section id="overview"><h2>三层证据怎么读</h2><div class="tiles"><div class="tile"><strong>53 / 7</strong>冻结状态 / 生态型<br><span class="muted">12父类型；没有新增或删减</span></div><div class="tile"><strong>42</strong>两队列均≥5细胞的状态<br><span class="muted">外部marker资格，不是验证通过数</span></div><div class="tile"><strong>34 / 47</strong>共同可输入状态中两队列marker均正差值<br><span class="muted">表达方向，不是独立准确率</span></div></div>
<p>结直肠：18供者、9,752输入细胞，父类型抽样上限100。回肠：30供者、4,818输入细胞，上限20；47个可输入状态，无Glial注释。设上限抽样不能代表整块黏膜比例，两队列也不是配对供者。结直肠HC为非肿瘤区域。</p>
<figure><img src="图1_生态型共关联与覆盖.png"><figcaption><b>图1｜成员共同变化及覆盖边界。</b>左为1,000次行置换标准化共关联Z，右为可计算/冻结成员。无生物学置信区间或未定义的显著性阈值；高Z不证明每个状态谱系正确。<a href="图1_生态型共关联与覆盖.pdf">PDF</a></figcaption></figure></section>
<section id="ecotypes"><h2>7生态型：支持与疑点放在一起</h2>{''.join(ce_cards)}
<figure><img src="图3_生态型成员关联结构.png"><figcaption><b>图3｜冻结成员关联结构。</b>同一53状态顺序和[-1,1]色标，灰色为NA。边界为发现阶段分组，不是新聚类。回肠胶质仅为缺失展示，原47状态计算不变。<a href="图3_生态型成员关联结构.pdf">PDF</a> · <a href="图3_矩阵行列顺序.tsv">逐状态行序</a></figcaption></figure></section>
<section id="states"><h2>53状态：展开具体依据</h2><p>图2先看覆盖与marker方向，再展开卡片查看基因、供者、文献、作者标签及疑点。短面板命中不自动成为细胞亚型标签。</p>
<figure class="tall"><img src="图2_53状态证据总览.png"><figcaption><b>图2｜全部状态证据。</b>左为分配细胞数，颜色log10(1+n)，数字为原始n；*表示不足5。右为自身marker减同父类型最佳其他状态得分，不是独立准确率；NA不补0。<a href="图2_53状态证据总览.pdf">PDF</a></figcaption></figure>
<div class="filters"><input id="search" placeholder="查状态、候选程序或基因，例如 CCL21"><select id="type"><option value="">全部父类型</option>{type_options}</select><select id="ce"><option value="">全部生态型</option>{ce_options}</select><button id="expand">展开可见状态</button><button id="collapse">收起</button></div><p id="shown" class="muted">显示53个状态</p>{''.join(cards)}</section>
<section id="methods"><h2>实际操作及未完成步骤</h2><ol>{method_html}</ol><div class="notice"><b>仍需补足：</b><ul>{pending_html}</ul></div>
<p>依据Cell 2021的人工程序对应和CE共关联方法；Cancer Cell 2021补充亚型比较思路。原文未详述的NA合并为已记录项目选择，不宣称完整数值复现。没有开展临床关联；模型缺Plasma/DC，不能称完整GIMATS。</p>
<p><a href="53状态_完整审阅表.tsv">53状态完整表</a> · <a href="7生态型_共关联与覆盖审阅表.tsv">CE统计与覆盖</a> · <a href="文献程序_逐项对应及范围.tsv">文献逐项证据</a> · <a href="../操作步骤与方法边界.md">原计算步骤</a> · <a href="../环境复现说明.md">环境记录</a></p></section>
<section id="sources"><h2>文献方法和数据库来源</h2><ul>{reference_html}</ul><p>各状态保留SourceID，索引说明实际使用及限制；下面URL可复制查阅。</p>{source_html}<p class="muted">来源索引M1已更新为当前已完成操作；其余来源沿用既有核查，原始来源表不覆盖。</p></section></main>
<script>
const cards=[...document.querySelectorAll('.state')];
function filter(){{const q=document.querySelector('#search').value.toLowerCase(),type=document.querySelector('#type').value,ce=document.querySelector('#ce').value;let n=0;cards.forEach(c=>{{c.hidden=!(c.dataset.search.includes(q)&&(!type||c.dataset.type===type)&&(!ce||c.dataset.ce===ce));if(!c.hidden)n++;}});document.querySelector('#shown').textContent='显示'+n+' / 53个状态';}}
['search','type','ce'].forEach(id=>document.querySelector('#'+id).addEventListener('input',filter));
document.querySelector('#expand').onclick=()=>cards.filter(c=>!c.hidden).forEach(c=>c.open=true);
document.querySelector('#collapse').onclick=()=>cards.forEach(c=>c.open=false);
function openAnchor(){{const target=document.getElementById(decodeURIComponent(location.hash.slice(1)));if(target?.classList.contains('state')){{document.querySelector('#search').value='';document.querySelector('#type').value='';document.querySelector('#ce').value='';filter();target.open=true;target.scrollIntoView();}}}}
window.addEventListener('hashchange',openAnchor);openAnchor();
</script></html>'''
(out / '结果审阅.html').write_text(page, encoding='utf8')

# Step6. 核查覆盖、源数据及下载原件均未变化
after = {str(p.relative_to(project)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
assert before == after
assert page.count('class="state"') == 53 and page.count('class="ecotype"') == 7
manifest = json.loads((root / 'download_manifest.json').read_text(encoding='utf8'))
for item in manifest:
    blob = (root / item['File']).read_bytes()
    assert len(blob) == item['Bytes'] and hashlib.sha256(blob).hexdigest() == item['SHA256'], item['File']
checks = {'Status': 'PASS', 'OriginalManifestFilesUnchanged': len(manifest), 'SourceFilesUnchanged': before, 'StateCards': 53, 'EcotypeCards': 7, 'MarkerEligibleBoth': int((states.MarkerEligible_Colon & states.MarkerEligible_Ileum).sum()), 'OwnMarkerPositiveBoth': int(both_positive.sum()), 'NoModelRefit': True, 'NoNewClinicalAnalysis': True}
(out / '展示数据核查.json').write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding='utf8')
print('报告完成；原始清单', len(manifest), '个文件及本轮来源表校验一致。')
