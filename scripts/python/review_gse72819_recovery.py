# 审阅既有 GSE72819 recovery：核对文件、状态分配、生态型分配和初步基因程序，不重跑或修改算法。
# 本机 cibersortx_review 环境运行；读取本地服务器副本与冻结discovery，输出审阅表和PNG。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Step1. 核对回传文件和73例输入注释；所有分配比例保留完整样本分母。
outdir = Path('D:/IBD_ecotype/IBD_ecotype_analysis/06_validation/GSE72819_recovery_review_20261009')
rawdir = outdir / 'raw_server'
snapshot = rawdir / 'audit_snapshot'
dataset = 'GSE72819_BP73_20261009_185536'
resultdir = rawdir / 'RecoveryOutput' / dataset
discovery_dir = Path('D:/IBD_ecotype/IBD_ecotype_analysis/ctspec_region_ibd')
definition_dir = Path('D:/IBD_ecotype/IBD_ecotype_analysis/06_validation/bayesprism_recovery_readiness_20261009')
manifest = json.loads((rawdir / 'audit_download_manifest.json').read_text(encoding='utf-8'))
for row in manifest:
    data = (rawdir / row['path']).read_bytes()
    assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'], row['path']
extra_file = snapshot / 'discovery_initial_ecotype_assignment.txt'
assert hashlib.sha256(extra_file.read_bytes()).hexdigest() == '3a9dfe9232ac6cb483620dbe43bce75f6a040a2bb070444acb225ef48bddb930'
annotation = pd.read_csv(snapshot / 'annotation.txt', sep='\t')
sample_ids = annotation['ID'].tolist()
assert len(sample_ids) == 73 and annotation['ID'].is_unique
assert annotation['Disease'].eq('UC').all() and annotation['Timepoint'].eq('Baseline').all()
ranks = pd.read_csv(snapshot / 'rank_data.txt', sep='\t')
ecotypes = pd.read_csv(snapshot / 'ecotypes.txt', sep='\t')
definitions = pd.read_csv(definition_dir / '53状态初步定义.tsv', sep='\t')
assert len(ranks) == 12 and len(ecotypes) == 53 and ecotypes['ID'].is_unique
assert set(definitions['CellState']) == set(ecotypes['ID'])
inventory = pd.DataFrame(manifest)
inventory['Group'] = np.where(inventory['path'].str.startswith('RecoveryOutput/'), '官方导出', np.where(inventory['path'].str.startswith('audit_snapshot/'), '内部结果及模型只读副本', '日志'))
inventory.to_csv(outdir / '输出文件清单.tsv', sep='\t', index=False)
print('Input: 73 UC baseline biopsies; verified files:', len(manifest), flush=True)

# Step2. 逐类型核对初始H、保留H和原始胜出分配；不重新归一化保留状态。
type_rows = []
state_rows = []
winner_rows = []
coefficient_tables = []
marker_rows = []
for _, rank_row in ranks.iterrows():
    cell_type = rank_row['CellType']
    rank = str(rank_row['Chosen_Rank'])
    internal_dir = snapshot / 'state_internal' / cell_type / rank
    model_dir = snapshot / 'frozen_model' / cell_type / rank
    W = pd.read_csv(model_dir / 'W.txt', sep='\t', index_col=0)
    mapping = pd.read_csv(model_dir / 'mapping_to_initial_states.txt', sep='\t')
    initial_info = pd.read_csv(model_dir / 'initial_gene_info.txt', sep='\t', index_col=0)
    gene_info = pd.read_csv(model_dir / 'gene_info.txt', sep='\t', index_col=0)
    raw_H = pd.read_csv(internal_dir / 'H.txt', sep='\t', index_col=0)
    initial_H = pd.read_csv(internal_dir / 'initial_state_abundances.txt', sep='\t', index_col=0)
    H = pd.read_csv(internal_dir / 'state_abundances.txt', sep='\t', index_col=0)
    exported_H = pd.read_csv(resultdir / cell_type / 'state_abundances.txt', sep='\t', index_col=0)
    assignment = pd.read_csv(internal_dir / 'state_assignment.txt', sep='\t')
    initial_assignment = pd.read_csv(internal_dir / 'initial_state_assignment.txt', sep='\t')
    scaled = pd.read_csv(internal_dir / 'expression_matrix_scaled.txt', sep='\t', index_col=0)
    assert list(H.columns) == sample_ids and list(initial_H.columns) == sample_ids
    assert np.isfinite(H).all().all() and H.ge(0).all().all() and np.isfinite(raw_H).all().all()
    assert np.allclose(initial_H.sum(axis=0), 1) and np.allclose(initial_H, raw_H.div(raw_H.sum(axis=0), axis=1))
    expected_H = initial_H.loc[mapping['InitialState']].copy()
    expected_H.index = mapping['State'].tolist()
    assert list(H.index) == mapping['State'].tolist() and np.allclose(H, expected_H, atol=1e-12)
    assert H.equals(exported_H) and H.sum(axis=0).le(1 + 1e-8).all()
    initial_winners = initial_H.idxmax(axis=0)
    final_winners = initial_winners.map(mapping.set_index('InitialState')['State'])
    expected_assignment = final_winners.dropna().sort_index()
    actual_assignment = assignment.set_index('ID')['State'].sort_index()
    assert expected_assignment.equals(actual_assignment.rename(expected_assignment.name))
    assert initial_assignment['ID'].is_unique and set(initial_assignment['ID']) == set(sample_ids)
    winner_rows.append(pd.DataFrame({'ID': sample_ids, 'CellType': cell_type, 'InitialWinner': initial_winners.reindex(sample_ids).values, 'RetainedWinner': final_winners.reindex(sample_ids).values}))

    # 只检查实际预测用到的正负展开顺序；这是结果核对，不修改W或再拟合H。
    training_genes = W.index[W.index.str.endswith('__pos')].str.removesuffix('__pos').tolist()
    expected_signed_names = [g + '__pos' for g in training_genes] + [g + '__neg' for g in training_genes]
    assert W.index.tolist() == expected_signed_names and set(training_genes) <= set(scaled.index)
    scaled = scaled.loc[training_genes, sample_ids]
    target = np.vstack([np.maximum(scaled.values, 0), np.maximum(-scaled.values, 0)])
    variable_rows = target.var(axis=1, ddof=1) > 0
    prediction = W.values[variable_rows] @ raw_H.loc[W.columns, sample_ids].values
    residual = prediction - target[variable_rows]
    relative_error = np.linalg.norm(residual) / np.linalg.norm(target[variable_rows])
    usable_weight = W.values[variable_rows].sum(axis=0) / W.values.sum(axis=0)
    retained_mass = H.sum(axis=0)
    discovery_H = pd.read_csv(discovery_dir / cell_type / 'state_abundances.txt', sep='\t', index_col=0)
    discovery_assignment = pd.read_csv(discovery_dir / cell_type / 'state_assignment.txt', sep='\t', index_col=0)
    type_rows.append({'CellType': cell_type, 'InitialRank': int(rank), 'RetainedStates': len(H), 'WinnerStates': assignment['State'].nunique(), 'AssignedSamples': len(assignment), 'UnassignedSamples': 73 - len(assignment), 'AssignedPct': len(assignment) / 73 * 100, 'MedianRetainedMass': retained_mass.median(), 'MinRetainedMass': retained_mass.min(), 'MaxRetainedMass': retained_mass.max(), 'MedianDiscoveryRetainedMass': discovery_H.sum(axis=0).median(), 'TrainingGenes': len(training_genes), 'VariableTrainingGenes': int(scaled.var(axis=1).gt(0).sum()), 'MinUsableWeightFraction': usable_weight.min(), 'RelativeProjectionError': relative_error})

    # MaxFC用于排序既有状态关联基因。bulk中的程序相关仅是内部一致性，不能独立证明谱系或准确率。
    for state in H.index:
        state_id = cell_type + '_' + state
        state_markers = gene_info.loc[gene_info['State'].eq(state)].sort_values('MaxFC', ascending=False)
        ordered_genes = state_markers['Gene'].tolist()
        present_genes = [g for g in ordered_genes if g in scaled.index]
        top_genes = present_genes[:20]
        marker_score = scaled.loc[top_genes].mean(axis=0)
        rho = spearmanr(H.loc[state], marker_score).statistic if H.loc[state].nunique() > 1 and marker_score.nunique() > 1 else np.nan
        n_assigned = int(assignment['State'].eq(state).sum())
        assigned_ids = assignment.loc[assignment['State'].eq(state), 'ID'].tolist()
        assigned_marker_mean = marker_score.loc[assigned_ids].mean() if assigned_ids else np.nan
        initial_state = mapping.set_index('State').loc[state, 'InitialState']
        state_rows.append({'CellState': state_id, 'CellType': cell_type, 'State': state, 'Ecotype': ecotypes.set_index('ID').loc[state_id, 'Ecotype'], 'InitialState': initial_state, 'WinnerSamples': n_assigned, 'WinnerPctAll73': n_assigned / 73 * 100, 'WinnerPctAssignedType': n_assigned / len(assignment) * 100, 'MedianAbundance': H.loc[state].median(), 'MaxAbundance': H.loc[state].max(), 'MeanAbundance': H.loc[state].mean(), 'AbundanceSD': H.loc[state].std(), 'DiscoveryWinnerSamples': int(discovery_assignment['State'].eq(state).sum()), 'DiscoveryMeanAbundance': discovery_H.loc[state].mean(), 'AssociatedGenes': len(ordered_genes), 'MarkersAvailableInScaledInput': len(present_genes), 'MarkerScoreSpearman_InternalOnly': rho, 'AssignedMarkerMeanZ_InternalOnly': assigned_marker_mean, 'Top10Genes': ','.join(ordered_genes[:10])})
        for gene in top_genes:
            marker_rows.append({'CellState': state_id, 'Gene': gene, 'MarkerRank': ordered_genes.index(gene) + 1, 'MeanZInWinnerSamples': scaled.loc[gene, assigned_ids].mean() if assigned_ids else np.nan, 'StateAbundanceSpearman_InternalOnly': spearmanr(H.loc[state], scaled.loc[gene]).statistic if scaled.loc[gene].nunique() > 1 else np.nan})
    renamed_H = H.copy()
    renamed_H.index = [cell_type + '_' + state for state in H.index]
    coefficient_tables.append(renamed_H)

type_audit = pd.DataFrame(type_rows)
state_audit = pd.DataFrame(state_rows)
state_audit = state_audit.merge(definitions[['CellState', 'InitialProgram', 'ProgramConfidence', 'ReviewReason', 'CellLineageValidated', 'TopPanelOverlaps']], on='CellState', validate='one_to_one')
state_audit['RecoveryDescription'] = np.where(state_audit['WinnerSamples'].eq(0), '本次无胜出样本，仍有软丰度', np.where(state_audit['WinnerSamples'].lt(3), '仅1–2例胜出，外部支持样本少', '存在胜出样本；未做独立显著性验证'))
all_H = pd.concat(coefficient_tables)
all_winners = pd.concat(winner_rows, ignore_index=True)
type_audit.to_csv(outdir / '12类状态恢复摘要.tsv', sep='\t', index=False)
state_audit.to_csv(outdir / '53状态注释与恢复审阅.tsv', sep='\t', index=False)
pd.DataFrame(marker_rows).to_csv(outdir / '状态marker内部一致性.tsv', sep='\t', index=False)
all_winners.to_csv(outdir / '73例逐类型初始与保留状态胜出分配.tsv', sep='\t', index=False)
all_H.to_csv(outdir / '53状态_73例恢复丰度.tsv', sep='\t', index_label='CellState')
print(type_audit[['CellType', 'AssignedSamples', 'WinnerStates', 'RetainedStates', 'MedianRetainedMass', 'RelativeProjectionError']].to_string(index=False), flush=True)

# Step3. 按官方规则核对生态型丰度和分配；阈值只做描述性敏感性计数，不改变正式结果。
eco_assignment = pd.read_csv(snapshot / 'ecotype_internal/initial_ecotype_assignment.txt', sep='\t', index_col=0).set_index('ID')
eco_assignment = eco_assignment.loc[sample_ids]
eco_H = pd.read_csv(resultdir / 'Ecotypes/ecotype_abundance.txt', sep='\t', index_col=0).loc[:, sample_ids]
eco_labels = ['E' + str(i) for i in range(1, 8)]
assert set(eco_H.index) == set(eco_labels) and np.isfinite(eco_H).all().all() and eco_H.ge(0).all().all()
assert np.allclose(eco_H.sum(axis=0), 1)
expected_eco = pd.DataFrame({eco: all_H.loc[ecotypes.loc[ecotypes['Ecotype'].eq(eco), 'ID']].mean(axis=0) for eco in eco_labels}).T
expected_eco = expected_eco.div(expected_eco.sum(axis=0), axis=1)
assert np.allclose(eco_H.loc[eco_labels], expected_eco.loc[eco_labels], atol=1e-12)
any_winner = all_winners.groupby('ID')['RetainedWinner'].apply(lambda x: x.notna().any()).reindex(sample_ids)
assert np.array_equal(eco_assignment['AssignedToEcotypeStates'].values, any_winner.values)
assert np.array_equal(eco_assignment['MaxEcotype'].values, eco_H.idxmax(axis=0).values)
expected_final = np.where(eco_assignment['AssignmentQ'].lt(0.25) & any_winner, eco_assignment['MaxEcotype'], 'Unassigned')
assert np.array_equal(eco_assignment['Ecotype'].values, expected_final)
initial_disc = pd.read_csv(extra_file, sep='\t', index_col=0).set_index('ID')
matched_disc = initial_disc.loc[initial_disc['Disease'].eq('UC') & initial_disc['Region3'].isin(['Colon', 'Rectum'])]
eco_rows = []
member_winner_rows = []
for eco in eco_labels + ['Unassigned']:
    members = ecotypes.loc[ecotypes['Ecotype'].eq(eco), 'ID'].tolist()
    assigned_ids = eco_assignment.index[eco_assignment['Ecotype'].eq(eco)].tolist()
    member_winner_counts = []
    for sample_id in assigned_ids:
        winners = all_winners.loc[all_winners['ID'].eq(sample_id)].dropna(subset=['RetainedWinner'])
        winner_names = (winners['CellType'] + '_' + winners['RetainedWinner']).tolist()
        member_winner_counts.append(len(set(winner_names) & set(members)))
        member_winner_rows.append({'ID': sample_id, 'Ecotype': eco, 'OwnMemberWinnerCount': len(set(winner_names) & set(members)), 'AnyRetainedWinnerCount': len(winner_names)})
    n = len(assigned_ids)
    matched_count = int(matched_disc['Ecotype'].eq(eco).sum())
    eco_rows.append({'Ecotype': eco, 'RecoverySamples': n, 'RecoveryPctAll73': n / 73 * 100, 'DiscoveryAllSamples': int(initial_disc['Ecotype'].eq(eco).sum()), 'MatchedUCColonRectumSamples': matched_count, 'MatchedUCColonRectumPct': matched_count / len(matched_disc) * 100, 'MemberStates': len(members), 'MedianMemberWinnerCount': np.median(member_winner_counts) if members and n else np.nan, 'MedianAssignmentQ': eco_assignment.loc[assigned_ids, 'AssignmentQ'].median(), 'MaxAssignmentQ': eco_assignment.loc[assigned_ids, 'AssignmentQ'].max(), 'MedianAbundance': eco_H.loc[eco, assigned_ids].median() if members and n else np.nan})
eco_audit = pd.DataFrame(eco_rows)
member_winner_audit = pd.DataFrame(member_winner_rows)
member_winner_audit.to_csv(outdir / '73例本生态型成员胜出计数.tsv', sep='\t', index=False)
eco_audit.to_csv(outdir / '7生态型恢复与匹配discovery分布.tsv', sep='\t', index=False)
eco_assignment.to_csv(outdir / '73例完整生态型分配含未分配.tsv', sep='\t', index_label='ID')
unassigned = eco_assignment.loc[eco_assignment['Ecotype'].eq('Unassigned')].copy()
unassigned['RetainedWinnerTypeCount'] = all_winners.groupby('ID')['RetainedWinner'].count().reindex(unassigned.index)
unassigned.to_csv(outdir / '8例未分配样本原因.tsv', sep='\t', index_label='ID')
q_counts = pd.DataFrame([{'AssignmentQThreshold': q, 'AssignedWithAnyRetainedWinner': int((eco_assignment['AssignmentQ'].lt(q) & any_winner).sum())} for q in [0.01, 0.05, 0.10, 0.25]])
q_counts.to_csv(outdir / '生态型Q阈值描述性计数.tsv', sep='\t', index=False)

# Step4. 既有基因程序和跨类型共变只作内部证据；用于发现共享程序及谱系待核对状态。
members_mean = all_H.T.join(eco_assignment['Ecotype']).groupby('Ecotype').mean().T.reindex(columns=eco_labels)
members_mean.to_csv(outdir / '各生态型组内平均状态丰度.tsv', sep='\t', index_label='CellState')
cross_type_state_rho = all_H.T.corr(method='spearman')
cross_type_state_rho.to_csv(outdir / '53状态跨样本相关矩阵_内部描述.tsv', sep='\t', index_label='CellState')
eco_program_rows = []
for eco in eco_labels:
    members = ecotypes.loc[ecotypes['Ecotype'].eq(eco), 'ID'].tolist()
    pair_rhos = []
    for i, first in enumerate(members):
        for second in members[i + 1:]:
            if state_audit.set_index('CellState').loc[first, 'CellType'] != state_audit.set_index('CellState').loc[second, 'CellType']:
                pair_rhos.append(cross_type_state_rho.loc[first, second])
    eco_program_rows.append({'Ecotype': eco, 'CrossTypeMemberPairs': len(pair_rhos), 'MedianCrossTypeMemberSpearman': np.median(pair_rhos) if pair_rhos else np.nan})
pd.DataFrame(eco_program_rows).to_csv(outdir / '生态型成员跨类型共变摘要.tsv', sep='\t', index=False)

# Step5. 绘制完整分母的数量摘要，以及各生态型组内状态丰度；不制造新的状态或分组。
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
figure_dir = outdir / 'figures'
figure_dir.mkdir(exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), gridspec_kw={'width_ratios': [1.6, 1]})
y = np.arange(len(type_audit))
axes[0].barh(y, type_audit['AssignedSamples'], color='#287d8e', label='Retained winning state')
axes[0].barh(y, type_audit['UnassignedSamples'], left=type_audit['AssignedSamples'], color='#d9d9d9', label='Filtered initial winner')
axes[0].set_yticks(y, type_audit['CellType'])
axes[0].invert_yaxis()
axes[0].set_xlim(0, 73)
axes[0].set_xlabel('Samples (all 73 biopsies)')
axes[0].set_title('Cell-state assignment by cell type')
for i, row in type_audit.iterrows():
    axes[0].text(row['AssignedSamples'] / 2, i, str(row['AssignedSamples']), ha='center', va='center', color='white')
fig.legend(*axes[0].get_legend_handles_labels(), loc='lower left', bbox_to_anchor=(0.08, 0.005), frameon=False, ncol=2)
colors = ['#d83b31', '#f3d74c', '#75b46a', '#df91bd', '#a87db5', '#a5a5a5', '#b8d5e4', '#ececec']
axes[1].barh(eco_audit['Ecotype'], eco_audit['RecoverySamples'], color=colors, edgecolor='#666666', linewidth=0.4)
axes[1].invert_yaxis()
axes[1].set_xlim(0, eco_audit['RecoverySamples'].max() + 5)
axes[1].set_xlabel('Samples')
axes[1].set_title('Ecotype assignment (65 / 73)')
for i, n in enumerate(eco_audit['RecoverySamples']):
    axes[1].text(n + 0.4, i, str(n), va='center')
for axis in axes:
    axis.spines[['top', 'right']].set_visible(False)
fig.suptitle('GSE72819: frozen-model recovery; 73 UC baseline biopsies', fontsize=13)
fig.tight_layout(rect=[0, 0.10, 1, 0.95])
fig.savefig(figure_dir / 'recovery_assignment_summary.png', dpi=200, bbox_inches='tight')
plt.close(fig)

state_order = ecotypes.sort_values(['Ecotype', 'CellType', 'State'])['ID'].tolist()
mean_matrix = members_mean.loc[state_order, eco_labels]
fig, ax = plt.subplots(figsize=(8.2, 13))
image = ax.imshow(mean_matrix.values, aspect='auto', cmap='viridis', vmin=0, vmax=0.4)
ax.set_xticks(np.arange(7), eco_labels)
ax.set_yticks(np.arange(53), state_order, fontsize=8)
ax.set_title('Mean recovered state abundance within assigned ecotype groups')
ax.set_xlabel('Assigned ecotype')
fig.colorbar(image, ax=ax, shrink=0.4, label='Retained-state abundance (not cell fraction)')
fig.tight_layout()
fig.savefig(figure_dir / 'state_mean_by_ecotype.png', dpi=180, bbox_inches='tight')
plt.close(fig)

# Step6. 保存可审阅摘要；无胜出不等于无表达，Q也不是独立细胞状态显著性。
summary = {
    'cohort': 'GSE72819 EUCALYPTUS UC baseline colonic biopsies; not RISK',
    'samples': 73, 'verified_manifest_files': len(manifest), 'additional_discovery_annotation_file': 1,
    'official_export_files': int(inventory['Group'].eq('官方导出').sum()),
    'cell_types': 12, 'retained_states': 53, 'states_with_at_least_one_winning_sample': int(state_audit['WinnerSamples'].gt(0).sum()),
    'states_with_fewer_than_three_winning_samples': int(state_audit['WinnerSamples'].lt(3).sum()),
    'zero_winner_states': state_audit.loc[state_audit['WinnerSamples'].eq(0), 'CellState'].tolist(),
    'ecotype_counts': dict(zip(eco_audit['Ecotype'], eco_audit['RecoverySamples'].astype(int))),
    'ecotype_assigned_samples': int(eco_assignment['Ecotype'].ne('Unassigned').sum()),
    'unassigned_with_no_retained_winner_in_any_type': int(unassigned['RetainedWinnerTypeCount'].eq(0).sum()),
    'matched_discovery_UC_colon_rectum_n': len(matched_disc),
    'final_lineage_validation_completed': False, 'recovery_significance_bootstrap_performed': False,
    'clinical_association_performed': False, 'recovery_rerun_or_model_modified': False,
    'program_confidence_counts': state_audit['ProgramConfidence'].value_counts().to_dict(),
    'median_internal_marker_rho': float(state_audit['MarkerScoreSpearman_InternalOnly'].median())
}
(outdir / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(eco_audit.to_string(index=False), flush=True)
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

# Step7. 展开53个状态的可阅读注释附录；程序名沿用既有初步定义，不把恢复当作谱系验证。
appendix = ['# GSE72819：53个状态的初步注释与恢复附录', '', '每个样本在每个类型内最多有一个胜出状态；同一样本可以在多个类型内胜出。表中的胜出例数以全部73例为分母。', '', '程序置信度仅表示既有基因程序命名的初步支持，全部状态的细胞谱系都尚未完成独立验证。关联基因与H的相关性是内部证据，不是准确率或独立显著性。', '']
for cell_type in ranks['CellType']:
    appendix.extend(['## ' + cell_type, '', '| 状态 | Ecotype | 胜出例数 | 程序初步名称 | 程序置信度 | 前10关联基因 |', '|---|---|---:|---|---|---|'])
    type_states = state_audit.loc[state_audit['CellType'].eq(cell_type)]
    for _, row in type_states.iterrows():
        appendix.append(f"| {row['CellState']} | {row['Ecotype']} | {row['WinnerSamples']} | {row['InitialProgram']} | {row['ProgramConfidence']} | {row['Top10Genes']} |")
    appendix.append('')
    for _, row in type_states.iterrows():
        appendix.append(f"- **{row['CellState']}**：{row['ReviewReason']}。恢复说明：{row['RecoveryDescription']}。")
    appendix.append('')
(outdir / '53状态注释与恢复附录.md').write_text('\n'.join(appendix) + '\n', encoding='utf-8')
