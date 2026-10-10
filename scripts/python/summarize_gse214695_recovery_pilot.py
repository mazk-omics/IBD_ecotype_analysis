# 汇总官方 GSE214695 单细胞 recovery 的状态分配、Z 分数、供者覆盖及冻结 marker 的表达。
# 在服务器 ibd_refaudit 环境运行；读取正式输出和冻结发现模型，写出便于本地绘图的小表。
# 不重新拟合，不删除失败状态；未分配细胞和缺失值保留在统计中。
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 读取正式输入记录和官方内部结果
root = Path('/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009')
assert (root / 'formal_exit_code.txt').read_text().strip() == '0'
output = root / 'audit'
output.mkdir(exist_ok=True)
annotation = pd.read_csv(root / 'input/GSE214695_Pilot18_20261009_annotation.tsv', sep='\t')
annotation = annotation.set_index('ID', drop=False)
model = Path('/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes')
recovery = 'GSE214695_Pilot18_20261009_CPM'
key = pd.read_csv(model / 'Analysis/rank_selection/rank_data.txt', sep='\t')
key = key[key.iloc[:, 0].isin(annotation.CellType)]
ecotypes = pd.read_csv(model / 'Ecotypes/discovery/ecotypes.txt', sep='\t')
ecotypes.to_csv(output / 'frozen_ecotype_members.tsv', sep='\t', index=False)
state_rows = []
cell_rows = []
donor_rows = []
program_rows = []
marker_rows = []

# Step2. 按类型逐一核对初始状态和保留状态；初始状态被过滤的细胞属于未分配
for cell_type, rank in key.iloc[:, :2].itertuples(index=False, name=None):
    discovery_dir = model / f'Cell_States/discovery/{cell_type}/{rank}'
    rec = model / f'Cell_States/recovery/scRNA/{recovery}/{cell_type}/{rank}'
    mapping = pd.read_csv(discovery_dir / 'mapping_to_initial_states.txt', sep='\t')
    initial = pd.read_csv(rec / 'initial_state_assignment.txt', sep='\t')
    initial = initial.rename(columns={'State': 'InitialState'})
    assignment = initial.merge(mapping, on='InitialState', how='left', validate='many_to_one')
    assignment['State'] = assignment.State.fillna('Unassigned')
    cells = annotation.loc[assignment.ID].reset_index(drop=True).merge(assignment, on='ID', validate='one_to_one')
    assert len(cells) == int(annotation.CellType.eq(cell_type).sum())
    cell_rows.append(cells)
    z = pd.read_csv(rec / 'recovery_z_scores.txt', sep='\t')
    z.to_csv(output / f'{cell_type}_official_z.tsv', sep='\t', index=False)
    scaled = pd.read_csv(rec / 'expression_matrix_scaled.txt', sep='\t', index_col=0)
    assert np.isfinite(scaled.to_numpy()).all()
    scaled = scaled.loc[:, cells.ID]
    log_expression = pd.read_csv(rec / 'expression_matrix.txt', sep='\t', index_col=0).loc[:, cells.ID]
    marker_info = pd.read_csv(discovery_dir / 'gene_info.txt', sep='\t', index_col=0)
    disc_expr = pd.read_csv(discovery_dir / 'heatmap_data.txt', sep='\t', index_col=0)
    disc_ann = pd.read_csv(discovery_dir / 'heatmap_top_ann.txt', sep='\t', index_col=0)
    disc_ann = disc_ann.loc[disc_expr.columns]
    disc_sd = disc_expr.std(axis=1, ddof=1).replace(0, np.nan)
    disc_scaled = disc_expr.sub(disc_expr.mean(axis=1), axis=0).div(disc_sd, axis=0).fillna(0)

    for state in mapping.State:
        state_id = f'{cell_type}_{state}'
        win = cells.State.eq(state)
        z_row = z[z.State.eq(state)].iloc[0].to_dict()
        marker_genes = marker_info.loc[marker_info.State.eq(state), 'Gene'].astype(str)
        marker_detected = log_expression.loc[marker_genes].gt(0)
        row = {'CellState': state_id, 'CellType': cell_type, 'State': state, 'InputCells': len(cells), 'InputDonors': cells.Donor.nunique(), 'AssignedCells': int(win.sum()), 'AssignedDonors': cells.loc[win, 'Donor'].nunique(), 'AssignedFraction': float(win.mean()), 'MarkerGenes': len(marker_genes), 'MarkerGenesDetectedInLineage': int(marker_detected.any(axis=1).sum()), 'MarkerMeanDetectionInAssignedCells': float(marker_detected.iloc[:, win.to_numpy()].to_numpy().mean()) if win.any() else np.nan, **z_row}
        state_rows.append(row)
        for donor, group in cells.groupby('Donor', sort=True):
            donor_rows.append({'CellState': state_id, 'CellType': cell_type, 'State': state, 'Donor': donor, 'Sample': group.Sample.iloc[0], 'Stratum': group.Stratum.iloc[0], 'TypeCells': len(group), 'AssignedCells': int(group.State.eq(state).sum()), 'AssignedFraction': float(group.State.eq(state).mean())})

        genes = marker_info.loc[marker_info.State.eq(state), 'Gene'].astype(str)
        genes = genes[genes.isin(scaled.index) & genes.isin(disc_scaled.index)]
        top_genes = marker_info[marker_info.State.eq(state)].sort_values('MaxFC', ascending=False).Gene.astype(str).head(5)
        for target in mapping.State:
            rec_mask = cells.State.eq(target).to_numpy()
            disc_mask = disc_ann.State.eq(target).to_numpy()
            rec_mean = float(scaled.loc[genes].iloc[:, rec_mask].to_numpy().mean()) if rec_mask.any() and len(genes) else np.nan
            disc_mean = float(disc_scaled.loc[genes].iloc[:, disc_mask].to_numpy().mean()) if disc_mask.any() and len(genes) else np.nan
            program_rows.append({'CellType': cell_type, 'MarkerState': state_id, 'TargetState': f'{cell_type}_{target}', 'MarkerGenes': len(genes), 'DiscoveryMeanZ': disc_mean, 'RecoveryMeanZ': rec_mean})
            for gene in top_genes:
                if gene not in scaled.index or gene not in disc_scaled.index:
                    continue
                marker_rows.append({'CellType': cell_type, 'MarkerState': state_id, 'TargetState': f'{cell_type}_{target}', 'Gene': gene, 'DiscoveryMeanZ': float(disc_scaled.loc[gene].iloc[disc_mask].mean()) if disc_mask.any() else np.nan, 'RecoveryMeanZ': float(scaled.loc[gene].iloc[rec_mask].mean()) if rec_mask.any() else np.nan})
    print(cell_type, 'input', len(cells), 'assigned', int(cells.State.ne('Unassigned').sum()), 'states', len(mapping), flush=True)

# Step3. 保存中间表和未经重新命名的官方生态型输出
states = pd.DataFrame(state_rows)
cells = pd.concat(cell_rows, ignore_index=True)
donors = pd.DataFrame(donor_rows)
states.to_csv(output / 'state_recovery_summary.tsv', sep='\t', index=False)
cells.to_csv(output / 'cell_assignments_all.tsv', sep='\t', index=False)
donors.to_csv(output / 'donor_state_assignments.tsv', sep='\t', index=False)
pd.DataFrame(program_rows).to_csv(output / 'state_program_expression.tsv', sep='\t', index=False)
pd.DataFrame(marker_rows).to_csv(output / 'top5_marker_expression.tsv', sep='\t', index=False)
cells.groupby(['CellType', 'State']).size().to_csv(output / 'assignment_counts_including_unassigned.tsv', sep='\t')
eco_dir = model / f'Ecotypes/recovery/{recovery}'
for filename in ['ecotype_assignment.txt', 'ecotype_abundance.txt', 'heatmap_all_samples.pdf']:
    shutil.copy2(eco_dir / filename, output / filename)
for filename in ['smoke.log', 'formal.log', 'source_before.sha256', 'source_after.sha256', 'run_started.txt', 'run_finished.txt']:
    shutil.copy2(root / filename, output / filename)
summary = {'cells': len(cells), 'donors': cells.Donor.nunique(), 'types': cells.CellType.nunique(), 'states': len(states), 'states_with_assigned_cells': int(states.AssignedCells.gt(0).sum()), 'states_with_finite_z': int(np.isfinite(states.Z).sum()), 'states_above_nominal_z_1_96': int(states.Z.gt(1.96).sum()), 'assigned_cells': int(cells.State.ne('Unassigned').sum()), 'official_shuffle_repeats': 30, 'source_unchanged': (root / 'source_before.sha256').read_bytes() == (root / 'source_after.sha256').read_bytes()}
(output / 'recovery_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(summary, flush=True)
