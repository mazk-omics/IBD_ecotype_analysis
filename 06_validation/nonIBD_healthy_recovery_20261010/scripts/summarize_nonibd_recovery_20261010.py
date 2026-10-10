# 审阅 nonIBD/健康 recovery 的官方输出，保留全部冻结状态，包括缺失谱系和未分配状态。
# 在服务器 ibd_refaudit 环境运行；读取官方结果及新输入，输出小表供本地阅读。
# 不重新命名状态，不把最高状态分配或内部 marker 一致性当成独立验证准确率。
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 检查任务结束和模型指纹，读取冻结状态及生态型
parser = argparse.ArgumentParser()
parser.add_argument('--root', required=True)
args = parser.parse_args()
root = Path(args.root)
assert (root / 'run_exit_code.txt').read_text().strip() == '0'
out = root / 'audit'
out.mkdir(exist_ok=True)
model = Path('/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes')
rank = pd.read_csv(model / 'Analysis/rank_selection/rank_data.txt', sep='\t')
members = pd.read_csv(model / 'Ecotypes/discovery/ecotypes.txt', sep='\t').rename(columns={'ID': 'CellState'})
states = members[['CellState', 'CellType', 'State', 'Ecotype']].copy()
assert len(states) == 53 and states.CellState.is_unique
states.to_csv(out / 'frozen_states.tsv', sep='\t', index=False)
eco = members.groupby('Ecotype', sort=True).agg(MemberStates=('CellState', 'size'), Members=('CellState', ';'.join)).reset_index()
eco.to_csv(out / 'frozen_ecotypes.tsv', sep='\t', index=False)
manifest = pd.read_csv(root / 'input/run_manifest.tsv', sep='\t')
summary_rows = []
state_rows = []
eco_rows = []
output_issues = []

# Step2. 每套输入分别计算状态的分配量、供者数和内部 marker 方向
for task in manifest.itertuples(index=False):
    name = task.Dataset
    assert (root / (name + '_exit_code.txt')).read_text().strip() == '0'
    a = pd.read_csv(root / 'input' / (name + '_annotation.tsv'), sep='\t', dtype=str)
    a['DonorKey'] = a.Donor if task.Assay == 'scRNA' else a.PatientID if 'PatientID' in a else a.ID
    dst = out / name
    dst.mkdir(exist_ok=True)
    all_assignments = []
    for cell_type, n_states in rank.iloc[:, :2].itertuples(index=False, name=None):
        frozen = states.loc[states.CellType.eq(cell_type)]
        selected = a.loc[a.CellType.eq(cell_type)].copy() if task.Assay == 'scRNA' else a.copy()
        rec = model / 'Cell_States/recovery' / 'scRNA' / name / cell_type / str(n_states) if task.Assay == 'scRNA' else model / 'Cell_States/recovery' / name / cell_type / str(n_states)
        if selected.empty:
            for state in frozen.itertuples(index=False):
                state_rows.append({'Dataset': name, 'Assay': task.Assay, 'CellState': state.CellState, 'CellType': cell_type, 'State': state.State, 'Ecotype': state.Ecotype, 'LineagePresent': False, 'InputUnits': 0, 'InputDonors': 0, 'AssignedUnits': np.nan, 'AssignedDonors': np.nan, 'Z': np.nan, 'OwnMinusBestOther_InternalOnly': np.nan})
            continue
        assignment_file = rec / 'state_assignment.txt'
        assert assignment_file.exists(), f'Missing official assignment: {rec}'
        assignment = pd.read_csv(assignment_file, sep='\t')[['ID', 'State']]
        assert assignment.ID.is_unique and set(assignment.ID) <= set(selected.ID)
        selected = selected.merge(assignment, on='ID', how='left', validate='one_to_one')
        selected['State'] = selected.State.fillna('Unassigned')
        selected['CellType'] = cell_type
        selected['CellState'] = np.where(selected.State.eq('Unassigned'), 'Unassigned', cell_type + '_' + selected.State)
        all_assignments.append(selected)
        scaled = pd.read_csv(rec / 'expression_matrix_scaled.txt', sep='\t', index_col=0).loc[:, selected.ID]
        marker_info = pd.read_csv(model / 'Cell_States/discovery' / cell_type / str(n_states) / 'gene_info.txt', sep='\t')
        z_file = rec / 'recovery_z_scores.txt'
        z_table = pd.read_csv(z_file, sep='\t') if task.Assay == 'scRNA' and z_file.exists() else pd.DataFrame(columns=['State', 'Z'])
        formal_type = root / 'formal' / name / cell_type
        required = ['state_assignment.txt', 'state_assignment_heatmap.pdf'] + (['recovery_z_scores.txt'] if task.Assay == 'scRNA' else ['state_abundances.txt'])
        for filename in required:
            if not (formal_type / filename).exists():
                output_issues.append({'Dataset': name, 'CellType': cell_type, 'MissingFormalOutput': filename})
        for state in frozen.itertuples(index=False):
            win = selected.State.eq(state.State)
            genes = marker_info.loc[marker_info.State.eq(state.State), 'Gene']
            genes = genes[genes.isin(scaled.index)]
            scores = scaled.loc[genes].mean(axis=0) if len(genes) else pd.Series(np.nan, index=scaled.columns)
            grouped = scores.groupby(selected.set_index('ID').State.reindex(scores.index)).mean()
            other = grouped.drop([state.State, 'Unassigned'], errors='ignore')
            own = grouped.get(state.State, np.nan)
            z = z_table.loc[z_table.State.eq(state.State), 'Z']
            state_rows.append({'Dataset': name, 'Assay': task.Assay, 'CellState': state.CellState, 'CellType': cell_type, 'State': state.State, 'Ecotype': state.Ecotype, 'LineagePresent': True, 'InputUnits': len(selected), 'InputDonors': selected.DonorKey.nunique(), 'AssignedUnits': int(win.sum()), 'AssignedDonors': selected.loc[win, 'DonorKey'].nunique(), 'Z': z.iloc[0] if len(z) else np.nan, 'OwnMinusBestOther_InternalOnly': own - other.max() if len(other) else np.nan})
    assignments = pd.concat(all_assignments, ignore_index=True)
    assignments.to_csv(dst / 'unit_assignments_all.tsv', sep='\t', index=False)
    if task.Assay == 'scRNA':
        assignments.rename(columns={'DonorKey': 'RecoveryDonor'}).to_csv(dst / 'cell_assignments_all.tsv', sep='\t', index=False)
    donor_states = assignments.groupby(['CellType', 'State', 'DonorKey']).size().rename('AssignedUnits').reset_index()
    donor_states.to_csv(dst / 'donor_state_assignments.tsv', sep='\t', index=False)

    # Step3. 生态型官方分配与最高分配分别保存，未列出的样本保留未分配
    rec_eco = model / 'Ecotypes/recovery' / name
    official = pd.read_csv(rec_eco / 'ecotype_assignment.txt', sep='\t', index_col=0)
    assert official.ID.is_unique
    units = a[['Sample', 'DonorKey']].drop_duplicates().rename(columns={'Sample': 'ID'}) if task.Assay == 'scRNA' else a[['ID', 'DonorKey']].copy()
    assert units.ID.is_unique
    units = units.merge(official, on='ID', how='left', validate='one_to_one')
    units['Ecotype'] = units.Ecotype.fillna('Unassigned')
    abundance = pd.read_csv(rec_eco / 'ecotype_abundance.txt', sep='\t', index_col=0)
    assert set(units.ID) <= set(abundance.columns)
    units['HighestEcotypeAllInputs'] = abundance.loc[:, units.ID].idxmax(axis=0).reindex(units.ID).to_numpy()
    units.to_csv(dst / 'ecotype_assignments_complete.tsv', sep='\t', index=False)
    for e in list(eco.Ecotype) + ['Unassigned']:
        chosen = units.Ecotype.eq(e)
        eco_rows.append({'Dataset': name, 'Assay': task.Assay, 'Ecotype': e, 'InputBiopsies': len(units), 'AssignedBiopsies': int(chosen.sum()), 'AssignedDonors': units.loc[chosen, 'DonorKey'].nunique(), 'HighestEcotypeBiopsies': int(units.HighestEcotypeAllInputs.eq(e).sum()), 'AssignedWithQBelow025': int((chosen & units.AssignmentQ.lt(0.25)).sum())})
    current = pd.DataFrame(state_rows).query('Dataset == @name')
    summary_rows.append({'Dataset': name, 'Assay': task.Assay, 'InputUnits': len(a), 'InputDonors': a.DonorKey.nunique(), 'LineagesPresent': int(current.loc[current.LineagePresent, 'CellType'].nunique()), 'FrozenStates': len(current), 'StatesWithAssignments': int(current.AssignedUnits.gt(0).sum()), 'StatesWithAtLeast5Units': int(current.AssignedUnits.ge(5).sum()), 'StatesWithPositiveInternalMarkerDelta': int(current.OwnMinusBestOther_InternalOnly.gt(0).sum()), 'StatesWithFiniteOfficialZ': int(np.isfinite(current.Z).sum()), 'StatesAboveNominalZ196': int(current.Z.gt(1.96).sum()) if task.Assay == 'scRNA' else np.nan, 'EcotypeAssignedBiopsies': int(units.Ecotype.ne('Unassigned').sum()), 'AssignedUnitsAcrossTypes': int(assignments.State.ne('Unassigned').sum())})
    for f in ['ecotype_assignment.txt', 'ecotype_abundance.txt', 'heatmap_all_samples.pdf']:
        if (rec_eco / f).exists():
            shutil.copy2(rec_eco / f, dst / f)
    print(name, summary_rows[-1], flush=True)

# Step4. 保存汇总和校验；单细胞名义 Z 不等于多重检验后的显著性
pd.DataFrame(state_rows).to_csv(out / '53_states_all_datasets.tsv', sep='\t', index=False)
pd.DataFrame(eco_rows).to_csv(out / '7_ecotypes_all_datasets.tsv', sep='\t', index=False)
pd.DataFrame(summary_rows).to_csv(out / 'dataset_summary.tsv', sep='\t', index=False)
pd.DataFrame(output_issues, columns=['Dataset', 'CellType', 'MissingFormalOutput']).to_csv(out / 'formal_output_issues.tsv', sep='\t', index=False)
hashes = {}
for group in ['frozen_model_source_before', 'source_inputs_before']:
    before = json.loads((root / 'input' / (group + '.json')).read_text())
    changed = [p for p, digest in before.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest() != digest]
    assert not changed, changed
    hashes[group] = {'FilesChecked': len(before), 'Unchanged': True}
(out / 'integrity_check.json').write_text(json.dumps({'Status': 'PASS' if not output_issues else 'INCOMPLETE_FORMAL_OUTPUTS', 'Fingerprints': hashes, 'Datasets': len(manifest), 'FrozenStatesPerDataset': 53, 'NoClinicalAnalysis': True, 'NoStateRenaming': True, 'BulkOfficialStateSignificance': False, 'ScRNAOfficialShuffles': 30}, indent=2), encoding='utf-8')
