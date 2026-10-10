# 独立核对服务器R结果的矩阵分母、四套相关、置换Z和输入/回传校验。
# 本机cibersortx_review环境运行；只读本轮输出，写核查摘要，不再拟合模型。
# 用Python重新核对保存结果，防止行列错配和方法实现错误。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 核对回传清单、实际包版本与冻结记录数
project = Path(__file__).resolve().parents[3] / 'IBD_ecotype_analysis'
out = project / '06_validation/ecotyper_original_method_annotation_20261010'
manifest = json.loads((out / 'download_manifest.json').read_text(encoding='utf8'))
for item in manifest:
    blob = (out / item['File']).read_bytes()
    assert len(blob)==item['Bytes'] and hashlib.sha256(blob).hexdigest()==item['SHA256'], item['File']
versions = pd.read_csv(out / 'environment/smoke_actual_versions.tsv', sep='\t').set_index('Package').Version
assert versions['Seurat']=='3.1.3' and versions['fgsea']=='1.14.0'
assert (out / 'analysis.exit').read_text().strip()=='0'
assert (out / 'environment/numerical_smoke_completed.txt').read_text().strip()=='NUMERICAL_SMOKE_OK'
states = pd.read_csv(out / 'inputs/frozen_states.tsv', sep='\t')
ecotypes = pd.read_csv(out / 'inputs/frozen_ecotypes.tsv', sep='\t')
assert len(states)==53 and len(ecotypes)==7 and states.CellState.is_unique

# Step2. 从原始分配表重新核对F1/F2，而不是只核对文件互相一致
results = pd.read_csv(out / 'ecotype_coassociation_both_datasets.tsv', sep='\t')
checks = []
for dataset, folder in [('Colon_GSE214695', 'GSE214695_scRNA_recovery_20261009'), ('Ileum_IBDverse', 'IBDverse_scRNA_recovery_20261009')]:
    cells = pd.read_csv(project / '06_validation' / folder / 'raw_server/audit/cell_assignments_all.tsv', sep='\t')
    cells['CellState'] = cells.CellType + '_' + cells.State
    f1 = pd.read_csv(out / (dataset + '_F1_assigned_denominator.tsv'), sep='\t', index_col=0)
    f2 = pd.read_csv(out / (dataset + '_F2_all_cells_denominator.tsv'), sep='\t', index_col=0)
    for state in f1.index:
        cell_type = states.set_index('CellState').loc[state, 'CellType']
        parent = cells[cells.CellType.eq(cell_type)]
        numerator = parent[parent.CellState.eq(state)].groupby('Sample').size().reindex(f1.columns, fill_value=0)
        d1 = parent[parent.State.ne('Unassigned')].groupby('Sample').size().reindex(f1.columns, fill_value=0)
        d2 = parent.groupby('Sample').size().reindex(f1.columns, fill_value=0)
        expected1 = numerator / d1.replace(0, np.nan)
        expected2 = numerator / d2.replace(0, np.nan)
        np.testing.assert_allclose(f1.loc[state], expected1, atol=1e-12, equal_nan=True)
        np.testing.assert_allclose(f2.loc[state], expected2, atol=1e-12, equal_nan=True)

    # Step3. Python独立重算Spearman，与服务器R结果逐边比较
    matrices = [f1, f2, f1.replace(0, np.nan), f2.replace(0, np.nan)]
    recomputed = [x.T.corr(method='spearman', min_periods=2).to_numpy() for x in matrices]
    c = pd.read_csv(out / (dataset + '_C_four_matrix_mean.tsv'), sep='\t', index_col=0)
    expected = sum(recomputed) / 4
    np.fill_diagonal(expected, np.nan)
    np.testing.assert_allclose(c, expected, atol=1e-12, equal_nan=True)
    np.testing.assert_allclose(c, c.T, atol=1e-12, equal_nan=True)
    assert np.isnan(np.diag(c)).all()

    # Step4. 核对保存的零分布、Z和实际可计算成员
    null = pd.read_csv(out / (dataset + '_permutation_null.tsv'), sep='\t')
    assert null.shape==(1000, 7)
    for row in results[results.Dataset.eq(dataset)].itertuples(index=False):
        members = ecotypes.set_index('Ecotype').loc[row.Ecotype, 'Members'].split(';')
        inside = c.index.intersection(members)
        outside = c.columns.difference(inside)
        delta = c.loc[inside, inside].mean(axis=1) - c.loc[inside, outside].mean(axis=1)
        np.testing.assert_allclose(delta.mean(), row.AvgDif, atol=1e-12, equal_nan=True)
        assert delta.notna().sum()==row.DefinedReferenceMembers
        np.testing.assert_allclose((row.AvgDif-null[row.Ecotype].mean()) / null[row.Ecotype].std(ddof=1), row.Z, atol=1e-12, equal_nan=True)
    checks.append({'Dataset':dataset, 'InputCells':len(cells), 'Samples':cells.Sample.nunique(), 'StateRows':len(c), 'DefinedUndirectedEdges':int(np.isfinite(c.to_numpy()).sum()//2), 'DenominatorsChecked':True, 'R_PythonSpearmanAgreement':True, 'NullIterations':1000})

# Step5. 核对冻结输入校验，写入可审阅摘要
input_hashes = pd.read_csv(out / 'input_md5_unchanged.tsv', sep='\t')
for row in input_hashes.itertuples(index=False):
    path = out / 'inputs' / Path(row.File).name
    assert hashlib.md5(path.read_bytes()).hexdigest()==row.MD5
eligibility = pd.read_csv(out / 'marker_eligibility_both_datasets.tsv', sep='\t')
pivot = eligibility.pivot(index='CellState', columns='Dataset', values='PaperMarkerEligible')
summary = {'Checks':'PASS', 'VerifiedTransferredFiles':len(manifest), 'Seurat':versions['Seurat'], 'fgsea':versions['fgsea'], 'EligibleBoth':int(pivot.all(axis=1).sum()), 'EligibleOne':int(pivot.sum(axis=1).eq(1).sum()), 'EligibleNeither':int(pivot.sum(axis=1).eq(0).sum()), 'Datasets':checks}
(out / '独立结果核查.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
