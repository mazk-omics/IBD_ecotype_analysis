# 在服务器ibd_outcomes_20261015环境复核本轮Recovery与GEO原始值。
# 这里只比较投影结果、源文件与冻结W，不重新训练或更改分配。
from pathlib import Path
import hashlib, tarfile, gzip, io, shutil, json
import pandas as pd
import numpy as np

# Step1. 核对RISK逐个GEO原始样本的存档值，证明负数不是读取错列。
root = Path('/home/mazekai/IBD_EcoTyper')
out = root / '06_validation/external_clinical_outcomes_20261015'
risk = root / '06_validation/01_bulk_progression/GSE93624'
x = pd.read_csv(risk / 'prepared/06_ecotyper_recovery_input/GSE93624_EcoTyper_recovery_expression.txt', sep='\t', index_col=0)
mapping = pd.read_csv(risk / 'prepared/02_gene_harmonization/GSE93624_canonical_gene_included_features_FINAL.tsv', sep='\t')
print('RISK gene mapping columns', mapping.columns.tolist(), flush=True)
raw_rows = []
with tarfile.open(risk / 'raw/GSE93624_RAW.tar') as archive:
    for member in archive.getmembers():
        if not member.isfile():
            continue
        b = archive.extractfile(member).read()
        if member.name.endswith('.gz'):
            b = gzip.decompress(b)
        d = pd.read_csv(io.BytesIO(b), sep='\t', index_col=0)
        sample = d.columns[0]
        shared = x.index.intersection(d.index)
        diff = np.abs(x.loc[shared, sample].to_numpy() - d.loc[shared].iloc[:, 0].to_numpy())
        assert len(shared) > 10000 and diff.max() < 1e-10
        mapped = d.loc[mapping.SourceSymbol].iloc[:, 0].copy()
        mapped.index = mapping.FinalSymbol
        mapped_diff = np.abs(x.loc[mapped.index, sample].to_numpy() - mapped.to_numpy())
        assert mapping.FinalSymbol.is_unique and mapped_diff.max() < 1e-10 and len(mapped) == 13323
        raw_rows.append({'RawFile': member.name, 'SampleID': sample, 'RawGenes': len(d), 'UnchangedDirectSymbolGenes': len(shared), 'MaxAbsDifferenceDirectGenes': diff.max(), 'AllMappedGenes': len(mapped), 'MaxAbsDifferenceAllMappedGenes': mapped_diff.max(), 'RawNegativeN': int((d.iloc[:, 0] < 0).sum()), 'InputNegativeN': int((x[sample] < 0).sum())})
pd.DataFrame(raw_rows).to_csv(out / '02_RISK_disease_progression/deposited_sample_value_verification.tsv', sep='\t', index=False)
assert len(raw_rows) == 245

# Step2. 审计W实际使用的基因覆盖、保留H与初始H的关系及标准化结果。
model = root / 'software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes'
ranks = pd.read_csv(model / 'Analysis/rank_selection/rank_data.txt', sep='\t')
qc = []
for acc, dataset, source in [('GSE72819', 'GSE72819_BP73_20261009_185536', root / '06_validation/02_bulk_treatment/GSE72819/prepared/06_ecotyper_recovery_input/GSE72819_bulk_recovery_expression.txt'), ('GSE73661', 'GSE73661_RMA_W0_67_outcomes_20261015', out / '03_GSE73661_mucosal_healing/GSE73661_RMA_W0_67_outcomes_20261015.txt'), ('GSE93624', 'GSE93624_processed245_outcomes_20261015', out / '02_RISK_disease_progression/GSE93624_processed245_outcomes_20261015.txt')]:
    expression = pd.read_csv(source, sep='\t', index_col=0)
    if acc == 'GSE73661':
        metadata = pd.read_csv(root / '06_validation/02_bulk_treatment/GSE73661/prepared/04_canonical_metadata/GSE73661_canonical_metadata_178_samples.tsv', sep='\t')
        baseline_ids = metadata.loc[metadata.Timepoint.eq('W0'), 'GSM']
        assert expression.shape[1] == 67 and set(expression.columns) == set(baseline_ids)
    if expression.ge(0).all().all() and expression.max().max() > 50:
        expression = np.log2(expression + 1)
    for _, rank in ranks.iterrows():
        parent, k = rank.iloc[0], str(rank.iloc[1])
        frozen = model / 'Cell_States/discovery' / parent / k
        internal = out / 'recovery_runtime/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes/Cell_States/recovery' / dataset / parent / k
        if acc == 'GSE72819':
            internal = root / '06_validation/02_bulk_treatment/GSE72819/official_recovery_20261009_185536/audit_snapshot/state_internal' / parent / k
        W = pd.read_csv(frozen / 'W.txt', sep='\t', index_col=0)
        genes = W.index[W.index.str.endswith('__pos')].str.removesuffix('__pos')
        present = genes.isin(expression.index)
        scaled = pd.read_csv(internal / 'expression_matrix_scaled.txt', sep='\t', index_col=0)
        expected = expression.reindex(scaled.index).fillna(0)
        expected = expected.sub(expected.mean(axis=1), axis=0).div(expected.std(axis=1, ddof=1), axis=0).fillna(0)
        assert np.allclose(scaled.loc[expected.index, expected.columns], expected, atol=1e-10)
        initial = pd.read_csv(internal / 'initial_state_abundances.txt', sep='\t', index_col=0)
        retained = pd.read_csv(internal / 'state_abundances.txt', sep='\t', index_col=0)
        m = pd.read_csv(frozen / 'mapping_to_initial_states.txt', sep='\t')
        expected_H = initial.loc[m.InitialState].copy()
        expected_H.index = m.State
        assert np.allclose(initial.sum(axis=0), 1) and np.allclose(expected_H, retained)
        target = scaled.reindex(genes).fillna(0).to_numpy()
        signed = np.vstack([np.maximum(target, 0), np.maximum(-target, 0)])
        variable = signed.var(axis=1, ddof=1) > 0
        assert list(W.index) == [g + '__pos' for g in genes] + [g + '__neg' for g in genes]
        usable_weight = W.to_numpy()[variable].sum(axis=0) / W.sum(axis=0).to_numpy()
        qc.append({'Accession': acc, 'ParentLineage': parent, 'ModelGenes': len(genes), 'PresentGenes': int(present.sum()), 'CoverageFraction': float(present.mean()), 'VariableSignedRows': int(variable.sum()), 'MinUsableWWeightFraction': usable_weight.min(), 'RetainedStates': len(retained), 'MedianRetainedHMass': retained.sum(axis=0).median(), 'StandardizationMatchesOfficialCode': True, 'FrozenStateMappingMatches': True})
pd.DataFrame(qc).to_csv(out / 'audit/model_gene_coverage_and_projection_QC.tsv', sep='\t', index=False)

# Step3. 保存代码与输入完整性，所有冻结W与模型配置必须未变。
before = pd.read_csv(out / 'audit/immutable_input_hashes_before.tsv', sep='\t')
before['SHA256_after'] = [hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in before.Path]
before['Unchanged'] = before.SHA256.eq(before.SHA256_after)
assert before.Unchanged.all()
before.to_csv(out / 'audit/immutable_input_hashes_verified.tsv', sep='\t', index=False)
for name in ['state_recovery_bulk.R', 'ecotypes_recovery.R', 'lib/ecotyper.R']:
    f = out / 'recovery_runtime/pipeline' / name
    original = root / 'software/ecotyper/pipeline' / name
    assert f.read_bytes() == original.read_bytes()
shutil.copy2(risk / 'prepared/01_authoritative_matrix/01b_build_and_audit_authoritative_245_matrix.log', out / '02_RISK_disease_progression/prior_server_scale_audit.log')
shutil.copy2(risk / 'prepared/00_inventory/GSE93624_data_processing_unique.txt', out / '02_RISK_disease_progression/GEO_data_processing.txt')
print('QC_COMPLETE', pd.DataFrame(qc).groupby('Accession').CoverageFraction.agg(['min', 'median', 'max']).to_dict(), flush=True)
