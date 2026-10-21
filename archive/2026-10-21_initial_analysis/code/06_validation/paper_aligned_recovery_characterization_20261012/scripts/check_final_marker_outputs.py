# 独立检查两版n4的维度、七项等权rank、留出泄漏和原模型完整性；不重新拟合模型。
# 服务器ibd_refaudit或本机cibersortx_review运行，读取已完成排名及检查点，保存核查JSON。
# 用原始指标重建rank，用完整53状态×两队列×两版检查不可评价记录。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 核验训练数据集隔离和n4的可评价上限
out = Path(__file__).resolve().parents[1]
metrics = pd.read_csv(out / '02_external_markers/marker_n1_n7.tsv', sep='\t')
ranking = pd.read_csv(out / '02_external_markers/marker_composite_ranking.tsv', sep='\t')
loocv = pd.read_csv(out / '02_external_markers/marker_LOOCV_per_dataset.tsv', sep='\t')
assert (metrics.n2<=metrics.n1).all()
np.testing.assert_allclose(metrics.n3, metrics.n2/metrics.n1)
for variant in ['DiseaseGroup', 'Region3']:
    evaluated = metrics['n4_'+variant].notna()
    assert (metrics.loc[evaluated, 'n4_'+variant]<=metrics.loc[evaluated, 'Evaluable_'+variant]).all()
for held_out, rows in ranking.groupby('HeldOut'):
    if held_out!='All':
        assert rows.TrainingDatasets.map(lambda x: held_out not in x.split(';')).all()
assert len(loocv)==212 and not loocv[['HeldOut', 'Variant', 'CellState']].duplicated().any()
assert loocv.groupby(['HeldOut', 'Variant']).size().eq(53).all()
assert loocv.loc[loocv.Status.eq('EVALUABLE'), 'PresentMarkers'].eq(10).all()

# Step2. 从原始七指标重新求平均rank，与保存分数逐基因核对
columns = ['n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7']
for key, rows in ranking.groupby(['HeldOut', 'Variant', 'CellState']):
    assert rows[columns].notna().all(axis=None)
    expected = rows[columns].rank(method='average', ascending=True).mean(axis=1)
    np.testing.assert_allclose(expected, rows.CompositeScore, atol=1e-10)
    assert rows.WithinStateOrder.min()==1 and rows.WithinStateOrder.max()==len(rows)

# Step3. 服务器再次核对冻结NMF和成员文件；本机没有服务器原件时不冒称此项通过
source_check = 'NOT_RUN_LOCAL'
if Path('/home/mazekai/IBD_EcoTyper').exists():
    verified = 0
    for name in ['frozen_model_manifest.tsv', 'NMF_selected_model_manifest.tsv']:
        manifest = pd.read_csv(out / '00_audit' / name, sep='\t')
        for row in manifest.itertuples():
            assert hashlib.sha256(Path(row.File).read_bytes()).hexdigest()==row.SHA256, row.File
            verified += 1
    source_check = {'UnchangedFrozenModelFiles': verified}
summary = {'Status': 'PASS', 'MetricRows': len(metrics), 'RankingRows': len(ranking), 'LOOCVRows': len(loocv), 'SevenEqualRanksRecomputed': True, 'HeldOutDatasetExcludedFromTraining': True, 'n4WithinEvaluableCoverage': True, 'FrozenSourceIntegrity': source_check}
(out / '00_audit/final_marker_independent_QC.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
