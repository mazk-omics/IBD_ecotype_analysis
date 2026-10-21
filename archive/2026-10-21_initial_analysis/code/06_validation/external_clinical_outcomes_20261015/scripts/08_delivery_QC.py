# 用Python核对最终统计与所有图件，并保留字体、碰撞和源代码预检日志。
# 不修改统计定义；每次图件改变后重新执行本脚本。
from pathlib import Path
import subprocess, sys, json, hashlib
import pandas as pd
import numpy as np

# Step1. 检查完整检验家庭、患者分母及统计数值。
out = Path(__file__).resolve().parents[1]
d = pd.read_csv(out / '04_cross_cohort_summary/all_associations.tsv', sep='\t')
assert len(d) == 360 and not d.duplicated(['Endpoint', 'Family', 'Feature']).any()
assert d.groupby(['Endpoint', 'Family']).size().tolist().count(7) == 6
assert d.groupby(['Endpoint', 'Family']).size().tolist().count(53) == 6
assert d.Status.eq('OK').all()
assert (d.CI_low < d.OR).all() and (d.OR < d.CI_high).all()
assert d[['OR', 'CI_low', 'CI_high', 'P_Firth', 'Q_Firth', 'Q_Wilcoxon']].apply(np.isfinite).all().all()
baseline_scores = pd.read_csv(out / '03_GSE73661_mucosal_healing/continuous_ecotype_scores.tsv', sep='\t', index_col=0)
all_samples = pd.read_csv(out / '03_GSE73661_mucosal_healing/all_samples_longitudinal_mapping.tsv', sep='\t')
assert len(baseline_scores) == 67 and set(baseline_scores.index) == set(all_samples.loc[all_samples.Timepoint.eq('W0'), 'GSM'])
for folder in ['01_GSE72819_treatment_outcome', '02_RISK_disease_progression', '03_GSE73661_mucosal_healing']:
    manifest = pd.read_csv(out / folder / 'cohort_endpoint_manifest.tsv', sep='\t')
    assert manifest.RecoveryComplete.all()
    for endpoint, m in manifest.groupby('Endpoint'):
        t = pd.read_csv(out / folder / (endpoint + '_analysis_patients.tsv'), sep='\t')
        assert t.PatientID.is_unique and t.SampleID.is_unique
        assert len(t) == int(d.loc[d.Endpoint.eq(endpoint), 'N'].iloc[0])

# Step2. 对所有最终主图与队列图进行实际PDF几何审计。
qa = out / 'scripts/figure_qa'
source = subprocess.run([sys.executable, str(qa / 'validate_figure.py'), str(out / 'scripts/06_plot_outcomes.py')], capture_output=True, text=True)
(out / 'logs/plot_source_preflight.log').write_text(source.stdout + source.stderr)
assert source.returncode == 0
figures = sorted(p for p in out.rglob('*.pdf') if p.parent.name == 'figures' and '.collision' not in p.name)
assert len(figures) == 7
results = []
for pdf in figures:
    stem = pdf.with_suffix('')
    assert pdf.with_suffix('.png').exists() and pdf.with_suffix('.svg').exists()
    text = subprocess.run([sys.executable, str(qa / 'audit_pdf_text.py'), str(pdf), '--min-pt', '5'], capture_output=True, text=True)
    collision = subprocess.run([sys.executable, str(qa / 'audit_figure_collisions.py'), str(pdf), '--json-out', str(stem) + '.collision.json', '--overlay-pdf', str(stem) + '.collision.pdf'], capture_output=True, text=True)
    (out / 'logs' / (pdf.stem + '_final_QA.log')).write_text(text.stdout + text.stderr + collision.stdout + collision.stderr)
    assert text.returncode == 0 and collision.returncode == 0
    results.append({'Figure': str(pdf.relative_to(out)), 'PDFText': 'PASS', 'Collision': 'PASS', 'Alignment': json.loads(Path(str(stem) + '.alignment.json').read_text()).get('verdict', 'See alignment JSON')})
pd.DataFrame(results).to_csv(out / 'audit/final_figure_QC.tsv', sep='\t', index=False)
(out / 'progress.json').write_text(json.dumps({'Date': '2026-10-15', 'Cohorts': 3, 'Endpoints': 6, 'PrimaryModels': 360, 'RiskAgeSexModels': 7, 'FinalFigures': 7, 'AllPrimaryFits': 'OK', 'Phase1_InputAudit': 'Completed with unresolved exact GSE93624 transformation', 'Phase2_GSE72819': 'Complete', 'Phase3_RISK': 'Conditional processed-scale association complete; provenance limitation remains', 'Phase4_GSE73661': 'Complete for 4 distinct endpoints', 'Phase5_Integration': 'Complete', 'Phase6_Delivery': 'Numerical and rendered QA complete'}, indent=2))
print('DELIVERY_QC_PASS', len(d), len(figures), flush=True)
