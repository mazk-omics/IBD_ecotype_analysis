# 保存本轮Group研究实现、来源证据、环境和验收摘要；不提交Git、不改历史运行。
# 在本机项目科学Python环境运行，读取现有结果，写第九轮manifest及历史入口提示。
from pathlib import Path
import ast
import hashlib
import json
import platform
import pandas as pd

# Step1. 逐项确认研究来源仍一致，计算文件哈希使后续审计可追溯
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
provenance = json.loads((audit / 'full_target_benchmark_inputs/input_provenance.json').read_text(encoding='utf-8'))
for item in provenance['sources']:
    assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() == item['sha256']
scripts = ['scripts/R/reconstruct_group_expression_20261011.R', 'scripts/R/export_group_math_checkpoints_20261011.R', 'scripts/R/diagnose_group_legacy_residual_20261011.R', 'scripts/R/reconstruct_hires_auxiliary_refill_20261011.R', 'scripts/python/audit_group_reconstruction_20261011.py', 'scripts/python/verify_group_math_20261011.py', 'scripts/python/record_group_reconstruction_stage_20261011.py', 'scripts/bash/run_fixed_fraction_expression_research_20261011.sh']
script_records = []
for relative in scripts:
    path = project / relative
    if path.suffix == '.py':
        ast.parse(path.read_text(encoding='utf-8'))
    script_records.append({'path': relative, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
evidence_names = ['vendor_reference/R_function_fragments/hires_f0bfeaddfa2e_runCIBERSORTxGEP.R', 'vendor_reference/R_function_fragments/hires_f0bfeaddfa2e_res.R', 'vendor_reference/R_function_fragments/hires_f0bfeaddfa2e_QCplots.R', 'decompiled_numeric/CIBERSORTxHiRes/numeric_001279fa.c', 'paper_methods/CIBERSORTx_2019_Supplement_pages.txt']
evidence = [{'path': name, 'sha256': hashlib.sha256((audit / name).read_bytes()).hexdigest()} for name in evidence_names]

# Step2. 官方CV条件验证与独立估计分开记录，不能用前者替代后者
validation = json.loads((audit / 'group_independent_validation/group_validation_summary.json').read_text(encoding='utf-8'))
math = json.loads((audit / 'group_math_checks/group_math_summary.json').read_text(encoding='utf-8'))
assert math['all_formula_checks_passed'] and validation['conditional_filter_mismatch'] == 0
manifest = {'date': '2026-10-11', 'stage': 9, 'status': 'Independent Group expression and CV filtering research prototype, fixed official fractions upstream; not web equivalence', 'report': 'Group独立估计与版本疑点.md', 'scripts': script_records, 'static_evidence': evidence, 'protected_inputs': provenance['sources'], 'run': pd.read_csv(audit / 'group_independent_run1/group_summary.tsv', sep='\t').to_dict('records'), 'validation': validation, 'math': math, 'python': platform.python_version(), 'R_session': 'group_independent_run1/sessionInfo.txt', 'paper_source': 'https://www.nature.com/articles/s41587-019-0114-2', 'version_choices': ['Full matched residual, not cached R short-vector recycling.', 'Paper zero guard, cached legacy one guard differs; Job24 branch did not trigger.', 'Independent per-gene L\'Ecuyer-CMRG seed 20261015+raw row index.', 'Bootstrap SD and analytical regression SE exported separately.'], 'missing': ['Exact Group sampling/random-stream and web-version equivalence.', 'HiRes 128 missed variable and 82 extra variable pairs persist before auxiliary.', 'Signature construction, fractions and S-mode independent implementation.', 'Linux/server integrated execution, full-transcriptome HiRes and formal EcoTyper output validation.'], 'execution_limits': ['No vendor executable run or authentication modification.', 'No SSH connection, commit or push.', 'Bash syntax checked only; R stages run on local Windows.', 'Legacy recycled-residual diagnostic kept separate and not used as calculation input.']}
(audit / 'group_reconstruction_stage_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')

# Step3. 历史manifest保存历史SHA，新增最新入口而不冒充重新验收旧运行
old_path = audit / 'full_target_stage_manifest.json'
old = json.loads(old_path.read_text(encoding='utf-8'))
old['latest_stage_report'] = 'Group独立估计与版本疑点.md'
old['latest_stage_manifest'] = 'group_reconstruction_stage_manifest.json'
old['historical_snapshot_note'] = 'Stage8 script hashes describe that run. Stage9 auxiliary script changes only provenance text; new Group results are separate.'
old_path.write_text(json.dumps(old, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print(json.dumps({'stage': 9, 'scripts': len(scripts), 'input_hashes_unchanged': True, 'math_passed': True, 'conditional_filter_mismatch': 0, 'independent_Group_filter_mismatch': validation['Group_independent_filter_NA_mismatch'], 'independent_Group_HiRes_NA_mismatch': validation['hires_NA_missed'] + validation['hires_false_NA']}, ensure_ascii=False, indent=2))
