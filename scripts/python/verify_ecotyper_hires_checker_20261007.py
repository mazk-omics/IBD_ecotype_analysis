# 用小型人工数据检查 HiRes 接口校验能否识别样本错位和误用 fraction；不涉及真实队列。
# 本地 cibersortx_review 环境运行；测试文件与结果放在 construction_chain_20261007/QA，勿上传。
from pathlib import Path
import json
import subprocess
import sys
import pandas as pd

# Step1. 建立四样本、两类型的人工输入，仅测试数据契约
root = Path(__file__).resolve().parents[2]
out = root / 'CIBERSORTx_signature_audit/construction_chain_20261007/QA'
good = out / 'valid_expression'
bad = out / 'wrong_sample_order'
fraction = out / 'fraction_table'
for directory in [good, bad, fraction]:
    directory.mkdir(parents=True, exist_ok=True)
samples = ['sample1', 'sample2', 'sample3', 'sample4']
annotation = pd.DataFrame({'ID': samples, 'Disease': ['CD', 'CD', 'UC', 'UC'], 'Region3': ['Ileum', 'Ileum', 'Colon', 'Rectum']})
annotation.to_csv(out / 'annotation.tsv', sep='\t', index=False)
for cell_type in ['Goblet', 'Paneth']:
    matrix = pd.DataFrame([[10, 20, 5, 30], [1, 4, 2, 8]], index=['gene1', 'gene2'], columns=samples)
    matrix.index.name = 'GeneSymbol'
    matrix.to_csv(good / f'{cell_type}.txt', sep='\t')
    matrix.loc[:, samples[::-1]].to_csv(bad / f'{cell_type}.txt', sep='\t')
    pd.DataFrame({'Mixture': samples, 'Goblet': [.5] * 4, 'Paneth': [.5] * 4}).to_csv(fraction / f'{cell_type}.txt', sep='\t', index=False)

# Step2. 正常输入须通过，错位和 fraction 须被拒绝
checks = []
for name, directory, expected_pass in [('valid', good, True), ('wrong_order', bad, False), ('fractions', fraction, False)]:
    command = [sys.executable, str(root / 'scripts/python/check_cibersortx_hires_for_ecotyper.py'), '--expression-dir', str(directory), '--annotation', str(out / 'annotation.tsv'), '--outdir', str(out / f'{name}_check'), '--expected-samples', '4', '--expected-types', 'Goblet,Paneth']
    result = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', errors='replace')
    assert (result.returncode == 0) == expected_pass
    (out / f'{name}.log').write_text(result.stdout + result.stderr, encoding='utf-8')
    checks.append({'case': name, 'exit_code': result.returncode, 'expected_behavior_verified': True})
(out / 'checker_test_summary.json').write_text(json.dumps({'artificial_test_only': True, 'checks': checks, 'actual_MSCCR_HiRes_tested': False}, indent=2), encoding='utf-8')
print(checks)
