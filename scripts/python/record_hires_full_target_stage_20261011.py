# 登记本轮完整目标面板、种子修正、辅助回填和独立数学验收的版本与来源。
# 本机运行；只保存轻量清单，不打包官方软件、反编译底稿或私人研究矩阵。
from pathlib import Path
import hashlib
import json
import platform
import numpy
import pandas
import scipy

# Step1. 记录实际使用的研究脚本与未完成阶段
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
scripts = ['scripts/python/prepare_hires_full_target_benchmark_20261011.py', 'scripts/python/prepare_hires_seed_consistency_check_20261011.py', 'scripts/R/reconstruct_hires_entry_statistics_20261011.R', 'scripts/R/reconstruct_hires_worker_statistics_20261011.R', 'scripts/R/reconstruct_hires_refinement_20261011.R', 'scripts/R/reconstruct_hires_auxiliary_refill_20261011.R', 'scripts/R/check_hires_refinement_checkpoints_20261011.R', 'scripts/python/verify_hires_full_target_math_20261011.py', 'scripts/python/audit_hires_full_target_20261011.py', 'scripts/bash/run_hires_full_target_benchmark_20261011.sh', 'scripts/python/record_hires_full_target_stage_20261011.py']
files = [{'path': name, 'sha256': hashlib.sha256((project / name).read_bytes()).hexdigest()} for name in scripts]
benchmark = json.loads((audit / 'full_target_validation/full_target_validation_summary.json').read_text(encoding='utf-8'))
math_checks = json.loads((audit / 'full_target_math_checks/full_target_math_summary.json').read_text(encoding='utf-8'))
manifest = {'date': '2026-10-11', 'status': 'Complete Job24 target-gene numerical benchmark; full independent CIBERSORTx incomplete', 'platform': platform.platform(), 'python': platform.python_version(), 'numpy': numpy.__version__, 'pandas': pandas.__version__, 'scipy': scipy.__version__, 'r': 'R 4.4.1, nnls 1.6; sessionInfo.txt retained for each run', 'scripts': files, 'benchmark': benchmark, 'math_checks': math_checks, 'method_corrections': ['Explicit one-to-one gene-symbol mapping for dot/hyphen output aliases.', 'Canonical original-mixture row indices replace inherited Adjusted-export indices for gene-specific independent seeds.', 'Auxiliary refill all-one guard evaluated over all supplied target genes.'], 'missing_stages': ['Independent signature construction and cell fractions/S-mode', 'Independent Group estimation and adaptive reliability filtering', 'Official random stream and current-version detection/shrinkage parity', 'Remaining NA and constant/variable disagreements', 'Validated formal EcoTyper expression export'], 'source_boundary': 'Research numerical reconstruction only; no vendor program execution or authentication patching.', 'server': 'No SSH or Linux execution this round; Bash syntax checked only', 'formal_analysis': 'No benchmark matrix used for formal discovery/recovery.'}
(audit / 'full_target_stage_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
previous = json.loads((audit / 'implementation_stage_manifest.json').read_text(encoding='utf-8'))
previous['latest_stage_report'] = '完整目标基因验证与回填修正.md'
previous['latest_stage_manifest'] = 'full_target_stage_manifest.json'
previous['historical_manifest_note'] = 'Previous-stage results retained; current canonical seed convention is documented in the latest manifest.'
(audit / 'implementation_stage_manifest.json').write_text(json.dumps(previous, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print('完整目标阶段清单已登记：', len(files), '份脚本；完整独立流程尚未完成。')
