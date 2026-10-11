# 记录本轮独立实现的文件版本、运行环境和验收状态，不把局部统计阶段标成完整软件。
# 在本机科学Python环境运行；写出轻量清单，供后续服务器运行与代码审阅定位。
from pathlib import Path
import hashlib
import json
import platform
import numpy
import pandas
import scipy

# Step1. 清单只包含本轮独立脚本和统计证据，不打包官方二进制或研究矩阵
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
script_paths = ['scripts/python/audit_hires_input_space_20261011.py', 'scripts/R/reconstruct_hires_entry_statistics_20261011.R', 'scripts/R/reconstruct_hires_worker_statistics_20261011.R', 'scripts/R/benchmark_hires_numeric_core_20261011.R', 'scripts/python/check_hires_statistics_and_input_fix_20261011.py', 'scripts/python/inventory_public_cibersortx_repositories_20261011.py', 'scripts/python/inspect_ecotyper_public_interfaces_20261011.py', 'scripts/bash/run_hires_core_smoke_20261011.sh', 'scripts/bash/run_hires_statistics_smoke_20261011.sh', 'scripts/R/reconstruct_hires_refinement_20261011.R', 'scripts/R/check_hires_refinement_checkpoints_20261011.R', 'scripts/python/audit_hires_refinement_20261011.py', 'scripts/python/trace_hires_remaining_discrepancies_20261011.py']
files = [{'path': path, 'sha256': hashlib.sha256((project / path).read_bytes()).hexdigest()} for path in script_paths]
validation = json.loads((audit / 'statistics_validation/validation_summary.json').read_text(encoding='utf-8'))
manifest = {'date': '2026-10-11', 'status': 'Independent numerical stages validated; full workflow incomplete', 'platform': platform.platform(), 'python': platform.python_version(), 'numpy': numpy.__version__, 'pandas': pandas.__version__, 'scipy': scipy.__version__, 'r_environment': 'R 4.4.1, nnls 1.6; exact loaded packages in each R result sessionInfo.txt', 'files': files, 'validation': validation, 'refinement_validation': json.loads((audit / 'refinement_validation/refinement_validation_summary.json').read_text(encoding='utf-8')), 'remaining_discrepancies': json.loads((audit / 'refinement_validation/remaining_discrepancies_summary.json').read_text(encoding='utf-8')), 'outer_refill_context_diagnostic': json.loads((audit / 'refinement_validation/outer_refill_context_summary.json').read_text(encoding='utf-8')), 'server_status': 'No SSH or Linux execution this round', 'implemented_stages': ['fixed official fractions plus explicit original expression input', 'full gene-space entry bootstrap, analytic t and BH', 'panel global detection, two-segment NNLS, Z and randomized differential statistics', 'uncalibrated sliding core and independent numerical checks', 'residual Spearman AS89/rank-t and active-type refit', 'direction-specific squared shrinkage and weak-type residual NNLS', 'sliding smooth and OLS calibration, sample order and lower floor', 'same-seed repeatability and limited refinement seed sensitivity'], 'missing_stages': ['independent signature construction', 'independent fractions, permutation diagnostics and S-mode', 'full-context auxiliary refill and final NA/reliability rules', 'full real-data coverage of rare-type override and official random stream', 'validated EcoTyper expression export and reusable full CLI'], 'official_parity': 'Not established; old numerical evidence and paper methods remain distinct from current website version'}
(audit / 'implementation_stage_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print('阶段清单已保存：', len(files), '份脚本，完整版尚未完成。')
