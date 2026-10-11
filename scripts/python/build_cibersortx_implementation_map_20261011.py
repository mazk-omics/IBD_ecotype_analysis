# 为已取得的原始 R 与数值伪代码建立阅读目录，标出仍未展开的直接依赖。
# 本机 cibersortx_review 环境运行；只读审计结果，写索引和摘要，不能据此计算完整算法恢复比例。
from pathlib import Path
from collections import Counter
import csv
import json
import hashlib

# Step1. 函数目录以程序版本和地址区分；同名函数、模板和不同二进制不能合并计数
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011'
functions = []
remaining = []
for directory in sorted((outdir / 'decompiled_numeric').iterdir()):
    if not directory.is_dir():
        continue
    if directory.name == 'CIBERSORTxHiRes_typed':
        continue # 分析器补参假设单列在 HiRes 接口与输出审计，不混入自动导出函数计数
    with (directory / 'numerical_functions.tsv').open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    for row in rows:
        relative = 'decompiled_numeric/' + directory.name + '/' + row['PseudocodeFile'] if row['PseudocodeFile'] else ''
        functions.append({'Program': directory.name, **row, 'RelativeFile': relative, 'SHA256': hashlib.sha256((outdir / relative).read_bytes()).hexdigest() if relative else ''})
    with (directory / 'numerical_dependency_edges.tsv').open(encoding='utf-8', newline='') as stream:
        dependencies = list(csv.DictReader(stream, delimiter='\t'))
    for row in dependencies:
        if row['Selected'] == 'true' or any(name in row['Callee'] for name in ['std::', 'boost::', 'Eigen::', '__gnu', '<EXTERNAL>']) or row['Callee'].startswith(('_', 'FUN_')):
            continue
        remaining.append({'Program': directory.name, **row})

# Step2. 保存完整函数索引和直接依赖；“Selected”只表示已选择导出，不表示动态分支已覆盖
with (outdir / 'numerical_source_index.tsv').open('w', encoding='utf-8', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(functions[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(functions)
remaining = list({(row['Program'], row['Caller'], row['Callee'], row['Address']): row for row in remaining}.values())
with (outdir / 'unexpanded_application_dependency_edges.tsv').open('w', encoding='utf-8', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=['Program', 'Caller', 'Callee', 'Address', 'Selected'], delimiter='\t')
    writer.writeheader()
    writer.writerows(remaining)

# Step3. 人工阅读入口集中到科学算法；原始函数、编译器生成代码与未确认语义分别说明
key_names = ['CIBERSORTxHiRes', 'runHiRes', 'doBootstrap', 'utility::sample<int>', 'utility::sequence', 'utility::order<double>', 'nnls', 'subsetNNLS', 'nnlsFTWrapper', 'lh_nnls_', 'findBestRSSThres', 'LinearRegression::calcSlope', 'LinearRegression::calcIntercept', 'LinearRegression::predict', 'CoreAlg', 'doPerm', 'runrmbatchSmode', 'adjustSigMatrix', 'adjustMixtures', 'buildSignatureMatrix', 'utility::rmTokensAndIG', 'GETable::rmTokensAndIG', 'GETable::performOptimization']
lines = ['# 数值源码阅读目录', '', '这里只索引已经取得的静态材料。伪代码不是作者原始源码，尚未编译或与网页结果核对；函数数量包含模板、第三方求解器和版本重复。', '', '最新进度判断与补参证据见 [HiRes 接口与输出审计](HiRes接口与输出审计.md)。此前的科学步骤见 [算法步骤与可实现性](算法步骤与可实现性.md)。全部自动导出记录见 [函数索引 TSV](numerical_source_index.tsv)。', '', '| 程序版本 | 函数 | 状态 | 伪代码入口 |', '|---|---|---|']
for row in functions:
    if row['Name'] not in key_names:
        continue
    link = '[' + row['Address'] + '](' + row['RelativeFile'] + ')' if row['RelativeFile'] else '暂未导出'
    lines.append('| ' + row['Program'] + ' | `' + row['Name'] + '` | ' + row['Status'] + ' | ' + link + ' |')
lines += ['', '原始 R 的函数及行号来源见 [R 语法审计 JSON](R_function_expression_audit.json)，通过逐行核对的导出片段位于 `vendor_reference/R_function_fragments/`。', '', '尚未展开的应用直接依赖见 [依赖边 TSV](unexpanded_application_dependency_edges.tsv)。其中包含计时、路径、通用 R 命令封装和读写辅助；不等于全部是未知的数值算法。间接调用、全局状态、异常分支和数据格式仍需另行核对。', '', '本轮未运行供应商程序，也未修改认证机制。']
(outdir / '数值源码阅读目录.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
counts = dict(Counter(row['Program'] for row in functions if row['Status'] == 'pseudocode_exported'))
summary = {'program_exports': counts, 'total_export_records': sum(counts.values()), 'unexpanded_application_direct_edges': len(remaining), 'coverage_percentage_claimed': False, 'binary_executed': False, 'dynamic_equivalence_tested': False}
(outdir / 'implementation_map_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
print('Export records:', counts)
print('Unexpanded application direct edges:', len(remaining), '(including paths, timing, and R command wrappers)')
