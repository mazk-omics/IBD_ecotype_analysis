# 汇总完整层扫描、R 语法审计和 Ghidra 数值函数恢复，检查导出文件与原始来源对应。
# 本机 cibersortx_review 环境运行；只核对文件和摘要，不执行供应商算法，不将私有参考文件加入 Git。
from pathlib import Path
from collections import Counter
import hashlib
import json
import csv
import re
import importlib.metadata
import tarfile

# Step1. 对齐三类证据，避免将第三方库头文件当成项目原始源码
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011'
layers = json.loads((outdir / 'full_image_layers_checked.json').read_text())
images = json.loads((outdir / 'full_image_scan_summary.json').read_text())
applications = json.loads((outdir / 'full_image_application_versions.json').read_text())
R_functions = json.loads((outdir / 'R_function_expression_audit.json').read_text())
R_status = json.loads((outdir / 'R_parse_status.json').read_text())
assert len(images) == 3 and len(layers) == 63
assert all(row['status'] == 'parsed' for row in R_status)
project_source_names = ['CIBERSORTxHiRes', 'CIBERSORTxFractions', 'GETable', 'SVRWrapper', 'CellStatTable', 'LinearRegression', 'Stopwatch']
source_matches = []
for image in images:
    with (outdir / (image['Image'] + '_all_layer_files.tsv')).open(encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            filename = Path(row['Path']).name
            if filename.lower().endswith(('.cpp', '.cc', '.cxx', '.c', '.h', '.hpp')) and any(name.lower() in filename.lower() for name in project_source_names):
                source_matches.append({'Image': image['Image'], **row})
(outdir / 'project_source_filename_search.json').write_text(json.dumps(source_matches, indent=2) + '\n', encoding='utf-8')

# Step2. 每个导出的 R 片段与原始行号核对；换行统一不改变代码内容或方法
cache = project.parent / '.tools' / 'cibersortx_image_layers'
reference = outdir / 'vendor_reference'
fragments = []
for row in R_functions:
    if row['status'] != 'original_function_fragment_exported':
        continue
    version = next(item for item in applications if item['Image'] == row['image'] and item['Path'] == row['path'] and item['FileSHA256'] == row['source_sha256'])
    with tarfile.open(cache / ('sha256_' + version['LayerSHA256'] + '.tar.gz'), 'r:gz') as archive:
        source_bytes = archive.extractfile(row['path']).read()
    assert hashlib.sha256(source_bytes).hexdigest() == row['source_sha256']
    source_lines = source_bytes.decode('utf-8-sig').replace('\r\n', '\n').split('\n')
    expected = '\n'.join(source_lines[row['start_line'] - 1:row['end_line']])
    path = reference / 'R_function_fragments' / row['file']
    observed = '\n'.join(path.read_text(encoding='utf-8').splitlines()[1:])
    assert observed.rstrip('\n') == expected.rstrip('\n'), row['file']
    fragments.append({'File': str(path.relative_to(outdir)).replace('\\', '/'), 'Function': row['function_name'], 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest()})

# Step3. 检查反编译伪代码、认证筛选和占位类型；伪代码导出成功不等于算法验证通过
native = []
keys = []
typed_exports = []
markers = re.compile(r'(?i)auth|token|password|credential|username|login|license|https?://|http::|[a-z0-9._%+-]+@[a-z0-9.-]+\.[a-z]{2,}')
for directory in sorted((outdir / 'decompiled_numeric').iterdir()):
    if not directory.is_dir():
        continue
    if directory.name == 'CIBERSORTxHiRes_typed':
        # 补参导出是独立分析假设，不混入原始自动反编译的 687 条记录
        for path in sorted(directory.glob('*.c')):
            text = path.read_text(encoding='utf-8')
            assert not markers.search(text)
            assert text.startswith(('/* Static pseudocode with analyst-supplied parameter types;', '/* Analyst ABI interpretation of a numerical matrix helper;'))
            typed_exports.append({'File': str(path.relative_to(outdir)).replace('\\', '/'), 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest(), 'analyst_supplied_types': True})
        assert len(typed_exports) == 3
        continue
    with (directory / 'numerical_functions.tsv').open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    warnings = 0
    untyped = 0
    outputs = []
    for row in rows:
        if row['Status'] != 'pseudocode_exported':
            continue
        path = directory / row['PseudocodeFile']
        text = path.read_text(encoding='utf-8')
        assert not markers.search(text.replace('rmTokensAndIG', 'rmGeneLabelsAndIG')), 'Review withheld content before exporting any further code'
        assert text.startswith('/* Ghidra pseudocode, not original source; not compiled or validated. */')
        warnings += 'WARNING:' in text
        untyped += bool(re.search(r'\bundefined(?:[1248])?\b', text))
        outputs.append({'File': str(path.relative_to(outdir)).replace('\\', '/'), 'Name': row['Name'], 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest()})
        if row['Name'] in ['runHiRes', 'CIBERSORTxHiRes', 'CoreAlg', 'CIBERSORTxFractions', 'runrmbatchSmode', 'runrmbatchBmode', 'adjustSigMatrix', 'adjustMixtures', 'doPerm', 'lh_nnls_']:
            keys.append({'Program': directory.name, **row, 'RelativeFile': str(path.relative_to(outdir)).replace('\\', '/')})
    native.append({'Program': directory.name, 'SelectedFunctions': len(rows), 'Statuses': dict(Counter(row['Status'] for row in rows)), 'FilesWithDecompilerWarnings': warnings, 'FilesWithUndefinedTypes': untyped, 'Files': outputs})

# Step4. 保存覆盖范围和运行环境；不报告不存在依据的“源码恢复百分比”
R_versions = [row for row in applications if row['Kind'] == 'R']
ELF_versions = [row for row in applications if row['Kind'] == 'ELF']
summary = {'date': '2026-10-11', 'scope': 'all manifest layers of three pinned official latest images; algorithm-only static decompilation', 'server_connected': False, 'supplier_program_executed': False, 'authentication_modified': False, 'layers_scanned_with_duplicates': len(layers), 'unique_layer_digests': len({row['LayerSHA256'] for row in layers}), 'unique_compressed_bytes': sum({row['LayerSHA256']: row['CompressedBytes'] for row in layers}.values()), 'archive_members_with_duplicates': sum(row['ArchiveMembers'] for row in images), 'project_source_filename_matches': source_matches, 'R_source_versions': len(R_versions), 'R_unique_file_hashes': len({row['FileSHA256'] for row in R_versions}), 'R_function_definitions': len(R_functions), 'R_function_fragment_exports': len(fragments), 'R_functions_withheld': sum(row['status'] != 'original_function_fragment_exported' for row in R_functions), 'ELF_file_versions': len(ELF_versions), 'ELF_unique_file_hashes': len({row['FileSHA256'] for row in ELF_versions}), 'R_fragment_files': fragments, 'native_programs': native, 'numerical_equivalence_tested': False, 'pseudocode_compiled': False, 'python_packages': {name: importlib.metadata.version(name) for name in ['pyelftools', 'capstone', 'cxxfilt']}}
summary['screening_identifier_exception'] = {'identifier': 'rmTokensAndIG', 'reason': 'gene-label filtering method, not a credential literal; exact identifier exception only', 'diagnostic': 'signature_screening_kind.tsv', 'supplier_source_rewritten': False}
summary['analyst_typed_exports'] = typed_exports
(outdir / 'source_recovery_validation.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
(outdir / 'key_algorithm_functions.json').write_text(json.dumps(keys, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print('Scanned layers:', len(layers), 'unique:', summary['unique_layer_digests'], 'original C++ basename matches:', len(source_matches))
print('R source versions:', len(R_versions), 'unique:', summary['R_unique_file_hashes'], 'verified function fragments:', len(fragments))
for program in native:
    print(program['Program'], program['Statuses'], 'warnings:', program['FilesWithDecompilerWarnings'], 'undefined types:', program['FilesWithUndefinedTypes'])
print('Saved:', outdir / 'source_recovery_validation.json')
