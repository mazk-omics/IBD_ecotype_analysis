# 只读检查已下载的官方应用层：区分 R 文本与 ELF 编译文件，判断源码可读范围。
# 本机 cibersortx_review 环境运行；不执行程序、不提取认证内容、不修改镜像，输出静态结构审计。
from pathlib import Path
import tarfile
import hashlib
import struct
import re
import csv
import json

# Step1. 复用昨日的小型应用层，逐层核对摘要，不重新下载整个容器
project = Path(__file__).resolve().parents[2]
source = project / 'CIBERSORTx_server_plan_20261010'
outdir = project / 'CIBERSORTx_server_plan_20261011'
outdir.mkdir(exist_ok=True)
files = []
sections = []
source_structure = []
markers = ['RInside', 'Rcpp', 'parseEval', 'Rscript', 'run_SVR', 'run_optimization', 'run_combatba', 'run_adjust_w_pseudo', 'run_make_HiRes_heatmap', 'CIBERSORTxGEP', 'CIBERSORTxHiRes']

for layer in sorted(source.glob('sha256_*.tar.gz')):  # 容器由多个层构成；只读取四个已保存的应用层
    archive_bytes = layer.read_bytes()
    expected = layer.name.removeprefix('sha256_').removesuffix('.tar.gz')
    assert hashlib.sha256(archive_bytes).hexdigest() == expected, layer.name
    with tarfile.open(layer, 'r:gz') as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            stream = archive.extractfile(member)
            data = stream.read()
            item = {'LayerSHA256': expected, 'Path': member.name, 'Bytes': len(data), 'FileSHA256': hashlib.sha256(data).hexdigest()}

            # Step2. ELF 头与节表仅说明编译结构；不反汇编认证逻辑，不输出二进制字符串
            if data[:4] == b'\x7fELF':
                assert data[4] == 2 and data[5] == 1, 'Only ELF64 little-endian is supported in this audit'
                header = struct.unpack('<16sHHIQQQIHHHHHH', data[:64])
                section_offset, section_size, section_count, string_index = header[6], header[11], header[12], header[13]
                assert section_size == 64 and section_count > 0
                table = [struct.unpack('<IIQQQQIIQQ', data[section_offset + i * 64:section_offset + (i + 1) * 64]) for i in range(section_count)]
                string_header = table[string_index]
                names = data[string_header[4]:string_header[4] + string_header[5]]
                section_names = []
                for entry in table:
                    end = names.find(b'\0', entry[0])
                    name = names[entry[0]:end].decode('ascii', errors='replace')
                    assert re.fullmatch(r'[A-Za-z0-9_.-]{0,100}', name), 'Unexpected section name'
                    section_names.append(name)
                    sections.append({'Path': member.name, 'Section': name, 'Bytes': entry[5], 'Type': entry[1]})
                dependencies = sorted(set(match.decode('ascii') for match in re.findall(rb'lib[A-Za-z0-9_.+-]{1,80}\.so(?:\.\d+)*', data)))
                item.update({'Kind': 'ELF64 compiled executable', 'Machine': header[2], 'HasSymbolTable': '.symtab' in section_names, 'HasDebugInfo': '.debug_info' in section_names, 'HasDynamicSymbols': '.dynsym' in section_names, 'Dependencies': dependencies, 'AlgorithmMarkerPresence': {word: word.encode() in data for word in markers}})
                symbol_names = []
                source_names = []
                for section_name, entry in zip(section_names, table):
                    if section_name != '.symtab':
                        continue
                    string_header = table[entry[6]]
                    symbol_strings = data[string_header[4]:string_header[4] + string_header[5]]
                    assert entry[9] == 24
                    for offset in range(entry[4], entry[4] + entry[5], 24):
                        symbol = struct.unpack('<IBBHQQ', data[offset:offset + 24])
                        end = symbol_strings.find(b'\0', symbol[0])
                        name = symbol_strings[symbol[0]:end].decode('ascii', errors='replace')
                        if symbol[1] & 15 == 2:
                            symbol_names.append(name)
                        if symbol[1] & 15 == 4 and re.fullmatch(r'[A-Za-z0-9_.-]+\.(cpp|cc|c)', name):
                            source_names.append(name)
                groups = ['svm', 'optimization', 'window', 'deconvol', 'expression', 'imput', 'normalize', 'quantile', 'fraction', 'hires', 'median']
                item.update({'FunctionSymbolCount': len(symbol_names), 'CompilationUnitBasenames': source_names, 'AlgorithmSymbolCounts': {word: sum(word.lower() in name.lower() for name in symbol_names) for word in groups}})

            # Step3. R 文件仅保存包名、函数名和算法标记位置；不保存源码行、字面量或认证字段
            elif member.name.lower().endswith('.r'):
                text = data.decode('utf-8-sig')
                definitions = re.findall(r'(?m)^\s*([A-Za-z][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function\s*\(', text)
                algorithm_functions = [name for name in definitions if re.search(r'(?i)(SVR|optimi|adjust|combat|normaliz|heatmap|CIBERSORT|imput|geps|highres)', name) and not re.search(r'(?i)(auth|login|token|license|credential)', name)]
                packages = sorted(set(re.findall(r'\b(?:library|require)\s*\(\s*[\"\x27]?([A-Za-z][A-Za-z0-9.]*)', text)))
                operations = ['svm', 'nnls', 'optim', 'optimize', 'lm', 'normalize.quantiles', 'ComBat', 'combatba', 'predict', 'coef', 'cor', 'quantile', 'median', 'kmeans', 'solve', 'mclapply']
                numerical_calls = {name: len(re.findall(r'(?<![A-Za-z0-9_.])' + re.escape(name) + r'\s*\(', text)) for name in operations}
                item.update({'Kind': 'Readable R source', 'Lines': len(text.splitlines()), 'AlgorithmFunctions': algorithm_functions, 'Packages': packages, 'HighResTextMentions': len(re.findall(r'(?i)hires|high.?res', text)), 'NumericalCallOccurrences': numerical_calls})
                hits = {word: [i for i, line in enumerate(text.splitlines(), 1) if word in line] for word in markers}
                source_structure.append({'Path': member.name, 'Lines': item['Lines'], 'AlgorithmFunctions': algorithm_functions, 'Packages': packages, 'NumericalCallOccurrences': numerical_calls, 'AlgorithmMarkerLines': hits})
            else:
                item['Kind'] = 'Other application file'
            files.append(item)

# Step4. 保存可重复核对的结构证据；程序、原始 R 文本和研究数据均未执行或修改
record = {'date': '2026-10-11', 'scope': 'four existing public COPY application layers only', 'runtime_executed': False, 'authentication_modified': False, 'raw_source_exported': False, 'files': files}
(outdir / 'application_structure_audit.json').write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
(outdir / 'readable_R_algorithm_structure.json').write_text(json.dumps(source_structure, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
with (outdir / 'ELF_section_inventory.tsv').open('w', encoding='utf-8', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=['Path', 'Section', 'Bytes', 'Type'], delimiter='\t')
    writer.writeheader()
    writer.writerows(sections)
for item in files:
    print(item['Path'], item['Kind'], item['Bytes'])
    if item['Kind'].startswith('ELF'):
        print('  symbol table:', item['HasSymbolTable'], 'debug info:', item['HasDebugInfo'], 'dependencies:', item['Dependencies'])
        print('  function symbols:', item['FunctionSymbolCount'], 'source units:', item['CompilationUnitBasenames'], 'algorithm symbol counts:', item['AlgorithmSymbolCounts'])
        print('  known algorithm markers:', item['AlgorithmMarkerPresence'])
    if item['Kind'] == 'Readable R source':
        print('  algorithm functions:', item['AlgorithmFunctions'], 'packages:', item['Packages'])
        print('  known numerical calls:', {name: count for name, count in item['NumericalCallOccurrences'].items() if count})
print('Output:', outdir)
