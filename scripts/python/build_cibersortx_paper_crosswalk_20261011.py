# 对照原论文补充材料与固定 HiRes 数值函数，保存八步流程和实现差异的证据索引。
# 本机 cibersortx_review 环境运行；只读数值代码，写审计表，不运行供应商程序或修改数据。
from pathlib import Path
import csv
import hashlib
import json
import struct
from elftools.elf.elffile import ELFFile

# Step1. 对照表区分论文规则、旧镜像证据与尚未验证的完整行为
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011'
worker = outdir / 'decompiled_numeric/CIBERSORTxHiRes_refinement_typed/runHiRes_typed.c'
main = outdir / 'decompiled_numeric/CIBERSORTxHiRes/numeric_001144fc.c'
rows = [
    (1, '12.3.1', 'PDF 34–35', '全局表达与显著性合并', 'runHiRes', '0010cc9b', 'Table::colMedians', 'bootstrap 系数中位数、非正频率与 BH 辅助表合并；全基因 BH 在目标基因子集限制之前'),
    (2, '12.3.2', 'PDF 35', 'bulk 表达排序并重排比例行', 'runHiRes', '0010cc9b', 'utility::order<double>', '逐基因排序，不使用 CD/UC/炎症标签定义两段'),
    (3, '12.3.3', 'PDF 35', '两段 NNLS 选择最小 RSS 分割', 'runHiRes', '0010cc9b', 'findBestRSSThres', '旧代码在 N-2W<10 时调整搜索边界至 W/2；边界与论文迭代描述需分别保留'),
    (4, '12.3.4', 'PDF 35–36', '差异显著性检验', 'runHiRes', '0010cc9b', 'Table::Table((Table *)local_15f8,0x32', 'Z 检验与 50 次随机参考合并；旧代码使用带符号差、排序和最近位置构造尾概率，不能直接替换为绝对差超过次数'),
    (5, '12.3.5–12.3.6', 'PDF 36–37', '重拟合与置信度收缩', 'runHiRes', '0010cc9b', 'pow<double,int>((undefined *)local_1878', '平方权重的收缩公式已核对；旧代码备用 Spearman 检验、+1 倍数、弱变化集合与最小平均比例类型覆盖需逐项保留'),
    (6, '12.3.7', 'PDF 37', '滑窗 NNLS 和边缘对齐平均', 'runHiRes', '0010cc9b', 'operator+((undefined *)local_1818,local_18d8', '前后补 W-1 行，平均为 N×C；补充正文写 row c，实际代码取最后有效窗口行'),
    (7, '12.3.7', 'PDF 37', '按两段参考作线性幅度校准', 'runHiRes', '0010cc9b', 'LinearRegression::LinearRegression', '旧代码加入 sd=0.001 扰动；完整随机流尚未与官方输出验证'),
    (8, '12.3.8', 'PDF 37', '样本逆排序及按类型导出', 'runHiRes', '0010cc9b', 'utility::order<int>', 'worker 逆排列可读；主函数有额外标签重排、常量替换、下限和可选 range>0.1 行筛选')
]
lines = worker.read_text(encoding='utf-8').splitlines()
table = []
for step, section, pages, topic, function, address, needle, finding in rows:
    matches = [i + 1 for i, text in enumerate(lines) if needle in text]
    assert matches, '证据标记未找到：' + needle
    table.append({'step': step, 'paper_section': section, 'paper_pdf_pages': pages, 'topic': topic, 'function': function, 'address': address, 'code_file': str(worker.relative_to(outdir)), 'code_lines': ','.join(map(str, matches)), 'finding': finding, 'official_execution_validated': False})
with (outdir / 'paper_code_step_crosswalk.tsv').open('w', encoding='utf-8', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(table[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(table)

# Step2. 只取已经定位的浮点常量；固定 ELF 摘要避免把版本证据混在一起
binary = outdir / 'vendor_reference/hires/src/CIBERSORTxHiRes'
binary_hash = hashlib.sha256(binary.read_bytes()).hexdigest()
assert binary_hash == '87511e7db27a973757e12d401b5184d306a9e85d24ae0748ab7199a640ae5f89'
addresses = [(0x22a298, 'global_alpha'), (0x22a2b0, 'minus_log10_alpha'), (0x22a2b8, 'weak_fold_change'), (0x22a2c0, 'anchor_jitter_sd'), (0x22a2c8, 'high_significance_and_regression_difference'), (0x22a2e0, 'suppress_other_q'), (0x22a2f0, 'optional_output_range')]
constants = []
with binary.open('rb') as stream:
    elf = ELFFile(stream)
    for address, role in addresses:
        section = next(item for item in elf.iter_sections() if item['sh_addr'] <= address - 0x100000 < item['sh_addr'] + item['sh_size'])
        offset = address - 0x100000 - section['sh_addr']
        value = struct.unpack('<d', section.data()[offset:offset + 8])[0]
        constants.append({'address': hex(address), 'analyst_role': role, 'value': value})

# Step3. 新补参版本另存摘要，不覆写早期自动导出或宣称获得作者源码
typed_files = sorted(worker.parent.glob('*.c'))
typed_hashes = [{'path': str(path.relative_to(outdir)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in typed_files]
sources = json.loads((outdir / 'paper_supplement_provenance.json').read_text(encoding='utf-8'))
evidence = {'binary_sha256': binary_hash, 'official_program_executed': False, 'original_source_recovered': False, 'new_typed_export_count': len(typed_hashes), 'typed_exports': typed_hashes, 'publisher_supplement_sha256': sources['sha256'], 'constants': constants, 'crosswalk_steps': len(table), 'main_code_sha256': hashlib.sha256(main.read_bytes()).hexdigest(), 'interpretation_limit': 'Typed static evidence from one historical ELF; paper and this version can differ. Numerical parity with current website is untested.'}
(outdir / 'paper_code_crosswalk_evidence.json').write_text(json.dumps(evidence, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print('Paper steps:', len(table), '; new typed exports:', len(typed_hashes), '; numerical constants:', len(constants))
print('Output:', outdir / 'paper_code_step_crosswalk.tsv')
