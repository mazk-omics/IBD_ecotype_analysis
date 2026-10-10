# 将已审计的单张 recovery PDF 按汇报阅读顺序合并，原图及源数据继续分别保存。
# 在本地 .tools/cibersortx_review 环境运行；读取本轮 figures，输出 15 页汇总 PDF 和目视审计记录。
# 不重绘、不改变数值、颜色和文字；合并后检查页数，并重新核对服务器源文件哈希记录。
from pathlib import Path
import json
from pypdf import PdfReader, PdfWriter

# Step1. 明确汇报顺序，先看总览和预先关注的谱系，再看其余类型
project = Path(__file__).resolve().parents[3]
root = project / 'IBD_ecotype_analysis/06_validation/IBDverse_scRNA_recovery_20261009'
figures = root / 'figures'
first = ['Fig1_recovery_overview', 'Fig2_Mast_marker_recovery', 'Fig2_Goblet_marker_recovery', 'Fig2_Enteroendocrine_marker_recovery', 'Fig3_donor_state_support', 'Fig4_lineage_marker_checks', 'Fig5_ecotype_projection']
other_types = ['B', 'CD4_T', 'CD8_T', 'Endothelial', 'Fibroblast', 'ILC', 'Monocyte_Macrophage', 'Pericyte']
order = first + [f'Fig2_{t}_marker_recovery' for t in other_types]
assert len(order) == 15 and set(order) == {p.stem for p in figures.glob('*.pdf')}, '汇总必须包含全部最终图件，不能漏掉弱状态图'

# Step2. 按顺序合并可编辑 PDF，源文件保持不变
writer = PdfWriter()
for name in order:
    writer.append(figures / f'{name}.pdf', outline_item=name)
writer.add_metadata({'/Title': 'IBDverse scRNA recovery of frozen MSCCR Region3 CD+UC discovery', '/Subject': 'Exploratory recovery: 30 donors, 4818 cells, 11 lineages; official code unchanged'})
output = root / 'recovery图件汇总.pdf'
writer.write(output)
writer.close()
assert len(PdfReader(output).pages) == 15
print('Saved:', output, '; pages:', len(PdfReader(output).pages))

# Step3. 记录已完成的最终图件检查；这些检查不代表统计或谱系验证通过
qa = root / 'figure_qa'
checks = json.loads((qa / 'audit_exit_codes.json').read_text(encoding='utf-8'))
assert len(checks) == 15 and all(x['text_exit_code'] == 0 and x['collision_exit_code'] == 0 for x in checks)
before = (root / 'raw_server/audit/source_before.sha256').read_text(encoding='utf-8')
after = (root / 'raw_server/audit/source_after.sha256').read_text(encoding='utf-8')
assert before == after, '官方源文件运行前后不一致，需停止解释并核对'
review = {'figure_count': 15, 'packet_pages': 15, 'manual_review': '已逐张查看全部 15 张 PDF 的渲染预览；47 状态、NA、未分配分母和字体均保留，修正图 1 与图 3 的文字间距后重新审计', 'pdf_text_and_collision': '15/15 passed', 'export_warning': 'PNG 非 TIFF 提示保留；PDF/SVG 为交付主格式', 'official_source_hashes_unchanged': True, 'page_order': order}
(qa / 'final_review.json').write_text(json.dumps(review, ensure_ascii=False, indent=2), encoding='utf-8')
print('QA and source hashes verified')
