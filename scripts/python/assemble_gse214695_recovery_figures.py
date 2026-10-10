# 合并本轮全部17张已审计图件；图中数值和独立源文件保持原样。
# 在本机cibersortx_review环境运行；写出阅读用PDF，目视记录在逐图查看后另存。
from pathlib import Path
import json
from pypdf import PdfReader, PdfWriter

# Step1. 先看总体和四个重点谱系，再看供者、生态型与其余状态
project = Path(__file__).resolve().parents[3]
root = project / 'IBD_ecotype_analysis/06_validation/GSE214695_scRNA_recovery_20261009'
figures = root / 'figures'
order = ['Fig1_recovery_overview', 'Fig2_Mast_marker_recovery', 'Fig2_Goblet_marker_recovery', 'Fig2_Enteroendocrine_marker_recovery', 'Fig2_Glial_marker_recovery', 'Fig3_donor_state_support', 'Fig4_lineage_marker_checks', 'Fig5_ecotype_projection', 'Fig6_colon_ileum_comparison']
other_types = ['B', 'CD4_T', 'CD8_T', 'Endothelial', 'Fibroblast', 'ILC', 'Monocyte_Macrophage', 'Pericyte']
order += [f'Fig2_{t}_marker_recovery' for t in other_types]
assert len(order) == 17 and set(order) == {p.stem for p in figures.glob('*.pdf')}
checks = json.loads((root / 'figure_qa/audit_exit_codes.json').read_text())
assert len(checks) == 17 and all(x['text_exit_code'] == 0 and x['collision_exit_code'] == 0 for x in checks)

# Step2. 合并并检查页数；完整单图及源数据仍分开保存
writer = PdfWriter()
for name in order:
    writer.append(figures / f'{name}.pdf', outline_item=name)
writer.add_metadata({'/Title': 'Independent GSE214695 recovery of frozen MSCCR Region3 CD+UC model', '/Subject': 'Exploratory recovery; 18 donors, 9752 cells, 12 lineages; official code unchanged'})
output = root / 'recovery图件汇总.pdf'
writer.write(output)
writer.close()
assert len(PdfReader(output).pages) == 17
(root / 'figure_qa/packet_check.json').write_text(json.dumps({'figures': 17, 'pages': 17, 'automated_qa_passed': True, 'manual_review_record': 'See final_review.json when available', 'page_order': order}, indent=2), encoding='utf-8')
print('汇总17页：', output)
