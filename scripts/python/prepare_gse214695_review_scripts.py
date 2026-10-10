# 为第二轮 recovery 复制已检查的审阅和图件结构，明确替换队列、样本数和谱系覆盖。
# 在本机运行；只创建本轮新脚本，不修改 IBDverse 脚本或官方 EcoTyper。
from pathlib import Path

# Step1. 继承同一统计定义，替换实际输入与分组信息
project = Path(__file__).resolve().parents[3]
scripts = project / 'IBD_ecotype_analysis/scripts/python'
text = (scripts / 'plot_ibdverse_recovery_nature.py').read_text(encoding='utf-8')
text = text.replace('IBDverse', 'GSE214695').replace("'Fibroblast', 'Goblet'", "'Fibroblast', 'Glial', 'Goblet'")
text = text.replace("['Healthy_uninflamed', 'CD_uninflamed', 'CD_inflamed']", "['HC_nonIBD', 'CD_inflamed', 'UC_inflamed']")
text = text.replace('set_xlim(0, 31)', 'set_xlim(0, 19)').replace('set_xticks([0, 10, 20, 30])', 'set_xticks([0, 6, 12, 18])')
text = text.replace('Independent terminal ileum recovery', 'Independent colon/rectum recovery')
text = text.replace("'30 donors; 4,818 cells; 11 lineages; Glial unavailable. Dashed line: nominal Z = 1.96.'", "f'{cells.Donor.nunique()} donors; {len(cells):,} cells; {cells.CellType.nunique()} lineages. Dashed line: nominal Z = 1.96.'")
text = text.replace("['Healthy / non-inflamed', 'CD / non-inflamed', 'CD / inflamed']", "['HC / non-tumor', 'CD / inflamed', 'UC / inflamed']")
text = text.replace('i * 10', 'i * 6').replace('(i + 1) * 10', '(i + 1) * 6').replace('(i + 0.5) * 10', '(i + 0.5) * 6')
text = text.replace('At most 20 cells', 'At most 100 cells').replace('Ecotype projection | missing Glial lineage', 'Ecotype projection | all 12 model lineages')
text = text.replace("'COL1A1', 'DCN', 'PECAM1', 'RGS5']", "'COL1A1', 'DCN', 'PECAM1', 'RGS5', 'SOX10', 'S100B', 'MPZ']")
(scripts / 'plot_gse214695_recovery_nature.py').write_text(text, encoding='utf-8')

# Step2. 审阅口径不变；全 53 状态保留，12 类型和 18 样本不套用上一队列分母
text = (scripts / 'review_ibdverse_recovery_results.py').read_text(encoding='utf-8')
text = text.replace('IBDverse', 'GSE214695').replace('47状态', '53状态').replace('11类型', '12类型').replace('30活检', '18样本')
(scripts / 'review_gse214695_recovery_results.py').write_text(text, encoding='utf-8')
text = (scripts / 'audit_ibdverse_recovery_figures.py').read_text(encoding='utf-8').replace('IBDverse', 'GSE214695').replace('plot_ibdverse', 'plot_gse214695')
(scripts / 'audit_gse214695_recovery_figures.py').write_text(text, encoding='utf-8')
print('Prepared plotting, review and QA scripts')
