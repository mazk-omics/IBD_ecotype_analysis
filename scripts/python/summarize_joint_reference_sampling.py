# 展示已用于 GSE282122 联合 signature 的实际抽样组成，不重新抽样、不修改原矩阵。
# 在本机项目环境运行；读取 1,000 个细胞的冻结清单，保存部位、疾病、供者和作者亚型汇总。
# 这里的等额抽样用于建参考，不代表肠黏膜真实细胞丰度；验证供者单独列出。
import json
from pathlib import Path

import pandas as pd

# Step1. 读取实际抽样清单，按原参考中的类型顺序核对
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
source = project / 'IBDverse_reference/three_reference_pilots_20261006/run_20261006_183211/gse_joint'
out = project / 'CIBERSORTx_signature_audit/joint_reference_sampling_display_20261006'
out.mkdir(exist_ok=True)
summary = json.loads((source / 'summary.json').read_text(encoding='utf-8'))
types = summary['web_reference_types']
cells = pd.read_csv(source / 'selected_cell_manifest.tsv', sep='\t')
split = pd.read_csv(source / 'donor_split.tsv', sep='\t')
labels = pd.read_csv(source / 'selected_author_label_composition.tsv', sep='\t')
assert len(cells) == 1000 and cells.cell_id.is_unique and cells.split.eq('training').all()
assert cells.groupby('pilot_type').size().reindex(types).eq(50).all()
assert not set(cells.Patient) & set(split.loc[split.split.eq('holdout'), 'Patient'])
print('实际参考:', len(cells), '细胞;', cells.Patient.nunique(), '供者;', len(types), '类型')

# Step2. 保留五个原始部位，再另列便于阅读的三部位汇总
# 合并 Ascending/Descending/Sigmoid 只用于展示，没有改变抽样或任何表达值。
sites = pd.crosstab(cells.pilot_type, cells.Site).reindex(types)
diseases = pd.crosstab(cells.pilot_type, cells.Disease).reindex(types)
donors = cells.groupby('pilot_type').Patient.nunique().reindex(types)
display = pd.DataFrame(index=pd.Index(types, name='cell_type'))
display['cells'] = cells.groupby('pilot_type').size().reindex(types)
display['donors'] = donors
display['Terminal_Ileum'] = sites.Terminal_Ileum
display['Colon_combined'] = sites[['Ascending_Colon', 'Descending_Colon', 'Sigmoid']].sum(axis=1)
display['Rectum'] = sites.Rectum
display[['CD', 'UC', 'Healthy']] = diseases[['CD', 'UC', 'Healthy']]
display.to_csv(out / 'type_region_disease_donor_counts.tsv', sep='\t')
sites.to_csv(out / 'type_five_original_sites.tsv', sep='\t')
labels.to_csv(out / 'type_author_subtype_counts.tsv', sep='\t', index=False)
print(display.to_string())

# Step3. 保存逐供者贡献和 QC 分布，避免把供者平衡当作亚型平衡
donor_counts = cells.groupby(['Patient', 'Disease']).size().rename('cells').reset_index()
donor_counts.to_csv(out / 'donor_contributions.tsv', sep='\t', index=False)
donor_type = pd.crosstab(cells.Patient, cells.pilot_type).reindex(columns=types)
donor_type.to_csv(out / 'donor_by_type_counts.tsv', sep='\t')
qc = cells.groupby('pilot_type').agg(mt_median=('pct_counts_mt', 'median'), mt_max=('pct_counts_mt', 'max'), counts_median=('total_counts', 'median'), detected_genes_median=('n_genes_by_counts', 'median')).reindex(types)
qc.to_csv(out / 'selected_QC_by_type.tsv', sep='\t')
assert int(donor_type.to_numpy().max()) == 4 and cells.Patient.nunique() == 33
roles = pd.read_csv(project / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006/server_outputs/validation_joint20_20261006/donor_roles.tsv', sep='\t')
roles.to_csv(out / 'validation_donor_roles.tsv', sep='\t', index=False)

# Step4. 写阅读入口和作者亚型详情，数字都来自冻结的细胞清单
columns = display.reset_index().columns.tolist()
table = '| ' + ' | '.join(columns) + ' |\n| ' + ' | '.join(['---'] * len(columns)) + ' |\n'
for row in display.reset_index().itertuples(index=False, name=None):
    table += '| ' + ' | '.join(map(str, row)) + ' |\n'
subtypes = '| 类型 | 作者原始注释 | 细胞数 |\n|---|---|---:|\n'
for row in labels.itertuples(index=False, name=None):
    subtypes += '| ' + ' | '.join(map(str, row)) + ' |\n'
report = '# GSE282122 联合参考实际抽样\n\n2026-10-06。统计的是用于建联合 signature 的训练参考，不是留出模拟样本，也没有重新抽样。\n\n'
report += '20 类 × 每类 50 个 = 1,000 个细胞；20,612 基因；33 位训练供者，8 位供者留出。原始表达为非 log 整数 counts。网页最终 signature 为 5,816 × 20，不等同输入参考维度。\n\n'
report += '抽样种子 20261006。每类按随机供者顺序轮流分配名额，再在各供者内无放回抽取细胞；设定每供者每类最多 5 个，实际最多 4 个。未强制疾病、部位或作者亚型等额。每位供者实际总贡献 23–36 个细胞。\n\n'
report += '输入沿用已审计的 CONFIDENT 身份与 KEEP 类型；本轮额外纳入作者确认的 Glial/Tuft。IBD 使用治疗前细胞，Healthy 按健康来源纳入；排除预测 doublet 和无有效 counts/n_genes 的行。本轮没有新增线粒体比例筛选，不把原始作者 QC 自动等同为本项目充分 QC。\n\n'
report += '## 各类型组成\n\nColon_combined 仅在展示时合并升结肠、降结肠、乙状结肠；原始五部位表单独保存。\n\n' + table
report += '\n## 总体分布\n\n部位：Terminal_Ileum 212；Ascending_Colon 118；Descending_Colon 263；Sigmoid 117；Rectum 290。疾病：CD 415、UC 520、Healthy 65。炎症注释：Inflamed 639、Non_Inflamed 296、Healthy 65。训练供者：CD 13、UC 18、Healthy 2。\n\n'
report += '留出供者：CD9、CD10、CD12、H2、UC4、UC11、UC15、UC17。其中当前模拟测试使用 CD9、CD10、H2、UC4、UC11、UC15；CD12、UC17 继续保留。\n\n'
report += '## 解读时要注意\n\n1. 每类 50 个是参考平衡设计，不是组织真实占比。部位/亚型没有被强制平衡。\n2. Paneth 包括作者 Ileal paneth 8 个、Non ileal paneth 42 个；这是联合类，不能直接当作纯回肠 Paneth 参考。\n3. Tuft 是 Ileal tuft 7 个、Non ileal tuft 43 个；Glial 全部为作者 Glia 注释。\n4. BEST4 在此次选中细胞中对应 Non ileal BEST4 OTOP2 3 个，合并到 Colonic_absorptive，未单独设为一列；也没有单列 T_other。\n5. DC 包含 CD1Chi 32、LAMP3pos 8、LAMP3pos IL1Bpos 2、XCR1pos 5、pDC 3 个，展示了大类型内部并不等额。\n6. 回肠吸收细胞的线粒体比例中位数为 45.10%，较高；当前 50 细胞测试没有因此重抽。是否调整应结合 marker、文库质量及留出预测，不能把本清单当作已经完成最终 QC 优化。\n\n'
report += '## 全部作者亚型\n\n' + subtypes
report += '\n逐细胞原始依据：[selected_cell_manifest.tsv](' + (source / 'selected_cell_manifest.tsv').as_posix() + ')。新增文件均为汇总与展示，没有修改原 signature、输入矩阵或留出设计。\n'
(out / '联合参考实际抽样展示.md').write_text(report, encoding='utf-8')
print('已保存:', out / '联合参考实际抽样展示.md')
