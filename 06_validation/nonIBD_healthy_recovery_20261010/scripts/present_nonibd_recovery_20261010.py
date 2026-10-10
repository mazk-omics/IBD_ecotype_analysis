# 汇总健康/nonIBD recovery 的原始编号、覆盖和共同变化，不进行文献命名。
# 在本地 cibersortx_review Python 环境运行；读取本轮 audit，另存阅读报告和可核对的小表。
# bulk 的最高状态分配、单细胞 Z 和生态型共关联分别展示，不能互相替代。
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

# Step1. 读取本轮正式结果，核对全体冻结状态和任务状态
parser = argparse.ArgumentParser()
parser.add_argument('--root', required=True)
args = parser.parse_args()
root = Path(args.root)
audit = root / 'audit'
summary = pd.read_csv(audit / 'dataset_summary.tsv', sep='\t')
states = pd.read_csv(audit / '53_states_all_datasets.tsv', sep='\t')
ecotypes = pd.read_csv(audit / 'frozen_ecotypes.tsv', sep='\t')
counts = pd.read_csv(audit / '7_ecotypes_all_datasets.tsv', sep='\t')
coassociation = pd.read_csv(audit / 'healthy_ecotype_coassociation.tsv', sep='\t')
integrity = json.loads((audit / 'integrity_check.json').read_text())
assert len(summary) == 6 and len(states) == 318 and len(coassociation) == 14
assert integrity['Status'] == 'PASS'
print(summary[['Dataset', 'InputUnits', 'StatesWithAssignments', 'StatesAboveNominalZ196']].to_string(index=False))

# Step2. 补充成员实际出现数量；它是覆盖描述，不能把部分成员当作完整生态型恢复
support_rows = []
for task in summary.itertuples(index=False):
    assignments = pd.read_csv(audit / task.Dataset / 'unit_assignments_all.tsv', sep='\t')
    eco_assignment = pd.read_csv(audit / task.Dataset / 'ecotype_assignments_complete.tsv', sep='\t')
    assignments['RecoveryUnit'] = assignments.Sample if task.Assay == 'scRNA' else assignments.ID
    presence = pd.crosstab(assignments.RecoveryUnit, assignments.CellState).gt(0)
    current = states.loc[states.Dataset.eq(task.Dataset)].set_index('CellState')
    for eco in ecotypes.itertuples(index=False):
        members = eco.Members.split(';')
        selected = current.loc[members]
        member_presence = presence.reindex(columns=members, fill_value=False).sum(axis=1)
        chosen = eco_assignment.loc[eco_assignment.Ecotype.eq(eco.Ecotype), 'ID']
        chosen_presence = member_presence.reindex(chosen)
        support_rows.append({'Dataset': task.Dataset, 'Assay': task.Assay, 'Ecotype': eco.Ecotype, 'FrozenMembers': len(members), 'LineageCoveredMembers': int(selected.LineagePresent.sum()), 'MembersWithAnyAssignments': int(selected.AssignedUnits.gt(0).sum()), 'MembersWithAtLeast5Assignments': int(selected.AssignedUnits.ge(5).sum()), 'MembersWithNominalZAbove196': int(selected.Z.gt(1.96).sum()) if task.Assay == 'scRNA' else np.nan, 'OfficialAssignedUnits': len(chosen), 'MedianDetectedMembersInAssignedUnits': chosen_presence.median(), 'AllMembersDetectedAssignedUnits': int(chosen_presence.eq(len(members)).sum())})
support = pd.DataFrame(support_rows)
support.to_csv(audit / 'ecotype_member_presence.tsv', sep='\t', index=False)
for value in ['AssignedUnits', 'AssignedDonors', 'Z', 'OwnMinusBestOther_InternalOnly']:
    states.pivot(index=['Ecotype', 'CellType', 'CellState'], columns='Dataset', values=value).to_csv(audit / ('state_comparison_' + value + '.tsv'), sep='\t')

# Step3. 同一队列对照、外部 bulk 和既有单细胞健康子集分开阅读
labels = {'MSCCR_nonIBD_Ileum_20261010': 'MSCCR nonIBD 回肠', 'MSCCR_nonIBD_Colon_20261010': 'MSCCR nonIBD 结肠', 'MSCCR_nonIBD_Rectum_20261010': 'MSCCR nonIBD 直肠', 'GSE57945_nonIBD_Ileum_20261010': 'GSE57945 nonIBD 回肠', 'IBDverse_Healthy_Ileum_20261010': 'IBDverse 健康回肠', 'GSE214695_HC_Colon_20261010': 'GSE214695 结肠对照'}
overview = summary[['Dataset', 'Assay', 'InputUnits', 'InputDonors', 'LineagesPresent', 'StatesWithAssignments', 'StatesAboveNominalZ196', 'EcotypeAssignedBiopsies']].copy()
overview['Dataset'] = overview.Dataset.map(labels)
overview.columns = ['数据', '模态', '输入活检或细胞', '供者键数量', '父类型', '有分配的状态/53', '单细胞名义Z>1.96', '生态型分配活检数']
overview.to_csv(root / '恢复结果概览.tsv', sep='\t', index=False)
eco_count_table = counts.pivot(index='Ecotype', columns='Dataset', values='AssignedBiopsies').rename(columns=labels)
health_z = states.loc[states.Assay.eq('scRNA') & states.Z.gt(1.96), ['Dataset', 'Ecotype', 'CellType', 'CellState', 'AssignedUnits', 'AssignedDonors', 'Z', 'OwnMinusBestOther_InternalOnly']].copy()
health_z['Dataset'] = health_z.Dataset.map(labels)
health_z.to_csv(root / '单细胞_名义Z候选状态.tsv', sep='\t', index=False)
z_comparison = states.loc[states.Assay.eq('scRNA')].pivot(index='CellState', columns='Dataset', values='Z')
shared = z_comparison.gt(1.96).all(axis=1)
z_comparison.loc[shared].to_csv(root / '两个对照子集_共同名义Z候选.tsv', sep='\t')

# Step4. 保存操作过程与边界；不把投影分配解释成疾病特异性或健康丰度
steps = '''# 健康/nonIBD recovery 操作与阅读说明

模型固定为 MSCCR_CORE12_CTSPEC_REGION_IBD：CD+UC discovery、Region3 组内标准化、BayesPrism presort，12 个父类型、53 个状态、7 个生态型。原始输入、冻结模型和官方 recovery 源码没有改写；指纹结果见 audit/integrity_check.json。

## 1. 选择对照

MSCCR 从 canonical annotation 中选 Disease=Control、Inflammation=Non-inflamed，461 个活检、243 个 PatientID。与 CD+UC discovery 的样本 ID、PatientID 均无重叠。回肠 121、结肠 115、直肠 225 个活检；供者可跨部位，不能把分部位人数相加。Control 是项目对照定义，不自动等于完全无症状健康人。

GSE57945 从已有 canonical metadata 选 Non-IBD，42 个回肠样本；未提供可独立确认的患者 ID，以样本 ID 为供者键，不能据此宣称 42 个已核验独立供者。IBDverse 选既有 recovery pilot 的 Healthy 10 位供者、1,362 细胞；GSE214695 选 HC 6 位供者、3,426 细胞。后二者是已有数据的对照子集复算，不是新独立队列；计数受原 pilot 的父类型抽样上限影响，不代表组织真实细胞比例。

GSE214695 的非IBD手术对照存在肿瘤背景：论文说明取非肿瘤区组织。它可作 nonIBD 比较，但不能无条件标为严格健康。来源：[原论文](https://www.nature.com/articles/s41467-023-40156-6)。GSE57945 的研究背景见 [Haberman et al., JCI 2014](https://jci.org/articles/view/75436)。本轮数量均来自实际文件筛选，不从论文总人数推算。

## 2. 准备表达

MSCCR 使用原始 counts 全部基因行计算每个样本总文库量，再转线性 CPM；模型需要的基因使用已审计 HGNC 映射，同符号计数求和。仅导出 4,553 个冻结投影基因的并集，不在导出后重算 CPM。GSE57945 沿用已整理线性 RPKM；两个单细胞子集直接保留现有 CPM 数值。未另行归一化、聚类、细胞注释或筛选 marker。

各 MSCCR 类型的 1,000 个投影基因全部覆盖；GSE57945 各类型覆盖 940–982 个。IBDverse 无 Glial 注释，该谱系记不可评估；GSE214695 周细胞仅 7 个，稀有类型缺乏统计支撑。详见 input/model_gene_coverage.tsv 与单细胞 donor_type_counts.tsv。

## 3. 调用官方 recovery

使用既有 ecotyper conda 环境，tmux 中顺序运行六套数据，-t 4。bulk 调用 EcoTyper_recovery_bulk.R；MSCCR 按 Region3 分开，以官方行标准化执行。单细胞调用 EcoTyper_recovery_scRNA.R，-s -1 -z TRUE，保留子集所有输入细胞。未改变官方代码。

bulk 官方流程没有计算 cell-state recovery Z，最高状态分配只能作投影描述。单细胞使用当前官方默认 30 次状态打乱；并非 Cell 2021 的 100 次状态打乱完整复现。本轮 Z>1.96 仅列为名义探索结果，未做 106 个状态检验的多重校正。

## 4. 审计状态和生态型

保留每队列全部 53 状态，包括无分配和缺失谱系；内部 marker 差值是“本状态已分配对象的 marker 平均分−最佳其他状态”的同源一致性检查，不是独立准确率。状态候选表应同时看细胞数、供者数、Z 和 marker 方向。

对两个单细胞健康子集，沿用上一轮 Cell 2021 方法：F1 分母为已分配的父类型细胞，F2 分母为全部输入父类型细胞；分别计算 F1/F2、零改 NA 后的 F1/F2，共四套 Spearman 矩阵，四套均有限才取平均。生态型统计量为成员内关联减成员外关联，行置换 1,000 次，对角线 NA 移回原位。参数、F1/F2、边覆盖、零分布、Z 全部保存。

四套矩阵 NA 合并和有限边平均是已明确记录的实现选择，不能称为完整数值复现。仅 6/10 位供者、欠抽样和缺失成员使共关联不稳定；少数可计算成员的 Z 不能证明完整生态型恢复。依据：[Luca et al., Cell 2021 方法段](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/)。

## 5. 如何使用本轮结果

先看恢复结果概览.tsv；再看 audit/53_states_all_datasets.tsv 和 7_ecotypes_all_datasets.tsv；生态型另看 healthy_ecotype_coassociation.tsv、ecotype_member_presence.tsv 及 pair_N / correlation_coverage。formal 下本地保存官方状态分配和丰度；audit 各队列目录保存生态型热图。较大的官方图件和热图表达数据留在服务器，路径和校验值见 server_retained_large_outputs.tsv；完整压缩包也保留在服务器，不需要重新计算。

健康子集检出某状态，提示它未必仅出现在 IBD；未检出不等于健康中不存在。MSCCR 各健康部位内部标准化，不能把它与 discovery 的状态丰度直接作为绝对大小比较，也不能由不同部位的分配比例推断疾病风险。本轮没有差异检验、临床关联、功能富集或文献状态命名；下一步解释需用户先审阅文献。

环境锁定文件、脚本、运行日志和时间戳均随结果保存。Python 输入审计使用 ibd_refaudit；共关联使用此前独立 ecotyper_interpretation 环境；原 ecotyper 环境不安装新包。
'''
(root / '操作步骤与方法边界.md').write_text(steps, encoding='utf-8')
findings = f'''# 本轮初步结果：健康 / nonIBD recovery

六套官方 recovery 均正常结束。冻结模型仍为 MSCCR_CORE12_CTSPEC_REGION_IBD，未重训、未改写原始输入；93 个模型/源码文件、9 个源输入指纹保持不变。357 个回传文件逐项 SHA256 核验通过；本地另生成的展示表不在该传输快照清单中。

## 状态恢复

- MSCCR 的461个nonIBD非炎症活检、243位供者，按回肠/结肠/直肠分别恢复：121/115/225个活检中，有状态分配的冻结状态数分别为53/52/53。每个部位7个生态型都存在官方分配样本。bulk没有官方状态recovery Z，因此这些是投影覆盖，不是53个状态均通过独立验证。
- 外部GSE57945的42个nonIBD回肠样本：51个状态有分配，6个生态型有分配。仅样本ID被核对，供者ID未独立确认。
- IBDverse健康回肠：10位供者、1,362细胞；11个父类型，47个冻结状态有谱系输入，44个状态有细胞分配，40个名义Z>1.96。无胶质注释的6个状态不可评估。
- GSE214695结肠对照：6位供者、3,426细胞；12个父类型，50个状态有细胞分配，37个名义Z>1.96。周细胞只有7个；HC是nonIBD非肿瘤区对照，不能无条件称严格健康。
- 两个单细胞对照子集共同有 **{int(shared.sum())}个状态名义Z>1.96**。完整编号、细胞数、供者数、内部marker差值都保存。这是当前官方30次状态打乱的探索读出，未做多重检验校正，不能直接称31个状态验证成功。

## 生态型共同变化

采用四套F1/F2 Spearman矩阵、1,000次行置换的共关联检查。E3在两个对照子集的探索Z分别为2.38（回肠）和4.16（结肠）；可计算成员分别为10/12和11/12，尚未形成全部12成员的完整恢复证据。它是后续检查“共享肠道组织状态”时值得保留的对照对象，不在本轮赋予功能名称。

E6仅回肠出现名义正向共关联（Z=2.17、6/9成员可计算）；E5仅结肠出现较强读出（Z=3.72、4/4可计算）。这些对照仅10/6位供者，稀有状态、相关边缺失和原pilot抽样都会影响结果。E4回肠没有可计算成员内关联，不能把缺失统计量当作不恢复的证据。

## 目前能回答什么

部分冻结状态在健康/nonIBD背景中也可检出并获得恢复支持；不能仅依据它们来自CD+UC discovery而命名为“IBD专属”。生态型还需要成员共同变化和覆盖支持，最高分配不够。这里也没有证明健康与IBD无差异：两者并未在共同标度、供者层级的差异分析中比较。

本轮结果是既有BayesPrism presort模型的对照恢复，不是新的CIBERSORTx discovery。两个单细胞组是此前recovery队列的对照子集，不是新增独立队列。单细胞限量抽样数不代表组织细胞比例；MSCCR各健康部位内部标准化后的丰度不能直接与discovery绝对比较。

## 文件和执行记录

入口：健康nonIBD恢复审阅.html；过程：操作步骤与方法边界.md；新候选：候选健康数据清单.md。audit含318行全状态表、生态型分配、成员出现覆盖、共关联矩阵和零分布。formal本地保留官方分配/丰度小表；大状态热图和热图表达数据保留服务器，路径见server_retained_large_outputs.tsv。

服务器：/home/mazekai/IBD_EcoTyper/06_validation/nonIBD_healthy_recovery_20261010；tmux运行开始2026-10-10 20:34:31、结束21:08:05（北京时间）。ecotyper用于官方恢复，ibd_refaudit用于准备/审计，ecotyper_interpretation的R4.0.3用于共关联；环境、脚本、日志、退出码均保存。原始全量官方结果快照nonIBD_healthy_recovery_review.zip也保留服务器，无需重算。

文献命名、功能解析和临床关联按用户安排暂缓。样本来源和方法文献链接见操作说明及候选清单。
'''
(root / '本轮初步结果.md').write_text(findings, encoding='utf-8')
tables = [('六套恢复的输入和分配覆盖', overview), ('官方生态型分配数量（不是独立验证）', eco_count_table), ('单细胞名义 Z 候选状态', health_z), ('健康子集生态型共关联与可计算成员', coassociation[['Dataset', 'Ecotype', 'FrozenMembers', 'DefinedReferenceMembers', 'AvgDif', 'Z', 'ValidPermutations']]), ('生态型成员实际出现覆盖', support)]
body = ['<h1>健康 / nonIBD 恢复审阅</h1><p>冻结：CD+UC、Region3 归一化、BayesPrism presort｜53 状态｜7 生态型。保持原始编号，不在本轮命名。</p>', '<p class="note">先区分“被分配”“状态恢复支持”“生态型共同变化”。三者并不等价。仅有 6/10 位单细胞供者，Z 为探索结果；bulk 无官方状态 Z。对照中检出不等于已证明与 IBD 无差异。</p>', f'<p>两个单细胞对照子集共同有{int(shared.sum())}个状态达到名义Z&gt;1.96。E3在两组均有正向共关联读出（Z=2.38/4.16），但可计算成员只有10/12、11/12。完整初步判断见下面入口。</p>', '<p><a href="本轮初步结果.md">本轮初步结果</a> · <a href="操作步骤与方法边界.md">操作步骤与方法边界</a> · <a href="候选健康数据清单.md">候选数据及来源</a> · <a href="audit/53_states_all_datasets.tsv">完整 53 状态×6 数据表</a> · <a href="audit/integrity_check.json">数据与源码指纹核验</a></p>']
for title, table in tables:
    body.append('<h2>' + title + '</h2><div class="table">' + table.to_html(index=table is eco_count_table, na_rep='不可评估', float_format=lambda x: f'{x:.3f}') + '</div>')
body.append('<p>MSCCR：461 个对照活检、243 位供者，各部位分开。GSE57945：42 个 nonIBD 回肠样本。IBDverse 和 GSE214695：既有恢复数据的对照子集；回肠无胶质注释，结肠周细胞仅7个。健康细胞的抽样数量不是组织比例。详细限制和执行命令见说明与 scripts。</p>')
html = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>健康nonIBD recovery审阅</title><style>body{font:16px/1.6 Segoe UI,Microsoft YaHei,sans-serif;color:#203044;margin:30px auto;max-width:1500px;padding:0 24px}h1,h2{color:#17415a}.note{background:#fff3d9;padding:18px;border-left:4px solid #d39823}.table{overflow:auto}table{border-collapse:collapse;font-size:13px}th,td{padding:7px 10px;border:1px solid #d5dfe5;white-space:nowrap}th{background:#eef3f6}tr:nth-child(even){background:#fafcfd}a{color:#126b91}</style>' + ''.join(body) + '</html>'
(root / '健康nonIBD恢复审阅.html').write_text(html, encoding='utf-8')
print('Saved:', root / '健康nonIBD恢复审阅.html')
