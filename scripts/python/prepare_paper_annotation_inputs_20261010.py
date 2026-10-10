# 为EcoTyper原文解释步骤准备冻结定义和文献marker，原始数据只读。
# 本机cibersortx_review环境运行；输出单独目录，供服务器R脚本读取。
# 文献短面板是人工挑选的读出基因，不冒充完整signature或正式富集。
from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET
import pandas as pd

# Step1. 保存输入副本，保持53状态和7生态型的冻结定义
project = Path(__file__).resolve().parents[3] / 'IBD_ecotype_analysis'
validation = project / '06_validation'
out = validation / 'ecotyper_original_method_annotation_20261010'
inputs = out / 'inputs'
inputs.mkdir(exist_ok=True)
shutil.copy2(validation / 'recovery_interpretation_20261010/53状态_初步判断及来源.tsv', inputs / 'frozen_states.tsv')
shutil.copy2(validation / 'recovery_interpretation_20261010/7生态型_初步判断及来源.tsv', inputs / 'frozen_ecotypes.tsv')
shutil.copy2(validation / 'bayesprism_recovery_readiness_20261009/53状态全部关联基因.tsv', inputs / 'discovery_associated_genes.tsv')

# Step2. 完整提取Martin正文明确列出的五个T细胞程序
paragraphs = [p.findtext('text', '') for p in ET.parse(out / 'sources/PMC7060942_bioc.xml').findall('.//passage')]
programs = [('cell cycling', '细胞周期', 'Any_shared_program'), ('immunoregulation', '免疫调节T程序', 'CD4_T'), ('naive/central memory (CM)', '初始_中央记忆T程序', 'CD4_T;CD8_T'), ('CD8/cytotoxic', '细胞毒程序', 'CD8_T;ILC'), ('resident memory T cell', '组织驻留记忆T程序', 'CD4_T;CD8_T')]
rows = []
for original, label, lineage in programs:
    match = re.search(re.escape(original) + r' \(([^)]+)\)', paragraphs[61])
    assert match is not None, original
    genes = [g.strip() for g in match.group(1).split(',')]
    for gene in genes:
        rows.append([label, gene, lineage, 'Martin2019', 'PMC7060942', 61, '正文该程序完整基因列表', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7060942/'])

# Step3. 肠道文献短面板，用于有来源的程序解释而非自动亚型命名
panels = [
    ['驻留巨噬细胞读出', 'MAFB,CSF1R,DNASE1L3,C1QA,C1QB,C1QC,MRC1', 'Monocyte_Macrophage', 'Martin2019', 'PMC7060942', 13],
    ['炎症巨噬细胞读出', 'CXCL2,CXCL3,CXCL8,IL23A,IL6,IL1B,IL1A,TNF,OSM', 'Monocyte_Macrophage', 'Martin2019', 'PMC7060942', 14],
    ['活化DC读出', 'LAMP3,CD40,CD86,CD274,CCR7,CCL17,CCL22,CCL19', 'DC', 'Martin2019', 'PMC7060942', 14],
    ['活化内皮读出', 'ACKR1,SELE,SELP', 'Endothelial', 'Martin2019', 'PMC7060942', 16],
    ['活化成纤维细胞读出', 'THY1,PDPN,CTHRC1,CHI3L1,IL6,IL11,CXCL2,CXCL8,CXCL1,CXCL5,CCL2,CCL7', 'Fibroblast', 'Martin2019', 'PMC7060942', 16],
    ['周细胞身份读出', 'CSPG4,RGS5,MCAM', 'Pericyte', 'Martin2019', 'PMC7060942', 16],
    ['平滑肌程序读出', 'ACTG2,MYH11,DES', 'Pericyte', 'Martin2019', 'PMC7060942', 16],
    ['FRC样淋巴基质读出', 'CCL19,CCL21,TNFSF14,CD74,CLU,CD24,FDCSP', 'Fibroblast', 'Kinchen2018', 'PMC6176871', '59;64'],
    ['炎症相关成纤维细胞IAF读出', 'IL11,IL24,IL13RA2,TNFRSF11B', 'Fibroblast', 'Smillie2019', 'PMC6662628', '25;39'],
    ['BEST4上皮读出', 'BEST4,OTOP2,CA7', 'Epithelial', 'Smillie2019', 'PMC6662628', 16],
    ['上皮支持基质读出', 'WNT2B,WNT4,DKK3,RSPO3', 'Fibroblast', 'Smillie2019', 'PMC6662628', '17;18'],
    ['胶质身份读出', 'FOXD3,MPZ,CDH19,PLP1,SOX10,S100B,ERBB3', 'Glial', 'Elmentaite2021', 'PMC8426186', 85],
    ['T区网状基质读出', 'CCL19,GREM1,TNFSF13B', 'Fibroblast', 'Elmentaite2021', 'PMC8426186', 87],
    ['滤泡树突基质读出', 'CXCL13,CR1,CR2', 'Fibroblast', 'Elmentaite2021', 'PMC8426186', 87],
    ['肌成纤维细胞读出', 'ACTA2,TAGLN,DCN', 'Fibroblast', 'Elmentaite2021', 'PMC8426186', 87]
]
for label, genes, lineage, source, pmc, passage in panels:
    for gene in genes.split(','):
        rows.append([label, gene, lineage, source, pmc, passage, '正文人工选择的短读出面板，非完整signature', 'https://pmc.ncbi.nlm.nih.gov/articles/' + pmc + '/'])
reference = pd.DataFrame(rows, columns=['Program', 'Gene', 'ExpectedLineage', 'Publication', 'PMC', 'BioCPassageIndex0', 'PanelStatus', 'URL'])
reference['SymbolNote'] = '正文直接列出；短面板保留选择性质'
reference.loc[reference.Gene.isin(['C1QA', 'C1QB', 'C1QC']) & reference.Publication.eq('Martin2019'), 'SymbolNote'] = '正文C1Qs家族读出，人工展开成人基因符号，非正文逐项列表'
reference.loc[reference.Gene.eq('IL23A'), 'SymbolNote'] = '正文IL-23蛋白读出，人工选择p19编码基因IL23A；正文未逐项指定亚基'
reference.loc[reference.Gene.isin(['IL6', 'IL1B', 'IL1A', 'CD274']) & reference.Publication.eq('Martin2019'), 'SymbolNote'] = '正文IL-6/IL-1β/IL-1α/PD-L1名称转换为对应人基因符号'
reference.loc[reference.Gene.eq('GREM1'), 'SymbolNote'] = '正文GREM1hi按基因GREM1记录，hi是表达限定词'
assert not reference.duplicated(['Program', 'Gene']).any()
reference.to_csv(inputs / 'literature_marker_panels.tsv', sep='\t', index=False)
print('文献程序数、基因记录数：', reference.Program.nunique(), len(reference))
print(reference.groupby('Program').size().to_string())
