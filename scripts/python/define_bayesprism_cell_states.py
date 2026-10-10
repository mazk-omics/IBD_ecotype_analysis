# 读取冻结discovery的53个状态基因，附上初步程序定义与MSigDB数据库重叠证据。
# 本机cibersortx_review环境运行；不改状态标签、表达矩阵或discovery，也不计算缺乏真实背景的富集P值。
from pathlib import Path
from urllib.request import urlopen
import hashlib
import json
from datetime import datetime, timezone, timedelta
import pandas as pd

# Step1. 读取模型和人工初步程序，按原状态ID对应
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
discovery = project / 'ctspec_region_ibd'
outdir = project / '06_validation/bayesprism_recovery_readiness_20261009'
database = outdir / 'database'
database.mkdir(parents=True, exist_ok=True)
ecotypes = pd.read_csv(discovery / 'Ecotypes/ecotypes.txt', sep='\t')
draft = pd.read_csv(outdir / 'cell_state_program_draft.tsv', sep='\t')
assert len(draft) == 53 and draft.CellState.is_unique and set(draft.CellState) == set(ecotypes.ID)
print('模型：', ecotypes.CellType.nunique(), '类型；', len(draft), '状态；', ecotypes.Ecotype.nunique(), 'ecotype')

# Step2. 下载预先选定的16个Hallmark程序，保存官方原件和哈希
# 这是初步解释的定向数据库面板，不是全库富集；不以命中数最高者自动决定状态名称。
programs = ['INTERFERON_ALPHA_RESPONSE', 'INTERFERON_GAMMA_RESPONSE', 'TNFA_SIGNALING_VIA_NFKB', 'INFLAMMATORY_RESPONSE', 'E2F_TARGETS', 'G2M_CHECKPOINT', 'MITOTIC_SPINDLE', 'HYPOXIA', 'UNFOLDED_PROTEIN_RESPONSE', 'GLYCOLYSIS', 'OXIDATIVE_PHOSPHORYLATION', 'FATTY_ACID_METABOLISM', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'TGF_BETA_SIGNALING', 'COMPLEMENT', 'APOPTOSIS']
gene_sets = {}
provenance = []
for program in programs:
    name = 'HALLMARK_' + program
    url = 'https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?fileType=json&geneSetName=' + name
    local = database / (name + '.json')
    if not local.exists():
        with urlopen(url, timeout=30) as response:
            content = response.read()
        parsed = json.loads(content)
        assert name in parsed and len(parsed[name]['geneSymbols']) > 0
        local.write_bytes(content)
    parsed = json.loads(local.read_bytes())
    gene_sets[name] = set(parsed[name]['geneSymbols'])
    provenance.append({'GeneSet': name, 'Genes': len(gene_sets[name]), 'URL': url, 'SHA256': hashlib.sha256(local.read_bytes()).hexdigest(), 'RetrievedDate': '2026-10-09'})
    print('数据库：', name, len(gene_sets[name]))
pd.DataFrame(provenance).to_csv(database / 'MSigDB来源记录.tsv', sep='\t', index=False)

# Step3. 逐类型读全部状态关联基因，保留原MaxFC及排名方便复核
# MaxFC是此流程的状态差异分数，不当作原始倍数变化或独立差异表达显著性。
markers_all = []
definitions = []
overlaps = []
for cell_type in ecotypes.CellType.unique():
    gene_info = pd.read_csv(discovery / cell_type / 'gene_info.txt', sep='\t')
    for state, markers in gene_info.groupby('State', sort=True):
        markers = markers.sort_values('MaxFC', ascending=False).copy()
        marker_set = set(markers.Gene)
        state_id = cell_type + '_' + state
        markers['CellState'] = state_id
        markers['WithinStateRank'] = range(1, len(markers) + 1)
        markers_all.append(markers[['CellState', 'Gene', 'MaxFC', 'InitialState', 'WithinStateRank']])
        for program, gene_set in gene_sets.items():
            hits = sorted(marker_set & gene_set)
            overlaps.append({'CellState': state_id, 'GeneSet': program, 'StateAssociatedGenes': len(marker_set), 'PanelGeneSetSize': len(gene_set), 'OverlapCount': len(hits), 'OverlapGenes': ','.join(hits), 'StatisticalEnrichmentTested': False})
        definitions.append({'CellState': state_id, 'CellType': cell_type, 'State': state, 'Ecotype': ecotypes.loc[ecotypes.ID.eq(state_id), 'Ecotype'].iloc[0], 'MarkerCount': len(markers), 'Top20Genes': ','.join(markers.Gene.head(20)), 'CellLineageValidated': False, 'DefinitionStatus': '初步程序名称；非最终细胞亚型定义'})
pd.concat(markers_all, ignore_index=True).to_csv(outdir / '53状态全部关联基因.tsv', sep='\t', index=False)
overlap_table = pd.DataFrame(overlaps)
overlap_table.to_csv(outdir / '53状态_MSigDB定向面板重叠.tsv', sep='\t', index=False)
definition_table = pd.DataFrame(definitions).merge(draft, on='CellState', validate='one_to_one')

# Step4. 保存易读定义表，展示数据库实际命中而不宣称统计富集
panel_summary = []
for state_id, evidence in overlap_table.groupby('CellState', sort=False):
    best = evidence.loc[evidence.OverlapCount.gt(0)].sort_values('OverlapCount', ascending=False).head(3)
    panel_summary.append({'CellState': state_id, 'TopPanelOverlaps': '; '.join(row.GeneSet.replace('HALLMARK_', '') + ':' + str(row.OverlapCount) for row in best.itertuples())})
definition_table = definition_table.merge(pd.DataFrame(panel_summary), on='CellState', validate='one_to_one')
definition_table.to_csv(outdir / '53状态初步定义.tsv', sep='\t', index=False)
lines = ['# 53个cell state的初步程序定义', '', '读取本地冻结discovery及MSigDB官方16个Hallmark程序。名称是人工基因程序解释，不改变S编号，不证明谱系来源；定向重叠未计算富集P值。真实NMF基因背景、服务器W及外部覆盖尚未验收。', '', '| 状态 | Ecotype | 初步程序 | 程序判断 | 主要标志 |', '|---|---|---|---|---|']
for row in definition_table.itertuples():
    lines.append('| ' + row.CellState + ' | ' + row.Ecotype + ' | ' + row.InitialProgram + ' | ' + row.ProgramConfidence + ' | ' + ', '.join(row.Top20Genes.split(',')[:6]) + ' |')
lines += ['', '详细依据、数据库命中基因和待审原因见同目录TSV。MSigDB：[Human Hallmark](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)。谱系疑点辅助来源：[CSF3R/HPA](https://www.proteinatlas.org/ENSG00000119535-CSF3R/single%2Bcell)、[CUBN/HPA](https://www.proteinatlas.org/ENSG00000107611-CUBN/single%2Bcell)。', '', '低判断状态保留标志组合名称，不强制命名Th17、耗竭T细胞、M1/M2、TLS、成熟胶质等细胞亚型。']
(outdir / '53状态初步定义.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
print(definition_table[['CellState', 'InitialProgram', 'TopPanelOverlaps']].to_string(index=False))

# Step5. 复核本地保存的GSE72819注释；表达当前仍在服务器
annotation_path = project / '06_validation/recovery_inventory_20261007/server_snapshot/06_validation/02_bulk_treatment/GSE72819/prepared/06_ecotyper_recovery_input/GSE72819_bulk_recovery_annotation.txt'
annotation = pd.read_csv(annotation_path, sep='\t')
assert annotation.ID.is_unique and len(annotation) == 73 and annotation.Disease.eq('UC').all() and annotation.Timepoint.eq('Baseline').all()
summary = {'time_Asia_Shanghai': datetime.now(timezone(timedelta(hours=8))).isoformat(), 'states_defined': len(definition_table), 'database_panel_sets': len(gene_sets), 'formal_enrichment_tested': False, 'GSE72819_samples_local_annotation': len(annotation), 'GSE72819_unique_patient_ids_local_annotation': int(annotation.PatientID.nunique()), 'GSE72819_response_labels': annotation.ResponseWeek10.fillna('Missing').value_counts().to_dict(), 'GSE72819_source': 'EUCALYPTUS Phase II etrolizumab UC baseline colonic biopsies; not RISK', 'current_server_connection': 'single attempt timed out before login, exit 255', 'server_ecotyper_conda_verified_this_turn': False, 'recovery_started': False}
(outdir / '本次准备状态.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
