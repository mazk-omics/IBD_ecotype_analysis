# 复用冻结discovery，为原文Fig2A/B、4A/B、5A整理输入和可核查源表。
# 服务器ibd_refaudit环境运行；读原始导出与既有Hallmark结果，写本批次source_data。
# 不重跑NMF、不改状态或生态型；比例分支使用已有层级合并BayesPrism结果。
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 定位冻结输出，记录真正读取的输入SHA256
root = Path('/home/mazekai/IBD_EcoTyper')
out = root / '06_validation/paper_figures_20261012'
prior = root / '06_validation/paper_aligned_recovery_characterization_20261012'
model = root / 'software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes'
export = root / '05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd'
source = out / 'source_data'
states = pd.read_csv(prior / '00_audit/frozen_states.tsv', sep='\t').sort_values(['Ecotype','CellType','State'])
metadata = pd.read_csv(prior / '06_ecotype_features/ecotype_dominant_assignments.tsv', sep='\t').set_index('ID')
assert len(states)==53 and len(metadata)==2025 and metadata.index.is_unique
states.to_csv(source / 'frozen_states.tsv', sep='\t', index=False)
metadata.to_csv(source / 'all_2025_biopsies.tsv', sep='\t')
manifest = []
inputs = [prior/'00_audit/frozen_states.tsv', prior/'06_ecotype_features/ecotype_dominant_assignments.tsv']
heatmap_audit = []

# Step2. Fig2A沿用官方heatmap_data和heatmap_simple(scale_rows=TRUE)，保留所有已导出的基因和分配样本
for ct in states.CellType.unique():
    folder = export / ct
    data = pd.read_csv(folder/'heatmap_data.txt', sep='\t', index_col=0)
    genes = pd.read_csv(folder/'gene_info.txt', sep='\t')
    annotation = pd.read_csv(folder/'heatmap_top_ann.txt', sep='\t')
    assignment = pd.read_csv(folder/'state_assignment.txt', sep='\t')
    abundance = pd.read_csv(folder/'state_abundances.txt', sep='\t', index_col=0)
    assert set(data.columns)==set(assignment.ID) and set(data.index)==set(genes.Gene)
    annotation = annotation.set_index('ID').loc[data.columns]
    genes = genes.set_index('Gene').loc[data.index]
    # 官方热图以样本为列做每行均值0、样本标准差1；已有Region3归一化的来源不改动
    relative = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1, ddof=1), axis=0)
    assert np.isfinite(relative.to_numpy()).all()
    relative.to_csv(source/(ct+'_Fig2A_relative_expression.tsv.gz'), sep='\t', compression='gzip')
    genes.to_csv(source/(ct+'_Fig2A_gene_order.tsv'), sep='\t')
    annotation.to_csv(source/(ct+'_Fig2A_sample_order.tsv'), sep='\t')
    abundance.loc[:,assignment.ID].T.to_csv(source/(ct+'_Fig2B_Hstar_assigned.tsv'), sep='\t')
    heatmap_audit.append({'CellType':ct,'Genes':len(data),'AssignedBiopsies':len(assignment),'AllBiopsies':len(metadata),'ExcludedUnassignedBiopsies':len(metadata)-len(assignment)})
    inputs.extend(folder/x for x in ['heatmap_data.txt','gene_info.txt','heatmap_top_ann.txt','state_assignment.txt','state_abundances.txt'])
    print('Fig2 input:', ct, relative.shape, flush=True)
pd.DataFrame(heatmap_audit).to_csv(source/'Fig2_input_coverage.tsv',sep='\t',index=False)

# Step3. Fig4A只显示官方已分配到生态型的活检；排序按冻结Ecotype，保持现有组成
assigned = pd.read_csv(model/'Ecotypes/discovery/ecotype_assignment.txt', sep='\t')
assigned = assigned[assigned.Ecotype.isin(states.Ecotype)].copy()
assigned = assigned.sort_values(['Ecotype','ID']).set_index('ID')
assert assigned.index.is_unique and set(assigned.index)<=set(metadata.index)
combined = pd.read_csv(model/'Ecotypes/discovery/combined_state_abundances.txt',sep='\t',index_col=0)
combined = combined.loc[states.CellState,assigned.index]
combined.to_csv(source/'Fig4A_state_abundance.tsv',sep='\t')
assigned.to_csv(source/'Fig4A_assigned_sample_order.tsv',sep='\t')
jaccard = pd.read_csv(model/'Ecotypes/discovery/jaccard_matrix.txt',sep='\t',index_col=0)
jaccard = jaccard.loc[states.CellState,states.CellState]
jaccard.to_csv(source/'Fig4B_frozen_jaccard.tsv',sep='\t')
# 直接读取冻结矩阵，保留官方P屏蔽规则；不改成新统计检验，也不重新聚类
edges = []
for ecotype, member in states.groupby('Ecotype',sort=False):
    for i,a in enumerate(member.CellState):
        for b in member.CellState.to_list()[i+1:]:
            edges.append({'Ecotype':ecotype,'From':a,'To':b,'Jaccard':jaccard.loc[a,b]})
pd.DataFrame(edges).to_csv(source/'Fig4B_all_within_ecotype_edges.tsv',sep='\t',index=False)
inputs.extend(model/'Ecotypes/discovery'/x for x in ['ecotype_assignment.txt','combined_state_abundances.txt','jaccard_matrix.txt','ecotypes.txt'])
print('Fig4A:',combined.shape,'Not ecotype-assigned:',2025-len(assigned),flush=True)

# Step4. Fig5A细胞组成：全活检比例求dominant-ecotype组均值，不把12类重新归一化
fraction_path = root/'03_reference/combined_reference/output/22_hierarchical_fraction_merge/bayesprism_hierarchical_full20_final_fraction.tsv'
fractions = pd.read_csv(fraction_path,sep='\t').set_index('sample_id')
assert fractions.index.is_unique and set(metadata.index)<=set(fractions.index)
fractions = fractions.loc[metadata.index]
broad = pd.read_csv(fraction_path.parent/'bayesprism_stage1_broad14_final_fraction.tsv',sep='\t').set_index('sample_id').loc[metadata.index]
epithelial = pd.read_csv(fraction_path.parent/'bayesprism_stage2_epithelial7_whole_biopsy_fraction.tsv',sep='\t').set_index('sample_id').loc[metadata.index]
non_epi = broad.columns.drop('Epithelial')
np.testing.assert_allclose(fractions[non_epi],broad[non_epi],atol=1e-12)
np.testing.assert_allclose(fractions[epithelial.columns],epithelial,atol=1e-12)
np.testing.assert_allclose(fractions[epithelial.columns].sum(axis=1),broad.Epithelial,atol=1e-12)
assert np.isfinite(fractions.to_numpy()).all() and np.allclose(fractions.sum(axis=1),1,atol=1e-5)
fractions.to_csv(source/'Fig5A_BayesPrism_whole_biopsy_fractions20.tsv',sep='\t')
means = fractions.groupby(metadata.DominantEcotype).mean().T.reindex(columns=['E'+str(i) for i in range(1,8)])
means.to_csv(source/'Fig5A_fraction_mean_all20.tsv',sep='\t')
means12 = means.loc[states.CellType.unique()]
means12.to_csv(source/'Fig5A_fraction_mean_CORE12.tsv',sep='\t')
# 原文注明averaged and scaled而未给出尺度公式；显示采用每细胞类型0-1缩放，原始均值另存
display12 = means12.sub(means12.min(axis=1),axis=0).div(means12.max(axis=1)-means12.min(axis=1),axis=0)
display12.to_csv(source/'Fig5A_fraction_display_CORE12.tsv',sep='\t')
inputs.append(fraction_path)

# Step5. Fig5A复用既有原文式连续/离散检验；肠段按原文离散特征方式做one-versus-rest
continuous = pd.read_csv(prior/'06_ecotype_features/ecotype_continuous_features.tsv',sep='\t')
discrete = pd.read_csv(prior/'06_ecotype_features/ecotype_discrete_features.tsv',sep='\t')
region_rows=[]
for region in ['Ileum','Colon','Rectum']:
    for ecotype in ['E'+str(i) for i in range(1,8)]:
        a=metadata.loc[metadata.Region3.eq(region),ecotype]
        b=metadata.loc[~metadata.Region3.eq(region),ecotype]
        test=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
        region_rows.append({'Ecotype':ecotype,'Feature':'Region3_'+region,'GroupA':region,'GroupB':'Other_regions','N_A':len(a),'N_B':len(b),'MeanDifference':a.mean()-b.mean(),'P':test.pvalue,'Unit':'Biopsy'})
all_features = pd.concat([continuous[['Ecotype','Feature','MeanDifference','P']],discrete[['Ecotype','Feature','MeanDifference','P']],pd.DataFrame(region_rows)[['Ecotype','Feature','MeanDifference','P']]],ignore_index=True)
# 原文要求跨本次所有feature检验BH；新增肠段后不能继续沿用分开的旧Q
all_features['Q']=stats.false_discovery_control(all_features.P.to_numpy(),method='bh')
all_features['BHFamily']='55 features x 7 ecotypes = 385 comparisons'
all_features.to_csv(source/'Fig5A_all_385_feature_tests.tsv',sep='\t',index=False)
pd.DataFrame(region_rows).to_csv(source/'Fig5A_region_discrete_tests.tsv',sep='\t',index=False)
inputs.extend(prior/'06_ecotype_features'/x for x in ['ecotype_continuous_features.tsv','ecotype_discrete_features.tsv','hallmark_sample_scores.tsv'])

# Step6. 冻结输入、原文图注方法和展示合同留存，便于结束时检查原始输入未改变
for p in inputs:
    manifest.append({'Path':str(p),'Bytes':p.stat().st_size,'SHA256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(manifest).to_csv(source/'frozen_inputs_SHA256.tsv',sep='\t',index=False)
r=ET.parse(prior/'00_audit/PMC8526411_bioc.xml')
paragraphs=[]
for p in r.findall('.//passage'):
    t=p.findtext('text','')
    info={i.attrib.get('key'):i.text for i in p.findall('infon')}
    if (info.get('type')=='fig_caption' and info.get('id') in ['F2','F4','F5']) or 'To create the UMAP projections' in t or 'For the network schematic depicted in Figure 4B' in t or 'For the analysis shown in Figure 5A' in t:
        paragraphs.append(t)
(out/'reports/paper_method_excerpts.txt').write_text('\n\n'.join(paragraphs))
contract={'Backend':'Python rendering; R native umap and igraph numerical analysis','Claim':'Display frozen state heterogeneity, community membership and phenotype associations with the original paper panel semantics','Archetype':'Quantitative grids and aligned multi-block overview','FontFloorPt':5,'Export':['PDF','SVG','PNG'],'Input':{'States':53,'Ecotypes':7,'Biopsies':2025,'Fig4A_Assigned':len(assigned)},'Boundaries':['BayesPrism sensitivity, not CIBERSORTx','No survival or mutation feature panels','No NMF or ecotype refitting','No inferred cell communication'],'Fig5A_BH':385,'ScaleDisclosure':'Fig5A cell-type mean min-max scaling is a display choice; exact original formula not specified'}
(out/'reports/figure_contract.json').write_text(json.dumps(contract,indent=2))
print('INPUTS_READY',flush=True)
