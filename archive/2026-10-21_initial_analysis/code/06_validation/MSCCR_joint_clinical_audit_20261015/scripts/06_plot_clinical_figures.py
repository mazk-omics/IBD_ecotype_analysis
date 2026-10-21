# 用 Python 绘制 MSCCR 临床关联图；输入均为本次核验统计表和真实 BayesPrism 比例。
# conda 环境 msccr_clinical_audit_20261015；每图保存 PNG/PDF/SVG、源数据和对齐记录。
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0,str(Path('C:/Users/马泽凯/.codex/skills/nature-figure/scripts')))
from audit_panel_alignment import require_matplotlib_panel_alignment

# Step1. 固定一致尺度和显示标签，所有临床对比保留完整来源
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = root/'06_validation/MSCCR_joint_clinical_audit_20261015'
figdir = out/'figures'
source = out/'source_data'
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':9,'axes.titlesize':11,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.7,'legend.frameon':False,'savefig.dpi':300})
cmap = LinearSegmentedColormap.from_list('ClinicalEffect',['#2E5F87','#8BB1C5','#FAFAF8','#DDA392','#A63E3A'])
members = pd.read_csv(source/'frozen_state_membership.tsv',sep='\t')
metadata = pd.read_csv(source/'analysis_metadata_2025.tsv',sep='\t',index_col=0)
adj = pd.read_csv(out/'statistics/patient_cluster_adjusted_contrasts.tsv',sep='\t')
adj = adj.loc[adj.Weighting.eq('Biopsy_weighted')].copy()
spec = [('Inflammation','Inflamed vs Non-inflamed','Inflamed vs non-inflamed'),('Disease','CD vs UC','CD vs UC'),('Region3','Ileum vs Colon','Ileum vs colon'),('Region3','Rectum vs Colon','Rectum vs colon'),('Age','per 1 SD Age','Age'),('Sex','Male vs Female','Male vs female'),('CRP_log2','per 1 SD CRP_log2','CRP (log2)'),('Fecal_calprotectin_log2','per 1 SD Fecal_calprotectin_log2','Fecal calprotectin (log2)'),('Local_GHAS','per 1 SD Local_GHAS','Local GHAS [CD]'),('Local_Nancy','per 1 SD Local_Nancy','Local Nancy [UC]'),('SES_CD','per 1 SD SES_CD','SES-CD [CD]'),('Mayo_endo','per 1 SD Mayo_endo','Mayo endoscopy [UC]'),('HBI','per 1 SD HBI','HBI [CD]'),('SCCAI','per 1 SD SCCAI','SCCAI [UC]'),('Max_GHAS','per 1 SD Max_GHAS','Maximum GHAS [CD]'),('Max_Nancy','per 1 SD Max_Nancy','Maximum Nancy [UC]'),('Clinical_activity','Active vs Inactive','Clinical active vs inactive'),('Endoscopic_severity','Mild vs Inactive','Mild vs inactive endoscopy'),('Endoscopic_severity','Moderate vs Inactive','Moderate vs inactive endoscopy'),('Endoscopic_severity','Severe vs Inactive','Severe vs inactive endoscopy'),('Endoscopic_remission','1.0 vs 0.0','Endoscopic remission'),('Derived_histology_zero','1.0 vs 0.0','Derived zero histology')]
labels = [s[2] for s in spec]
adj['Display'] = adj.apply(lambda r: next(s[2] for s in spec if s[0]==r.Feature and s[1]==r.Contrast),axis=1)
scale = np.ceil(adj.Beta.abs().max()*10)/10
eco_order = ['E'+str(i) for i in range(1,8)]

# Step2. 53 状态和 7 ecotype 的完整患者调整热图，不省略非显著结果
for entity,prefix,height in [('state','Fig1_state_clinical_adjusted',14.5),('ecotype','Fig2_ecotype_clinical_adjusted',6.8)]:
    frame = adj.loc[adj.Entity.eq(entity)]
    order = members.sort_values(['Ecotype','CellType','State']).ID.tolist() if entity=='state' else eco_order
    beta = frame.pivot(index='Marker',columns='Display',values='Beta').reindex(index=order,columns=labels)
    q = frame.pivot(index='Marker',columns='Display',values='Q').reindex(index=order,columns=labels)
    beta.to_csv(source/(prefix+'_beta.tsv'),sep='\t')
    q.to_csv(source/(prefix+'_q.tsv'),sep='\t')
    frame.to_csv(source/(prefix+'_full_statistics.tsv'),sep='\t',index=False)
    fig = plt.figure(figsize=(12.8,height))
    ax = fig.add_axes([0.225,0.36 if entity=='ecotype' else 0.15,0.68,0.44 if entity=='ecotype' else 0.77])
    image = ax.imshow(beta.values,aspect='auto',cmap=cmap,vmin=-scale,vmax=scale,interpolation='none')
    ax.set_xticks(np.arange(len(labels)),labels=labels,rotation=90,ha='right',rotation_mode='anchor')
    yl = [m+'  | '+members.set_index('ID').loc[m,'Ecotype'] for m in order] if entity=='state' else order
    ax.set_yticks(np.arange(len(order)),labels=yl,fontsize=7.5 if entity=='state' else 10)
    ax.tick_params(length=0,pad=5)
    yy,xx = np.where(q.values<0.05)
    ax.scatter(xx,yy,s=4 if entity=='state' else 10,c='#212121',marker='o',linewidths=0)
    ax.set_xticks(np.arange(-0.5,len(labels),1),minor=True)
    ax.set_yticks(np.arange(-0.5,len(order),1),minor=True)
    ax.grid(which='minor',color='#FFFFFF',linewidth=0.35)
    ax.tick_params(which='minor',length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cax = fig.add_axes([0.935,0.61 if entity=='state' else 0.48,0.012,0.19])
    cb = fig.colorbar(image,cax=cax)
    cb.set_label('Adjusted beta (abundance SD)',fontsize=8)
    title = '53 cell states' if entity=='state' else '7 ecotypes'
    fig.text(0.225,0.96,title+' | associations with measured clinical features',fontsize=13,weight='bold')
    fig.text(0.225,0.936 if entity=='state' else 0.90,'Patient-clustered SE; covariates: disease, region, age and sex (excluding the exposure)',fontsize=8.5)
    fig.text(0.225,0.022,'Dot: BH q < 0.05 across '+str(len(frame))+' contrasts. Continuous exposure: per 1 SD; binary exposure: displayed contrast.\nAll 2,025 biopsies retained when clinical values are available; 313 officially Unassigned biopsies included.',fontsize=8)
    require_matplotlib_panel_alignment(fig,axes=[ax],panel_ids=['heatmap'],json_out=figdir/(prefix+'.alignment.json'),strict=True)
    fig.savefig(figdir/(prefix+'.png'),dpi=300)
    fig.savefig(figdir/(prefix+'.pdf'))
    fig.savefig(figdir/(prefix+'.svg'))
    plt.close(fig)

# Step3. Cell 2021 Figure 5A 的 IBD 适配：覆盖、真实父细胞组成、全部连续临床特征
fractions = pd.read_csv(root/'06_validation/paper_figures_20261012/source_data/Fig5A_BayesPrism_whole_biopsy_fractions20.tsv',sep='\t',index_col=0)
assert fractions.index.is_unique and set(metadata.index)==set(fractions.index)
fractions = fractions.reindex(metadata.index)
assert np.allclose(fractions.sum(axis=1),1,atol=1e-5)
parents = members.CellType.drop_duplicates().tolist()
means = fractions.groupby(metadata.DominantEcotype).mean().T.reindex(index=parents,columns=eco_order)
display = means.sub(means.min(axis=1),axis=0).div(means.max(axis=1)-means.min(axis=1),axis=0)
means.to_csv(source/'Fig3_CORE12_whole_biopsy_mean_fraction.tsv',sep='\t')
display.to_csv(source/'Fig3_CORE12_row_minmax_display.tsv',sep='\t')
dom = pd.read_csv(out/'statistics/dominant_ecotype_continuous_features.tsv',sep='\t')
features = ['Age','CRP_log2','Fecal_calprotectin_log2','Local_GHAS','Local_Nancy','SES_CD','Mayo_endo','HBI','SCCAI','Max_GHAS','Max_Nancy']
dom_beta = dom.pivot(index='Feature',columns='Ecotype',values='RankBiserial').reindex(index=features,columns=eco_order)
dom_q = dom.pivot(index='Feature',columns='Ecotype',values='Q').reindex(index=features,columns=eco_order)
dom.to_csv(source/'Fig3_dominant_CE_clinical_statistics.tsv',sep='\t',index=False)
coverage = pd.crosstab(metadata.DominantEcotype,metadata.AssignedEcotype.ne('Unassigned')).reindex(eco_order)
coverage.to_csv(source/'Fig3_assignment_coverage.tsv',sep='\t')
fig = plt.figure(figsize=(9,10.8))
ax_a = fig.add_axes([0.29,0.78,0.57,0.13])
ax_b = fig.add_axes([0.29,0.445,0.57,0.255])
ax_c = fig.add_axes([0.29,0.11,0.57,0.255])
x = np.arange(7)
ax_a.bar(x,coverage[True],color='#547B8C',width=0.7,label='Assigned')
ax_a.bar(x,coverage[False],bottom=coverage[True],color='#CDD3D4',width=0.7,label='Unassigned')
ax_a.set_xlim(-0.5,6.5)
ax_a.set_xticks(x,labels=eco_order)
ax_a.set_ylabel('Biopsies')
ax_a.set_ylim(0,coverage.sum(axis=1).max()*1.18)
ax_a.legend(loc='upper left',bbox_to_anchor=(1.01,1.0),ncol=1,fontsize=8)
ax_a.set_title('a   Highest-abundance ecotype and assignment coverage',loc='left',pad=10,fontsize=10)
image_b = ax_b.imshow(display.values,aspect='auto',cmap='Blues',vmin=0,vmax=1,interpolation='none')
ax_b.set_yticks(np.arange(len(parents)),labels=[p.replace('_',' ') for p in parents],fontsize=8)
ax_b.set_xticks(x,labels=eco_order)
ax_b.tick_params(length=0,pad=5)
ax_b.set_title('b   Whole-biopsy parental cell proportions',loc='left',pad=16,fontsize=10)
cb = fig.colorbar(image_b,cax=fig.add_axes([0.895,0.48,0.014,0.17]))
cb.set_label('Row min-max scaled mean',fontsize=8)
image_c = ax_c.imshow(dom_beta.values,aspect='auto',cmap=cmap,vmin=-1,vmax=1,interpolation='none')
feature_labels = [next(s[2] for s in spec if s[0]==f) for f in features]
ax_c.set_yticks(np.arange(len(features)),labels=feature_labels,fontsize=8)
ax_c.set_xticks(x,labels=eco_order)
ax_c.tick_params(length=0,pad=5)
ax_c.set_title('c   Continuous clinical features | dominant CE versus rest',loc='left',pad=16,fontsize=10)
yy,xx = np.where(dom_q.values<.05)
ax_c.scatter(xx,yy,s=8,c='#202020',linewidths=0)
cb = fig.colorbar(image_c,cax=fig.add_axes([0.895,0.14,0.014,0.17]))
cb.set_label('Rank-biserial effect',fontsize=8)
for axis in [ax_b,ax_c]:
    for spine in axis.spines.values():
        spine.set_visible(False)
fig.text(0.29,0.974,'MSCCR ecotype clinical characterization',fontsize=13,weight='bold')
fig.text(0.29,0.953,'Cell 2021 Figure 5A adaptation | 2,025 biopsies, 919 patients',fontsize=9)
fig.text(0.09,0.032,'Panel b: CORE12 selected from real whole-biopsy fractions; no renormalization across the 12 types.\nPanel c: two-sided rank-sum, BH across 77 tests; dot q < 0.05. Biopsy-level descriptive analysis.\nPatient-adjusted evidence is shown separately in Figures 1, 2 and 4.',fontsize=8)
require_matplotlib_panel_alignment(fig,axes=[ax_a,ax_b,ax_c],panel_ids=['a','b','c'],row_groups=[],column_groups=[['a','b','c']],json_out=figdir/'Fig3_Cell2021_Fig5A_adaptation.alignment.json',strict=True)
for ext in ['png','pdf','svg']:
    fig.savefig(figdir/('Fig3_Cell2021_Fig5A_adaptation.'+ext),dpi=300)
plt.close(fig)

# Step4. 代表性 E1 效应与 95% 区间，调整前后使用相同完整病例和标准化
effects = pd.read_csv(source/'representative_effects_before_after_adjustment.tsv',sep='\t')
order = effects.Feature.drop_duplicates().tolist()
models = ['Unadjusted_cluster_SE','Adjusted','Adjusted_equal_patient_weight']
colors = ['#939B9D','#A63E3A','#2E5F87']
names = ['Unadjusted, patient-clustered SE','Adjusted, patient-clustered SE','Adjusted, equal patient weights']
fig,ax = plt.subplots(figsize=(9.7,6.7))
fig.subplots_adjust(left=0.34,right=0.96,bottom=0.24,top=0.85)
for i,(model,color,name) in enumerate(zip(models,colors,names)):
    frame = effects.loc[effects.Model.eq(model)].set_index('Feature').reindex(order)
    y = np.arange(len(order))+(i-1)*0.2
    ax.errorbar(frame.Beta,y,xerr=np.vstack([frame.Beta-frame.CI_low,frame.CI_high-frame.Beta]),fmt='o',color=color,markersize=4,capsize=2,linewidth=1.1,label=name)
yl = []
for f in order:
    row = effects.loc[effects.Feature.eq(f)].iloc[0]
    label = next(s[2] for s in spec if s[0]==f)
    yl.append(label+'\n'+str(int(row.N))+' biopsies / '+str(int(row.Patients))+' patients')
ax.set_yticks(np.arange(len(order)),labels=yl,fontsize=9)
ax.invert_yaxis()
ax.axvline(0,color='#5C6264',linewidth=0.7,linestyle='--')
ax.set_xlabel('E1 abundance difference (SD); continuous clinical exposure per 1 SD')
ax.legend(loc='upper center',bbox_to_anchor=(0.43,-0.15),fontsize=9,ncol=1)
fig.text(0.08,0.95,'E1 associations persist after patient-aware adjustment',fontsize=13,weight='bold')
fig.text(0.08,0.90,'Point estimates and 95% t intervals; disease-specific scores analyzed within CD or UC',fontsize=9)
require_matplotlib_panel_alignment(fig,axes=[ax],panel_ids=['forest'],json_out=figdir/'Fig4_E1_effects_adjustment.alignment.json',strict=True)
for ext in ['png','pdf','svg']:
    fig.savefig(figdir/('Fig4_E1_effects_adjustment.'+ext),dpi=300)
plt.close(fig)

# Step5. 严格同患者同详细部位的全部 30 对活检，展示三种主要炎症方向
abundance = pd.read_csv(source/'ecotype_abundance_2025x7.tsv',sep='\t',index_col=0)
merged = metadata[['PatientID','RegionDetailed','Inflammation']].join(abundance)
medians = merged.groupby(['PatientID','RegionDetailed','Inflammation'])[abundance.columns].median()
inf = medians.xs('Inflamed',level='Inflammation')
non = medians.xs('Non-inflamed',level='Inflammation')
both = inf.index.intersection(non.index)
inf = inf.loc[both].groupby('PatientID').mean()
non = non.loc[both].groupby('PatientID').mean()
pair_source = pd.concat({'Non-inflamed':non,'Inflamed':inf},names=['Inflammation','PatientID']).reset_index()
pair_source.to_csv(source/'Fig5_all_exact_region_patient_pairs.tsv',sep='\t',index=False)
paired = pd.read_csv(out/'statistics/patient_detailed_region_matched_inflammation.tsv',sep='\t').set_index('Marker')
fig,axes = plt.subplots(1,3,figsize=(10.6,4.5),sharey=True)
fig.subplots_adjust(left=0.08,right=0.96,bottom=0.23,top=0.75,wspace=0.26)
for ax,marker in zip(axes,['E1','E5','E6']):
    for patient in inf.index:
        ax.plot([0,1],[non.loc[patient,marker],inf.loc[patient,marker]],color='#AAB2B5',alpha=0.5,linewidth=0.7,zorder=1)
    ax.scatter(np.zeros(len(non)),non[marker],s=15,color='#2E5F87',zorder=2)
    ax.scatter(np.ones(len(inf)),inf[marker],s=15,color='#A63E3A',zorder=2)
    ax.set_xticks([0,1],labels=['Non-inflamed','Inflamed'],fontsize=9)
    ax.set_xlim(-0.2,1.2)
    ax.set_ylim(-0.035,1.035)
    row = paired.loc[marker]
    ax.set_title(marker+' | paired q = '+format(row.Q,'.3g'),fontsize=10,pad=14)
axes[0].set_ylabel('Ecotype abundance')
fig.text(0.08,0.94,'Within-patient inflammation differences at the same detailed site',fontsize=12,weight='bold')
fig.text(0.08,0.87,'30 patients / 30 exact-site pairs; all eligible patients displayed',fontsize=9)
fig.text(0.08,0.06,'Patient-site median abundance; gray lines connect each patient. Paired two-sided Wilcoxon; BH across all 7 ecotypes.\nThese are current-biopsy associations; sampling dates and future treatment response are unavailable.',fontsize=8)
require_matplotlib_panel_alignment(fig,axes=list(axes),panel_ids=['E1','E5','E6'],row_groups=[['E1','E5','E6']],column_groups=[],json_out=figdir/'Fig5_exact_region_paired.alignment.json',strict=True)
for ext in ['png','pdf','svg']:
    fig.savefig(figdir/('Fig5_exact_region_paired.'+ext),dpi=300)
plt.close(fig)
print('FIGURES_COMPLETE',len(list(figdir.glob('*.png'))),flush=True)
