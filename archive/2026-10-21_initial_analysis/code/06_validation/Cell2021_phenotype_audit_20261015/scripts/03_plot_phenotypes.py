# 用已完成的统计表制作16:9汇报图；固定状态ID，灰格表示不可评价。
# 本机msccr_clinical_audit_20261015环境运行，输出PNG/PDF/SVG与布局审计。
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
sys.path.insert(0,str(Path.home()/'.codex/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment

# Step1. 固定证据链与字号；主图颜色为Z，星号仅表示本项目经验P的BH q<0.05
out = Path(__file__).resolve().parents[1]
(out/'figures').mkdir(exist_ok=True)
members = pd.read_csv(out/'audit/frozen_states.tsv',sep='\t').sort_values(['CellType','State'])
discovery = pd.read_csv(out/'statistics/discovery_state_phenotype_enrichment.tsv',sep='\t')
concordance = pd.read_csv(out/'statistics/discovery_external_concordance_all_states.tsv',sep='\t')
distribution = pd.read_csv(out/'statistics/ecotype_phenotype_distribution.tsv',sep='\t')
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':12,'axes.titlesize':14,'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':10.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
parents = members.CellType.unique()
groups = [members[members.CellType.isin(parents[:6])],members[members.CellType.isin(parents[6:])]]
max_rows = max(len(g) for g in groups)
contracts = []

# Step2. 两种53状态主表型热图；不改变色阶追求显著效果
for mode,stem,title in [('BiopsyNominal','P1_53states_biopsy_phenotype','Frozen cell states: biopsy-level phenotype enrichment'),('OneBiopsyPerPatientRestricted','P2_53states_patient_restricted','Frozen cell states: patient and stratum sensitivity')]:
    data = discovery[discovery.Scope.eq('All')&discovery.Mode.eq(mode)].copy()
    data['Column'] = data.Feature+':'+data.Target
    columns = ['Disease:CD','Disease:UC','Inflammation:Inflamed','Inflammation:Non-inflamed','Region3:Colon','Region3:Ileum','Region3:Rectum']
    z = data.pivot(index='CellState',columns='Column',values='Z').reindex(columns=columns)
    q = data.pivot(index='CellState',columns='Column',values='EmpiricalQ').reindex(columns=columns)
    fig,axes = plt.subplots(1,2,figsize=(16,9),gridspec_kw={'left':.19,'right':.96,'bottom':.21,'top':.86,'wspace':.78})
    cmap = mpl.colormaps['RdBu_r'].copy()
    cmap.set_bad('#d9d9d9')
    for panel,(ax,g) in enumerate(zip(axes,groups)):
        values = np.full((max_rows,7),np.nan)
        values[:len(g)] = z.reindex(g.CellState).to_numpy()
        image = ax.imshow(values,cmap=cmap,vmin=-5,vmax=5,aspect='auto',interpolation='none')
        ax.set_xticks(range(7),['CD','UC','Infl.','Non-I','Colon','Ileum','Rectum'],rotation=35,ha='right',rotation_mode='anchor')
        ax.set_yticks(range(len(g)),g.CellState+' ['+g.Ecotype+']')
        ax.tick_params(length=0,pad=5)
        ax.set_ylim(max_rows-.5,-.5)
        ax.text(-.02,1.01,chr(97+panel),transform=ax.transAxes,fontweight='bold',fontsize=16)
        for row,state in enumerate(g.CellState):
            for col in range(7):
                if q.loc[state,columns[col]]<.05:
                    ax.text(col,row,'*',ha='center',va='center',fontsize=13,color='white' if abs(values[row,col])>2.8 else 'black')
                elif not np.isfinite(values[row,col]):
                    ax.text(col,row,'x',ha='center',va='center',fontsize=9,color='#555555')
        ends = np.flatnonzero(g.CellType.to_numpy()[:-1]!=g.CellType.to_numpy()[1:])+.5
        for boundary in ends:
            ax.axhline(boundary,color='white',linewidth=1.5)
    cax = fig.add_axes([.42,.075,.23,.018])
    fig.colorbar(image,cax=cax,orientation='horizontal',label='Permutation Z (display clipped at +/-5)',ticks=[-5,0,5])
    fig.suptitle(title,fontsize=19,y=.97)
    fig.text(.5,.91,'53 states | 12 parent types | 1,000 permutations | * empirical BH q < 0.05 | x not evaluable',ha='center',fontsize=12)
    require_matplotlib_panel_alignment(fig,json_out=str(out/'figures'/f'{stem}.alignment.json'),exclude_axes=[cax],strict=True)
    for ext in ['.png','.pdf','.svg']:
        fig.savefig(out/'figures'/f'{stem}{ext}',dpi=300)
    plt.close(fig)
    contracts.append(dict(Figure=stem,Claim=title,Archetype='quantitative grid',ReplicateUnit=mode,Interval='None: permutation standardized count; full null mean and SD supplied',Source='discovery_state_phenotype_enrichment.tsv',Stars='Empirical BH q<0.05 across all prespecified tests in this mode'))

# Step3. 同一状态的跨队列方向与可评价性；按真实可比表型分别显示
columns = []
values = []
nominal = []
for dataset,target in [('Colon_GSE214695','CD'),('Ileum_IBDverse','Inflamed')]:
    for rule,label in [('BiopsyNominal vs CellNominal','Nominal'),('OneBiopsyPerPatient vs DonorEqualFraction','Patient/donor')]:
        d = concordance[concordance.Dataset_External.eq(dataset)&concordance.Target.eq(target)&concordance.Rule.eq(rule)].set_index('CellState')
        for origin in ['Discovery','External']:
            columns.append(('Colon CD' if dataset.startswith('Colon') else 'Ileum CD infl.')+'\n'+label+'\n'+origin[:4])
            values.append(d['Z_'+origin].reindex(members.CellState).to_numpy())
            nominal.append(d['PaperNominal_'+origin].reindex(members.CellState).fillna(False).to_numpy())
matrix = pd.DataFrame(np.array(values).T,index=members.CellState,columns=range(8))
flags = pd.DataFrame(np.array(nominal).T,index=members.CellState,columns=range(8))
fig,axes = plt.subplots(1,2,figsize=(16,9),gridspec_kw={'left':.19,'right':.96,'bottom':.24,'top':.85,'wspace':.78})
cmap = mpl.colormaps['RdBu_r'].copy()
cmap.set_bad('#d9d9d9')
for panel,(ax,g) in enumerate(zip(axes,groups)):
    data = np.full((max_rows,8),np.nan)
    data[:len(g)] = matrix.reindex(g.CellState).to_numpy()
    image = ax.imshow(data,cmap=cmap,vmin=-4,vmax=4,aspect='auto',interpolation='none')
    ax.set_yticks(range(len(g)),g.CellState+' ['+g.Ecotype+']')
    ax.set_xticks(range(8),['C-D','C-E','C-P','C-S','I-D','I-E','I-P','I-S'],fontsize=10)
    ax.tick_params(length=0,pad=5)
    ax.text(-.02,1.01,chr(97+panel),transform=ax.transAxes,fontweight='bold',fontsize=16)
    for row,state in enumerate(g.CellState):
        for col in range(8):
            if flags.loc[state,col]:
                ax.text(col,row,'o',ha='center',va='center',fontsize=8,color='white' if abs(data[row,col])>2.5 else 'black')
            elif not np.isfinite(data[row,col]):
                ax.text(col,row,'x',ha='center',va='center',fontsize=9,color='#555555')
    ax.axvline(3.5,color='white',linewidth=3)
cax = fig.add_axes([.42,.075,.23,.018])
fig.colorbar(image,cax=cax,orientation='horizontal',label='Phenotype Z; o nominal Z > 1.96; x not evaluable',ticks=[-4,0,4])
fig.suptitle('Phenotype concordance is limited in external single-cell cohorts',fontsize=19,y=.97)
fig.text(.5,.91,'Colon: inflamed CD vs UC, 12 donors | Ileum: CD inflamed vs non-inflamed, 20 donors',ha='center',fontsize=12)
fig.text(.5,.16,'C colon, I ileum | D discovery biopsies, E external cells | P one biopsy/patient, S equal donor fractions',ha='center',fontsize=12)
stem = 'P3_discovery_external_phenotype_concordance'
require_matplotlib_panel_alignment(fig,json_out=str(out/'figures'/f'{stem}.alignment.json'),exclude_axes=[cax],strict=True)
for ext in ['.png','.pdf','.svg']:
    fig.savefig(out/'figures'/f'{stem}{ext}',dpi=300)
plt.close(fig)
contracts.append(dict(Figure=stem,Claim='External nominal enrichment does not establish donor-level global concordance',Archetype='quantitative grid',ReplicateUnit='Biopsy vs cell; independent patient vs donor sensitivity',Interval='None: standardized permutation statistics',Source='discovery_external_concordance_all_states.tsv',Stars='o paper nominal Z>1.96, not BH significance'))

# Step4. 7生态型加Unassigned，展示各表型内部的官方显著分配比例
fig,ax = plt.subplots(figsize=(16,9))
fig.subplots_adjust(left=.09,right=.84,bottom=.15,top=.84)
ecotypes = ['E1','E2','E3','E4','E5','E6','E7','Unassigned']
official_colors = pd.read_csv(out.parents[1]/'06_validation/paper_figures_20261012/source_data/official_annotation_colors.tsv',sep='\t')
official_colors = official_colors[official_colors.Kind.eq('Ecotype')].set_index('Label').Color
colors = [official_colors[e] for e in ecotypes[:-1]]+['#d9d9d9']
labels = []
x = []
base = np.zeros(7)
target_order = [('Disease','CD'),('Disease','UC'),('Inflammation','Inflamed'),('Inflammation','Non-inflamed'),('Region3','Colon'),('Region3','Ileum'),('Region3','Rectum')]
positions = [0,1,2.6,3.6,5.2,6.2,7.2]
for feature,target in target_order:
    g = distribution[distribution.Feature.eq(feature)&distribution.Target.eq(target)]
    labels.append(target+'\nn='+str(g.PhenotypeBiopsies.iloc[0])+'; p='+str(g.PhenotypePatients.iloc[0]))
for ecotype,color in zip(ecotypes,colors):
    height = np.array([distribution.loc[distribution.Feature.eq(f)&distribution.Target.eq(t)&distribution.Assignment.eq(ecotype),'PercentWithinPhenotype'].iloc[0] for f,t in target_order])
    ax.bar(positions,height,bottom=base,width=.68,color=color,label=ecotype,edgecolor='white',linewidth=.5)
    base += height
np.testing.assert_allclose(base,100)
ax.set_xticks(positions,labels,fontsize=12)
ax.set_ylabel('Frozen official assignment within phenotype (%)')
ax.set_ylim(0,100)
ax.legend(bbox_to_anchor=(1.02,1),loc='upper left',fontsize=13,frameon=False)
ax.set_title('Disease                     Inflammation                       Intestinal region',pad=17,fontsize=15)
fig.suptitle('Seven ecotypes: phenotype distribution with unassigned biopsies retained',fontsize=19,y=.96)
fig.text(.5,.04,'2,025 biopsies | 919 patients | 1,712 significant assignments | 313 unassigned | n biopsies; p patients',ha='center',fontsize=13)
stem = 'P4_ecotype_official_assignment_phenotype_distribution'
require_matplotlib_panel_alignment(fig,json_out=str(out/'figures'/f'{stem}.alignment.json'),strict=True)
for ext in ['.png','.pdf','.svg']:
    fig.savefig(out/'figures'/f'{stem}{ext}',dpi=300)
plt.close(fig)
contracts.append(dict(Figure=stem,Claim='Ecotype distribution depends on phenotype and assignment coverage',Archetype='single quantitative panel',ReplicateUnit='Biopsy; patients can appear in several groups',Interval='None: descriptive census proportions, no independence claim',Source='ecotype_phenotype_distribution.tsv',Stars='None; no significance test'))
(out/'figures/figure_contracts.json').write_text(json.dumps(contracts,indent=2))
print('PLOT_COMPLETE',len(contracts),flush=True)
