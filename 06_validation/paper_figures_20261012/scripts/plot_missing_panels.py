# 仅补原文Fig2B、4B、5A对应缺失图；Fig2A和4A直接复用官方原件，不重新绘制。
# 服务器ibd_refaudit环境运行；R计算UMAP/FR，本脚本仅用Python绘图，输出PDF/SVG/PNG及版面审计。
# 点为活检、节点为cell state、边为共现；BayesPrism敏感性结果不等同于官方CIBERSORTx结果。
from pathlib import Path
import json
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).parent/'figure_qa'))
from audit_panel_alignment import require_matplotlib_panel_alignment

# Step1. 固定原模型的编号，所有面板使用一致颜色；字号下限5pt，输出文字可编辑
out=Path('/home/mazekai/IBD_EcoTyper/06_validation/paper_figures_20261012')
src=out/'source_data'
figdir=out/'figures'
states=pd.read_csv(src/'frozen_states.tsv',sep='\t')
cts=['B','CD4_T','CD8_T','ILC','Monocyte_Macrophage','Mast','Fibroblast','Pericyte','Endothelial','Goblet','Enteroendocrine','Glial']
ctlabels={'B':'B cells','CD4_T':'CD4 T cells','CD8_T':'CD8 T cells','ILC':'ILCs','Monocyte_Macrophage':'Monocytes / macrophages','Mast':'Mast cells','Fibroblast':'Fibroblasts','Pericyte':'Pericytes','Endothelial':'Endothelial cells','Goblet':'Goblet cells','Enteroendocrine':'Enteroendocrine cells','Glial':'Glial cells'}
ecos=['E'+str(i) for i in range(1,8)]
eco_colors=dict(zip(ecos,['#d73027','#f5d94e','#71aa57','#d889ba','#9270a9','#999999','#a6cadd']))
state_colors=dict(zip(['S%02d'%i for i in range(1,8)],['#ce397c','#3c99ac','#62aa51','#8055a0','#ef8b3c','#d3ad31','#a86644']))
ct_colors=dict(zip(cts,['#d3a020','#dc4b92','#7465ab','#9bc871','#3199a8','#9bc3d8','#b0a1c7','#e777af','#d56c20','#eaaa7c','#219779','#c4bcae']))
meta_colors={'Disease':{'CD':'#d73027','UC':'#f5d94e'},'Inflammation':{'Inflamed':'#d73027','Non-inflamed':'#f5d94e'},'Region3':{'Ileum':'#f5d94e','Colon':'#d73027','Rectum':'#71aa57'}}
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':7,'axes.titlesize':7,'axes.labelsize':6,'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':0.5,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42,'svg.fonttype':'none'})
figures=[]

colors=pd.read_csv(src/'official_annotation_colors.tsv',sep='\t')
state_colors=colors[colors.Kind.eq('State')].set_index('Label').Color.to_dict()
ct_colors=colors[colors.Kind.eq('CellType')].set_index('Label').Color.to_dict()
eco_colors=colors[colors.Kind.eq('Ecotype')].set_index('Label').Color.to_dict()

# Step3. Fig2B输入是H*距离，不是表达矩阵；逐父类型全部通过QC的活检，坐标由论文版本R包计算
fig,axes=plt.subplots(3,4,figsize=(7.2,6.9))
fig.subplots_adjust(left=0.06,right=0.97,top=0.90,bottom=0.13,wspace=0.32,hspace=0.43)
for ax,ct in zip(axes.ravel(),cts):
    coords=pd.read_csv(src/(ct+'_Fig2B_UMAP.tsv'),sep='\t')
    ax.scatter(coords.UMAP1,coords.UMAP2,c=coords.State.map(state_colors),s=2.0,linewidths=0,alpha=0.75,rasterized=True)
    ax.set_title(ctlabels[ct]+' | n='+str(len(coords)),fontsize=6.8,pad=5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('UMAP 1',labelpad=2)
    ax.set_ylabel('UMAP 2',labelpad=2)
fig.suptitle('Cell state heterogeneity | UMAP of H* Euclidean distances',fontsize=9,y=0.96)
handles=[mpl.lines.Line2D([],[],marker='o',linestyle='none',color=c,markersize=4,label=s) for s,c in state_colors.items()]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(0.5,0.045),ncol=7,fontsize=6,handletextpad=0.4,columnspacing=1.3)
fig.text(0.06,0.017,'Fig. 2B analogue. R umap 0.2.0.0 defaults; seed 20261012. Each point is an assigned biopsy.',fontsize=5.8)
require_matplotlib_panel_alignment(fig,axes=list(axes.ravel()),json_out=figdir/'Fig2B.alignment.json',strict=True)
figures.append(('Fig2B_cell_state_UMAP',fig))

# Step5. Fig4B节点=状态，颜色=父类型；边宽严格正比于冻结Jaccard，不增加通讯边
edges=pd.read_csv(src/'Fig4B_all_within_ecotype_edges.tsv',sep='\t')
fig,axes=plt.subplots(2,4,figsize=(7.2,5.2))
fig.subplots_adjust(left=0.03,right=0.97,top=0.85,bottom=0.23,wspace=0.16,hspace=0.40)
for ax,ecotype in zip(axes.ravel(),ecos):
    nodes=pd.read_csv(src/(ecotype+'_Fig4B_FR_nodes.tsv'),sep='\t').set_index('CellState')
    xy=nodes[['X','Y']].to_numpy()
    # 仅统一显示坐标范围；不会改变边权或网络成员
    xy=(xy-xy.min(axis=0))/(xy.max(axis=0)-xy.min(axis=0)).clip(min=1e-12)
    nodes[['PlotX','PlotY']]=xy
    nodes.to_csv(src/(ecotype+'_Fig4B_display_nodes.tsv'),sep='\t')
    own=edges[edges.Ecotype.eq(ecotype)&edges.Jaccard.gt(0)]
    for e in own.itertuples():
        a=nodes.loc[e.From]
        b=nodes.loc[e.To]
        ax.plot([a.PlotX,b.PlotX],[a.PlotY,b.PlotY],color='#aaaaaa',linewidth=1.3*e.Jaccard,zorder=1,rasterized=True)
    for n in nodes.itertuples():
        ax.scatter(n.PlotX,n.PlotY,s=82,color=ct_colors[n.CellType],edgecolors='white',linewidths=0.45,zorder=2)
        ax.text(n.PlotX,n.PlotY,n.State[1:].lstrip('0'),fontsize=5.3,ha='center',va='center',zorder=3)
    ax.set_xlim(-0.17,1.17)
    ax.set_ylim(-0.17,1.17)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.set_title(ecotype+' | '+str(len(nodes))+' states',fontsize=8,color='black',fontweight='bold',pad=7)
axes.ravel()[-1].set_axis_off()
handles=[mpl.lines.Line2D([],[],marker='o',linestyle='none',color=ct_colors[c],markersize=4,label=ctlabels[c]) for c in cts]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(0.5,0.04),ncol=4,fontsize=5.8,columnspacing=1.2,handletextpad=0.4)
fig.suptitle('Ecotype composition | frozen cell state co-occurrence networks',fontsize=9,y=0.955)
fig.text(0.04,0.018,'Fig. 4B analogue. Fruchterman-Reingold; edge width = 1.3 x Jaccard. Node number = S##; no communication inference.',fontsize=5.6)
require_matplotlib_panel_alignment(fig,axes=list(axes.ravel()),json_out=figdir/'Fig4B.alignment.json',strict=True)
figures.append(('Fig4B_ecotype_member_networks',fig))

# Step6. Fig5A按原文组合均值比例和特征富集；无生存/突变数据则不建立相关面板
features=pd.read_csv(src/'Fig5A_all_385_feature_tests.tsv',sep='\t')
feature_order=['Disease','Inflammation','Region3_Ileum','Region3_Colon','Region3_Rectum']+sorted(features.loc[features.Feature.str.startswith('HALLMARK_'),'Feature'].unique())
effect=features.pivot(index='Feature',columns='Ecotype',values='MeanDifference').reindex(index=feature_order,columns=ecos)
q=features.pivot(index='Feature',columns='Ecotype',values='Q').reindex(index=feature_order,columns=ecos)
means=pd.read_csv(src/'Fig5A_fraction_display_CORE12.tsv',sep='\t',index_col=0).reindex(index=cts,columns=ecos)
counts=pd.read_csv(src/'all_2025_biopsies.tsv',sep='\t').DominantEcotype.value_counts().reindex(ecos)
fig=plt.figure(figsize=(7.2,10.9))
left,width=0.45,0.42
ax_header=fig.add_axes([left,0.900,width,0.026])
ax_header.imshow(np.arange(7)[None,:],aspect='auto',cmap=ListedColormap(list(eco_colors.values())),vmin=-0.5,vmax=6.5)
for i,e in enumerate(ecos): ax_header.text(i,0,e,ha='center',va='center',fontsize=7)
ax_header.set_axis_off()
fig.text(left,0.935,'Ecotype',fontsize=7)
for i,e in enumerate(ecos): fig.text(left+width*(i+.5)/7,0.889,'n='+str(counts[e]),ha='center',fontsize=5.8)
ax_fraction=fig.add_axes([left,0.704,width,0.165])
image=ax_fraction.imshow(means.to_numpy(),aspect='auto',cmap=LinearSegmentedColormap.from_list('Fraction',['#f4f0fa','#9eb2c7','#265f54']),vmin=0,vmax=1,interpolation='none')
ax_fraction.set_yticks(np.arange(12),labels=[ctlabels[c] for c in cts],fontsize=6.0)
ax_fraction.set_xticks([])
ax_fraction.tick_params(length=0,pad=6)
fig.text(left,0.875,'Relative cell type abundance (averaged and scaled)',fontsize=6.4)
cbax=fig.add_axes([0.902,0.73,0.012,0.11])
cb=fig.colorbar(image,cax=cbax,ticks=[0,1])
cb.set_label('Row-scaled mean',fontsize=6)
ax_features=fig.add_axes([left,0.102,width,0.570])
for i in range(len(effect)):
    for j in range(7):
        if q.iloc[i,j]<0.05:
            ax_features.scatter(j,i,s=180*abs(effect.iloc[i,j]),color='black' if effect.iloc[i,j]>0 else '#b3b3b3',linewidths=0,zorder=2)
ax_features.set_xlim(-.5,6.5)
ax_features.set_ylim(len(effect)-.5,-.5)
ax_features.set_xticks([])
labels=[{'Disease':'CD vs UC','Inflammation':'Inflamed vs non-inflamed','Region3_Ileum':'Ileum vs other regions','Region3_Colon':'Colon vs other regions','Region3_Rectum':'Rectum vs other regions'}.get(f,f.replace('HALLMARK_','').replace('_',' ').capitalize()) for f in feature_order]
ax_features.set_yticks(np.arange(len(labels)),labels=labels,fontsize=5.8)
ax_features.tick_params(length=0,pad=6)
ax_features.axhline(4.5,color='#888888',linewidth=0.5)
for old,new in {'Dna':'DNA','Il2':'IL2','Il6':'IL6','Pi3k':'PI3K','Mtorc1':'mTORC1','Myc':'MYC','Kras':'KRAS','Wnt':'WNT','Nfkb':'NF-kB','Stat3':'STAT3','Stat5':'STAT5','Tgf beta':'TGF-beta','Tnfa':'TNF-alpha','Uv':'UV'}.items():
    labels=[x.replace(old,new) for x in labels]
ax_features.set_yticklabels(labels,fontsize=5.8)
for y in [0,1,2,3,4]: ax_features.axhspan(y-.5,y+.5,color='#f5f5f5',zorder=0)
for spine in ax_features.spines.values(): spine.set_visible(True)
fig.text(left,0.680,'Clinical and biological features | BH q < 0.05',fontsize=6.4)
fig.suptitle('Clinical and biological features of IBD ecotypes',fontsize=9,y=0.975)
handles=[mpl.lines.Line2D([],[],marker='o',linestyle='none',color='black',markersize=5,label='Enriched (q < 0.05)'),mpl.lines.Line2D([],[],marker='o',linestyle='none',color='#b3b3b3',markersize=5,label='Depleted (q < 0.05)')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(0.62,0.052),ncol=1,fontsize=6)
fig.text(0.035,0.025,'Fig. 5A analogue. BayesPrism sensitivity; all 2,025 CD/UC biopsies. Circle area proportional to absolute mean difference.\nTwo-sided Wilcoxon; BH across 385 tests. Repeated biopsies are dependent. No survival or mutation data analyzed.',fontsize=5.7)
require_matplotlib_panel_alignment(fig,axes=[ax_header,ax_fraction,ax_features],panel_ids=['ecotypes','fractions','features'],row_groups=[],column_groups=[['ecotypes','fractions','features']],json_out=figdir/'Fig5A.alignment.json',strict=True)
figures.append(('Fig5A_ecotype_feature_overview',fig))

# Step7. 同版多格式导出，保留文字矢量；大矩阵及散点栅格化仅降低体积，不减少观测
for name,fig in figures:
    fig.savefig(figdir/(name+'.pdf'),dpi=600)
    fig.savefig(figdir/(name+'.svg'),dpi=600)
    fig.savefig(figdir/(name+'.png'),dpi=350)
    plt.close(fig)
    print('EXPORTED',name,flush=True)
(out/'reports/plot_environment.json').write_text(json.dumps({'Python':sys.version,'Matplotlib':mpl.__version__,'Pandas':pd.__version__,'Numpy':np.__version__},indent=2))
print('THREE_MISSING_PANELS_RENDERED',flush=True)
