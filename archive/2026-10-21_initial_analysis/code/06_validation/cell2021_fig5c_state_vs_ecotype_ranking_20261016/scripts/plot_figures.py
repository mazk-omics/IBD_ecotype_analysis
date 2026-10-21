"""Python-only scientific figures; inputs are the independently recalculated tables."""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42})
plt.rcParams.update({'font.size':7,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'legend.frameon':False})
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
qa_source=Path(__file__).resolve().parent/'figure_qa'
if not qa_source.exists():qa_source=Path.home()/'.codex/skills/nature-figure/scripts'
sys.path.insert(0,str(qa_source))
from audit_panel_alignment import require_matplotlib_panel_alignment

OUT=Path(__file__).resolve().parents[1]
FIG=OUT/'figures'
EXT=pd.read_csv(OUT/'statistics/external_all60_rankings.tsv',sep='\t')
DISC=pd.read_csv(OUT/'statistics/discovery_one_biopsy_all60_rankings.tsv',sep='\t')
BOOT=pd.read_csv(OUT/'statistics/discovery_E1_member_effect_bootstrap.tsv',sep='\t')
SENS=pd.read_csv(OUT/'statistics/discovery_selection_rank_sensitivity_summary.tsv',sep='\t')
COLORS={'Ecotype':'#B45E26','State':'#7B8794','E1':'#B45E26','Member':'#346B86'}
SHORT={'ETR_W10_clinical_remission':'ETR W10','IFX_W4_6_endoscopic_healing':'IFX W4–6','VDZ_W6_endoscopic_healing':'VDZ W6','RISK_3yr_B2_B3':'RISK 3y*','VDZ_W12_endoscopic_healing':'VDZ W12','VDZ_W52_endoscopic_healing':'VDZ W52','MSCCR_Inflammation_one_biopsy':'Discovery inflammation'}
PRIMARY=list(SHORT)[:4]
qa=[]

def name(f):
    return f.replace('Monocyte_Macrophage','Mono/Mac').replace('Enteroendocrine','Enteroendo').replace('_S',' (S')+')' if '_S' in f else f

def export(fig,base,axes=None,labels=None,exclude=None):
    fig.canvas.draw()
    kwargs={'axes':axes,'panel_ids':labels} if axes is not None else {}
    if exclude:kwargs['exclude_axes']=exclude
    require_matplotlib_panel_alignment(fig,json_out=str(FIG/(base+'.alignment.json')),overlay_svg=str(FIG/(base+'.alignment.svg')),strict=True,**kwargs)
    fig.savefig(FIG/(base+'.svg'),facecolor='white')
    fig.savefig(FIG/(base+'.pdf'),facecolor='white')
    fig.savefig(FIG/(base+'.png'),dpi=600,facecolor='white')
    plt.close(fig)

def heatmap(table,endpoints,base,title,order,ordering_note,footer):
    # The same feature occupies the same row in every endpoint; no pooled meta-Z.
    rows=len(order)
    fig=plt.figure(figsize=(7.2,9.2))
    ax=fig.add_axes([.12,.16,.32,.64])
    z=table.pivot(index='Feature',columns='Endpoint',values='SignedZ_favorable').loc[order,endpoints]
    image=ax.imshow(z.to_numpy(),aspect='auto',cmap='viridis',vmin=-3,vmax=3,interpolation='nearest')
    ax.set_xticks(np.arange(len(endpoints)))
    ax.set_xticklabels([SHORT[e] for e in endpoints],rotation=60,ha='left',rotation_mode='anchor',fontsize=7)
    ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=22)
    for j,e in enumerate(endpoints):
        r=table[table.Endpoint.eq(e)].iloc[0]
        ax.text(j,-1.5,f'{int(r.N)}/{int(r.Events)}',ha='center',va='bottom',fontsize=6.5,clip_on=False)
    ax.set_yticks([])
    ax.set_ylabel('Candidate features (n = 60)\n'+ordering_note,fontsize=8,labelpad=24)
    ax.set_xticks(np.arange(-.5,len(endpoints),1),minor=True)
    ax.set_yticks(np.arange(-.5,rows,1),minor=True)
    ax.grid(which='minor',color='white',lw=.18)
    ax.tick_params(which='minor',length=0)
    for spine in ax.spines.values():spine.set_visible(True)
    # Nominal P<.05 is a small central mark; BH60 significance is a ring.
    tt=table.set_index(['Feature','Endpoint'])
    for i,f in enumerate(order):
        for j,e in enumerate(endpoints):
            r=tt.loc[(f,e)]
            if r.P_two_sided<.05:ax.plot(j,i,'o',ms=1.9,color='white',mec='none')
            if r.Q_BH60<.05:ax.plot(j,i,'o',ms=4.0,mfc='none',mec='white',mew=.7)
    # Key feature labels use the reference's leader-line convention.
    selected=set(table.loc[table.Endpoint.isin(endpoints)&table.DisplayOrder_absZ.le(3),'Feature'])
    selected.update(['E'+str(i) for i in range(1,8)])
    selected.update(table.loc[table.E1_member,'Feature'])
    inds=[i for i,f in enumerate(order) if f in selected]
    positions=np.array(inds,dtype=float)
    gap=1.65
    for i in range(1,len(positions)):positions[i]=max(positions[i],positions[i-1]+gap)
    if positions[-1]>59:
        positions-=positions[-1]-59
        for i in range(len(positions)-2,-1,-1):positions[i]=min(positions[i],positions[i+1]-gap)
    for i,yy in zip(inds,positions):
        f=order[i];eco=f in ['E'+str(k) for k in range(1,8)]
        ax.plot([len(endpoints)-.43,len(endpoints)-.1,len(endpoints)+.12],[i,yy,yy],color='#444444',lw=.45,clip_on=False)
        ax.text(len(endpoints)+.21,yy,name(f)+('  [Ecotype]' if eco else ''),va='center',fontsize=6.4,color=COLORS['Ecotype'] if eco else '#222222',fontweight='bold' if eco else 'normal',clip_on=False)
    ax.set_xlim(-.5,len(endpoints)-.5)
    ax.set_ylim(rows-.5,-.5)
    fig.text(.12,.95,title,fontsize=13,ha='left',va='top')
    fig.text(.12,.895,'Baseline scores | two-sided Wilcoxon | N/events below endpoint',fontsize=8)
    cax=fig.add_axes([.79,.79,.025,.10])
    cb=fig.colorbar(ScalarMappable(norm=Normalize(-3,3),cmap='viridis'),cax=cax,ticks=[-3,0,3])
    cb.ax.tick_params(labelsize=6,length=2)
    fig.text(.76,.917,'Favorable',fontsize=7)
    fig.text(.76,.76,'Unfavorable',fontsize=7)
    fig.text(.75,.716,'Signed Z\nclipped at ±3',fontsize=6.5,linespacing=1.4)
    fig.text(.12,.10,'White dot: nominal P < 0.05   Ring: BH60 q < 0.05\nPositive Z: higher scores associate with favorable outcome; RISK sign reversed.\nRows are shared across columns; full 60-feature labels and ranks are in the TSV.',fontsize=6.5,linespacing=1.6)
    fig.text(.12,.027,footer,fontsize=6.4,va='bottom',linespacing=1.5)
    export(fig,base,axes=[ax],labels=['a'],exclude=[cax])
    z.reset_index().to_csv(FIG/(base+'_source.tsv'),sep='\t',index=False)
    qa.append(dict(Figure=base,Panel='a',Claim='Feature-specific directions across distinct endpoints',Unit='One baseline biopsy per patient',Uncertainty='Z association evidence; no pooled estimator',Limit=footer))

order=EXT[EXT.Endpoint.eq(PRIMARY[0])].sort_values('DisplayOrder_absZ').Feature.tolist()
heatmap(EXT,PRIMARY,'Fig1_reference_style_external_heatmap','Clinical outcome association: Ecotypes and states',order,'ordered by ETR W10 |Z|',
        '*RISK: deposited expression scale unresolved; age/sex-adjusted E2 BH7 q=0.0613.\nSmall event groups and unadjusted comparisons; association ranking only.')
heatmap(EXT,list(SHORT)[4:6],'FigS1_later_VDZ_heatmap','Later VDZ outcomes: supplementary comparisons',order,'same ETR row order as main figure',
        'W12: only 3 events. W52 overlaps W6/W12 patients; time points are not independent.\nNo combined Liptak score; supplementary associations only.')

# Rank distributions use all 7 and all 53 features, with no feature-independence test.
fig,axes=plt.subplots(1,2,figsize=(7.2,4.6),gridspec_kw={'left':.20,'right':.98,'bottom':.24,'top':.80,'wspace':.38})
ys=np.arange(8)
for i,e in enumerate(PRIMARY):
    for j,family in enumerate(['Ecotype','State']):
        g=EXT[EXT.Endpoint.eq(e)&EXT.Family.eq(family)]
        y=2*i+j
        offsets=np.linspace(-.16,.16,len(g))
        for ax,value in [(axes[0],g.Rank_absZ.to_numpy()),(axes[1],g.RankBiserial.abs().to_numpy())]:
            ax.scatter(value,y+offsets,s=9,color=COLORS[family],alpha=.7,linewidths=0)
            ax.plot(np.median(value),y,'|',color='black',ms=9,mew=1.2)
axes[0].set_yticks(ys)
axes[0].set_yticklabels([SHORT[e]+'  '+f for e in PRIMARY for f in ['E (7)','S (53)']],fontsize=6.6)
axes[1].set_yticks(ys);axes[1].set_yticklabels([])
for i,ax in enumerate(axes):
    ax.invert_yaxis();ax.set_ylim(7.5,-.6)
    for y in [1.5,3.5,5.5]:ax.axhline(y,color='#DDDDDD',lw=.5)
    ax.text(-.10,1.05,chr(97+i),transform=ax.transAxes,weight='bold',fontsize=8)
axes[0].set_xlim(0,61);axes[0].set_xlabel('Association rank among 60 (smaller is stronger)')
axes[1].set_xlim(0,1);axes[1].set_xlabel('Absolute rank-biserial effect')
fig.text(.08,.95,'Ecotype ranks show no consistent external advantage',fontsize=12)
fig.text(.08,.89,'All features shown; black tick marks the within-family median',fontsize=7.5)
fig.text(.08,.055,'N/events: ETR 70/12; IFX 23/8; VDZ W6 27/6; RISK 210/27.\nFeature distributions are descriptive: features are correlated, and family sizes differ.\n*RISK scale unresolved; rank ties are retained. Unadjusted baseline associations.',fontsize=6.5,linespacing=1.5)
export(fig,'Fig2_rank_and_effect_distributions',axes=axes.tolist(),labels=['a','b'])
qa.extend([dict(Figure='Fig2_rank_and_effect_distributions',Panel=p,Claim=c,Unit='Feature; 7 Ecotypes / 53 states',Uncertainty='Full distribution, median; no between-feature test',Limit='Correlated features; heterogeneous endpoints') for p,c in [('a','Relative ranking'),('b','Magnitude of score-group separation')]])

# Discovery: point estimates and common-subject bootstrap; E1 is shown alongside every member.
dd=DISC.set_index('Feature').loc[BOOT.Feature].reset_index().merge(BOOT,on='Feature',suffixes=('','_boot'))
dd=dd.sort_values('Rank_absZ').reset_index(drop=True)
n=int(DISC.N.iloc[0]);events=int(DISC.Events.iloc[0])
fig,axes=plt.subplots(1,2,figsize=(7.2,4.8),gridspec_kw={'left':.27,'right':.97,'bottom':.33,'top':.78,'wspace':.45})
y=np.arange(len(dd))
for i,r in dd.iterrows():
    color=COLORS['E1'] if r.Feature=='E1' else COLORS['Member']
    axes[0].errorbar(r.RankBiserial,i,xerr=[[r.RankBiserial-r.RB_CI_low],[r.RB_CI_high-r.RankBiserial]],fmt='o',color=color,ms=4,lw=1)
    rr=SENS[SENS.Feature.eq(r.Feature)].iloc[0]
    axes[1].errorbar(rr.RankMedian,i,xerr=[[rr.RankMedian-rr.RankP10],[rr.RankP90-rr.RankMedian]],fmt='o',color=color,ms=4,lw=1.4)
axes[0].set_yticks(y);axes[0].set_yticklabels([name(r.Feature)+'  #'+f'{r.Rank_absZ:g}' for r in dd.itertuples()],fontsize=6.7)
axes[1].set_yticks(y);axes[1].set_yticklabels([])
for i,ax in enumerate(axes):
    ax.set_ylim(len(dd)-.4,-.6)
    ax.text(-.12,1.06,chr(97+i),transform=ax.transAxes,weight='bold',fontsize=8)
axes[0].set_xlim(.45,.75);axes[0].set_xlabel('Rank-biserial effect\n95% bootstrap interval')
axes[1].set_xlim(0,15);axes[1].set_xticks([1,5,10,15]);axes[1].set_xlabel('Rank in 100 biopsy selections\nMedian and 10th–90th percentiles')
fig.text(.08,.95,'E1 and its members carry closely related inflammation signals',fontsize=11)
fig.text(.08,.89,f'Discovery only | one biopsy/patient | N={n}; inflamed={events}; non-inflamed={n-events}',fontsize=7.5)
fig.text(.08,.07,'a: Fixed seed 20261016; one biopsy per patient selected without outcome/score information.\nPaired-feature, outcome-stratified patient bootstrap (2,000 draws; descriptive intervals).\nb: 100 new one-biopsy selections; same 919 patients and all 60 features each time.\nInternal Discovery association; no external validation or incremental prediction claim.',fontsize=6.5,linespacing=1.5)
export(fig,'Fig3_discovery_E1_and_members',axes=axes.tolist(),labels=['a','b'])
dd.to_csv(FIG/'Fig3_discovery_E1_and_members_source.tsv',sep='\t',index=False)
qa.extend([dict(Figure='Fig3_discovery_E1_and_members',Panel=p,Claim=c,Unit='919 patients; one biopsy/patient',Uncertainty=u,Limit='Discovery internal association') for p,c,u in [('a','E1 versus all ten member effects','Paired 2000 stratified patient bootstrap percentile 95% CI'),('b','Biopsy selection rank sensitivity','100 selections; median and P10-P90')]])

# A complete, labelled ranked strip for each endpoint prevents selective callout display.
def ranked_strip(table,endpoint,base):
    d=table[table.Endpoint.eq(endpoint)].sort_values('DisplayOrder_absZ')
    fig=plt.figure(figsize=(5.1,9.2));ax=fig.add_axes([.13,.14,.18,.68])
    ax.imshow(d.SignedZ_event.to_numpy()[:,None],aspect='auto',cmap='viridis',vmin=-3,vmax=3,interpolation='nearest')
    ax.set_xticks([]);ax.set_yticks([])
    for i,r in enumerate(d.itertuples()):
        ax.text(1.0,i,f'{r.Rank_absZ:g}  '+name(r.Feature),fontsize=6.1,va='center',color=COLORS[r.Family],fontweight='bold' if r.Family=='Ecotype' else 'normal',clip_on=False)
    fig.text(.10,.95,SHORT[endpoint]+': complete association ranking',fontsize=11)
    first=d.iloc[0]
    fig.text(.10,.89,f'N={int(first.N)}; event=1: {int(first.Events)}; event=0: {int(first.Nonevents)}\nRank: descending |Z|; equal ranks retained',fontsize=7,linespacing=1.5)
    fig.text(.10,.035,'Color: signed Z for event=1 (yellow=higher; purple=lower), clipped at ±3.\nEcotypes in brown; states in gray. Full-precision statistics in companion TSV.\n'+('Conditional deposited scale; age/sex-adjusted E2 BH7 q=0.0613.' if endpoint.startswith('RISK') else 'Discovery internal association.' if endpoint.startswith('MSCCR') else 'Small event groups; unadjusted baseline association.'),fontsize=6.2,linespacing=1.5)
    export(fig,base,axes=[ax],labels=['a'])
for k,e in enumerate(PRIMARY+list(SHORT)[4:6]):ranked_strip(EXT,e,'FigS2_'+str(k+1)+'_'+e)
ranked_strip(DISC,'MSCCR_Inflammation_one_biopsy','FigS3_discovery_all60_ranked')
pd.DataFrame(qa).to_csv(OUT/'audit/figure_panel_contract.tsv',sep='\t',index=False)
print('Rendered',len(list(FIG.glob('*.pdf'))),'PDF/SVG/PNG figure triplets')
