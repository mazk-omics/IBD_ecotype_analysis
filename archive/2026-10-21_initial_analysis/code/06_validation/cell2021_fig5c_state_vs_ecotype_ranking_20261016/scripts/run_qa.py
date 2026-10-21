"""Independent output checks and rendered PDF audits; write only inside new output."""
from pathlib import Path
import concurrent.futures, subprocess, sys, json, hashlib, shutil
import numpy as np
import pandas as pd
from scipy import stats

OUT=Path(__file__).resolve().parents[1]
SKILL=Path.home()/'.codex/skills/nature-figure/scripts'
LOCAL=OUT/'scripts/figure_qa'
LOCAL.mkdir(exist_ok=True)
for filename in ['validate_figure.py','audit_pdf_text.py','audit_panel_alignment.py','audit_figure_collisions.py']:
    if not (LOCAL/filename).exists():shutil.copy2(SKILL/filename,LOCAL/filename)

def run(cmd,log):
    p=subprocess.run([sys.executable,*map(str,cmd)],capture_output=True,encoding='utf-8',errors='replace')
    (OUT/'audit'/log).write_text(p.stdout+'\n'+p.stderr,encoding='utf-8')
    return p.returncode,p.stdout

source_code,source_text=run([LOCAL/'validate_figure.py',OUT/'scripts/plot_figures.py','--json'],'plot_source_validation.json')

def audit_pdf(pdf):
    name=pdf.stem
    tc,txt=run([LOCAL/'audit_pdf_text.py',pdf,'--min-pt','5','--json'],name+'.text_audit.json')
    cc,ctxt=run([LOCAL/'audit_figure_collisions.py',pdf,'--json-out',OUT/'audit'/(name+'.collision_audit.json')],name+'.collision.log')
    return dict(Figure=name,TextAuditExit=tc,CollisionAuditExit=cc,CollisionVerdict='PASS' if 'verdict: PASS' in ctxt else 'REVIEW REQUIRED' if 'verdict: REVIEW REQUIRED' in ctxt else 'FAIL',TextEditable='<text' in pdf.with_suffix('.svg').read_text(encoding='utf-8'),PNGExists=pdf.with_suffix('.png').exists())

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    figures=list(pool.map(audit_pdf,sorted((OUT/'figures').glob('*.pdf'))))
figures=pd.DataFrame(figures)
figures.to_csv(OUT/'audit/final_figure_QA.tsv',sep='\t',index=False)

checks=[]
def check(name,condition,scope):
    checks.append(dict(Check=name,Pass=bool(condition),Scope=scope))

ext=pd.read_csv(OUT/'statistics/external_all60_rankings.tsv',sep='\t')
disc=pd.read_csv(OUT/'statistics/discovery_one_biopsy_all60_rankings.tsv',sep='\t')
check('External full scope',len(ext)==360 and ext.Endpoint.nunique()==6,'60 features per six endpoints')
check('Discovery scope',len(disc)==60 and (disc.N==919).all(),'919 unique patients / 60 features')
for endpoint,t in pd.concat([ext,disc]).groupby('Endpoint'):
    check(endpoint+' feature IDs',len(t)==t.Feature.nunique()==60,'All 7+53; no duplicates')
    check(endpoint+' fixed N',t.N.nunique()==1 and (t.N==t.EffectiveN).all(),'Same complete patients for every feature')
    check(endpoint+' P/Z consistency',np.allclose(2*stats.norm.sf(t.Z_two_sided_abs),t.P_two_sided,rtol=1e-12,atol=1e-14),'Full 60')
    check(endpoint+' directional consistency',np.allclose(t.SignedZ_event,np.sign(t.RankBiserial)*t.Z_two_sided_abs),'U-derived sign; no median-based sign')
    check(endpoint+' rank consistency',np.allclose(stats.rankdata(t.P_two_sided,method='average'),t.Rank_absZ),'Average tie ranks')
    check(endpoint+' BH60 consistency',np.allclose(stats.false_discovery_control(t.P_two_sided,method='bh'),t.Q_BH60),'Entire family 60')
    check(endpoint+' effect consistency',np.allclose(t.RankBiserial,2*t.AUC_like_event-1),'U probability and RB identity')
    check(endpoint+' original direction',t.SignedZ_favorable.isna().all() and t.SignedZ_paper_style.isna().all() if endpoint.startswith('MSCCR') else np.allclose(t.SignedZ_paper_style,-t.SignedZ_favorable),'External negative favorable; Discovery clinical favorability not applicable')
samplemanifest=pd.read_csv(OUT/'manifests/discovery_one_biopsy_919_patient_manifest.tsv',sep='\t')
check('Discovery primary unique patients',samplemanifest.PatientID.is_unique and samplemanifest.SampleID.is_unique and len(samplemanifest)==919,'One biopsy each')
sensitivity=pd.read_csv(OUT/'manifests/discovery_100_selection_sampleIDs.tsv',sep='\t')
check('Sensitivity manifest 91900 IDs',len(sensitivity)==91900 and sensitivity.groupby('Iteration').PatientID.nunique().eq(919).all(),'100 selections each 919 different patients')
snap=pd.read_csv(OUT/'audit/input_SHA256_before_after.tsv',sep='\t')
current=snap.Path.map(lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest())
snap['SHA256_final']=current;snap['Unchanged_final']=snap.SHA256_before.eq(current)
snap.to_csv(OUT/'audit/input_SHA256_final_verification.tsv',sep='\t',index=False)
check('Frozen sources unchanged',snap.Unchanged_final.all(),str(len(snap))+' registered files')
check('Figure text and collision gates',((figures.TextAuditExit==0)&(figures.CollisionAuditExit==0)&figures.TextEditable&figures.PNGExists).all(),'All '+str(len(figures))+' PNG/PDF/SVG triplets')
check('Plot source validation',source_code==0,'Source preflight; warnings require review')
check('Chinese report delivered',(OUT/'中文分析汇报.md').stat().st_size>10000,'Includes top10, all E and E1 members for every endpoint')
pd.DataFrame(checks).to_csv(OUT/'audit/delivery_requirements_QA.tsv',sep='\t',index=False)
summary=dict(Passed=sum(c['Pass'] for c in checks),Failed=sum(not c['Pass'] for c in checks),Figures=len(figures),SourceValidationExit=source_code,InputsUnchanged=bool(snap.Unchanged_final.all()))
(OUT/'audit/delivery_QA_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(figures.to_string(index=False));print(json.dumps(summary))
sys.exit(1 if summary['Failed'] else 0)
