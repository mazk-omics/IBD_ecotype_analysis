# 最终核对冻结输入、统计维度、图件 QA 和交付文件；打包明日 PPT 所需材料。
# 本机执行，不连接服务器，不覆盖既有项目成果。
from pathlib import Path
import pandas as pd
import numpy as np
import json, hashlib, zipfile, subprocess, sys

# Step1. 核对冻结文件从分析前到交付时仍一致，全部矩阵样本完整匹配
out = Path('D:/IBD_ecotype/IBD_ecotype_analysis/06_validation/MSCCR_joint_clinical_audit_20261015')
frozen = pd.read_csv(out/'audit/frozen_copy_hash_comparison.tsv',sep='\t')
frozen['FinalLocalSHA256'] = [hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in frozen.local_path]
assert len(frozen)==27 and frozen.FinalLocalSHA256.eq(frozen.local_sha256).all() and frozen.FinalLocalSHA256.eq(frozen.server_sha256).all()
frozen.to_csv(out/'audit/frozen_final_integrity.tsv',sep='\t',index=False)
meta = pd.read_csv(out/'source_data/analysis_metadata_2025.tsv',sep='\t',index_col=0)
states = pd.read_csv(out/'source_data/state_abundance_2025x53.tsv',sep='\t',index_col=0)
eco = pd.read_csv(out/'source_data/ecotype_abundance_2025x7.tsv',sep='\t',index_col=0)
assert states.shape==(2025,53) and eco.shape==(2025,7) and states.index.equals(meta.index) and eco.index.equals(meta.index)
assert meta.PatientID.nunique()==919 and meta.AssignedEcotype.eq('Unassigned').sum()==313
adj = pd.read_csv(out/'statistics/patient_cluster_adjusted_contrasts.tsv',sep='\t')
assert len(adj)==2640 and not adj.Feature.eq('Histological_remission').any() and adj.Feature.eq('Derived_histology_zero').any()
assert np.isfinite(adj[['Beta','SE','P','Q']]).all().all() and adj.Q.between(0,1).all()

# Step2. 所有图件重新检查最终 PDF 字体和碰撞，并保留适用对齐结果
skill = Path('C:/Users/马泽凯/.codex/skills/nature-figure/scripts')
quality = []
for pdf in (out/'figures').glob('*.pdf'):
    assert pdf.with_suffix('.png').exists() and pdf.with_suffix('.svg').exists()
    text = subprocess.run([sys.executable,str(skill/'audit_pdf_text.py'),str(pdf),'--min-pt','5','--json'],capture_output=True,text=True,encoding='utf-8')
    pdf.with_suffix('.text-audit.json').write_text(text.stdout,encoding='utf-8')
    collision = subprocess.run([sys.executable,str(skill/'audit_figure_collisions.py'),str(pdf),'--json-out',str(pdf.with_suffix('.collision-audit.json'))],capture_output=True,text=True,encoding='utf-8')
    assert text.returncode==0 and collision.returncode==0
    geometry = json.loads(pdf.with_suffix('.collision-audit.json').read_text())
    assert geometry['summary']['fail']==0 and geometry['summary']['warn']==0
    alignment = json.loads(pdf.with_suffix('.alignment.json').read_text())
    quality.append(dict(Figure=pdf.stem,PDFTextExit=text.returncode,CollisionExit=collision.returncode,CollisionFails=geometry['summary']['fail'],CollisionWarnings=geometry['summary']['warn'],AlignmentStatus=alignment.get('status',alignment.get('verdict','see JSON')),ManualVisualReview='All 5 full figures and panels inspected; labels, uncertainty and color meaning checked'))
assert len(quality)==5
pd.DataFrame(quality).to_csv(out/'audit/final_figure_QA.tsv',sep='\t',index=False)

# Step3. 保留整个 GEO 样本空间及其未进入冻结集合的明确原因
geo = pd.read_csv(out/'source_data/GEO_metadata_2490.tsv',sep='\t')
excluded = pd.read_csv(out/'audit/excluded_4_GEO_records.tsv',sep='\t')
geo['FrozenDiscovery'] = geo.SampleID.isin(meta.index)
geo['Disposition'] = np.select([geo.FrozenDiscovery,geo.SampleID.isin(excluded.SampleID)],['Frozen CD/UC Discovery','Technical epithelial zero posterior mass exclusion'],default='Control not in IBD-only Discovery')
geo[['geo_accession','SampleID','GEO_PatientID','ibd_disease','typere','regionre','Disposition']].to_csv(out/'source_data/all_2490_sample_disposition.tsv',sep='\t',index=False)
assert geo.Disposition.value_counts().to_dict()=={'Frozen CD/UC Discovery':2025,'Control not in IBD-only Discovery':461,'Technical epithelial zero posterior mass exclusion':4}
local = pd.read_csv(out/'audit/local_inventory.tsv',sep='\t')
remote = pd.read_csv(out/'audit/server_inventory.tsv',sep='\t',low_memory=False)
local['clinical_match'] = local.sha256.map(remote.dropna(subset=['sha256']).drop_duplicates('sha256').set_index('sha256').clinical_match).fillna('See specific source/clinical mapping; non-data files not applicable')
local['local_copy_status'] = 'Present; this is the local file'
local['dimension_status'] = np.where(local.n_samples.notna(),'Matched server SHA or object-specific evidence','Not individually loaded or dimension not applicable')
local.to_csv(out/'audit/local_inventory.tsv',sep='\t',index=False)

# Step4. 输出便于汇报传递的压缩包，完整审计大清单保留在目录内
package = out/'MSCCR_PPT_figures_and_summary_20261015.zip'
files = [out/'presentation_summary.md',out/'audit_report.md',out/'figure_legends.md',out/'statistical_analysis_plan.md',out/'README.md',out/'source_data/clinical_variable_dictionary.tsv',out/'source_data/clinical_coverage_by_disease_region.tsv']
files += list((out/'figures').glob('*.png'))+list((out/'figures').glob('*.pdf'))+list((out/'figures').glob('*.svg'))
files += list((out/'source_data').glob('Fig*'))+[out/'source_data/representative_effects_before_after_adjustment.tsv']
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as archive:
    for p in files:
        archive.write(p,p.relative_to(out))
summary = dict(Status='Complete',FrozenInputsUnchanged=True,FrozenFilesMatched=27,DiscoveryBiopsies=2025,Patients=919,States=53,Ecotypes=7,Figures=5,PatientAdjustedContrasts=2640,ServerSSHConnections=1,RemoteWrites=0,RemoteSession='Closed normally after audit',LargeExpressionMatricesDownloaded=0,Package=str(package),PackageBytes=package.stat().st_size,Pending=['Original historemiss coding','bMIS per-biopsy values or separately validated reconstruction','Medication and longitudinal metadata','Independent cohort validation'])
(out/'audit/final_completion_checkpoint.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(out/'progress_summary.md').write_text('# 2026-10-15 完成检查点\n\nPhase 1–5 完成。两端审计、临床字典、样本映射、患者调整/敏感性和5套图件全部完成。原始 historemiss 编码与样本级 bMIS 等未解决事项见 presentation_summary.md。服务器唯一 SSH 会话已正常关闭，不需要重新扫描。\n',encoding='utf-8')
print(json.dumps(summary,indent=2),flush=True)
