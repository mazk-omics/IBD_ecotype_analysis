# 保存原排名副本、逐状态状态清单和实际进展，复核未提交状态。
# 所有写入仅在独立目录；不改动benchmark或既有结果。
from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib
import json
import shutil
import subprocess
import pandas as pd

# Step1. 保存历史全部排名与LOOCV，不改原件
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
base = root / '06_validation/paper_aligned_recovery_characterization_20261012'
for name in ['marker_n1_n7.tsv', 'marker_composite_ranking.tsv', 'marker_LOOCV_per_dataset.tsv', 'marker_LOOCV_per_state.tsv', 'n4_paired_LOOCV_comparison.tsv', 'DE_stratum_eligibility_parallel4.tsv']:
    shutil.copy2(base / '02_external_markers' / name, out / 'audit' / ('FROZEN_REUSE_' + name))
master = pd.read_csv(out / 'CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t').fillna('')
state_status = master[['CellState', 'CellType', 'Ecotype', 'CandidateAnnotation', 'Confidence', 'BiologicalRationale', 'Uncertainty']].copy()
state_status['ReviewStatus'] = state_status.Confidence.map({'High confidence': 'PASS_CANDIDATE_PROGRAM', 'Moderate confidence': 'PASS_PROVISIONAL_PROGRAM', 'Low confidence': 'PARTIAL_WEAK_SUPPORT', 'Unresolved': 'BLOCKED_IDENTITY_OR_FUNCTION'})
state_status['Scope'] = 'Review completeness is not identity accuracy; frozen ID retained'
state_status.to_csv(out / 'STATE_PASS_PARTIAL_BLOCKED.tsv', sep='\t', index=False)

# Step2. 置信度标准为人工证据审阅，非数据拟合分数
criteria = '''# 置信度判定标准

这是人工、多证据候选程序审阅，不是计算出的分类概率，也不以Ecotype标签或表型关联作身份依据。

- High confidence：多个一致程序marker；父谱系锚点与候选程序可兼容；两外部队列都有≥3供者内匹配支持，marker排名/方向一致，且无足以推翻身份的严重矛盾。当前只用在CD8_T_S01的IFN程序与CD8_T_S02的TCF7/LEF1低分化程序，不代表细亚型已验证。
- Moderate confidence：多个相关marker支持可解释程序，有明确外部marker/锚点/至少一队列供者支持；第二背景或具体细亚型证据仍有限。可写谨慎候选描述，不定名经过验证的功能亚型。
- Low confidence：有多基因程序线索，但共享/泛化信号、稀少细胞、供者内对比不足或缺第二可评价谱系降低置信度。
- Unresolved：严重父谱系矛盾，或仅基础核内/加工信号不足以识别具体功能，或没有可评价状态间证据且身份不明。保留原ID；描述性程序不能代替正式命名。

任何Hallmark显著性、恢复Z、LOOCV成功或单基因均不能单独提升等级。ORA不显著不自动降低等级。细胞身份需未来同细胞父谱系+程序+阴性谱系证据；Mast_S01/Mast_S07/Glial_S02优先验证。

文献表保存全部面板比较；只有父谱系兼容且Discovery Top20至少2个面板基因命中，才列为Top20多基因候选对应。全关联基因弱重叠、单基因和跨谱系共享另标，不能用于强制命名。该门槛是本次透明审阅规则，不是Luca原文精确复现；没有宣称signature显著匹配或亚型等价。
'''
(out / 'audit/CONFIDENCE_AND_MATCHING_CRITERIA.md').write_text(criteria, encoding='utf-8')

# Step3. Git只读比较；历史改动仍保留，新文件只在隔离目录
current = subprocess.run(['git', '-C', str(root), '-c', 'core.quotepath=false', 'status', '--porcelain=v1'], capture_output=True, text=True, encoding='utf-8', check=True).stdout
(out / 'audit/git_status_after.txt').write_text(current, encoding='utf-8')
before = set((out / 'audit/git_status_before.txt').read_text(encoding='utf-8-sig').splitlines())
after = set(current.splitlines())
git_audit = {'ExistingStatusEntriesRemoved': sorted(before - after), 'NewStatusEntries': sorted(after - before), 'OnlyNewIsolatedDirectory': after - before == {'?? analysis_reports/'}, 'NoGitMutationCommands': True}
assert not git_audit['ExistingStatusEntriesRemoved']
assert git_audit['OnlyNewIsolatedDirectory']
(out / 'audit/GIT_STATUS_COMPARISON.json').write_text(json.dumps(git_audit, ensure_ascii=False, indent=2), encoding='utf-8')
tracked_modified = [line[3:] for line in before if line[:2].strip() == 'M']
pd.DataFrame([{'Path': str(root / name), 'SHA256AtDelivery': hashlib.sha256((root / name).read_bytes()).hexdigest(), 'Scope': 'delivery fingerprint; no byte-level baseline was taken at first git inspection'} for name in tracked_modified]).to_csv(out / 'audit/EXISTING_MODIFIED_FILES_DELIVERY_FINGERPRINT.tsv', sep='\t', index=False)
benchmark_scripts = [root / '06_validation/deconv_fraction_sensitivity_joint20_20261013/scripts/run_fraction_sensitivity_joint20_20261013.R', root / '06_validation/deconv_fraction_sensitivity_joint20_20261013/scripts/run_external_pseudobulk_joint20_20261013.R']
pd.DataFrame([{'Path': str(path), 'Bytes': path.stat().st_size, 'SHA256AtDelivery': hashlib.sha256(path.read_bytes()).hexdigest(), 'AccessByThisTask': 'read-only configuration review; not a start-of-task integrity baseline'} for path in benchmark_scripts]).to_csv(out / 'audit/BENCHMARK_READ_ONLY_CONFIGURATION_FINGERPRINT.tsv', sep='\t', index=False)

# Step4. 进展记录按实际落盘阶段记录时间，不冒称绘图完成
zone = timezone(timedelta(hours=8))
progress = ['# 实际执行进展（Asia/Shanghai）', '', '所有阶段均在本地独立目录执行；未打开新SSH或附接已有benchmark连接。', '']
stages = [('P0', 'audit/git_status_before.txt', '已只读检查目录/Git/本机进程与资源，创建版本化隔离目录'), ('P1', 'logs/01_inspect_sources.log', '已审计旧版结果小表结构，确认两版n4与既有谱系审计'), ('P1/P2/P3', 'logs/02_integrate_legacy_annotation.log', '已保存53状态基因证据、2650 Hallmark ORA、1325文献面板比较和159疾病-部位支持行'), ('P2数值QA', 'logs/03_verify_ora.log', 'R phyper与BH独立复核通过，最大P/Q差<1e-10'), ('P4报告QA', 'logs/04_write_reports_and_validate_final.log', '8项核心文件落盘，E1–E7逐成员解释，保存功能和marker图源；仅修正文献对应层级/表型单位/报告显示，未重读全DE或重跑模型'), ('P4图件', 'logs/05_plot_source_preflight_final.log', '绘图源静态检查19 PASS、2 WARN、0 FAIL；两个WARN为当前组会宽版布局与仅PNG预览。matplotlib缺失，渲染/最终图QA未执行'), ('交付', 'STATE_PASS_PARTIAL_BLOCKED.tsv', '保存53状态逐项PASS/PARTIAL/BLOCKED及历史全部排名副本')]
for stage, name, detail in stages:
    date = datetime.fromtimestamp((out / name).stat().st_mtime, zone).isoformat(timespec='seconds')
    progress.append(f'- {date} | {stage} | {detail}。产物：`{name}`。')
(out / 'PROGRESS_LOG.md').write_text('\n'.join(progress), encoding='utf-8')
counts = master.Confidence.value_counts().to_dict()
readme = f'''# 旧版BayesPrism注释交付入口

冻结53状态、E1–E7全部成员，2026-10-13。独立于正在运行的Joint20 benchmark；仅本地小表整合和单线程ORA。

| 产物 | 内容 |
|---|---|
| [CELL_STATE_ANNOTATION_MASTER.tsv](CELL_STATE_ANNOTATION_MASTER.tsv) | 53状态候选注释、置信度、marker/DE/外部/功能/表型和验证候选 |
| [CELL_STATE_MARKER_EVIDENCE.tsv](CELL_STATE_MARKER_EVIDENCE.tsv) | 38,233行marker、全部fold两版Top50排名和已有DE方向/效应 |
| [CELL_STATE_FUNCTIONAL_ENRICHMENT.tsv](CELL_STATE_FUNCTIONAL_ENRICHMENT.tsv) | 2650项Hallmark条件ORA，背景、命中、全局/状态内BH |
| [CELL_STATE_LITERATURE_MATCHING.tsv](CELL_STATE_LITERATURE_MATCHING.tsv) | 1325项状态×真实文献面板比较，含零/弱/跨谱系匹配 |
| [CELL_STATE_EXTERNAL_SUPPORT.tsv](CELL_STATE_EXTERNAL_SUPPORT.tsv) | 159项状态×疾病–部位背景支持，明确可评价范围 |
| [ECOTYPE_FUNCTIONAL_COMPOSITION.tsv](ECOTYPE_FUNCTIONAL_COMPOSITION.tsv) | 53冻结成员与E1–E7功能组成，不更改成员 |
| [ANNOTATION_METHODS_AND_AUDIT.md](ANNOTATION_METHODS_AND_AUDIT.md) | 方法对应、资源/来源审计、补充方法和缺口 |
| [CELL_STATE_AND_ECOTYPE_INTERPRETATION.md](CELL_STATE_AND_ECOTYPE_INTERPRETATION.md) | 53状态和7生态型完整证据卡及解释 |

候选程序置信度：{json.dumps(counts, ensure_ascii=False)}。High仅指CD8_T_S01和CD8_T_S02的功能程序，不是最终细亚型身份。Mast_S01/Mast_S07/Glial_S02均Unresolved。

已有保存全队列DE：41状态双队列、9单队列，Enteroendocrine_S01/Glial_S03/Glial_S02无可复用对比。三个联合背景不是旧n4重定义；实际保存n4为DiseaseGroup与Region3两版。UC疾病DE混合Colon/Rectum，不当作精确UC_Colon marker验证。

Hallmark ORA有17状态的63项全局FDR<0.05，仅针对冻结NMF候选空间。GO/Reactome与未预筛GSEA尚缺；文献面板对应不是正式亚型分类。

状态列表见[STATE_PASS_PARTIAL_BLOCKED.tsv](STATE_PASS_PARTIAL_BLOCKED.tsv)，模块状态见[PASS_PARTIAL_BLOCKED.tsv](PASS_PARTIAL_BLOCKED.tsv)，实际阶段日志见[PROGRESS_LOG.md](PROGRESS_LOG.md)。图源和Python脚本已保存；matplotlib缺失，新图尚未渲染，无最终图QA通过声明。置信度标准见[audit/CONFIDENCE_AND_MATCHING_CRITERIA.md](audit/CONFIDENCE_AND_MATCHING_CRITERIA.md)。

Git status前后原有条目完整保留，唯一新条目为隔离analysis_reports目录。34来源哈希一致，R数值复核通过。远端实时资源/benchmark状态未通过新登录确认；本任务未修改或控制benchmark。修改文件交付指纹与benchmark配置指纹是交付时记录，不冒充首次检查时的字节级基线。

脚本执行顺序：01检查→02整合→03 R复核→04报告/验证→06交付审计；05是绘图环境就绪后的待执行入口。所有方法与结果均为旧版探索性补充，不限制未来CIBERSORTx模型。
'''
(out / 'README.md').write_text(readme, encoding='utf-8')
delivery = []
for path in sorted(out.rglob('*')):
    if path.is_file() and path.name != 'DELIVERY_MANIFEST.tsv':
        delivery.append({'RelativePath': str(path.relative_to(out)), 'Bytes': path.stat().st_size, 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(delivery).to_csv(out / 'DELIVERY_MANIFEST.tsv', sep='\t', index=False)
print('Finalized', len(delivery), 'files; Git additions only isolated analysis_reports/')
print(state_status.ReviewStatus.value_counts().to_string())
