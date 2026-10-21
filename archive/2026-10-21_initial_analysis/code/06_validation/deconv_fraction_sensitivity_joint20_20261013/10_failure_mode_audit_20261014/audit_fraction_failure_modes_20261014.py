# 定位外部验证暴露的参考选择、signature特征和基因过滤问题；不改变任何模型或参数。
# 本机科学Python运行；仅保存已有矩阵的诊断表，不调用服务器，不以外部表现筛选新marker。
from pathlib import Path
import shutil
import pandas as pd

# Step1. 把原始图谱、实际训练子集和signature分开记录
project = Path(__file__).resolve().parents[2]
root = project / '06_validation/deconv_fraction_sensitivity_joint20_20261013'
out = root / '10_failure_mode_audit_20261014'
out.mkdir(exist_ok=False)
capacity = pd.read_csv(root / '01_reference/reference_type_capacity.tsv', sep='\t')
capacity.to_csv(out / 'actual_reference_type_capacity.tsv', sep='\t', index=False)
signature = pd.read_csv(project / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions/GJ20v2_Signature20.txt', sep='\t', index_col=0)
markers = ['IGKC', 'IGHA1', 'IGLC2', 'IGLC3', 'JCHAIN', 'IGHM', 'MZB1', 'IL7R', 'KIT', 'NCR2', 'CD3D', 'TRAC', 'NKG7', 'OLFM4', 'LGR5', 'MKI67', 'POU2F3', 'TRPM5', 'SH2D6']
columns = ['Plasma', 'ILC', 'CD4_T', 'CD8_T', 'NK', 'Stem_TA_progenitor', 'Tuft']
signature.reindex(markers)[columns].to_csv(out / 'frozen_signature_target_marker_expression.tsv', sep='\t', index_label='gene')

# Step2. 核实marker在BP实际基因集合中是否保留，不能把分母差异误判为缺失细胞类型
tables = []
for dataset in ['GSE214695', 'IBDverse']:
    audit = pd.read_csv(root / '08_external_pseudobulk' / dataset / '03_bayesprism/gene_filter_audit.tsv', sep='\t')
    selected = audit.loc[audit.gene.isin(markers)].copy()
    selected.insert(0, 'dataset', dataset)
    selected['present_in_official_signature'] = selected.gene.isin(signature.index)
    tables.append(selected)
pd.concat(tables, ignore_index=True).to_csv(out / 'BP_target_marker_retention.tsv', sep='\t', index=False)

# Step3. 结论限于可核验现象；下一轮定位实验不等于新的独立外部验证
text = '''# 外部反卷积误差定位：2026-10-14

不因两套外部验证的偏差直接放弃GSE282122。原图谱覆盖广，但当前实际参考是9294个选中训练细胞、32位供者、212个汇总profiles；CIBERSORTx再提取4474基因signature。原图谱规模不等于实际参考对每个亚型和状态均有代表性。

## 已核验的事实

- Plasma与Stem/TA各600个训练细胞；ILC有402个、来自10位供者；Tuft有221个、来自6位供者。不存在这四类完全未进入参考的问题。这些数量不直接证明类型内状态覆盖充足，也不直接证明不足。
- 冻结signature包含Tuft的POU2F3、TRPM5、SH2D6；ILC的IL7R、KIT、NCR2；Stem/TA的OLFM4、LGR5、MKI67。因此不能把失败自动解释为未选到任何对应marker。
- 两个外部BP运行中，IGKC、IGHA1、IGLC2、IGLC3、JCHAIN均在clean reference中，但不在new.prism最终保留基因集合。日志分别记录14与15个outlier基因被过滤。过滤参数沿用默认outlier.cut=0.01、outlier.fraction=0.1；这不是为本次比较事后调参。
- GSE214695 Plasma的实际保留基因RNA量仅占其全基因RNA量17.79%，不同细胞类型保留比例不同。两种truth并非统一缩放；CIBERSORTx相对cell-count truth的Plasma高估不能与BP相对保留基因truth的误差直接混为同一估计量。
- IBDverse作者ILC标签的IL7R/KIT表达与ILC身份相符，是支持注释合理性的初步marker检查，不能替代细胞级复核；Job26全零预测更需定位特征分离、平台差异和混合求解，不是删除该真值类型的理由。

## 推荐的处理顺序

1. 保留GSE282122、Joint20 v2及现有无校正正式结果。保持53 states、7 ecotypes冻结。
2. 优先审查基因分母与Ig基因过滤。若试验保留这些基因，另立明确参数分支，重新计算对应RNA truth；不能只挑性能最好的一组替换主结果，也不能默认关闭过滤就更正确。
3. 仅对Plasma、ILC、Stem/TA、Tuft做小规模诊断：用开发/训练来源的独立holdout细胞构建单类型与已知配比混合物，检查同源条件下能否恢复；再区分外部平台/标签差异。训练来源控制只作算法定位，不能称独立外部accuracy。
4. S-mode可作探索性敏感性，但现有1060单细胞子集不等于Job13实际212列汇总构建输入；不能保证修复上述问题，也不能作为唯一修复手段。
5. 若确需修改signature，在开发集完成并另存v3，然后冻结后用新的未参与调参的供者/数据验证。当前两套外部结果已经看过，不能反复据其表现调参后仍称未触碰的独立验证。暂不增加新数据集或重跑discovery。

当前定位仍不提供ILC全零或Tuft高估的唯一因果解释，不做未经验证的RNA比例到细胞比例换算。

方法依据：[BayesPrism官方FAQ](https://github.com/Danko-Lab/BayesPrism/blob/main/README.md)：new.prism默认过滤阈值以及fractions对应所用基因集合reads贡献；[反卷积benchmark原始研究](https://www.nature.com/articles/s41467-020-19015-1)：参考矩阵、marker和归一化对性能的影响。
'''
(out / 'FAILURE_MODE_AUDIT_AND_NEXT_STEPS.md').write_text(text, encoding='utf-8')
shutil.copyfile(Path(__file__), out / Path(__file__).name)
print('Read-only diagnosis saved:', out)
print('No signature, BP parameters, frozen states or ecotypes modified')
