# 用真实MSCCR bulk计算50个Hallmark的平均log2单位方差表达，再比较dominant ecotype。
# 服务器 ibd_refaudit 环境运行；读取既有gene-symbol harmonized counts及2025样本定义，原件只读。
# 这是论文算法的IBD适配：CPM基于harmonized完整输入库，MSigDB7.4版本明确但原文版本尚未确认。
from pathlib import Path
import csv
import gzip
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 读取固定的七生态型和2025个CD/UC样本，不按统计结果筛选样本
root = Path('/home/mazekai/IBD_EcoTyper')
out = root / '06_validation/paper_aligned_recovery_characterization_20261012'
annotation = pd.read_csv(out / '06_ecotype_features/ecotype_dominant_assignments.tsv', sep='\t', index_col=0)
assert len(annotation) == 2025 and annotation.index.is_unique
gmt_path = out / '00_audit/h.all.v7.4.symbols.gmt'
gene_sets = {}
for line in gmt_path.read_text().splitlines():
    fields = line.split('\t')
    gene_sets[fields[0]] = fields[2:]
assert len(gene_sets)==50
hallmark_genes = set().union(*gene_sets.values())
source = root / '03_reference/combined_reference/output/09_bayesprism_input_harmonization/GSE193677_MSCCR_Biopsy_counts_GRCh37_Ensembl75_gene_symbol_unique.txt.gz'
source_hash_before = hashlib.sha256(source.read_bytes()).hexdigest()

# Step2. 逐行累计完整harmonized库计数，只保留Hallmark所需基因；CPM不以通路基因作分母
selected = []
symbols = []
seen = set()
with gzip.open(source, 'rt') as handle:
    header = next(csv.reader([handle.readline().strip()], delimiter=' ', skipinitialspace=True))
    assert len(set(header)) == len(header)
    positions = pd.Index(header).get_indexer(annotation.index)
    assert (positions>=0).all()
    library = np.zeros(len(header))
    for line in handle:
        name, values = line.split(None, 1)
        gene = name.strip('"')
        values = np.fromstring(values, sep=' ')
        assert len(values)==len(header) and gene not in seen and np.isfinite(values).all() and (values>=0).all()
        seen.add(gene)
        library += values
        if gene in hallmark_genes:
            selected.append(values[positions])
            symbols.append(gene)
assert (library>0).all()
counts = np.vstack(selected)
cpm = counts * (1e6/library[positions])
log2_expression = np.log2(cpm+1)
sd = log2_expression.std(axis=1, ddof=1)
variable = sd>0
scaled = (log2_expression[variable]-log2_expression[variable].mean(axis=1, keepdims=True))/sd[variable, None]
expression = pd.DataFrame(scaled, index=np.array(symbols)[variable], columns=annotation.index)
print('Bulk:', len(seen), 'harmonized genes;', len(annotation), 'biopsies;', len(expression), 'variable Hallmark genes', flush=True)

# Step3. 每条通路以实际可用且有方差的全部基因取均值，不使用GSEA/ssGSEA替换算法
scores = []
coverage = []
for name, genes in gene_sets.items():
    available = expression.index.intersection(genes)
    assert len(available)>0
    score = expression.loc[available].mean(axis=0)
    score.name = name
    scores.append(score)
    coverage.append({'Feature': name, 'GMTGenes': len(genes), 'AvailableVariableGenes': len(available), 'UsedGenes': ';'.join(available), 'MissingOrZeroVarianceGenes': ';'.join(sorted(set(genes)-set(available))), 'MSigDBVersion': '7.4'})
scores = pd.DataFrame(scores).T
scores.to_csv(out / '06_ecotype_features/hallmark_sample_scores.tsv', sep='\t')
pd.DataFrame(coverage).to_csv(out / '06_ecotype_features/hallmark_gene_coverage.tsv', sep='\t', index=False)
results = []
for ecotype in ['E'+str(i) for i in range(1, 8)]:
    own = annotation.DominantEcotype.eq(ecotype)
    other = annotation.DominantEcotype.notna() & ~own
    assert own.sum()>0 and other.sum()>0
    for feature in scores.columns:
        a, b = scores.loc[own, feature], scores.loc[other, feature]
        p = stats.mannwhitneyu(a, b, alternative='two-sided', method='asymptotic').pvalue
        results.append({'Dataset': 'MSCCR_discovery', 'Ecotype': ecotype, 'Feature': feature, 'N_Ecotype': len(a), 'N_OtherEcotypes': len(b), 'MeanDifference': a.mean()-b.mean(), 'P': p, 'Method': 'PAPER_ADAPTED_IBD', 'BHFamily': 'All 50 Hallmarks x 7 ecotypes; 350 hypotheses'})
results = pd.DataFrame(results)
results['Q'] = stats.false_discovery_control(results.P.to_numpy(), method='bh')
results.to_csv(out / '06_ecotype_features/ecotype_continuous_features.tsv', sep='\t', index=False)
results.pivot(index='Feature', columns='Ecotype', values='MeanDifference').to_csv(out / '06_ecotype_features/ecotype_feature_matrix.tsv', sep='\t')
assert hashlib.sha256(source.read_bytes()).hexdigest()==source_hash_before
(out / '06_ecotype_features/hallmark_input_audit.json').write_text(json.dumps({'Input': str(source), 'SHA256': source_hash_before, 'InputUnchanged': True, 'Samples': len(annotation), 'RawHarmonizedGeneRows': len(seen), 'HallmarkVariableGenes': len(expression), 'Normalization': 'CPM denominator all rows of harmonized source, then log2(CPM+1), per-gene centered unit SD ddof1', 'MSigDBVersion': '7.4', 'OriginalMSigDBVersionConfirmed': False, 'SampleUnit': 'Biopsy; repeated patients present', 'PaperParity': 'PAPER_ADAPTED_IBD'}, indent=2))
print('HALLMARK_FEATURES_COMPLETE:', int(results.Q.lt(0.05).sum()), 'of', len(results), 'q<0.05', flush=True)
