# 有限核对MSCCR原始counts中的免疫球蛋白/浆细胞相关基因，与当前网页CPM比较；不重处理整个队列。
# 在服务器ibd_refaudit环境运行；原始文件只读，读入明确gene ID清单与120样本名单，输出counts/CPM差异小表。
# CPM仍用全部源基因库总量作分母；不改gene symbol、补零、删marker或更改signature。
from pathlib import Path
import argparse
import csv
import gzip
import json
import numpy as np
import pandas as pd

# Step1. 明确原件和输出目录，参数只是文件位置
parser = argparse.ArgumentParser()
parser.add_argument('--raw-counts', required=True)
parser.add_argument('--gene-list', required=True)
parser.add_argument('--sample-list', required=True)
parser.add_argument('--web-mixture', required=True)
parser.add_argument('--outdir', required=True)
args = parser.parse_args()
source = Path(args.raw_counts)
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=False)
gene_list = pd.read_csv(args.gene_list, sep='\t', dtype=str)
sample_list = pd.read_csv(args.sample_list, sep='\t', dtype=str)
assert gene_list.msccr_ensembl_id.is_unique and gene_list.ensembl75_gene_name.is_unique and sample_list.ID.is_unique
targets = set(gene_list.msccr_ensembl_id)
source_stat = source.stat()

# Step2. 扫描原始文件，分母累计全部源基因，保留命中基因的原counts；循环是大文件读取所需
library = np.zeros(len(sample_list))
matched = {}
row_count = 0
with gzip.open(source, 'rt', encoding='utf-8-sig', newline='') as handle:
    # 本次已核对GEO原件：R write.table式空格分隔、名称有引号，标题只列2490个样本，数据另有首列gene ID。
    reader = csv.reader(handle, delimiter=' ', skipinitialspace=True)
    sample_header = next(reader)
    assert len(sample_header) == len(set(sample_header))
    header = ['gene_id'] + sample_header
    columns = [header.index(sample) for sample in sample_list.ID]
    assert len(set(columns)) == len(sample_list)
    for row in reader:
        assert len(row) == len(header)
        values = np.array([float(row[i]) for i in columns])
        assert np.isfinite(values).all() and (values >= 0).all()
        library += values
        gene_id = row[0].split('.')[0]
        if gene_id in targets:
            assert gene_id not in matched, '同一基因ID有多行，先保存来源诊断，不任取一行'
            matched[gene_id] = values
        row_count += 1
assert (library > 0).all()
counts = pd.DataFrame.from_dict(matched, orient='index', columns=sample_list.ID)
counts.to_csv(out / 'selected_raw_counts_by_Ensembl.tsv', sep='\t', index_label='msccr_ensembl_id')
pd.DataFrame({'sample': sample_list.ID, 'full_source_library_counts': library}).to_csv(out / 'full_source_library_counts.tsv', sep='\t', index=False)
missing = gene_list.loc[~gene_list.msccr_ensembl_id.isin(counts.index)]
missing.to_csv(out / 'missing_source_gene_IDs.tsv', sep='\t', index=False)
print('Read source gene rows:', row_count, '; matched targeted genes:', len(counts), '/', len(targets), flush=True)

# Step3. 同一分母重算少数基因CPM，与已上传的120样本输入比较，不把原件全表重导出
names = gene_list.set_index('msccr_ensembl_id').ensembl75_gene_name
raw_cpm = counts.div(library, axis=1) * 1e6
raw_cpm.index = raw_cpm.index.map(names)
web_parts = []
web_header = pd.read_csv(args.web_mixture, sep='\t', nrows=0).columns
for block in pd.read_csv(args.web_mixture, sep='\t', usecols=[web_header[0]] + sample_list.ID.tolist(), index_col=0, chunksize=500):
    hits = block.index.intersection(raw_cpm.index)
    if len(hits):
        web_parts.append(block.loc[hits, sample_list.ID])
web_cpm = pd.concat(web_parts)
common = raw_cpm.index.intersection(web_cpm.index)
differences = (raw_cpm.loc[common] - web_cpm.loc[common]).abs()
comparison = pd.DataFrame({'raw_median_CPM': raw_cpm.loc[common].median(axis=1), 'web_median_CPM': web_cpm.loc[common].median(axis=1), 'raw_mean_CPM': raw_cpm.loc[common].mean(axis=1), 'web_mean_CPM': web_cpm.loc[common].mean(axis=1), 'maximum_absolute_CPM_difference': differences.max(axis=1)})
comparison['agrees_with_8_significant_digit_export'] = [np.allclose(raw_cpm.loc[gene], web_cpm.loc[gene], rtol=1e-7, atol=1e-8) for gene in common]
comparison.to_csv(out / 'raw_vs_web_CPM.tsv', sep='\t', index_label='gene')
raw_cpm.to_csv(out / 'selected_raw_CPM_by_symbol.tsv', sep='\t', index_label='gene')
web_cpm.loc[common].to_csv(out / 'selected_web_CPM_by_symbol.tsv', sep='\t', index_label='gene')
assert source.stat().st_size == source_stat.st_size and source.stat().st_mtime_ns == source_stat.st_mtime_ns
summary = {'source_opened_read_only': True, 'source_size_and_mtime_unchanged': True, 'raw_format': 'space-delimited quoted sample header without gene-ID header', 'source_gene_rows': row_count, 'source_samples': len(sample_header), 'samples': len(sample_list), 'target_genes': len(targets), 'matched_source_genes': len(counts), 'common_web_genes': len(common), 'all_compared_values_agree_with_export_rounding': bool(comparison.agrees_with_8_significant_digit_export.all()), 'signature_modified': False, 'data_export_modified': False}
(out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(comparison.to_string(), flush=True)
