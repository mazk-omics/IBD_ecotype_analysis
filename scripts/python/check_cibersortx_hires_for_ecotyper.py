# 检查 HiRes 每类一个 gene×sample 非log表达矩阵能否接入 EcoTyper presorted；不改变或修补输入。
# 使用项目已有 numpy/pandas 的 Python 环境运行；读取指定目录及注释，写出 QC 表与来源哈希。
# 本脚本检查结构和动态范围，不把格式合格解释为生物学验证通过。
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

# Step1. 指定实际 HiRes 目录；不能将 signature、fractions 或 group-mode 表填入此处
parser = argparse.ArgumentParser()
parser.add_argument('--expression-dir', required=True)
parser.add_argument('--annotation', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--expected-samples', type=int, default=2025)
parser.add_argument('--expected-types', default='B,Plasma,CD4_T,CD8_T,NK,ILC,Monocyte_Macrophage,DC,Mast,Fibroblast,Pericyte,Endothelial,Ileal_absorptive,Colonic_absorptive,Goblet,Stem_TA_progenitor,Paneth,Enteroendocrine,Glial,Tuft')
args = parser.parse_args()
expression_dir = Path(args.expression_dir).resolve()
out = Path(args.outdir).resolve()
assert out != expression_dir and expression_dir not in out.parents, 'QC 输出必须放在表达输入目录外，EcoTyper 会读取目录中的每个文件'
out.mkdir(parents=True, exist_ok=True)
(out / 'input_check_status.json').write_text(json.dumps({'status': 'checking_not_yet_passed'}, indent=2), encoding='utf-8')
annotation_path = Path(args.annotation).resolve()
annotation = pd.read_csv(annotation_path, sep='\t', dtype=str, keep_default_na=False)
assert {'ID', 'Disease', 'Region3'}.issubset(annotation.columns)
assert annotation.ID.is_unique and annotation.ID.ne('').all()
assert len(annotation) == args.expected_samples
assert annotation.Disease.isin(['CD', 'UC']).all() and annotation.Region3.isin(['Colon', 'Ileum', 'Rectum']).all()
sample_ids = annotation.ID.tolist()
expected_types = args.expected_types.split(',')
files = sorted(expression_dir.iterdir())
assert len(files) == len(expected_types) and all(path.is_file() and path.suffix == '.txt' for path in files), '输入目录只放每类一个 .txt，说明文件和结果汇总放在外面'
assert {path.stem for path in files} == set(expected_types), '类型文件名须与选定 reference 的类型完全一致'
records = []
gene_sets = []
hashes = []

# Step2. 分块扫描全部表达；顺序必须与注释一致，禁止默默排序或填补 NaN
for path in files:
    with path.open(encoding='utf-8-sig', newline='') as handle:
        header = next(csv.reader(handle, delimiter='\t'))
    assert header[1:] == sample_ids, f'{path.name}: 样本顺序/数量与注释不一致'
    assert len(header) == len(set(header))
    with path.open('rb') as handle:
        before = hashlib.file_digest(handle, 'sha256').hexdigest()
    genes = []
    total_by_sample = np.zeros(len(sample_ids))
    expression_min = float('inf')
    expression_max = 0.0
    n_variable = 0
    n_all_zero = 0
    region_variable = {region: 0 for region in annotation.Region3.unique()}
    marker_max = {}
    for block in pd.read_csv(path, sep='\t', chunksize=500, keep_default_na=False):
        block_genes = block.iloc[:, 0].astype(str)
        assert block_genes.ne('').all() and not block_genes.str.contains(r'[\t\r\n]').any()
        values = block.iloc[:, 1:].to_numpy(dtype=float)
        assert np.isfinite(values).all() and (values >= 0).all(), f'{path.name}: 存在 NA/Inf/负值，请追溯 HiRes 输出，不直接填零'
        genes.extend(block_genes.tolist())
        total_by_sample += values.sum(axis=0)
        expression_min = min(expression_min, float(values.min()))
        expression_max = max(expression_max, float(values.max()))
        log_values = np.log2(values + 1)
        n_variable += int((np.ptp(log_values, axis=1) > 1e-12).sum())
        n_all_zero += int((values == 0).all(axis=1).sum())
        for region in region_variable:
            columns = annotation.Region3.eq(region).to_numpy()
            region_variable[region] += int((np.ptp(log_values[:, columns], axis=1) > 1e-12).sum())
        for gene in ['DEFA5', 'DEFA6', 'MUC2', 'OLFM4', 'LGR5', 'CHGA', 'POU2F3', 'SOX10', 'PLP1', 'TPSAB1']:
            if gene in block_genes.tolist():
                marker_max[gene] = float(values[block_genes.eq(gene).to_numpy()].max())
    assert len(genes) == len(set(genes)) and n_variable > 0
    assert total_by_sample.sum() > 0, f'{path.name}: 整个类型没有表达信息'
    with path.open('rb') as handle:
        after = hashlib.file_digest(handle, 'sha256').hexdigest()
    assert before == after
    gene_sets.append(set(genes))
    hashes.append({'file': str(path), 'sha256': after, 'unchanged': True})
    records.append({'cell_type': path.stem, 'n_genes': len(genes), 'n_samples': len(sample_ids), 'min': expression_min, 'max': expression_max,
        'genes_variable_across_samples': n_variable, 'all_zero_genes': n_all_zero,
        'variable_genes_by_region': json.dumps(region_variable), 'column_sum_min': float(total_by_sample.min()), 'column_sum_max': float(total_by_sample.max()),
        'n_all_zero_profiles_needing_review': int((total_by_sample == 0).sum()),
        'diagnostic_marker_max': json.dumps(marker_max), 'linear_expression_provenance_still_required': True})
    print(path.stem, 'genes/samples/variable:', len(genes), len(sample_ids), n_variable, flush=True)

# Step3. 保存通过范围；仍需人工核对 HiRes provenance、谱系表达及 Region3 修改源码
common_genes = set.intersection(*gene_sets)
assert common_genes, '没有跨类型共同基因，无法开展当前类型特异性筛选'
pd.DataFrame(records).to_csv(out / 'hires_input_QC.tsv', sep='\t', index=False)
pd.DataFrame(hashes).to_csv(out / 'input_hashes.tsv', sep='\t', index=False)
summary = {'status': 'structural_checks_passed_expression_provenance_and_biology_pending', 'expression_dir': str(expression_dir),
    'annotation': str(annotation_path), 'samples': len(sample_ids), 'cell_types': len(files), 'common_genes': len(common_genes),
    'identical_gene_sets': all(genes == gene_sets[0] for genes in gene_sets), 'fractions_are_not_expression_input': True,
    'Region3_modified_server_source_not_verified_by_this_script': True, 'NMF_not_started': True}
(out / 'input_check_status.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('结构检查通过；尚未确认表达尺度来源、谱系准确性或服务器 Region3 源码。')
