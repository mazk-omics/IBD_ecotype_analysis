# 用同一 NNLS 求解器和同一 signature 基因集，比较原始训练均值与网页 signature 对开发留出混合物的拟合。
# 在本地 cibersortx_review 环境运行；读取已审计的 Joint20、28 个自然和 8 个等RNA混合物，保存诊断。
# 这是开发集上的表达谱敏感性分析，不是新的反卷积方法验证，不读取最终留出 CD12/UC17。
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls

# Step1. 固定基因、类型和已有开发输入，三个模型均不调参
root = Path(__file__).resolve().parents[2]
chain = root / 'CIBERSORTx_signature_audit/construction_chain_20261007'
validation = root / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006/server_outputs/validation_joint20_20261006'
signature_path = next((root / 'CIBERSORTx_signature_audit/three_candidate_comparison_20261006/gse_joint').glob('*_signature_gene_symbols_restored.txt'))
signature = pd.read_csv(signature_path, sep='\t', index_col=0)
source_mean = pd.read_csv(chain / 'gse_joint/source_equal_cell_mean_CPM.tsv', sep='\t', index_col=0).loc[signature.index, signature.columns]
source_pool = pd.read_csv(chain / 'gse_joint/source_pooled_counts_CPM.tsv', sep='\t', index_col=0).loc[signature.index, signature.columns]
models = [('web_signature', signature), ('training_equal_cell_mean_CPM', source_mean), ('training_pooled_counts_CPM', source_pool)]
records = []
bias_records = []
for construct, mixture_file, truth_file in [('natural', 'GJ20_Holdout28_Mixture_CPM.txt', 'ground_truth_export_UMI_fractions.tsv'), ('equalRNA', 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt', 'ground_truth_equalRNA_fractions.tsv')]:
    mixture = pd.read_csv(validation / mixture_file, sep='\t', index_col=0).loc[signature.index]
    truth = pd.read_csv(validation / truth_file, sep='\t', index_col=0).loc[mixture.columns, signature.columns]
    assert np.allclose(truth.sum(axis=1), 1)
    regions = mixture.columns.str.split('_').str[0]
    # Step2. 在线性空间拟合，全部模型使用不加权 NNLS、同样的系数和归一化规则
    for model_name, profiles in models:
        matrix = profiles.to_numpy(float)
        scale = matrix.max()
        predictions = pd.DataFrame(index=mixture.columns, columns=signature.columns, dtype=float)
        coefficient_sums = []
        for sample in mixture.columns:
            coefficient, residual = nnls(matrix / scale, mixture[sample].to_numpy() / scale)
            assert coefficient.sum() > 0
            predictions.loc[sample] = coefficient / coefficient.sum()
            coefficient_sums.append(float(coefficient.sum()))
        predictions.index.name = 'sample'
        predictions.to_csv(chain / f'NNLS_development_{construct}_{model_name}.tsv', sep='\t')
        error = (predictions - truth).abs()
        for region in regions.unique():
            selected = regions == region
            for cell_type in signature.columns:
                bias_records.append({'construct': construct, 'model': model_name, 'region': region, 'cell_type': cell_type,
                    'truth_percent': float(truth.loc[selected, cell_type].mean() * 100), 'prediction_percent': float(predictions.loc[selected, cell_type].mean() * 100),
                    'MAE_pp': float(error.loc[selected, cell_type].mean() * 100)})
        for sample, coefficient_sum in zip(mixture.columns, coefficient_sums):
            records.append({'construct': construct, 'model': model_name, 'sample': sample, 'region': sample.split('_')[0],
                'MAE20_pp': float(error.loc[sample].mean() * 100), 'coefficient_sum_before_relative_normalization': coefficient_sum})

# Step3. 保存结果；模型差异仅定位参考谱/求解的敏感性，不能单独归因于网页筛选
metrics = pd.DataFrame(records)
metrics.to_csv(chain / 'NNLS_development_sample_metrics.tsv', sep='\t', index=False)
pd.DataFrame(bias_records).to_csv(chain / 'NNLS_development_region_type_metrics.tsv', sep='\t', index=False)
summary = metrics.groupby(['construct', 'model', 'region'], observed=True).MAE20_pp.mean().reset_index()
summary.to_csv(chain / 'NNLS_development_region_summary.tsv', sep='\t', index=False)
(chain / 'NNLS_diagnostic_limits.json').write_text(json.dumps({'same_genes_and_solver': True, 'models_do_not_use_holdout_to_fit_parameters': True,
    'development_holdout_results_are_used_for_diagnosis': True, 'final_holdout_used': False, 'not_CIBERSORTx': True,
    'not_validation_of_alternative_method': True, 'differences_include_profile_aggregation_and_web_processing': True}, indent=2), encoding='utf-8')
print(summary.to_string(index=False))
