# 用独立PyFirth Newton算法复核本次SciPy惩罚似然拟合。
# 不替换主分析；比较全部360个单特征模型的系数及似然比P值。
from pathlib import Path
import pandas as pd
import numpy as np
import json
from pyfirth.PyFirth import PyFirth
from scipy.special import logit
from scipy.stats import chi2

# Step1. 从已保存的真实完整病例与连续评分重建相同设计矩阵。
out = Path(__file__).resolve().parents[1]
stats = pd.read_csv(out / '04_cross_cohort_summary/all_associations.tsv', sep='\t')
rows = []
folders = {'GSE72819': '01_GSE72819_treatment_outcome', 'GSE93624': '02_RISK_disease_progression', 'GSE73661': '03_GSE73661_mucosal_healing'}
for (accession, endpoint), d in stats.groupby(['Accession', 'Endpoint']):
    folder = out / folders[accession]
    patients = pd.read_csv(folder / (endpoint + '_analysis_patients.tsv'), sep='\t')
    for family, score_file in [('Ecotype', 'continuous_ecotype_scores.tsv'), ('State', 'continuous_state_scores.tsv')]:
        scores = pd.read_csv(folder / score_file, sep='\t', index_col=0)
        for _, row in d.loc[d.Family.eq(family) & d.Status.eq('OK')].iterrows():
            x = scores.loc[patients.SampleID, row.Feature].to_numpy()
            check = pd.DataFrame({'z': (x - x.mean()) / x.std(ddof=1), 'y': patients.Outcome.to_numpy()})
            independent = PyFirth(check, ['z'], 'y')
            fit = independent._newton_rapheson(step_limit=100, num_iters=1000, convergence_limit=1e-6, max_step=5.0)
            beta = fit['ParamTable'].loc['z', 'BETA']
            # 在斜率=0的两参数模型中，完整Jeffreys惩罚给出可解析截距p=(事件+1)/(N+2)。
            null_beta = np.array([logit((check.y.sum() + 1) / (len(check) + 2)), 0.0])
            independent_p = chi2.sf(2 * (fit['LogLike'] - independent.firth_likelihood(null_beta)), 1)
            coef = fit['ParamTable'].BETA.to_numpy()
            p = independent.model.predict(coef)
            w = p * (1 - p)
            info = independent.X.T @ (w[:, None] * independent.X)
            h = w * np.sum((independent.X @ np.linalg.inv(info)) * independent.X, axis=1)
            adjusted_score = independent.X.T @ (check.y.to_numpy() - p + h * (.5 - p))
            rows.append({'Endpoint': endpoint, 'Feature': row.Feature, 'Beta_primary': row.Beta, 'Beta_PyFirth': beta, 'AbsBetaDifference': abs(row.Beta - beta), 'P_primary': row.P_Firth, 'P_independent_full_penalty_null': independent_p, 'AbsPDifference': abs(row.P_Firth - independent_p), 'PyFirthFlag': fit['Flag'], 'IndependentMaxAdjustedScore': abs(adjusted_score).max()})
validation = pd.DataFrame(rows)
validation.to_csv(out / 'audit/independent_PyFirth_validation.tsv', sep='\t', index=False)
assert len(validation) == 360
assert validation.AbsBetaDifference.max() < 1e-4 and validation.AbsPDifference.max() < 1e-5
assert validation.IndependentMaxAdjustedScore.max() < 1e-5
(out / 'audit/independent_Firth_summary.json').write_text(json.dumps({'Models': len(validation), 'MaxBetaDifference': float(validation.AbsBetaDifference.max()), 'MaxPDifference': float(validation.AbsPDifference.max()), 'NewtonWarningN': int(validation.PyFirthFlag.ne(0).sum()), 'MaxAdjustedScore': float(validation.IndependentMaxAdjustedScore.max()), 'NullMethod': 'Closed-form nuisance optimum retaining full two-parameter Jeffreys penalty', 'Verdict': 'PASS'}, indent=2))
print('INDEPENDENT_FIRTH_VALIDATION_PASS', len(validation), validation.AbsBetaDifference.max(), validation.AbsPDifference.max(), flush=True)
