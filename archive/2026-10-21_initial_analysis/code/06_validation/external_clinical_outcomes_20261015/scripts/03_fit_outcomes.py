# 在服务器 ibd_outcomes_20261015 环境执行冻结连续评分的单特征结局关联。
# 按终点分别保存7 ecotypes与53 states的完整结果；Firth算法循环仅用于拟合和剖面区间。
from pathlib import Path
import json, hashlib, warnings, sys
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize, brentq, check_grad
from scipy.stats import mannwhitneyu, chi2

# Step1. Firth Jeffreys惩罚似然；无截断、无正则化超参数、无特征筛选。
def firth_objective(beta, X, y):
    eta = X @ beta
    p = expit(eta)
    w = p * (1 - p)
    info = X.T @ (w[:, None] * X)
    sign, logdet = np.linalg.slogdet(info)
    if sign <= 0:
        return np.inf, np.full(len(beta), np.nan)
    leverage = w * np.sum((X @ np.linalg.inv(info)) * X, axis=1)
    score = X.T @ (y - p + leverage * (0.5 - p))
    loglik = np.sum(y * eta - np.logaddexp(0, eta)) + 0.5 * logdet
    return -loglik, -score

def firth_fit(X, y, interest=1):
    fit = minimize(firth_objective, np.zeros(X.shape[1]), args=(X, y), jac=True, method='BFGS', options={'gtol': 1e-8, 'maxiter': 1000})
    gradient = np.max(np.abs(firth_objective(fit.x, X, y)[1]))
    if not np.isfinite(fit.fun) or gradient > 1e-5:
        raise ValueError('Firth non-convergence; max adjusted score=' + str(gradient))
    free = np.arange(X.shape[1]) != interest
    # 固定目标系数后仍保留完整模型的Jeffreys惩罚，不能改成低维惩罚。
    def profile(b):
        def objective(nuisance):
            coef = fit.x.copy()
            coef[interest] = b
            coef[free] = nuisance
            val, grad = firth_objective(coef, X, y)
            return val, grad[free]
        restricted = minimize(objective, fit.x[free], jac=True, method='BFGS', options={'gtol': 1e-8, 'maxiter': 1000})
        g = np.max(np.abs(objective(restricted.x)[1]))
        if not np.isfinite(restricted.fun) or g > 1e-5:
            raise ValueError('Profile non-convergence; max score=' + str(g))
        return restricted.fun
    lr = max(0.0, 2 * (profile(0) - fit.fun))
    cutoff = chi2.ppf(0.95, 1)
    beta = fit.x[interest]
    bounds = []
    for direction in [-1, 1]:
        step = 0.5
        edge = beta + direction * step
        while 2 * (profile(edge) - fit.fun) < cutoff and step < 64:
            step *= 2
            edge = beta + direction * step
        if step >= 64 and 2 * (profile(edge) - fit.fun) < cutoff:
            raise ValueError('Unbounded profile CI')
        bounds.append(brentq(lambda b: 2 * (profile(b) - fit.fun) - cutoff, min(beta, edge), max(beta, edge), xtol=1e-7))
    return {'Beta': beta, 'OR': np.exp(beta), 'CI_low': np.exp(bounds[0]), 'CI_high': np.exp(bounds[1]), 'P_Firth': chi2.sf(lr, 1), 'PenalizedLR': lr, 'MaxAdjustedScore': gradient, 'Iterations': fit.nit}

def bh_family(values, family_size):
    p = np.asarray(values, dtype=float)
    q = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    adjusted = p[order] * family_size / np.arange(1, len(order) + 1)
    q[order] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return q

# Step2. 用可解析的完全分离2x2表和变换不变性验证Firth实现。
root = Path('/home/mazekai/IBD_EcoTyper')
out = root / '06_validation/external_clinical_outcomes_20261015'
X = np.column_stack([np.ones(20), np.repeat([0, 1], 10)])
y = np.repeat([0, 1], 10)
test = firth_fit(X, y)
expected = np.log((10.5 / 0.5) / (0.5 / 10.5))
assert abs(test['Beta'] - expected) < 1e-5
reverse = firth_fit(X, 1 - y)
assert abs(test['Beta'] + reverse['Beta']) < 1e-5
assert abs(test['P_Firth'] - reverse['P_Firth']) < 1e-7
scaled = X.copy()
scaled[:, 1] *= 2
assert abs(firth_fit(scaled, y)['Beta'] * 2 - test['Beta']) < 1e-5
gradient_error = check_grad(lambda b: firth_objective(b, X, y)[0], lambda b: firth_objective(b, X, y)[1], np.array([-.4, .7]))
assert gradient_error < 1e-5
(out / 'audit/Firth_analytic_validation.json').write_text(json.dumps({'CompleteSeparation': test, 'ExpectedLogOR': expected, 'OutcomeReversal': 'PASS', 'CovariateRescaling': 'PASS', 'AnalyticGradientError': gradient_error}, indent=2))

# Step3. 按实际输出读取连续H和官方ecotype abundance，独立核对定义。
membership = pd.read_csv(out / 'audit/frozen_state_ecotype_membership.tsv', sep='\t')
cohorts = [('GSE72819', '01_GSE72819_treatment_outcome', root / '06_validation/02_bulk_treatment/GSE72819/official_recovery_20261009_185536/RecoveryOutput/GSE72819_BP73_20261009_185536'), ('GSE73661', '03_GSE73661_mucosal_healing', out / '03_GSE73661_mucosal_healing/RecoveryOutput/GSE73661_RMA_W0_67_outcomes_20261015')]
risk_result = out / '02_RISK_disease_progression/RecoveryOutput/GSE93624_processed245_outcomes_20261015'
if (risk_result / 'Ecotypes/ecotype_abundance.txt').exists():
    cohorts.append(('GSE93624', '02_RISK_disease_progression', risk_result))
if len(sys.argv) > 1:
    cohorts = [c for c in cohorts if c[0] in sys.argv[1:]]
all_results = []
all_denominators = []
for acc, folder, recovered in cohorts:
    state_tables = []
    for parent in membership.CellType.unique():
        H = pd.read_csv(recovered / parent / 'state_abundances.txt', sep='\t', index_col=0)
        H.index = parent + '_' + H.index
        assert H.index.is_unique and H.columns.is_unique
        assert np.isfinite(H.to_numpy()).all() and H.ge(0).all().all()
        state_tables.append(H)
    states = pd.concat(state_tables).loc[membership.ID]
    eco = pd.read_csv(recovered / 'Ecotypes/ecotype_abundance.txt', sep='\t', index_col=0).loc[['E' + str(i) for i in range(1, 8)], states.columns]
    expected_eco = states.groupby(membership.set_index('ID').Ecotype).mean()
    expected_eco = expected_eco.div(expected_eco.sum(axis=0), axis=1)
    assert np.allclose(eco, expected_eco.loc[eco.index, eco.columns], atol=1e-10)
    assert len(states) == 53 and len(eco) == 7 and np.isfinite(eco.to_numpy()).all()
    states.T.rename_axis('SampleID').to_csv(out / folder / 'continuous_state_scores.tsv', sep='\t')
    eco.T.rename_axis('SampleID').to_csv(out / folder / 'continuous_ecotype_scores.tsv', sep='\t')
    assignment = pd.read_csv(recovered / 'Ecotypes/ecotype_assignment.txt', sep='\t')
    state_qc = pd.DataFrame({'Feature': states.index, 'N': states.shape[1], 'Min': states.min(axis=1), 'Max': states.max(axis=1), 'SD': states.std(axis=1), 'ExactZeroN': states.eq(0).sum(axis=1), 'NearZeroN_1e8': states.le(1e-8).sum(axis=1)})
    state_qc.to_csv(out / folder / 'state_recovery_QC.tsv', sep='\t', index=False)
    pd.DataFrame([{'AllContinuousSamples': states.shape[1], 'States': 53, 'Ecotypes': 7, 'DiscreteEcotypeAssignedN': len(assignment), 'EcotypeFormulaMaxError': np.abs(eco - expected_eco.loc[eco.index, eco.columns]).to_numpy().max(), 'ScoreMeaning': 'State retained H normalized over all initial states within parent; ecotype is normalized mean constituent H; not absolute cell fractions or recovery Z'}]).to_csv(out / folder / 'recovery_QC.tsv', sep='\t', index=False)
    manifest = pd.read_csv(out / folder / 'cohort_endpoint_manifest.tsv', sep='\t')
    manifest['RecoveryComplete'] = manifest.SampleID.isin(eco.columns)
    manifest['DiscreteEcotypeAssigned'] = manifest.SampleID.isin(assignment.ID)
    manifest['ExclusionReason'] = manifest.ExclusionReason.fillna('')
    if acc == 'GSE93624':
        manifest['ExclusionReason'] = np.where(manifest.Disease.eq('CD'), '', 'Non-IBD control')
    manifest.to_csv(out / folder / 'cohort_endpoint_manifest.tsv', sep='\t', index=False)
    cohort_results = []
    for endpoint, full in manifest.groupby('Endpoint', sort=False):
        d = full.loc[full.ExclusionReason.eq('') & full.Outcome.notna() & full.RecoveryComplete].copy()
        assert d.PatientID.is_unique and d.SampleID.is_unique
        assert d.Timepoint.isin(['Baseline', 'W0', 'Diagnosis; treatment naive']).all()
        y = d.Outcome.to_numpy(dtype=float)
        n1 = int(y.sum())
        n0 = len(y) - n1
        all_denominators.append({'Accession': acc, 'Endpoint': endpoint, 'N': len(d), 'Events': n1, 'Nonevents': n0, 'BaselineCandidates': len(full), 'RecoveryN': int(full.RecoveryComplete.sum()), 'ExcludedN': len(full) - len(d), 'Interpretation': 'Processed-scale conditional association; transformation provenance unresolved' if acc == 'GSE93624' else 'Baseline association in biopsy-outcome complete cases'})
        d.to_csv(out / folder / (endpoint + '_analysis_patients.tsv'), sep='\t', index=False)
        for family, scores, size in [('Ecotype', eco, 7), ('State', states, 53)]:
            results = []
            for feature in scores.index:
                x = scores.loc[feature, d.SampleID].to_numpy(dtype=float)
                row = {'Accession': acc, 'Endpoint': endpoint, 'Family': family, 'Feature': feature, 'ParentLineage': membership.set_index('ID').CellType.get(feature, ''), 'Ecotype': membership.set_index('ID').Ecotype.get(feature, feature), 'N': len(y), 'Events': n1, 'Nonevents': n0, 'FamilySize': size, 'Adjustment': 'None', 'ORUnit': 'Per 1 SD among endpoint-evaluable patients', 'Status': 'OK', 'InputStatus': 'Conditional_processed_scale' if acc == 'GSE93624' else 'Supported_input_scale'}
                row.update({'ScoreSD': np.std(x, ddof=1) if len(x) > 1 else np.nan, 'ScoreIQR': np.subtract(*np.quantile(x, [.75, .25])) if len(x) else np.nan, 'ExactZeroN': int(np.sum(x == 0)), 'NearZeroN_1e8': int(np.sum(x <= 1e-8)), 'UniqueScores': len(np.unique(x)), 'Median_event': np.median(x[y == 1]) if n1 else np.nan, 'Median_nonevent': np.median(x[y == 0]) if n0 else np.nan})
                row['Median_difference'] = row['Median_event'] - row['Median_nonevent']
                if n1 < 3 or n0 < 3:
                    row['Status'] = 'Descriptive_only: fewer than 3 patients in one outcome group'
                elif not np.isfinite(x).all() or np.std(x, ddof=1) <= 1e-12:
                    row['Status'] = 'Not_evaluable: nonfinite or constant score'
                else:
                    u = mannwhitneyu(x[y == 1], x[y == 0], alternative='two-sided', method='asymptotic', use_continuity=True)
                    row['P_Wilcoxon'] = u.pvalue
                    row['Rank_biserial'] = 2 * u.statistic / (n1 * n0) - 1
                    z = (x - np.mean(x)) / np.std(x, ddof=1)
                    X = np.column_stack([np.ones(len(x)), z])
                    try:
                        row.update(firth_fit(X, y))
                    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as e:
                        row['Status'] = 'Firth_failed: ' + str(e)
                results.append(row)
            table = pd.DataFrame(results)
            for p, q in [('P_Wilcoxon', 'Q_Wilcoxon'), ('P_Firth', 'Q_Firth')]:
                if p not in table:
                    table[p] = np.nan
                table[q] = bh_family(table[p], size)
            cohort_results.append(table)
            print(endpoint, family, 'N', len(y), 'events', n1, 'fit', table.Status.value_counts().to_dict(), 'q<0.05', int(table.Q_Firth.lt(.05).sum()), flush=True)
    cohort_table = pd.concat(cohort_results, ignore_index=True)
    cohort_table.loc[cohort_table.Family.eq('Ecotype')].to_csv(out / folder / 'ecotype_associations_complete.tsv', sep='\t', index=False)
    cohort_table.loc[cohort_table.Family.eq('State')].to_csv(out / folder / 'cell_state_associations_complete.tsv', sep='\t', index=False)
    all_results.append(cohort_table)

    # RISK的27事件支持预设的年龄/性别敏感性；不按主分析显著性筛选生态型。
    if acc == 'GSE93624':
        adjusted = []
        d = manifest.loc[manifest.Disease.eq('CD') & manifest.Outcome.notna() & manifest.RecoveryComplete].dropna(subset=['AgeAtDiagnosis', 'SexBinary']).copy()
        age = (d.AgeAtDiagnosis.to_numpy() - d.AgeAtDiagnosis.mean()) / d.AgeAtDiagnosis.std(ddof=1)
        for feature in eco.index:
            x = eco.loc[feature, d.SampleID].to_numpy()
            z = (x - x.mean()) / x.std(ddof=1)
            X = np.column_stack([np.ones(len(d)), z, age, d.SexBinary.to_numpy()])
            row = {'Endpoint': 'RISK_3yr_B2_B3', 'Feature': feature, 'N': len(d), 'Events': int(d.Outcome.sum()), 'Adjustment': 'Baseline age and sex; prespecified all 7 ecotypes', 'InputStatus': 'Conditional_processed_scale', 'Status': 'OK'}
            try:
                row.update(firth_fit(X, d.Outcome.to_numpy()))
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as e:
                row['Status'] = str(e)
            adjusted.append(row)
        adjusted = pd.DataFrame(adjusted)
        adjusted['Q_Firth'] = bh_family(adjusted.P_Firth, 7)
        adjusted.to_csv(out / folder / 'ecotype_age_sex_sensitivity.tsv', sep='\t', index=False)

# Step4. 完整保留所有终点，禁止合并不同药物和时间窗的P值。
combined = pd.concat([pd.read_csv(out / folder / name, sep='\t') for folder in ['01_GSE72819_treatment_outcome', '02_RISK_disease_progression', '03_GSE73661_mucosal_healing'] for name in ['ecotype_associations_complete.tsv', 'cell_state_associations_complete.tsv'] if (out / folder / name).exists()], ignore_index=True)
combined.to_csv(out / '04_cross_cohort_summary/all_associations.tsv', sep='\t', index=False)
denominator_path = out / '04_cross_cohort_summary/endpoint_denominators.tsv'
denominators = pd.DataFrame(all_denominators)
if denominator_path.exists():
    previous = pd.read_csv(denominator_path, sep='\t')
    denominators = pd.concat([previous.loc[~previous.Accession.isin(denominators.Accession)], denominators], ignore_index=True)
denominators.to_csv(denominator_path, sep='\t', index=False)
(out / 'logs/STATISTICS_COMPLETED').write_text('Completed all available cohorts; inspect InputStatus and Status before interpretation.\n')
