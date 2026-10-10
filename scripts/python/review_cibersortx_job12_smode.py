# 审阅 Job12 S-mode，用同一批 28 个开发集混合样本与 Job6 未校正结果配对比较。
# 在本机 cibersortx_review 环境运行；读网页输出和冻结真值，写 Job12_Smode_review，不修改参考或留出集。
# 主真值为混合输入的 UMI/RNA 贡献比例；软件输出的表达 Correlation/RMSE 不作为比例准确率。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pypdf import PdfReader
import pypdfium2 as pdfium

# Step1. 保存原始交付文件，并按样本与类型名称匹配
base = Path(r"D:\IBD_ecotype\IBD_ecotype_analysis\CIBERSORTx_signature_audit\main_candidate_validation_20261006")
out = base / "Job12_Smode_review"
original = out / "original_outputs"
qa = out / "QA"
original.mkdir(parents=True, exist_ok=True)
qa.mkdir(exist_ok=True)
downloads = Path(r"C:\Users\马泽凯\Downloads")
source_files = [downloads / "CIBERSORTx_Job12_Adjusted.txt", downloads / "CIBERSORTx_Job12_Results.pdf", downloads / "CIBERSORTx_Job12_stackedChart.pdf", Path(r"C:\Users\马泽凯\.codex\attachments\3b1ed6a6-2e80-4322-a177-1a437af8d6ff\已粘贴的文本.txt")]
integrity = []
for path in source_files:
    name = "CIBERSORTx_Job12_runtime_log.txt" if path.suffix == ".txt" and path.name != "CIBERSORTx_Job12_Adjusted.txt" else path.name
    destination = original / name
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if destination.exists():
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, "已有归档内容不同，不覆盖"
    else:
        shutil.copyfile(path, destination)
    integrity.append({"source": str(path), "archived_file": name, "sha256": digest, "bytes": path.stat().st_size})
pd.DataFrame(integrity).to_csv(out / "original_output_integrity.tsv", sep="\t", index=False)

truth_dir = base / "server_outputs" / "validation_joint20_20261006"
truth = pd.read_csv(truth_dir / "ground_truth_export_UMI_fractions.tsv", sep="\t", index_col=0)
cell_truth = pd.read_csv(truth_dir / "ground_truth_cell_fractions.tsv", sep="\t", index_col=0).loc[truth.index, truth.columns]
design = pd.read_csv(truth_dir / "mixture_design.tsv", sep="\t").set_index("sample").loc[truth.index]
baseline = pd.read_csv(base / "Job6_review" / "original_outputs" / "CIBERSORTx_Job6_Results.txt", sep="\t", index_col=0)
adjusted = pd.read_csv(original / "CIBERSORTx_Job12_Adjusted.txt", sep="\t", index_col=0)
assert baseline.index.is_unique and adjusted.index.is_unique and truth.index.is_unique
assert set(baseline.index) == set(adjusted.index) == set(truth.index)
assert set(baseline.columns) == set(adjusted.columns) == set(truth.columns) | {"P-value", "Correlation", "RMSE"}
before = baseline.loc[truth.index, truth.columns]
after = adjusted.loc[truth.index, truth.columns]
assert before.shape == after.shape == truth.shape == (28, 20)
assert np.isfinite(after.to_numpy()).all() and after.ge(0).all().all() and after.le(1).all().all()
assert np.allclose(before.sum(axis=1), 1) and np.allclose(after.sum(axis=1), 1) and np.allclose(truth.sum(axis=1), 1)
print("Aligned: 28 development samples x 20 types; Job6 versus Job12", flush=True)

# Step2. 同样本配对误差；百分点以 RNA/UMI 贡献真值为主，细胞数量真值另列
before_error = before - truth
after_error = after - truth
sample_metrics = design[["region", "scenario", "absent_type", "n_cells", "n_donors"]].copy()
sample_metrics["Job6_RNA_MAE_pp"] = before_error.abs().mean(axis=1) * 100
sample_metrics["Job12_RNA_MAE_pp"] = after_error.abs().mean(axis=1) * 100
sample_metrics["MAE_change_pp"] = sample_metrics.Job12_RNA_MAE_pp - sample_metrics.Job6_RNA_MAE_pp
sample_metrics["Job6_RNA_TV"] = before_error.abs().sum(axis=1) / 2
sample_metrics["Job12_RNA_TV"] = after_error.abs().sum(axis=1) / 2
sample_metrics["Job6_cell_MAE_pp"] = (before - cell_truth).abs().mean(axis=1) * 100
sample_metrics["Job12_cell_MAE_pp"] = (after - cell_truth).abs().mean(axis=1) * 100
sample_metrics["Job6_fraction_pearson_RNA"] = before.corrwith(truth, axis=1)
sample_metrics["Job12_fraction_pearson_RNA"] = after.corrwith(truth, axis=1)
sample_metrics["Job6_expression_correlation"] = baseline.loc[truth.index, "Correlation"]
sample_metrics["Job12_expression_correlation"] = adjusted.loc[truth.index, "Correlation"]
sample_metrics.to_csv(out / "paired_sample_metrics.tsv", sep="\t")

group_metrics = []
type_metrics = []
groups = {"All28": pd.Series(True, index=truth.index), "Ileum": design.region.eq("Ileum"), "ColonRectum": design.region.eq("ColonRectum"), "Joint": design.region.eq("Joint")}
for group, mask in groups.items():
    before_group = before.loc[mask]
    after_group = after.loc[mask]
    truth_group = truth.loc[mask]
    mae_before = float((before_group - truth_group).abs().to_numpy().mean() * 100)
    mae_after = float((after_group - truth_group).abs().to_numpy().mean() * 100)
    group_metrics.append({"group": group, "samples": int(mask.sum()), "cell_types_in_MAE": 20, "Job6_RNA_MAE_pp": mae_before, "Job12_RNA_MAE_pp": mae_after, "MAE_change_pp": mae_after - mae_before, "relative_MAE_reduction_percent": (mae_before - mae_after) / mae_before * 100, "samples_improved": int(sample_metrics.loc[mask, "MAE_change_pp"].lt(0).sum())})
    for cell_type in truth.columns:
        true = truth_group[cell_type]
        old = before_group[cell_type]
        new = after_group[cell_type]
        type_metrics.append({"group": group, "cell_type": cell_type, "samples": int(mask.sum()), "truth_mean_percent": true.mean() * 100, "Job6_mean_percent": old.mean() * 100, "Job12_mean_percent": new.mean() * 100, "Job6_MAE_pp": (old - true).abs().mean() * 100, "Job12_MAE_pp": (new - true).abs().mean() * 100, "Job6_zero_when_present": int((old.eq(0) & true.gt(0)).sum()), "Job12_zero_when_present": int((new.eq(0) & true.gt(0)).sum()), "true_present_samples": int(true.gt(0).sum())})
group_metrics = pd.DataFrame(group_metrics)
type_metrics = pd.DataFrame(type_metrics)
type_metrics["MAE_change_pp"] = type_metrics.Job12_MAE_pp - type_metrics.Job6_MAE_pp
group_metrics.to_csv(out / "paired_group_metrics.tsv", sep="\t", index=False)
type_metrics.to_csv(out / "paired_group_type_metrics.tsv", sep="\t", index=False)
print(group_metrics.to_string(index=False), flush=True)

# 对 560 个比例保留完整对照，缺失类型单列，不以总均值掩盖大类型偏差
comparisons = []
for sample in truth.index:
    for cell_type in truth.columns:
        comparisons.append({"sample": sample, "region": design.loc[sample, "region"], "scenario": design.loc[sample, "scenario"], "cell_type": cell_type, "RNA_truth": truth.loc[sample, cell_type], "cell_truth": cell_truth.loc[sample, cell_type], "Job6": before.loc[sample, cell_type], "Job12": after.loc[sample, cell_type]})
pd.DataFrame(comparisons).to_csv(out / "all_560_paired_fraction_comparisons.tsv", sep="\t", index=False)
absent = []
for sample, row in design.loc[design.absent_type.notna()].iterrows():
    cell_type = row.absent_type
    assert truth.loc[sample, cell_type] == 0
    absent.append({"sample": sample, "absent_type": cell_type, "truth_percent": 0, "Job6_percent": before.loc[sample, cell_type] * 100, "Job12_percent": after.loc[sample, cell_type] * 100})
pd.DataFrame(absent).to_csv(out / "designed_absent_type_predictions.tsv", sep="\t", index=False)

# Step3. 绘制同真值对照；区域均值只作描述，不将模拟重复当作独立患者
focus = ["Paneth", "Stem_TA_progenitor", "Goblet", "Enteroendocrine", "Ileal_absorptive", "Colonic_absorptive"]
fig, axes = plt.subplots(1, 2, figsize=(13, 5.3), sharex=True)
for axis, region in zip(axes, ["Ileum", "ColonRectum"]):
    table = type_metrics.loc[type_metrics.group.eq(region)].set_index("cell_type").loc[focus]
    position = np.arange(len(focus))
    axis.barh(position - 0.24, table.truth_mean_percent, height=0.23, label="Known RNA contribution", color="#333333")
    axis.barh(position, table.Job6_mean_percent, height=0.23, label="Job6: no correction", color="#B6BDC8")
    axis.barh(position + 0.24, table.Job12_mean_percent, height=0.23, label="Job12: S-mode", color="#008A86")
    axis.set_yticks(position, [name.replace("_", " ") for name in focus], fontsize=9)
    axis.invert_yaxis()
    axis.set_title(f"{region}: 12 matched mixtures")
    axis.set_xlabel("Mean fraction (%)")
    axis.set_xlim(0, 70)
    axis.grid(axis="x", alpha=0.2)
    axis.set_axisbelow(True)
axes[0].legend(fontsize=8, loc="lower right")
fig.suptitle("Joint20 development validation: S-mode retains major ileal bias")
fig.tight_layout()
fig.savefig(out / "Job12_vs_Job6_key_types.png", dpi=180)
plt.close(fig)

# Step4. 原 PDF 只读提取与渲染，保留完整页面用于查看图例与样本顺序
for path in [original / "CIBERSORTx_Job12_Results.pdf", original / "CIBERSORTx_Job12_stackedChart.pdf"]:
    reader = PdfReader(path)
    extracted = [f"PAGE {i + 1}\n{page.extract_text() or ''}" for i, page in enumerate(reader.pages)]
    (qa / f"{path.stem}_text.txt").write_text("\n\n".join(extracted), encoding="utf-8")
    document = pdfium.PdfDocument(str(path))
    for i in range(len(document)):
        page = document[i]
        bitmap = page.render(scale=1.6)
        bitmap.to_pil().save(qa / f"{path.stem}_page{i + 1}.png")
        bitmap.close()
        page.close()
    document.close()

# Step5. 保存解释边界：成功完成不等于已验证；本轮没有 MSCCR 或 HiRes 输出
log = (original / "CIBERSORTx_Job12_runtime_log.txt").read_text(encoding="utf-8-sig")
assert "rmbatchSmode: TRUE" in log and "All done." in log
summary = {"status": "review_completed", "job": "Job12", "baseline": "Job6", "samples": 28, "cell_types": 20, "truth": "export UMI contribution proxy for RNA contribution, not cell counts", "comparison": group_metrics.to_dict(orient="records"), "Job12_zero_reported_P_values": int(adjusted['P-value'].eq(0).sum()), "Job12_expression_correlation_median": float(adjusted.Correlation.median()), "Job6_expression_correlation_median": float(baseline.Correlation.median()), "explicit_error_or_warning_in_pasted_log": any(term in log.lower() for term in ["error:", "warning:"]), "adjusted_signature_file_provided": False, "output_is_fraction_table_not_HiRes_G": True, "final_holdout_used": False, "MSCCR_computed_by_this_review": False, "decision": "small overall improvement, persistent large ileal bias; not accepted as validated Joint20 reference; S-mode effect on real MSCCR remains untested"}
(out / "review_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print("Saved:", out, flush=True)
