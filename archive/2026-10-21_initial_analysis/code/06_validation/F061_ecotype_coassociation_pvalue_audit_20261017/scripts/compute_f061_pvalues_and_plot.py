# 补算F061外部Ecotype共现的单侧经验P值、正态近似P值和14检验BH q值。
# 在Windows项目本机运行；只读既有统计量与零分布，写入新的审计目录，不修改冻结模型或原始产物。
# 输入为两队列的14行汇总、各1000次零分布、成员覆盖及原F061图；输出TSV、中文报告和增强图。
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path.home() / ".codex/skills/nature-figure/scripts"))
from audit_panel_alignment import require_matplotlib_panel_alignment


# Step1. 参数和输出目录
parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(r"D:\IBD_ecotype\IBD_ecotype_analysis"))
parser.add_argument("--output-dir", type=Path, default=None)
args = parser.parse_args()

root = args.project_root.resolve()
source_root = root / "06_validation/paper_aligned_recovery_characterization_20261012"
coassoc_dir = source_root / "05_ecotype_coassociation"
original_figure_dir = source_root / "08_figures"
out = args.output_dir.resolve() if args.output_dir else root / "06_validation/F061_ecotype_coassociation_pvalue_audit_20261017"
statistics_dir = out / "statistics"
figures_dir = out / "figures"
reports_dir = out / "reports"
scripts_dir = out / "scripts"
source_snapshot_dir = out / "source_snapshot"
qa_dir = out / "qa"
for path in [statistics_dir, figures_dir, reports_dir, scripts_dir, source_snapshot_dir, qa_dir]:
    path.mkdir(parents=True, exist_ok=True)

summary_path = coassoc_dir / "ecotype_coassociation_1000permutations.tsv"
coverage_path = coassoc_dir / "ecotype_member_coverage.tsv"
null_paths = {
    "Colon_GSE214695": coassoc_dir / "Colon_GSE214695_permutation_null.tsv",
    "Ileum_IBDverse": coassoc_dir / "Ileum_IBDverse_permutation_null.tsv",
}
matrix_paths = {
    dataset: {
        "mean": coassoc_dir / f"{dataset}_C_four_matrix_mean.tsv",
        "F1": coassoc_dir / f"{dataset}_C_F1.tsv",
        "F2": coassoc_dir / f"{dataset}_C_F2.tsv",
        "F1_nozero": coassoc_dir / f"{dataset}_C_F1_nozero.tsv",
        "F2_nozero": coassoc_dir / f"{dataset}_C_F2_nozero.tsv",
        "fractions": coassoc_dir / f"{dataset}_F1_assigned_denominator.tsv",
    }
    for dataset in ["Colon_GSE214695", "Ileum_IBDverse"]
}
frozen_states_path = source_root / "00_audit/frozen_states.tsv"
original_paths = [
    original_figure_dir / "Figure_E_source.tsv",
    original_figure_dir / "Figure_E_ecotype_Z.png",
    original_figure_dir / "Figure_E_ecotype_Z.pdf",
    original_figure_dir / "Figure_E_ecotype_Z.svg",
    original_figure_dir / "Figure_D_Colon_GSE214695_source.tsv",
    original_figure_dir / "Figure_D_Colon_GSE214695.png",
    original_figure_dir / "Figure_D_Colon_GSE214695.pdf",
    original_figure_dir / "Figure_D_Colon_GSE214695.svg",
    original_figure_dir / "Figure_D_Ileum_IBDverse_source.tsv",
    original_figure_dir / "Figure_D_Ileum_IBDverse.png",
    original_figure_dir / "Figure_D_Ileum_IBDverse.pdf",
    original_figure_dir / "Figure_D_Ileum_IBDverse.svg",
]
matrix_input_paths = [path for dataset_paths in matrix_paths.values() for path in dataset_paths.values()]
input_paths = [summary_path, coverage_path, frozen_states_path, *null_paths.values(), *matrix_input_paths, *original_paths]
for path in input_paths:
    if not path.exists():
        raise FileNotFoundError(path)


# Step2. 记录输入哈希并复制只读快照
manifest_rows = []
for path in input_paths:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_rows.append({"SourcePath": str(path), "Bytes": path.stat().st_size, "SHA256": digest})
    shutil.copy2(path, source_snapshot_dir / path.name)
input_manifest = pd.DataFrame(manifest_rows)
input_manifest.to_csv(reports_dir / "input_manifest_sha256.tsv", sep="\t", index=False)


# Step3. 读入14行汇总与成员覆盖，逐行从实际零分布计数
summary = pd.read_csv(summary_path, sep="\t")
coverage = pd.read_csv(coverage_path, sep="\t")
datasets = ["Colon_GSE214695", "Ileum_IBDverse"]
ecotypes = [f"E{i}" for i in range(1, 8)]
assert summary.shape[0] == 14
assert set(summary["Dataset"]) == set(datasets)
assert set(summary["Ecotype"]) == set(ecotypes)
assert not summary.duplicated(["Dataset", "Ecotype"]).any()
assert not coverage.duplicated(["Dataset", "Ecotype"]).any()

result_rows = []
normality_rows = []
for dataset in datasets:
    null_table = pd.read_csv(null_paths[dataset], sep="\t")
    assert null_table.shape == (1000, 7)
    assert list(null_table.columns) == ecotypes
    for ecotype in ecotypes:
        row = summary.loc[summary["Dataset"].eq(dataset) & summary["Ecotype"].eq(ecotype)].iloc[0]
        cov = coverage.loc[coverage["Dataset"].eq(dataset) & coverage["Ecotype"].eq(ecotype)].iloc[0]
        null_series = pd.to_numeric(null_table[ecotype], errors="coerce")
        assert null_series.notna().all()
        null = null_series.to_numpy(dtype=float)
        assert len(null) == int(row["ValidPermutations"]) == 1000
        recomputed_mean = float(np.mean(null))
        recomputed_sd = float(np.std(null, ddof=1))
        recomputed_z = float((row["AvgDif"] - recomputed_mean) / recomputed_sd)
        np.testing.assert_allclose(recomputed_mean, row["NullMean"], rtol=0, atol=1e-12)
        np.testing.assert_allclose(recomputed_sd, row["NullSD"], rtol=0, atol=1e-12)
        np.testing.assert_allclose(recomputed_z, row["Z"], rtol=0, atol=1e-12)

        exceedances = int(np.sum(null >= float(row["AvgDif"])))
        p_empirical = (1 + exceedances) / (len(null) + 1)
        p_normal = float(stats.norm.sf(float(row["Z"])))
        missing_members = "" if pd.isna(cov["MissingMembers"]) else str(cov["MissingMembers"])
        missing_count = 0 if missing_members == "" else len(missing_members.split(";"))
        evidence_direction = "positive co-association" if row["AvgDif"] > row["NullMean"] else "opposite/no positive support"
        result_rows.append({
            "Dataset": dataset,
            "Ecotype": ecotype,
            "FrozenMembers": int(row["FrozenMembers"]),
            "InputLineageMembers": int(row["InputLineageMembers"]),
            "DefinedReferenceMembers": int(row["DefinedReferenceMembers"]),
            "MissingMemberCount": missing_count,
            "MissingMembers": missing_members,
            "AvgDif": float(row["AvgDif"]),
            "NullMean": recomputed_mean,
            "NullSD": recomputed_sd,
            "Z": recomputed_z,
            "ValidPermutations": len(null),
            "ExceedancesGreaterOrEqual": exceedances,
            "P_empirical_greater": p_empirical,
            "P_normal_approx_greater": p_normal,
            "Tail": "greater; stronger within-ecotype co-association",
            "EvidenceDirection": evidence_direction,
            "EmpiricalPFormula": "(1 + count(null >= AvgDif)) / (1000 + 1)",
            "EmpiricalPResolution": "1/1001 = 0.000999000999",
            "BHFamily": "2 external cohorts x 7 ecotypes = 14 tests",
            "PaperSpecificCorrection": "none stated for CE-specific P values",
            "GlobalTest": "not computed; distinct joint 1,000,000-permutation criterion",
        })
        normality_rows.append({
            "Dataset": dataset,
            "Ecotype": ecotype,
            "N": len(null),
            "Skewness": float(stats.skew(null, bias=False)),
            "ExcessKurtosis": float(stats.kurtosis(null, fisher=True, bias=False)),
            "ShapiroW": float(stats.shapiro(null).statistic),
            "ShapiroP": float(stats.shapiro(null).pvalue),
            "NormalApproximationUse": "sensitivity only; empirical P is primary",
        })

results = pd.DataFrame(result_rows)
results["Q_BH14_empirical"] = stats.false_discovery_control(results["P_empirical_greater"].to_numpy(), method="bh")
results["Q_BH14_normal_approx"] = stats.false_discovery_control(results["P_normal_approx_greater"].to_numpy(), method="bh")
column_order = [
    "Dataset", "Ecotype", "FrozenMembers", "InputLineageMembers", "DefinedReferenceMembers",
    "MissingMemberCount", "MissingMembers", "AvgDif", "NullMean", "NullSD", "Z",
    "ValidPermutations", "ExceedancesGreaterOrEqual", "P_empirical_greater", "Q_BH14_empirical",
    "P_normal_approx_greater", "Q_BH14_normal_approx", "Tail", "EvidenceDirection",
    "EmpiricalPFormula", "EmpiricalPResolution", "BHFamily", "PaperSpecificCorrection", "GlobalTest",
]
results = results[column_order]
results.to_csv(statistics_dir / "F061_ecotype_coassociation_pvalues_14tests.tsv", sep="\t", index=False, float_format="%.12g")
pd.DataFrame(normality_rows).to_csv(statistics_dir / "F061_null_normality_diagnostics.tsv", sep="\t", index=False, float_format="%.12g")
results.to_csv(figures_dir / "F061_enhanced_source.tsv", sep="\t", index=False, float_format="%.12g")
print(results[["Dataset", "Ecotype", "Z", "ExceedancesGreaterOrEqual", "P_empirical_greater", "Q_BH14_empirical"]].to_string(index=False))


# Step4. 绘制F061增强图；填充标记表示14检验BH q<0.05
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})
colors = {"Colon_GSE214695": "#2F6F8F", "Ileum_IBDverse": "#B06F40"}
labels = {"Colon_GSE214695": "Colon / GSE214695", "Ileum_IBDverse": "Ileum / IBDverse"}
offsets = {"Colon_GSE214695": -0.12, "Ileum_IBDverse": 0.12}

fig, ax = plt.subplots(figsize=(7.2, 3.9))
fig.subplots_adjust(left=0.095, right=0.975, bottom=0.22, top=0.81)
ax.axvspan(-5.0, 8.6, color="#F7F7F7", zorder=0)
ax.axvline(0, color="#888888", linewidth=0.8, zorder=1)
ax.axvline(8.9, color="#D0D0D0", linewidth=0.7, zorder=1)
for dataset in datasets:
    data = results.loc[results["Dataset"].eq(dataset)].set_index("Ecotype").reindex(ecotypes)
    y = np.arange(7) + offsets[dataset]
    significant = data["Q_BH14_empirical"].to_numpy() < 0.05
    ax.scatter(data.loc[~significant, "Z"], y[~significant], s=34, marker="o", facecolors="white", edgecolors=colors[dataset], linewidths=1.2, label=labels[dataset], zorder=3)
    ax.scatter(data.loc[significant, "Z"], y[significant], s=34, marker="o", facecolors=colors[dataset], edgecolors=colors[dataset], linewidths=1.0, zorder=4)

for i, ecotype in enumerate(ecotypes):
    colon = results.loc[results["Dataset"].eq("Colon_GSE214695") & results["Ecotype"].eq(ecotype)].iloc[0]
    ileum = results.loc[results["Dataset"].eq("Ileum_IBDverse") & results["Ecotype"].eq(ecotype)].iloc[0]
    colon_text = f"{colon.P_empirical_greater:.4f} / {colon.Q_BH14_empirical:.4f} / {colon.DefinedReferenceMembers}/{colon.FrozenMembers}"
    ileum_text = f"{ileum.P_empirical_greater:.4f} / {ileum.Q_BH14_empirical:.4f} / {ileum.DefinedReferenceMembers}/{ileum.FrozenMembers}"
    ax.text(9.25, i - 0.12, colon_text, color=colors["Colon_GSE214695"], fontsize=6.25, va="center")
    ax.text(13.05, i + 0.12, ileum_text, color=colors["Ileum_IBDverse"], fontsize=6.25, va="center")

ax.set_xlim(-5.0, 17.4)
ax.set_ylim(6.55, -0.95)
ax.set_xticks([-4, 0, 4, 8])
ax.set_yticks(np.arange(7), labels=ecotypes)
ax.set_xlabel("Ecotype co-association Z (AvgDif vs 1,000 row-permutation nulls)")
ax.tick_params(axis="y", length=0)
ax.set_title("F061 | External ecotype co-association with recalculated P and q", fontsize=9, loc="left", pad=28)
ax.text(9.25, -0.62, "Colon\nPemp / qBH14 / coverage", fontsize=6.2, fontweight="bold", va="center", linespacing=1.15)
ax.text(13.05, -0.62, "Ileum\nPemp / qBH14 / coverage", fontsize=6.2, fontweight="bold", va="center", linespacing=1.15)
ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=6.4, handletextpad=0.35, columnspacing=1.2)
fig.text(0.095, 0.085, "Greater-tail empirical P = (1 + null >= observed)/(1,001); q = BH across 14 tests. Filled markers: q < 0.05.", fontsize=6.2)
fig.text(0.095, 0.045, "Minimum resolvable empirical P = 0.000999. Normal-approximation P values are retained in the source table/report only.", fontsize=6.2, color="#555555")

figure_prefix = figures_dir / "F061_enhanced_ecotype_coassociation_P_and_q"
require_matplotlib_panel_alignment(fig, axes=[ax], json_out=str(figure_prefix) + ".alignment.json", overlay_svg=str(figure_prefix) + ".alignment.svg", tolerance_pt=1.5, strict=True)
fig.savefig(str(figure_prefix) + ".pdf")
fig.savefig(str(figure_prefix) + ".svg")
fig.savefig(str(figure_prefix) + ".png", dpi=600)
plt.close(fig)


# Step4b. 按Cell 2021 Figure 4D结构绘制两队列C矩阵；仅借用版式，不改变本项目统计定义
frozen_states = pd.read_csv(frozen_states_path, sep="\t")
assert frozen_states.shape[0] == 53
assert frozen_states["Ecotype"].tolist() == sorted(frozen_states["Ecotype"].tolist(), key=lambda x: int(x[1:]))
state_to_ecotype = frozen_states.set_index("CellState")["Ecotype"].to_dict()
ecotype_colors = {
    "E1": "#9E1B32",
    "E2": "#D9A21B",
    "E3": "#7FBF7B",
    "E4": "#5AAEAD",
    "E5": "#6B8EBD",
    "E6": "#8B6AAE",
    "E7": "#D67A2E",
}
dataset_titles = {
    "Colon_GSE214695": "Colon / GSE214695",
    "Ileum_IBDverse": "Ileum / IBDverse",
}

paper_style_prefix = figures_dir / "F061_paper_style_external_ecotype_cooccurrence"
annotation_rows = []
fig = plt.figure(figsize=(8.4, 4.65))
grid = fig.add_gridspec(1, 2, left=0.075, right=0.805, bottom=0.205, top=0.79, wspace=0.18)
heatmap_axes = []
image = None
for panel_index, dataset in enumerate(datasets):
    ax = fig.add_subplot(grid[0, panel_index])
    heatmap_axes.append(ax)
    signed_c = pd.read_csv(matrix_paths[dataset]["mean"], sep="\t", index_col=0)
    component_matrices = [
        pd.read_csv(matrix_paths[dataset][name], sep="\t", index_col=0)
        for name in ["F1", "F2", "F1_nozero", "F2_nozero"]
    ]
    for component in component_matrices:
        assert signed_c.index.equals(component.index)
        assert signed_c.columns.equals(component.columns)
    pre_na_c = (component_matrices[0] + component_matrices[1] + component_matrices[2] + component_matrices[3]) / 4
    finite = np.isfinite(signed_c.to_numpy())
    assert np.allclose(signed_c.to_numpy()[finite], pre_na_c.to_numpy()[finite], atol=1e-12)
    assert np.allclose(signed_c.to_numpy(), signed_c.to_numpy().T, equal_nan=True, atol=1e-12)

    present = [state for state in frozen_states["CellState"] if state in signed_c.index]
    assert set(present) == set(signed_c.index)
    signed_ordered = signed_c.loc[present, present].copy()
    pre_na_ordered = pre_na_c.loc[present, present]
    np.fill_diagonal(signed_ordered.values, np.diag(pre_na_ordered.to_numpy()))
    display_ordered = signed_ordered.clip(lower=0, upper=1)
    signed_ordered.to_csv(figures_dir / f"F061_paper_style_{dataset}_C_signed_with_display_diagonal.tsv", sep="\t", float_format="%.12g")
    display_ordered.to_csv(figures_dir / f"F061_paper_style_{dataset}_C_display_0_to_1.tsv", sep="\t", float_format="%.12g")

    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#D9D9D9")
    image = ax.imshow(display_ordered.to_numpy(), cmap=cmap, vmin=0, vmax=1, interpolation="none", aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])
    sample_count = pd.read_csv(matrix_paths[dataset]["fractions"], sep="\t", index_col=0).shape[1]
    state_count = len(present)
    ax.set_xlabel(f"{state_count} cell states", fontsize=8, labelpad=4)
    ax.set_ylabel(f"{state_count} cell states", fontsize=8, labelpad=4)
    ax.set_title(f"{dataset_titles[dataset]}\n(n = {sample_count} samples)", fontsize=9, pad=28)

    start = 0
    dataset_results = results.loc[results["Dataset"].eq(dataset)].set_index("Ecotype")
    for ecotype in ecotypes:
        block_states = [state for state in present if state_to_ecotype[state] == ecotype]
        width = len(block_states)
        assert width > 0
        center = start + (width - 1) / 2
        row = dataset_results.loc[ecotype]
        p_value = float(row["P_empirical_greater"])
        q_value = float(row["Q_BH14_empirical"])
        symbol = "*" if p_value < 0.05 else "+" if p_value < 0.10 else ""
        ax.add_patch(Rectangle((start - 0.5, start - 0.5), width, width, fill=False, edgecolor="black", linewidth=0.85, zorder=4))
        ax.add_patch(Rectangle((start - 0.5, -2.10), width, 0.72, facecolor=ecotype_colors[ecotype], edgecolor="none", clip_on=False, zorder=5))
        ax.text(center, -2.75, ecotype, ha="center", va="bottom", fontsize=6.3, fontweight="bold" if q_value < 0.05 else "normal", clip_on=False)
        if symbol:
            ax.text(center, -4.00, symbol, ha="center", va="bottom", fontsize=9, fontweight="bold", clip_on=False)
        annotation_rows.append({
            "Dataset": dataset,
            "Ecotype": ecotype,
            "MatrixStatesInBlock": width,
            "DefinedReferenceMembers": int(row["DefinedReferenceMembers"]),
            "FrozenMembers": int(row["FrozenMembers"]),
            "P_empirical_greater": p_value,
            "Q_BH14_empirical": q_value,
            "PaperStyleSymbol": symbol,
            "EcotypeLabelBold": bool(q_value < 0.05),
        })
        start += width
    assert start == state_count
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color("black")

pd.DataFrame(annotation_rows).to_csv(figures_dir / "F061_paper_style_annotation_source.tsv", sep="\t", index=False, float_format="%.12g")

fig.text(0.018, 0.945, "F061", fontsize=16, fontweight="bold", va="top")
fig.text(0.425, 0.960, "External ecotype co-occurrence in scRNA-seq", fontsize=11, ha="center", va="top")
legend_x = 0.835
fig.text(legend_x, 0.765, "Ecotype significance", fontsize=7.2, fontweight="bold")
fig.text(legend_x, 0.716, "*  Pemp < 0.05", fontsize=7)
fig.text(legend_x, 0.674, "+  0.05 <= Pemp < 0.10", fontsize=7)
fig.text(legend_x, 0.612, "Bold E label", fontsize=7, fontweight="bold")
fig.text(legend_x, 0.572, "qBH14 < 0.05", fontsize=7)
fig.text(legend_x, 0.486, "Defined / frozen", fontsize=7.2, fontweight="bold")
coverage_heading_y = {"Colon_GSE214695": 0.443, "Ileum_IBDverse": 0.330}
coverage_text_y = {"Colon_GSE214695": 0.405, "Ileum_IBDverse": 0.292}
for dataset in datasets:
    data = results.loc[results["Dataset"].eq(dataset)].set_index("Ecotype").reindex(ecotypes)
    coverage_parts = [f"{e} {int(data.loc[e, 'DefinedReferenceMembers'])}/{int(data.loc[e, 'FrozenMembers'])}" for e in ecotypes]
    coverage_text = "  ".join(coverage_parts[:4]) + "\n" + "  ".join(coverage_parts[4:])
    fig.text(legend_x, coverage_heading_y[dataset], dataset_titles[dataset], fontsize=6.6, fontweight="bold")
    fig.text(legend_x, coverage_text_y[dataset], coverage_text, fontsize=5.8, linespacing=1.35)

color_ax = fig.add_axes([0.318, 0.105, 0.245, 0.025])
colorbar = fig.colorbar(image, cax=color_ax, orientation="horizontal")
colorbar.set_ticks([0, 1], labels=["low", "high"])
colorbar.ax.tick_params(labelsize=6.5, length=0)
fig.text(0.275, 0.116, "Co-occurrence index", fontsize=7.2, ha="right", va="center")
fig.text(0.075, 0.053, "Paper-style display: the four Spearman matrices are averaged into C; C <= 0 is mapped to the low color and the pre-test diagonal is shown.", fontsize=6.0)
fig.text(0.075, 0.025, "Statistics use diagonal = NA. Gray = undefined. Symbols use unadjusted empirical P; bold E labels show the additional BH correction across 14 tests.", fontsize=6.0, color="#444444")
require_matplotlib_panel_alignment(fig, axes=heatmap_axes, json_out=str(paper_style_prefix) + ".alignment.json", overlay_svg=str(paper_style_prefix) + ".alignment.svg", tolerance_pt=1.5, strict=True)
fig.savefig(str(paper_style_prefix) + ".pdf")
fig.savefig(str(paper_style_prefix) + ".svg")
fig.savefig(str(paper_style_prefix) + ".png", dpi=600)
plt.close(fig)


# Step5. 写中文审计报告、运行元数据与原文件未变化验证
sig = results.loc[results["Q_BH14_empirical"] < 0.05, ["Dataset", "Ecotype", "P_empirical_greater", "Q_BH14_empirical"]]
borderline = results.loc[(results["P_empirical_greater"] < 0.05) & (results["Q_BH14_empirical"] >= 0.05), ["Dataset", "Ecotype", "P_empirical_greater", "Q_BH14_empirical"]]
max_abs_skew = pd.DataFrame(normality_rows)["Skewness"].abs().max()
max_abs_kurt = pd.DataFrame(normality_rows)["ExcessKurtosis"].abs().max()

report_lines = [
    "# F061 外部 Ecotype 共现 P 值补算审计报告",
    "",
    "## 结论",
    "",
    "本次从两队列各自保存的 1,000 次逐行置换零分布直接计数，完成 2 队列 x 7 Ecotypes 共 14 个单独检验。主结果采用与研究问题一致的单侧 greater-tail 经验 P；BH-FDR 校正族在计算前固定为全部 14 个检验。原始汇总、零分布、F061 源表和原图均未改写。",
    "",
    f"经验 P 经 14 检验 BH 后共有 {len(sig)} 项 q<0.05：" + "；".join(f"{r.Dataset}/{r.Ecotype} (P={r.P_empirical_greater:.6g}, q={r.Q_BH14_empirical:.6g})" for r in sig.itertuples(index=False)) + "。",
    "",
    ("名义 P<0.05 但未通过 BH14：" + "；".join(f"{r.Dataset}/{r.Ecotype} (P={r.P_empirical_greater:.6g}, q={r.Q_BH14_empirical:.6g})" for r in borderline.itertuples(index=False)) + "。") if len(borderline) else "没有名义 P<0.05 但未通过 BH14 的项目。",
    "",
    "Ileum/IBDverse 的 E4 为负向 AvgDif 和负 Z，greater-tail 经验 P=1，不能解释为共现恢复。覆盖核对为 E1 回肠 9/10、E3 回肠 10/12；缺成员和不可定义成员只在实际存在的父类型与可计算边上解释。",
    "",
    "## Figure 4D 结构复核与新图",
    "",
    "用户提供的第二张参考图对应 Luca et al. Cell 2021 Figure 4D 的共现矩阵结构，而不是 Figure S5H 的 sparsity/depletion 气泡图。新图据此使用两队列实际保存的 C 矩阵：Colon/GSE214695 为 53 x 53，Ileum/IBDverse 为 47 x 47；状态按冻结 E1-E7 顺序排列，黑框标示各 Ecotype 对角块，顶部色带标示成员归属。",
    "",
    "为接近原文的 low-high 顺序色标，图中将 C<=0 映射到最低颜色；带符号的完整 C 数值另存于图源表，未改写。统计检验仍使用对角线为 NA 的 C；热图对角仅恢复四个原始 Spearman 矩阵平均后的检验前数值，用于显示自相关。灰色表示四矩阵严格平均后仍不可定义。",
    "",
    "顶部 * 表示未校正经验单侧 P<0.05，+ 表示 0.05<=P<0.10，沿用原文单个 CE 的展示口径；加粗的 E 标签另表示本项目预先定义的 14 检验 BH q<0.05。因此 Colon/E4 虽有 *，但标签不加粗（P=0.042957，q=0.100233），避免把名义显著误写为 FDR 显著。",
    "",
    "## 方法核对",
    "",
    "- 原文先按每个样本、每个父类型计算状态比例：F1 的分母为可被分配到 EcoTyper 状态的细胞，F2 的分母为该父类型全部输入细胞。",
    "- 对 F1、F2 直接计算两套 Spearman 相关矩阵，再将 0 替换为 NA 后计算两套相关矩阵；四矩阵平均得到 C。",
    "- 对每个 Ecotype 成员状态，计算成员内平均共现减去成员外平均共现，再在成员间平均得到 AvgDif_ref。",
    "- 每次置换独立打乱 C 的每一行，并把原对角 NA 与新对角位置交换；共 1,000 次，得到该 Ecotype 的 AvgDif_null。",
    "- 原文研究问题明确为成员状态是否比随机情况更强地共现，因此本次检验尾为 greater；经验 P=(1+count(null>=AvgDif_ref))/(1000+1)。没有从 Z 反推极端计数。",
    "- 原文给出的单个 CE 统计量为 Z=(AvgDif_ref-mean(null))/sd(null)，并写明 Z 转 P，但在可核对的 STAR Methods 中没有明写 pnorm 的单双尾代码；官方当前公开仓库也未包含 Figure 4D 的四矩阵/置换显著性实现。因此 normal-approximation greater-tail P 仅作为原文口径敏感性，不作为经验置换 P。",
    "- 原文单个 CE 段未声明多重校正；本表同时保留未校正 P，并额外给出预先固定的 BH14 q。",
    "",
    "## 零分布与精度",
    "",
    f"14 个零分布均有 1,000 个有限值。偏度绝对值最大为 {max_abs_skew:.3f}，超额峰度绝对值最大为 {max_abs_kurt:.3f}；主体形态与正态近似大体相容，但 1,000 次置换不能校准极端正态尾。经验 P 的最小可分辨值为 1/1001=0.000999001，零极端计数不表示 P=0。",
    "",
    "## 全局 CE 检验边界",
    "",
    "原文的联合全局检测不是把单个 CE 的 1,000 次 P 值合并。它对每个 CE 的成员内均值使用共同的逐行置换，按每个 CE 的全部成员状态同时满足观测值不小于置换值的条件计数，并运行 1,000,000 次以评估所有 CE 同时检出的概率。当前 7-Ecotype 项目只保存了各 Ecotype 的 1,000 次 AvgDif 零分布，缺少原文联合判定所需的每状态成员内统计量、共同置换轨迹和 1,000,000 次运行，因此本次不计算也不声称全局 P。",
    "",
    "## 14 项结果",
    "",
    "| 队列 | Ecotype | coverage | Z | exceedances | 经验单侧 P | BH14 q | 正态近似单侧 P |",
    "|---|---:|---:|---:|---:|---:|---:|---:|",
]
for row in results.itertuples(index=False):
    report_lines.append(f"| {row.Dataset} | {row.Ecotype} | {row.DefinedReferenceMembers}/{row.FrozenMembers} | {row.Z:.3f} | {row.ExceedancesGreaterOrEqual} | {row.P_empirical_greater:.6g} | {row.Q_BH14_empirical:.6g} | {row.P_normal_approx_greater:.6g} |")
report_lines.extend([
    "",
    "## 运行代码与路径",
    "",
    f"- 脚本：`{scripts_dir / 'compute_f061_pvalues_and_plot.py'}`",
    f"- 14 行主表：`{statistics_dir / 'F061_ecotype_coassociation_pvalues_14tests.tsv'}`",
    f"- 零分布正态性诊断：`{statistics_dir / 'F061_null_normality_diagnostics.tsv'}`",
    f"- 增强图：`{figure_prefix}.png`、`.pdf`、`.svg`",
    f"- Figure 4D 参考结构热图：`{paper_style_prefix}.png`、`.pdf`、`.svg`",
    f"- 热图显著性与覆盖注释源表：`{figures_dir / 'F061_paper_style_annotation_source.tsv'}`",
    f"- 原始输入快照：`{source_snapshot_dir}`",
    f"- 输入哈希：`{reports_dir / 'input_manifest_sha256.tsv'}`",
    "- 运行命令：`D:\\anaconda3\\python.exe <脚本路径> --project-root D:\\IBD_ecotype\\IBD_ecotype_analysis`",
    "",
    "## 来源",
    "",
    "- Luca BA, Steen CB, Matusiak M, et al. Cell 2021. https://doi.org/10.1016/j.cell.2021.09.014",
    "- 作者官方 EcoTyper 仓库（核对快照 commit 11408e71900e15829cf4f33e1331ce47c053472d）：https://github.com/digitalcytometry/ecotyper",
])
(reports_dir / "F061_P值补算审计报告_20261017.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

run_metadata = {
    "Python": sys.version,
    "pandas": pd.__version__,
    "numpy": np.__version__,
    "scipy": stats.__version__ if hasattr(stats, "__version__") else __import__("scipy").__version__,
    "matplotlib": mpl.__version__,
    "Tail": "greater",
    "BHFamily": 14,
    "PermutationCountPerTest": 1000,
    "EmpiricalPMinimum": 1 / 1001,
    "FigureContract": {
        "CoreConclusion": "Archived row-permutation nulls support reproducible positive co-association for a subset of frozen ecotypes across two external atlases.",
        "Role": "external validation",
        "Archetype": "quantitative single-panel",
        "Backend": "Python/matplotlib",
        "Export": ["PNG", "PDF", "SVG"],
        "Correction": "BH across 14 tests",
        "PaperStyleFigure": "two external C-matrix heatmaps ordered by frozen ecotype; style-only adaptation of Cell 2021 Figure 4D",
    },
}
(reports_dir / "run_metadata.json").write_text(json.dumps(run_metadata, ensure_ascii=False, indent=2), encoding="utf-8")

unchanged_rows = []
for row in manifest_rows:
    path = Path(row["SourcePath"])
    digest_after = hashlib.sha256(path.read_bytes()).hexdigest()
    unchanged_rows.append({"SourcePath": str(path), "SHA256Before": row["SHA256"], "SHA256After": digest_after, "Unchanged": digest_after == row["SHA256"]})
unchanged = pd.DataFrame(unchanged_rows)
assert unchanged["Unchanged"].all()
unchanged.to_csv(reports_dir / "source_unchanged_verification.tsv", sep="\t", index=False)

shutil.copy2(Path(__file__), scripts_dir / Path(__file__).name)
print(f"OUTPUT_DIR={out}")
print(f"SIGNIFICANT_BH14={len(sig)}")
print("SOURCE_FILES_UNCHANGED=True")
