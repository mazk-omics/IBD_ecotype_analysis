"""Summarize the five server audit reports copied into the local project."""
import csv
import hashlib
import json
import shutil
from pathlib import Path

root = Path(__file__).resolve().parents[2]
out = root / "IBDverse_reference" / "audit_20261006_stage1"
out.mkdir(parents=True, exist_ok=True)
files = ["REPORT.txt", "h5ad_structure.json", "obs_columns.tsv", "obs_value_counts.tsv", "gene_columns.tsv"]
manifest = []
for name in files:
    source = root / name
    shutil.copy2(source, out / name)
    manifest.append({"file": name, "sha256": hashlib.sha256(source.read_bytes()).hexdigest()})
summary = json.loads((root / "h5ad_structure.json").read_text(encoding="utf-8"))
with (root / "obs_value_counts.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
labels = {r["value"]: int(r["n_cells"]) for r in rows if r["column"] == "label"}
categories = {r["value"]: int(r["n_cells"]) for r in rows if r["column"] == "category"}
assert sum(labels.values()) == sum(categories.values()) == summary["n_cells"]
groups = {
    "Enterocyte_excluding_BEST4": [k for k in labels if k.startswith("Enterocytes") and "BEST4" not in k],
    "BEST4": [k for k in labels if "BEST4" in k],
    "Stem_proliferative": [k for k in labels if k.startswith("Stem cells")],
    "Goblet": [k for k in labels if k.startswith("Goblet")],
    "Tuft": [k for k in labels if k.startswith("Tuft")],
    "Enteroendocrine": ["Enteroendocrine cells"],
    "Enterochromaffin": ["Enterochromaffin cells"],
}
with (out / "epithelial_candidate_groups.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["candidate_group", "n_cells", "author_labels", "status"])
    for group, members in groups.items():
        writer.writerow([group, sum(labels[k] for k in members), " | ".join(members), "proposal_only_not_final_identity_mapping"])
epithelial = sum(categories[k] for k in ["Enterocyte", "Secretory", "Stem cell"])
assert sum(sum(labels[k] for k in members) for members in groups.values()) == epithelial
report = f"""# IBDverse 回肠参考：第一批报告审阅

依据：服务器只读检查结果及本地复制的五份报告。未读取服务器表达值。

- 对象：`{summary['source']}`。
- {summary['n_cells']:,} 细胞 × {summary['n_genes']:,} 基因；343 个细胞注释中的样本 ID。
- 57 个 `label`，8 个 `category`；两套完整类型计数均与总细胞数一致。
- 上皮相关三大类合计 {epithelial:,}（{epithelial / summary['n_cells']:.2%}）。
- 运行环境为 `ibd_refaudit`；注释解码错误为 0。

## 上皮网页小试验候选

以下是按作者标签名称作出的暂定合并，不是已完成的细胞身份判定。需用 label×category×sample 表及 marker 进一步确认。

| 候选组 | 细胞数 |
|---|---:|
"""
for group, members in groups.items():
    report += f"| {group} | {sum(labels[k] for k in members):,} |\n"
report += """
吸收性 crypt progenitor 与成熟吸收细胞是否分开，Stem/MKI67 群是否单列，以及 Enterochromaffin 是否并入 Enteroendocrine，暂不冻结。

## 需要补齐的证据

1. `obs` 只有样本 ID，尚无 individual_id、disease_status 或炎症字段。必须用服务器现有三张 sidecar 元数据表进行键匹配，不能把样本 ID 直接当作供者 ID。
2. 57 个标签中没有明确命名 Paneth、Glial 或 Neutrophil 的类别。仅表示注释未单列，不能由此判断这些细胞不存在。上皮 pilot 不强行添加 Paneth 列；完整黏膜参考另查覆盖。
3. `counts`、`X` 和 `log1p_cp10k` 都是同形 CSR/float32。尚未验证数值，不能凭层名或 dtype 确认 raw count 或输入尺度。下一步只检查有限稀疏行。
4. 33,538 个唯一 Ensembl ID 对应 33,514 个唯一 gene symbol，有 24 个超过唯一符号数的额外行。不能把这等同于“24 个重复基因”；需导出真实重复列表后制定映射规则。
5. `conf_score` 是分类编码形式、约 118 万个不同值；存储为 categorical 不代表生物分组。尚未核实该分数的语义或设置阈值。
6. 第一轮 value_counts 对样本只导出了前 200 个，不能用该文件推断完整样本/供者平衡；下一轮应生成完整交叉覆盖表。

## 下一步

服务器读取 label、category、sanger_sample_id 三个注释字段，汇总后再连接 sample→individual→疾病/炎症元数据；只抽查少量 CSR 行的表达值和归一化关系。同时导出重复基因符号列表。完成这些检查后再确定身份合并、供者抽样和网页输入尺度。暂不生成 signature，不运行 EcoTyper。
"""
(out / "第一批报告审阅.md").write_text(report, encoding="utf-8")
(out / "source_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
print(report)
