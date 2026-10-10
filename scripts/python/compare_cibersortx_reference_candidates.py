"""Compare audited signatures on native and shared axes, without deconvolution."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


def condition(frame, unit=False):
    a = frame.to_numpy(float)
    if unit:
        a = a / np.linalg.norm(a, axis=0)
    return float(np.linalg.cond(a))


def run(root):
    root = Path(root)
    modes = ["gse_joint", "gse_ileum", "scp_rectum"]
    signatures = {m: pd.read_csv(next((root/m).glob("*_signature_gene_symbols_restored.txt")), sep="\t", index_col=0) for m in modes}
    summaries = {m: json.loads((root/m/"audit_summary.json").read_text(encoding="utf-8")) for m in modes}
    shared_types = sorted(set.intersection(*(set(s.columns) for s in signatures.values())))
    shared_genes = sorted(set.intersection(*(set(s.index) for s in signatures.values())))
    assert len(shared_types) == 17 and len(shared_genes) > len(shared_types)
    all_types = sorted(set.union(*(set(s.columns) for s in signatures.values())))
    stats, coverage, markers, categories, pairs, gene_sets, profile = [], [], [], [], [], [], []
    mt_sensitivity = []
    for mode, s in signatures.items():
        meta = summaries[mode]
        common = s.loc[shared_genes, shared_types]
        mt = [g for g in s.index if g.startswith("MT-")]
        stats.append({"candidate": mode, "signature_genes": len(s), "types": len(s.columns),
                      "native_condition_number": condition(s), "native_column_unit_condition_number": condition(s, True),
                      "shared_17_types_2227_genes_unit_condition_number": condition(common, True),
                      "n_mitochondrial_signature_genes": len(mt), "gene_names_restored": meta["signature_gene_symbols_restored"],
                      "rank": int(np.linalg.matrix_rank(s.to_numpy())),
                      "zero_rows_in_inferred_reference": meta["inferred_reference_zero_rows"]})
        for t in all_types:
            coverage.append({"candidate": mode, "pilot_type": t, "included": t in s.columns})
        p = pd.read_csv(root/mode/"marker_panel_audit.tsv", sep="\t")
        p.insert(0, "candidate", mode); markers.append(p)
        q = pd.read_csv(root/mode/"gene_category_signal.tsv", sep="\t")
        q.insert(0, "candidate", mode); categories.append(q)
        for t in s.columns:
            assert s[t].sum() > 0
            native_peak = s[t].nlargest(10)
            for rank, (gene, expression) in enumerate(native_peak.items(), 1):
                profile.append({"candidate": mode, "pilot_type": t, "rank": rank, "gene": gene, "expression": expression})
        for x, y in itertools.combinations(shared_types, 2):
            pairs.append({"candidate": mode, "type1": x, "type2": y,
                          "native_raw_pearson": float(s[x].corr(s[y])),
                          "native_log2p1_pearson": float(np.log2(s[x]+1).corr(np.log2(s[y]+1))),
                          "shared_gene_raw_pearson": float(common[x].corr(common[y])),
                          "shared_gene_log2p1_pearson": float(np.log2(common[x]+1).corr(np.log2(common[y]+1)))})
        if mt:
            without = s.drop(index=mt)
            for x, y in [("Stem_TA_progenitor", "Colonic_absorptive"), ("Stem_TA_progenitor", "Tuft"),
                         ("Colonic_absorptive", "Tuft"), ("CD8_T", "NK"), ("DC", "Monocyte_Macrophage")]:
                mt_sensitivity.append({"candidate": mode, "type1": x, "type2": y,
                    "original_raw_pearson": float(s[x].corr(s[y])), "without_mito_raw_pearson": float(without[x].corr(without[y])),
                    "original_log2p1_pearson": float(np.log2(s[x]+1).corr(np.log2(s[y]+1))),
                    "without_mito_log2p1_pearson": float(np.log2(without[x]+1).corr(np.log2(without[y]+1)))})
            stats[-1]["without_mito_native_condition_number"] = condition(without)
            stats[-1]["without_mito_unit_condition_number"] = condition(without, True)
    for left, right in itertools.combinations(modes, 2):
        a, b = set(signatures[left].index), set(signatures[right].index)
        gene_sets.append({"candidate1": left, "candidate2": right, "shared_genes": len(a&b),
                          "jaccard": len(a&b)/len(a|b), "fraction_covered_in_candidate1": len(a&b)/len(a),
                          "fraction_covered_in_candidate2": len(a&b)/len(b)})
    outframes = {
        "candidate_summary.tsv": pd.DataFrame(stats), "cell_type_coverage.tsv": pd.DataFrame(coverage),
        "marker_comparison.tsv": pd.concat(markers, ignore_index=True), "gene_category_comparison.tsv": pd.concat(categories, ignore_index=True),
        "pair_similarity_comparison.tsv": pd.DataFrame(pairs), "signature_gene_set_overlap.tsv": pd.DataFrame(gene_sets),
        "top_expression_genes_comparison.tsv": pd.DataFrame(profile), "mitochondrial_sensitivity.tsv": pd.DataFrame(mt_sensitivity),
        "common_signature_genes.tsv": pd.DataFrame({"gene": shared_genes}), "common_cell_types.tsv": pd.DataFrame({"pilot_type": shared_types})}
    for filename, df in outframes.items():
        df.to_csv(root/filename, sep="\t", index=False)
    refs = {}
    for mode in modes:
        r = pd.read_csv(next((root/mode).glob("*_inferred_refsample_gene_symbols_restored.txt")), sep="\t", index_col=0)
        types = [x.rsplit(".", 1)[0] if x.rsplit(".", 1)[-1].isdigit() else x for x in r.columns]
        refs[mode] = pd.DataFrame({t: r.loc[:, np.asarray(types)==t].mean(axis=1) for t in signatures[mode].columns})
    common_reference_genes = sorted(set.intersection(*(set(r.index) for r in refs.values())))
    consistency = []
    for left, right in itertools.combinations(modes, 2):
        for t in shared_types:
            a, b = refs[left].loc[common_reference_genes, t], refs[right].loc[common_reference_genes, t]
            consistency.append({"candidate1": left, "candidate2": right, "pilot_type": t,
                                "n_common_reference_genes": len(common_reference_genes),
                                "raw_pearson": float(a.corr(b)), "log2p1_pearson": float(np.log2(a+1).corr(np.log2(b+1)))})
    pd.DataFrame(consistency).to_csv(root/"cross_candidate_reference_profile_consistency.tsv", sep="\t", index=False)
    data = {"status": "completed", "candidate_statistics": stats, "common_types": shared_types,
            "common_signature_gene_count": len(shared_genes), "common_reference_gene_count": len(common_reference_genes),
            "gene_set_overlap": gene_sets, "mitochondrial_sensitivity": mt_sensitivity,
            "recommendation_for_whole_MSCCR": "gse_joint as provisional primary; gse_ileum and scp_rectum as regional sensitivity references",
            "real_bulk_gene_coverage_tested": False, "heldout_pseudobulk_accuracy_tested": False,
            "source_files_modified": False,
            "diagnostic_normalization_only": "L2 unit columns and log2p1 were used only for comparison; output signatures unchanged",
            "cautions": ["Shared-axis conditioning omits unique informative genes and types; not an accuracy ranking",
                         "Marker presence depends on reference gene availability and sampled subtypes",
                         "Mitochondrial deletion is a numerical sensitivity analysis, not a validated production filter",
                         "No biological or clinical association test was performed"]}
    (root/"comparison_summary.json").write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-root", required=True)
    run(parser.parse_args().audit_root)
