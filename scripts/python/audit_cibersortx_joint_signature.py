"""Read-only audit of a CIBERSORTx signature and its inferred references."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


PANELS = {
    "B": ["MS4A1", "CD79A", "CD79B", "CD19", "CD74"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1", "XBP1", "DERL3"],
    "CD4_T": ["CD3D", "CD3E", "TRAC", "CD4", "IL7R", "CCR7", "CD40LG"],
    "CD8_T": ["CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "CCL5"],
    "NK": ["NKG7", "GNLY", "KLRD1", "PRF1", "FCGR3A"],
    "ILC": ["IL7R", "KLRB1", "KIT", "RORC", "IL23R", "AHR"],
    "Monocyte_Macrophage": ["LST1", "TYROBP", "FCER1G", "CSF1R", "LYZ", "CD68"],
    "DC": ["FCER1A", "CD1C", "CLEC9A", "CLEC10A", "HLA-DRA"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"],
    "Pericyte": ["RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "ACTA2"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "KDR", "EMCN"],
    "Ileal_absorptive": ["SI", "ALPI", "FABP1", "FABP2", "APOA4", "SLC5A1"],
    "Colonic_absorptive": ["CA1", "CA2", "AQP8", "SLC26A3", "BEST4", "OTOP2"],
    "Goblet": ["MUC2", "CLCA1", "SPINK4", "FCGBP", "ZG16", "TFF3"],
    "Stem_TA_progenitor": ["OLFM4", "LGR5", "MKI67", "TOP2A", "PCNA", "SOX9"],
    "Paneth": ["DEFA5", "DEFA6", "LYZ", "PLA2G2A", "ITLN2"],
    "Enteroendocrine": ["CHGA", "CHGB", "TPH1", "NEUROD1", "PCSK1"],
    "Glial": ["SOX10", "PLP1", "S100B", "GFAP", "SOX2"],
    "Tuft": ["POU2F3", "DCLK1", "TRPM5", "AVIL", "SH2D6", "IL25", "GFI1B"],
}


def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def write_table(df, out, name):
    df.to_csv(out / name, sep="\t", index=False)


def main(args):
    root, out = Path(args.input_dir), Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    prefix = f"CIBERSORTx_{args.job}_"
    sig_path = next(root.glob(prefix + "*.bm.K999.txt"))
    ref_path = next(root.glob(prefix + "*_inferred_refsample.txt"))
    cls_path = next(root.glob(prefix + "*_inferred_phenoclasses.txt"))
    pdf_path = next(root.glob(prefix + "*.pdf"))
    sources = [sig_path, ref_path, cls_path, pdf_path]
    before = {p.name: sha256(p) for p in sources}
    source = Path(args.source_manifest)
    manifest = pd.read_csv(source, sep="\t")
    input_meta = json.loads((source.parent / "summary.json").read_text(encoding="utf-8"))
    input_web = source.parent / input_meta["web_input_filename"]
    with input_web.open(encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        assert header[1:] == manifest.web_column.tolist()
        input_gene_order = [line.split("\t", 1)[0] for line in f]
    sig = pd.read_csv(sig_path, sep="\t", index_col=0)
    ref = pd.read_csv(ref_path, sep="\t", index_col=0)
    cls = pd.read_csv(cls_path, sep="\t", header=None, index_col=0)
    # Restore only observed, one-to-one transformations verified against source.
    assert len(input_gene_order) == len(ref)
    expected_output = [g.replace("-", ".") for g in input_gene_order]
    expected_output = ["X" + g if re.match(r"^(?:[0-9]|\.[0-9])", g) else g for g in expected_output]
    assert ref.index.tolist() == expected_output
    backmap = dict(zip(ref.index, input_gene_order))
    assert len(backmap) == len(input_gene_order) == len(set(input_gene_order))
    changes = [{"web_output_gene": output, "source_gene_symbol": original,
                "in_signature": output in sig.index}
               for output, original in backmap.items() if output != original]
    assert sig.index.isin(ref.index).all()
    write_table(pd.DataFrame(changes), out, "gene_symbol_transformations.tsv")
    output_label = args.output_label
    for original, restored_name in [(sig_path, output_label + "_signature_gene_symbols_restored.txt"),
                                    (ref_path, output_label + "_inferred_refsample_gene_symbols_restored.txt")]:
        with original.open(encoding="utf-8", newline="") as f, (out/restored_name).open("w", encoding="utf-8", newline="") as w:
            w.write(f.readline())
            for line in f:
                symbol, values = line.split("\t", 1)
                w.write(backmap[symbol] + "\t" + values)
        with original.open(encoding="utf-8", newline="") as f, (out/restored_name).open(encoding="utf-8", newline="") as g:
            assert f.readline() == g.readline()
            for left, right in zip(f, g, strict=True):
                assert left.split("\t", 1)[1] == right.split("\t", 1)[1]
    sig.index = sig.index.map(backmap)
    ref.index = ref.index.map(backmap)
    a, r, c = sig.to_numpy(float), ref.to_numpy(float), cls.to_numpy(int)
    assert sig.index.is_unique and ref.index.is_unique and cls.index.is_unique
    assert sig.columns.is_unique and ref.columns.is_unique
    assert np.isfinite(a).all() and np.isfinite(r).all()
    assert (a >= 0).all() and (r >= 0).all() and (a.sum(0) > 0).all() and (r.sum(0) > 0).all()
    assert set(sig.columns) == set(cls.index) == set(manifest.pilot_type)
    assert set(sig.columns).issubset(PANELS)
    assert c.shape == (len(sig.columns), len(ref.columns))
    assert set(np.unique(c)).issubset({0, 1, 2}) and ((c == 1).sum(0) == 1).all()
    assert ((c == 1).sum(1) == 5).all()
    assert sig.index.isin(ref.index).all()
    labels = [str(x).rsplit(".", 1)[0] if str(x).rsplit(".", 1)[-1].isdigit() else str(x) for x in ref.columns]
    for j, label in enumerate(labels):
        assert cls.index[np.flatnonzero(c[:, j] == 1)[0]] == label
    sv = np.linalg.svd(a, compute_uv=False)
    raw_corr = sig.corr()
    log_corr = pd.DataFrame(np.log2(a + 1), index=sig.index, columns=sig.columns).corr()
    for matrix, name in [(raw_corr, "signature_raw_pearson.tsv"), (log_corr, "signature_log2p1_pearson.tsv")]:
        matrix.to_csv(out / name, sep="\t")
    pairs = []
    for i, x in enumerate(sig.columns):
        for y in sig.columns[i+1:]:
            pairs.append({"type1": x, "type2": y, "raw_pearson": raw_corr.loc[x, y], "log2p1_pearson": log_corr.loc[x, y]})
    pair_df = pd.DataFrame(pairs).sort_values("log2p1_pearson", ascending=False)
    write_table(pair_df, out, "signature_pair_similarity.tsv")
    top, diagnostics, markers, coherence = [], [], [], []
    genes = sig.index.to_numpy(); best = a.argmax(1)
    sorted_a = np.sort(a, axis=1)
    for i, name in enumerate(sig.columns):
        target = a[:, i]; other = np.delete(a, i, axis=1).max(1)
        fold = (target + 1) / (other + 1)
        for rank, ix in enumerate(np.argsort(fold)[::-1][:30], 1):
            top.append({"pilot_type": name, "rank": rank, "gene": genes[ix], "target_expression": target[ix], "max_other_expression": other[ix], "target_vs_max_other_ratio_plus1": fold[ix]})
        diagnostics.append({"pilot_type": name, "n_genes_highest_in_type": int((best == i).sum()),
            "n_genes_ratio_plus1_ge2_and_expression_ge10": int(((best == i) & (fold >= 2) & (target >= 10)).sum()),
            "column_sum_signature_subset": float(target.sum()), "nearest_raw_type": raw_corr[name].drop(name).idxmax(),
            "nearest_raw_pearson": float(raw_corr[name].drop(name).max()), "nearest_log_type": log_corr[name].drop(name).idxmax(),
            "nearest_log2p1_pearson": float(log_corr[name].drop(name).max())})
        use = np.asarray(labels) == name
        full_rep = r[:, use]
        mu = full_rep.mean(1)
        cpm = (full_rep / full_rep.sum(0) * 1e6).mean(1)
        loc = ref.index.get_indexer(sig.index)
        signature_target = target
        coherence.append({"pilot_type": name, "n_internal_replicates": int(use.sum()),
            "reference_column_sum_min": float(full_rep.sum(0).min()), "reference_column_sum_max": float(full_rep.sum(0).max()),
            "signature_vs_mean_ref_pearson": float(np.corrcoef(signature_target, mu[loc])[0,1]),
            "signature_vs_mean_ref_log2p1_pearson": float(np.corrcoef(np.log2(signature_target+1), np.log2(mu[loc]+1))[0,1]),
            "signature_vs_mean_ref_cpm_pearson": float(np.corrcoef(signature_target, cpm[loc])[0,1])})
        for gene in PANELS[name]:
            in_sig, in_ref = gene in sig.index, gene in ref.index
            rec = {"pilot_type": name, "gene": gene, "in_signature": in_sig, "in_inferred_reference": in_ref}
            if in_sig:
                z = sig.loc[gene]
                rec.update(target_signature=float(z[name]), highest_signature_type=z.idxmax(), highest_signature_value=float(z.max()),
                           target_vs_max_other_ratio_plus1=float((z[name]+1)/(z.drop(name).max()+1)))
            if in_ref:
                z = ref.loc[gene].to_numpy(float)
                type_means = {t: float(z[np.asarray(labels)==t].mean()) for t in sig.columns}
                rec.update(target_inferred_reference_mean=type_means[name], highest_reference_type=max(type_means, key=type_means.get))
            markers.append(rec)
    for frame, name in [(pd.DataFrame(top), "top_diagnostic_genes_per_type.tsv"), (pd.DataFrame(diagnostics), "cell_type_diagnostics.tsv"),
                        (pd.DataFrame(markers), "marker_panel_audit.tsv"), (pd.DataFrame(coherence), "signature_reference_coherence.tsv")]:
        write_table(frame, out, name)
    categories = {
        "mitochondrial": [g for g in genes if g.startswith("MT-")],
        "ribosomal_prefix": [g for g in genes if g.startswith(("RPL", "RPS"))],
        "immunoglobulin_prefix": [g for g in genes if g.startswith(("IGH", "IGK", "IGL"))],
        "HLA_prefix": [g for g in genes if g.startswith("HLA-")],
    }
    category_report = []
    for category, gg in categories.items():
        for name in sig.columns:
            category_report.append({"category": category, "pilot_type": name, "n_signature_genes": len(gg),
                                    "fraction_column_signal_signature_subset": float(sig.loc[gg, name].sum()/sig[name].sum()) if gg else 0})
    write_table(pd.DataFrame(category_report), out, "gene_category_signal.tsv")
    write_table(pd.DataFrame({"gene": genes, "highest_expression_type": sig.columns[best],
        "highest_expression": sorted_a[:, -1], "second_highest_expression": sorted_a[:, -2],
        "highest_vs_second_ratio_plus1": (sorted_a[:, -1]+1)/(sorted_a[:, -2]+1)}), out, "signature_gene_diagnostic_assignment.tsv")
    assert set(manifest.pilot_type) == set(sig.columns)
    input_genes = set(input_gene_order)
    assert set(ref.index) == input_genes
    selected = {
        name: {"n_cells": int((manifest.pilot_type==name).sum()),
               "n_donors": int(manifest.loc[manifest.pilot_type==name, "individual_id"].nunique()),
               "sites": manifest.loc[manifest.pilot_type==name, "Site"].value_counts().to_dict(),
               "author_labels": manifest.loc[manifest.pilot_type==name, "author_label"].value_counts().to_dict()}
        for name in sig.columns
    }
    summary = {
        "status": "structural_checks_passed", "job": args.job, "signature_shape": list(a.shape), "reference_shape": list(r.shape),
        "phenotype_shape_excluding_label": list(c.shape), "cell_types": list(sig.columns),
        "signature_genes_subset_of_input": True, "reference_gene_set_matches_input": True,
        "source_verified_gene_symbol_transformations": len(changes),
        "signature_gene_symbols_restored": sum(x["in_signature"] for x in changes),
        "gene_restoration_does_not_change_numeric_text": True,
        "nonfinite_or_negative_values": 0, "signature_min": float(a.min()), "signature_max": float(a.max()),
        "signature_exactly_one_fraction": float((a==1).mean()),
        "signature_rank": int(np.linalg.matrix_rank(a)), "signature_condition_number_2norm": float(sv[0]/sv[-1]),
        "signature_singular_max": float(sv[0]), "signature_singular_min": float(sv[-1]),
        "internal_reference_replicates_per_type": 5, "internal_replicates_are_not_independent_donors": True,
        "inferred_reference_zero_rows": int((r.sum(1) == 0).sum()),
        "inferred_reference_column_sum_min": float(r.sum(0).min()),
        "inferred_reference_column_sum_max": float(r.sum(0).max()),
        "signature_reference_raw_correlation_min": float(min(x["signature_vs_mean_ref_pearson"] for x in coherence)),
        "signature_reference_raw_correlation_max": float(max(x["signature_vs_mean_ref_pearson"] for x in coherence)),
        "categories": {k: len(v) for k,v in categories.items()}, "top_log_correlated_pairs": pair_df.head(8).to_dict("records"),
        "source_sampling": selected, "source_files": [{"path": str(p), "bytes": p.stat().st_size, "sha256": before[p.name]} for p in sources],
        "limitations": ["Marker panel is diagnostic, not an exhaustive identity test", "Argmax/ratio gene assignments are audit summaries, not tool per-phenotype DEG membership",
                        "Column correlation is a similarity warning, not observed deconvolution error", "No heldout-pseudobulk or MSCCR deconvolution has been run"]
    }
    assert {p.name: sha256(p) for p in sources} == before
    (out / "audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k not in {"source_sampling", "source_files"}}, ensure_ascii=False, indent=2))
    print("MARKER_COUNTS", pd.DataFrame(markers).groupby("pilot_type").in_signature.sum().to_dict())
    print("DIAGNOSTIC_COUNTS", diagnostics)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--job", default="Job4")
    parser.add_argument("--output-label", default="GSE282122_Joint20")
    main(parser.parse_args())
