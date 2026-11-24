# Cell-state annotation evidence method

## 1. Scope

This layer annotates the 85 frozen EcoTyper states after state discovery and validation. It does not redefine states, rerun NMF, change ecotype membership, or alter the frozen marker lists. A label is an interpretation of a transcriptional state, not automatically a conventional cell subtype.

## 2. Evidence hierarchy

Evidence is reviewed in the following order.

1. Frozen discovery markers and enriched pathways define the observed program.
2. Marker recovery in the four external single-cell cohorts establishes reproducibility outside discovery.
3. Leave-one-disease-out, leave-one-region-out and HiRes audits assess robustness and expression-purification support.
4. CellMarker 2.0, PanglaoDB and HPA v26 provide independent marker or expression cross-references.
5. Primary intestinal and IBD single-cell studies provide tissue, region, disease and inflammation context.

Database support alone cannot establish a novel state name. General immune or epithelial databases are supporting evidence only when the tissue and cell-type context agrees.

## 3. Label classes

- `CANONICAL_STATE`: the state matches a published intestinal or immune state with compatible lineage and at least two informative positive markers.
- `FUNCTIONAL_PROGRAM`: markers support an activity program, but not a distinct published subtype.
- `CYCLING_OR_STRESS_PROGRAM`: the dominant signal is cell cycle, interferon, hypoxia, mitochondrial, ribosomal or other state-like variation.
- `CROSS_LINEAGE_OR_CONTAMINATION`: markers conflict with the assigned mother cell type or are dominated by another lineage.
- `UNRESOLVED`: available markers do not support a reproducible biological label.

Labels must remain conservative. A marker name is not promoted to a canonical subtype when the remaining marker set, external recovery or tissue context disagrees.

## 4. Literature evidence grade

- `A`: exact or near-exact published intestinal/IBD state match; at least two informative positive markers; compatible external recovery; no major lineage conflict.
- `B`: compatible intestinal/IBD program supported by primary literature and database evidence, but the exact state identity is not established.
- `C`: biologically plausible functional program supported mainly by pathway or broad marker evidence; no exact intestinal state match.
- `D`: weak, contradictory or cross-lineage evidence; retain only a descriptive label or mark unresolved.

The evidence grade is separate from EcoTyper recovery significance and from the existing annotation-confidence column. Statistical recovery does not by itself prove the biological name.

## 5. Required per-state record

Every state must record the frozen label, proposed manuscript label, label class, positive markers, missing or contradictory markers, database-matched markers, database source IDs, literature-matched markers, literature source IDs, reference state, tissue and region, disease and inflammation context, evidence direction, grade, conflicts, reviewer decision and reviewer note.

The detailed database matches are stored separately so that each gene, reference cell name, tissue, PMID or database identifier, and HPA expression rank remains auditable.

## 6. Database interpretation

- CellMarker 2.0 is split into broad matches and preferred intestinal or immune-context matches. Entries may differ in assay and evidence type.
- PanglaoDB is a broad single-cell marker cross-reference. Its 2020-03-27 marker file is not intestine-specific.
- HPA v26 is used as a quantitative expression cross-reference. `preferred` means target cell-type rank at most 10 and target maximum nCPM at least 1. This transparent descriptive threshold is not an annotation acceptance threshold.
- Human Gut Cell Atlas v1.0 is registered but remains `TO_INTEGRATE` until a reproducible mapping and versioned evidence table are generated.

## 7. Manuscript wording

The Methods should describe state annotation as marker-, pathway-, validation-, database- and literature-supported interpretation. Results should use `program` when evidence does not establish a canonical subtype. Supplementary tables should report all evidence grades and conflicts, including low-confidence and unresolved states. Discussion should not equate a recovered transcriptional program with a stable cell lineage or subtype.

For presentation, grades A-B are literature-anchored named states, grade C is a marker-led functional program, and grade D is unresolved or cross-lineage. Grades A-C are biologically interpretable under the Cell 2021 convention, but grade C must be described as a program rather than a canonical subtype. Grade D is restricted to supplementary reporting unless a sensitivity analysis specifically concerns state uncertainty. This tier controls wording and placement, not whether the state remains in the frozen EcoTyper model.

## 8. Reproducibility

Database files, mapping rules, input hashes, detailed gene-level matches, state-level summaries and completion records are frozen under `inputs/reference_databases/20261123` and `output/10_state_annotation/literature_database_review`. The original 85-state table is read-only in this process.
