# Cell 2021 phenotype enrichment and Figure S2G/H contract

Core conclusion: IBD and nonIBD share most recoverable states, while cell-level phenotype enrichment is widespread and dominant-state over-representation is concentrated in a smaller, dataset-dependent subset.

Results-level question: Which frozen Signature-D states are enriched in IBD phenotypes, and which enrichments remain visible when each biological sample is reduced to its dominant state per parent type?

Figure archetype: Validation envelope with full-state heatmaps and a compact method summary.

Target journal/output: Supplementary working figures; vector PDF/SVG, 600-dpi TIFF and PNG preview.

Backend: R.

Final size: 183 mm wide; 285 mm high for full-state heatmaps and 105 mm high for summary figures.

Panel map:

- Known phenotype enrichment: permutation Z for every evaluable state, phenotype level and dataset.
- Figure S2G counterpart: percentage of eligible biological samples assigned to each dominant state.
- Figure S2H counterpart: one-sided Fisher enrichment with BH correction for the same state-by-phenotype grid.
- Summary: method-level significant fractions by dataset and phenotype axis.

Evidence hierarchy:

- hero evidence: complete state-by-phenotype heatmaps without selecting favorable states.
- validation evidence: dominant-state frequency and Fisher-BH significance shown side by side.
- controls/robustness: dataset stratification, explicit non-evaluable cells and the primary threshold of at least 20 parent-type cells.

Statistics needed: 1,000 state-label permutations and Z > 1.96 following Cell 2021; BH sensitivity; one-sided Fisher exact tests on dominant-state assignments with BH correction.

Source data needed: all permutation rows, all within-dataset Fisher rows, dominant-state percentages and method-level summaries.

Image-integrity notes: All panels are generated directly from tabular results; display caps do not alter source values.

Reviewer risk: cell-level permutations do not provide independent-sample inference; phenotype axes differ in evaluability across datasets; shared nonIBD donors in SCP1884 and SCP259 preclude treating those cohorts as fully independent replication.
