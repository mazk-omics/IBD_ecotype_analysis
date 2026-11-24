# IBD/nonIBD cell-state phenotype validation figure contract

Core conclusion: Cell-level phenotype enrichment is common, but only a small subset of state–dataset effects remains concordant at both sample and donor levels, and CD8 T S04 reverses direction between the CD-oriented SCP1884 cohort and the UC-oriented SCP259 cohort.

Results-level question: Which IBD/nonIBD cell-state differences survive a change from cell-level enrichment to biologically independent sample- and donor-level testing?

Figure archetype: Quantitative validation envelope.

Target journal/output: Manuscript working figure; vector PDF/SVG plus 600-dpi TIFF and PNG preview.

Backend: R.

Final size: 183 mm × 155 mm.

Panel map:

- a: Fraction of evaluable states called significant by each validation method in each dataset.
- b: Sample- and donor-level median-fraction differences for effects significant in both units.
- c: Orthogonal method support for the robust state–dataset effects.

Evidence hierarchy:

- hero evidence: sample/donor effect-size concordance and the cross-cohort CD8 T S04 reversal.
- validation evidence: Cell 2021-style permutation and dominant-state Fisher support.
- controls/robustness: method-wise denominator disclosure, dataset-stratified testing and the ≥20 parent-cell threshold.

Statistics needed: Cell 2021-style 1,000 label shuffles with Z > 1.96 and BH sensitivity; one-sided dominant-state Fisher with BH correction; within-dataset sample/donor Wilcoxon tests with BH correction and median-fraction effect sizes.

Source data needed: Method-level significant/evaluable state counts, paired sample/donor effects and long-format method-support table.

Image-integrity notes: No microscopy or raster manipulation; all panels are generated directly from tabular results.

Reviewer risk: Cell-level tests use cells as units and can overstate evidence; SCP1884 and SCP259 share 12 nonIBD donors; cohort-specific CD-versus-UC composition prevents interpreting the opposite CD8 T S04 directions as a universal IBD effect.
