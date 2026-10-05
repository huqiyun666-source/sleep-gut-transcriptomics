# GSE289089 cluster annotation evidence reconstruction notes

Update date: 2026-10-04  
Scope: provenance/figure reconstruction from the frozen Seurat object using the original analysis parameters. Cluster identities and working annotations were not changed.

## Provenance distinction

The original saved `FindAllMarkers` result and original DotPlot were not present in the archived project files. The files generated in this update are therefore labelled **RECONSTRUCTED_WITH_ORIGINAL_PARAMETERS**, not original saved outputs.

The contemporaneous analysis record recovered the original marker call:

```r
FindAllMarkers(
  combined.clean,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
```

The same invocation was applied once to the current frozen object, `processed/GSE289089_clustered.rds`, after verifying 27,671 cells and unchanged identities 0–17. No normalization, clustering, identity assignment, or downstream manuscript analysis was rerun.

## Reconstructed evidence

- Full marker table: `GSE289089_FindAllMarkers_reconstructed.csv` (26,565 rows).
- Annotation summary: `GSE289089_cluster_annotation_evidence_reconstructed.csv` (18 clusters).
- Canonical marker DotPlot: `GSE289089_annotation_DotPlot.png` and `GSE289089_annotation_DotPlot.pdf`.
- All 32 requested canonical-panel genes were present in the frozen RNA assay.

Representative top markers in the annotation summary are the first 10 unique genes after ordering each reconstructed cluster result by decreasing average log2 fold change, then adjusted P value and gene symbol. This is a presentation summary of the reconstructed marker table, not a recovered original top-marker selection rule.

## Annotation review

The reconstructed marker table and canonical panel support the retained working annotations:

- cluster 0 shows stress/immediate-early markers including `Areg`, `Dusp1`, `Atf3`, and `Fos`;
- cluster 1 includes interferon-response markers `Ddx60` and `Irf7`;
- clusters 2, 7, 10, 11, and 13 show proliferative/cell-cycle marker programs;
- cluster 3 shows secretory/Goblet evidence;
- clusters 4 and 15 show enterocyte markers;
- cluster 5 shows Paneth/defensin markers;
- cluster 6 shows ISC markers `Lgr5`, `Olfm4`, `Ascl2`, and `Smoc2`;
- cluster 8 shows `Dclk1`, `Pou2f3`, and `Trpm5`;
- clusters 9, 12, and 14 show EEC/EC-lineage evidence, including `Tph1` in cluster 12 and `Neurog3` in cluster 14;
- cluster 16 combines Paneth/defensin and ISC-lineage evidence;
- cluster 17 shows immune markers `Ptprc` and `Lst1`.

No clear discrepancy requiring a change to the frozen working annotations was identified. The absence of `Cd3d` from cluster 17's reconstructed positive-marker list at the original thresholds does not conflict with the deliberately broad annotation “immune”; it only means the reconstructed evidence does not support a more specific T-cell label.

## Files retained from the earlier archive audit

`GSE289089_cluster_annotation_evidence.csv` remains the record of what could be recovered before reconstruction. It must not be cited as the reconstructed marker table. The newly generated files above carry the explicit reconstructed provenance.

