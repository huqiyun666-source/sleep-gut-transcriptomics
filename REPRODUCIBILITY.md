# Reproducibility guide

## Scope

This package preserves the code and compact evidence needed to audit the frozen manuscript results. It does not recompute results during packaging and does not distribute raw sequencing data or large binary intermediates.

## Public inputs

- GSE289089 / PRJNA1221288: obtain the four single-cell libraries and available processed matrices from GEO/SRA.
- PRJNA1167170: obtain the eight bulk-RNA-seq runs from SRA/ENA.
- PRJNA862187: obtain only the twelve prespecified colon RNA-seq runs listed in the frozen plan; fecal metagenomic runs are outside scope.

See `accession_manifest.csv` and `dataset_search/dataset_candidate_registry.csv` for archive links and eligibility provenance.

## Reference resource

Bulk quantification used GENCODE Mouse release M39 on GRCm39. The transcriptome and Salmon index are not distributed. Download the matching GENCODE M39 transcript FASTA, retain its release identity, and build the index with Salmon 2.7.0, for example:

```bash
salmon index \
  -t reference/gencode_M39/gencode.vM39.transcripts.fa.gz \
  -i reference/salmon_mouse_GRCm39
```

Do not substitute a different genome build or annotation release when reproducing the frozen analysis.

## Execution order and expected outputs

1. **Dataset-search audit**
   - Script: `scripts/dataset_search/build_dataset_search_outputs.py`
   - Input: `dataset_search/dataset_candidate_registry.csv`
   - Expected outputs: the supplementary search table and selection-flow exports.
2. **GSE289089 single-cell discovery**
   - Scripts: `scripts/single_cell/`
   - Required large inputs: reconstructed/frozen Seurat objects under `processed/`.
   - Expected compact outputs: cell-type OXPHOS tables, sensitivity summaries, and annotation reconstruction records.
   - Important gap: the original raw-to-Seurat preprocessing script was not archived. Use `provenance/GSE289089_QC_provenance.md`, `provenance/GSE289089_scDblFinder_provenance.md`, and `provenance/Methods_provenance_table.csv`; do not infer unavailable parameters.
3. **PRJNA1167170 validation**
   - Scripts: `scripts/PRJNA1167170/`
   - Expected outputs: Salmon QC, tximport counts, edgeR results, primary OXPHOS validation, and secondary mitochondrial pathway tables.
4. **PRJNA862187 factorial validation**
   - Plan: `frozen_plans/PRJNA862187_analysis_plan_frozen.md`
   - Scripts: `scripts/PRJNA862187/`
   - Expected outputs: primary average sleep-deprivation effect, diet-specific contrasts, exploratory interaction results, and QC summaries.
5. **Cross-dataset validation and sensitivities**
   - Scripts: `scripts/three_dataset_validation/` and the relevant single-cell sensitivity scripts.
   - Expected outputs: concordance tables, pathway summary, and downsampling summary in `results_tables/`.
6. **Figures**
   - Scripts: `scripts/figures/`
   - Source mapping: `provenance/figure_source_map_v2.csv`
   - Expected outputs: manuscript figure PDF/PNG files. Final main-figure PDFs are supplied in `docs/final_figures/` for verification.

## Frozen-result mapping

- Manuscript headline values: `results_tables/manuscript_key_results.csv`
- Cell-type results: `results_tables/manuscript_celltype_results.csv` and `results_tables/celltype_OXPHOS_logFC_ranking.csv`
- PRJNA1167170: files prefixed `PRJNA1167170_` in `results_tables/`
- PRJNA862187: files prefixed `PRJNA862187_` in `results_tables/`
- Three-dataset concordance: `results_tables/three_dataset_*`
- Sensitivity and exploratory regeneration summaries: `results_tables/downsampled_OXPHOS_sensitivity.csv` and `results_tables/ISC_TA_regeneration_pathways.csv`
- Numeric audit: `results_tables/numeric_audit_v3.csv`

## Portable versus original scripts

Files in `scripts/` are publication copies with path-only portability changes. Files in `provenance/original_scripts/` are unmodified archived source copies. Path edits do not change accession lists, parameters, statistical models, or frozen numerical outputs.

## Exclusions

The package excludes FASTQ/SRA data, Salmon indices and per-sample quantification directories, reference files, `.rds`/`.RData` objects, temporary logs, downloaded publications, and historical result versions. These exclusions are intentional and listed in `.gitignore`.
