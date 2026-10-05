# Sleep deprivation is associated with cell-type-heterogeneous suppression of intestinal oxidative-phosphorylation transcriptional programs

## Study overview

This repository accompanies a transcriptomic study of intestinal responses to sleep deprivation. It preserves the frozen analysis scripts, manuscript-critical result tables, public-dataset search audit, provenance records, and final main-figure PDFs. It does not contain raw sequencing reads, reference resources, large intermediate objects, or unpublished credentials.

## Datasets and roles

| Dataset | Role | Design summary |
|---|---|---|
| GSE289089 / PRJNA1221288 | Single-cell discovery | Mouse small-intestinal crypt epithelium; control and sleep deprivation, n=2 biological samples per condition |
| PRJNA1167170 | First proximal-colon validation | Mouse colon bulk RNA-seq; control and sleep deprivation, n=4 per condition |
| PRJNA862187 | Second factorial colon validation | Mouse colon bulk RNA-seq; sleep condition by diet, n=3 per factorial cell |

PRJNA1395688 was recorded as potentially eligible but was not analysed because its public library metadata were inconsistent at the search cutoff.

## Analysis overview

The frozen workflow comprises scRNA-seq QC and annotation provenance, sample-level pseudobulk aggregation, edgeR differential-expression modelling, fgsea pathway testing, two independent bulk-RNA-seq validations, a documented public-dataset search, cross-dataset OXPHOS concordance, and prespecified sensitivity analyses. No fourth dataset was analysed.

## Repository contents

- `scripts/`: portable copies of scripts used for the frozen analyses and figures.
- `provenance/original_scripts/`: unmodified source-script copies retained for traceability.
- `frozen_plans/`: the genuine pre-analysis freeze available for PRJNA862187.
- `results_tables/`: manuscript-critical final tables only.
- `dataset_search/`: candidate registry, search log, selection table, flow record, and completeness audit.
- `figure_sources/`: selected frozen tables required by the plotting scripts.
- `docs/`: the submission-ready manuscript snapshot and final main-figure PDFs.
- `file_manifest.csv` and `checksums_sha256.txt`: package inventory and integrity records.

## Reproduction order

Run commands from the repository root and supply or download the inputs described in `REPRODUCIBILITY.md`.

1. **01 dataset search:** `python3 scripts/dataset_search/build_dataset_search_outputs.py --registry dataset_search/dataset_candidate_registry.csv --output-dir dataset_search`
2. **02 discovery single-cell:** reconstruct the archived Seurat inputs as documented in `provenance/GSE289089_QC_provenance.md`, then run the scripts in `scripts/single_cell/` with the repository root as their project-root argument where required.
3. **03 PRJNA1167170:** run `scripts/PRJNA1167170/PRJNA1167170_bulk_pipeline.sh`; the associated R analysis is `scripts/PRJNA1167170/PRJNA1167170_tximport_edgeR_fgsea.R`.
4. **04 PRJNA862187:** follow `frozen_plans/PRJNA862187_analysis_plan_frozen.md`, run the pilot and full acquisition modes in `scripts/PRJNA862187/PRJNA862187_pipeline.sh`, then run `scripts/PRJNA862187/PRJNA862187_edgeR_fgsea.R` with the repository root.
5. **05 cross-dataset:** run the scripts in `scripts/three_dataset_validation/` in the documented dependency order.
6. **06 figures:** run `Rscript scripts/figures/generate_all.R` after restoring the required frozen inputs listed in `provenance/figure_source_map_v2.csv`.

## Software requirements

Confirmed versions are recorded in `software_versions.csv`. Core requirements include R, Seurat, edgeR, fgsea, msigdbr, tximport, data.table, ggplot2, Salmon, FastQC, MultiQC, Python 3, Pillow, and ReportLab. The exact scDblFinder version is unavailable and is explicitly labelled as such.

## Public data access

Public accessions and source links are listed in `accession_manifest.csv`. Raw reads must be retrieved from the public archives; they are not redistributed here.

## Known limitations

- GSE289089 has n=2 biological samples per condition.
- No functional mitochondrial assays were available.
- The PRJNA862187 OXPHOS estimate was negative but not statistically significant.
- PRJNA1395688 had inconsistent public metadata and was not analysed.
- Transcriptomic associations do not establish causality.
- The original end-to-end GSE289089 preprocessing script was not archived. Its recovered QC and doublet-removal provenance are documented, but the exact `PercentageFeatureSet` call and scDblFinder package version remain unavailable.
- Large Seurat, tximport, edgeR, raw-read, and reference-index objects are excluded. Re-rendering every figure therefore requires restoring the frozen inputs named in the source map.
- Numeric and reference audits are included as final records; no separate generator script for those audits was present in the archived project.

## Citation

Please cite the associated manuscript and the archived repository release. Bibliographic details and the Zenodo DOI must be added after publication/deposition; see `CITATION.cff`.

## Contact

Corresponding author: Hu Zhang  
Email: 20234028@zcmu.edu.cn  
ORCID: none / not available

## Licensing

Code is licensed under the MIT License. Documentation and processed result tables are licensed under the Creative Commons Attribution 4.0 International License. See `LICENSE`, `LICENSE-CODE-MIT.txt`, and `LICENSE-DOCUMENTATION-CC-BY-4.0.md` for scope and terms. Third-party source data and metadata remain subject to their original terms.
