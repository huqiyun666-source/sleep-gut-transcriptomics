# PRJNA862187 external-validation analysis plan — FROZEN

Freeze timestamp: **2026-10-04 (Asia/Shanghai), before FASTQ download or inspection of expression results**

Status: **FROZEN BEFORE DATA ACQUISITION**

This analysis was initiated after a systematic public-dataset audit identified PRJNA862187 as an additional eligible intestinal sleep-deprivation transcriptomic study. The primary pathway and all contrasts below were specified from the previously completed and frozen GSE289089 discovery analysis and PRJNA1167170 independent validation. They must not be changed in response to PRJNA862187 results.

## Study design to be verified before download

Expected design: a 2 × 2 factorial experiment with 12 independent mouse colon/large-intestine RNA-seq samples, three biological replicates per cell:

| Group | Sleep condition | Diet | Expected n |
|---|---|---|---:|
| EC_SCD | exercise control | standard chow diet | 3 |
| SD_SCD | sleep deprivation | standard chow diet | 3 |
| EC_HFD | exercise control | high-fat diet | 3 |
| SD_HFD | sleep deprivation | high-fat diet | 3 |

Metadata must be rechecked against current BioProject, SRA/ENA, BioSample, and publication records. A serious conflict in tissue, assay, independence, replication, group labels, or read structure is a stop condition.

## Primary validation hypothesis

Sleep deprivation is associated with negative enrichment of `HALLMARK_OXIDATIVE_PHOSPHORYLATION` across dietary backgrounds.

- Prespecified direction: **NES < 0**.
- Statistical unit: mouse/sample, never individual reads or cells.
- Gene-set release: the same native-mouse MSigDB release used in the frozen manuscript analysis, expected `MSigDB 2026.1.Mm`; the resolved release must be recorded.
- Ranking: complete tested gene list ranked by edgeR logFC for the prespecified contrast.
- Enrichment method: `fgseaMultilevel`; do not substitute a different primary pathway if the result is weak or opposite.

## Prespecified group factor and contrasts

Use all 12 biological samples and `design <- model.matrix(~0 + group)`, with `group` levels `EC_SCD`, `SD_SCD`, `EC_HFD`, and `SD_HFD`.

### Primary contrast

Average sleep-deprivation effect across both diet backgrounds:

`0.5 * SD_SCD + 0.5 * SD_HFD - 0.5 * EC_SCD - 0.5 * EC_HFD`

Diet interaction must not be mixed into the primary conclusion.

### Secondary contrasts

1. Standard-chow sleep-deprivation effect: `SD_SCD - EC_SCD`
2. High-fat-diet sleep-deprivation effect: `SD_HFD - EC_HFD`

### Exploratory interaction

`(SD_HFD - EC_HFD) - (SD_SCD - EC_SCD)`

The interaction is exploratory and must not redefine the primary endpoint.

## Prespecified differential-expression workflow

1. Import Salmon transcript estimates with the GENCODE Mouse M39 transcript-to-gene map and `tximport(countsFromAbundance = "no")`.
2. Fit edgeR quasi-likelihood models using `DGEList`, `filterByExpr`, `calcNormFactors`, `estimateDisp`, `glmQLFit`, and `glmQLFTest`.
3. Retain and report every sample unless a technical failure makes quantification impossible. Flag outliers; do not automatically remove them.
4. Save complete gene-level results for the primary, both diet-specific, and interaction contrasts.

## Prespecified pathway tests

### Primary

- `HALLMARK_OXIDATIVE_PHOSPHORYLATION`
- Applied to the primary average sleep-deprivation contrast.

### Five secondary mitochondrial pathways

1. `REACTOME_RESPIRATORY_ELECTRON_TRANSPORT`
2. `REACTOME_MITOCHONDRIAL_TRANSLATION`
3. `REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING`
4. `REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE`
5. `REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION`

These five pathways will be evaluated for the primary contrast and both diet-specific contrasts. The primary pathway plus these five pathways will also be evaluated for the exploratory interaction. No additional pathway discovery is permitted.

## Prespecified interpretation

- Primary OXPHOS NES < 0 and FDR < 0.05: second independent dataset supports reproducible SD-associated OXPHOS transcriptional suppression.
- NES < 0 and FDR ≥ 0.05: directionally consistent but statistically weaker external replication.
- NES > 0: the third eligible dataset does not reproduce the direction; report without concealment and reassess the manuscript conclusion.
- Opposite SCD and HFD directions or clear interaction evidence: dietary context may modify the transcriptional response.
- A non-significant interaction will be described as “no clear evidence of interaction in this dataset,” not as proof that diet has no effect.

Transcriptional suppression must not be described as mitochondrial dysfunction. Interaction does not establish mechanism, association does not establish causation, and whole-colon RNA-seq is not ISC-specific replication.

## Technical stage gates

### Stage 1: metadata QC

Confirm exactly 12 independent colon RNA-seq samples, n=3 per factorial group, clear labels, paired/single-end structure, platform, read length, and public FASTQ size. Stop before download if the design materially conflicts with the study record.

### Stage 2: one-sample pilot

Download one `EC_SCD` run only. Inspect R1/R2 structure, read count and length, FastQC results, adapter/quality profiles, and quantify with Salmon 2.7.0 using `reference/salmon_mouse_GRCm39/`, `-l A`, and the same mapping options used for PRJNA1167170. Record mapping rate, inferred library type, processed fragments, orphan rate, and compatible-fragment ratio.

Stop if the pilot is not recognizable as mouse RNA-seq, has clearly abnormal mapping, cannot be parsed reliably, or conflicts with group metadata.

### Stage 3 onward

Only after the pilot passes, process the remaining 11 samples with identical tools, reference, index, and parameters; then perform the frozen tximport, edgeR, pathway, and cross-dataset analyses.

## Cross-dataset validation

For every gene in the complete native-mouse Hallmark OXPHOS gene set, retain logFC from:

1. GSE289089 ISC
2. PRJNA1167170 proximal/whole colon validation dataset
3. PRJNA862187 primary average sleep-deprivation contrast

Report all pairwise Spearman correlations, the number and proportion negative in all three datasets, and pathway-level directions. Do not restrict the comparison to statistically significant genes.

## Prohibited deviations

No new pathway discovery, WGCNA, hub-gene analysis, machine learning, random forest, molecular docking, CellChat, pseudotime, pathway fishing, arbitrary DEG-selected pathways, removal of unsupportive results, or post hoc primary-endpoint changes.

## Reference and software targets

- Genome/transcriptome: GRCm39 / GENCODE Mouse M39
- Existing Salmon index: `reference/salmon_mouse_GRCm39/`
- GTF/transcript FASTA: `reference/gencode_M39/`
- Salmon target: 2.7.0
- Library inference: `-l A`
- Gene-level inference: tximport with `countsFromAbundance = "no"`
- Differential expression: edgeR quasi-likelihood framework
- Enrichment: fgseaMultilevel with native-mouse MSigDB 2026.1.Mm if confirmed available

This file is the immutable pre-analysis specification for PRJNA862187. Any unavoidable deviation must be documented prospectively before viewing the affected result.
