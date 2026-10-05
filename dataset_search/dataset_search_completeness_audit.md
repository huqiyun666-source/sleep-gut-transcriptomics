# Dataset-search completeness audit

Search cutoff: **2026-10-04**  
Audit performed: **2026-10-05**  
Scope: verification and completeness checking only; no sequence or expression data were downloaded or analysed.

## Outcome

The focused audit added one unique candidate study but did not add an eligible dataset. The final registry contains **23 unique candidates: 3 ELIGIBLE, 1 POTENTIALLY_ELIGIBLE, and 19 EXCLUDED**.

## PMID 41094491 / PRJNA1395688

Wang et al., *Nuclear receptor Nr1d1 links sleep deprivation to intestinal homeostasis via microbiota-derived taurine* (PMID 41094491; DOI 10.1186/s12967-025-07089-8), was already present in the registry through **PRJNA1395688**. It is therefore a duplicate publication/accession representation of an existing study-level candidate, not a newly identified candidate.

The NCBI BioProject record identifies nine SRA experiments and nine BioSamples with raw sequence reads. The repository record is titled as colon RNA-seq, and the publication describes intestinal epithelial RNA sequencing. However, the linked SRA library fields are internally inconsistent (`AMPLICON / PCR / GENOMIC`). Because the frozen rules prohibit resolving contradictory metadata by inference, the study remains **POTENTIALLY_ELIGIBLE**. Public raw-read availability alone was not treated as sufficient evidence that the archived libraries can be validly processed as RNA-seq.

Sources:

- https://pubmed.ncbi.nlm.nih.gov/41094491/
- https://link.springer.com/article/10.1186/s12967-025-07089-8
- https://www.ncbi.nlm.nih.gov/bioproject/1395688

## Focused query audit

The following exact query concepts were rerun across article and repository search interfaces, with exact-title, DOI, accession, Data Availability, supplement, backward-citation, forward-citation, and related-article checks for plausible records:

1. `"sleep deprivation" "intestinal epithelial" RNA-seq`
2. `"sleep deprivation" "intestinal epithelium" transcriptome`
3. `"sleep deprivation" colon transcriptomics mouse`
4. `"sleep deprivation" intestinal epithelial cells transcriptome`
5. `"chronic sleep deprivation" intestinal RNA sequencing`
6. `"sleep deficiency" intestinal RNA-seq`
7. `"sleep loss" intestinal epithelial transcriptome`
8. `"sleep deprivation" gut "RNA sequencing"`

Resources checked were GEO, SRA, BioProject, PubMed, Europe PMC, ENA, and BioStudies/ArrayExpress. Overlapping repository, publication, and preprint records were merged at the study level.

## Newly added candidate

| Study | Transcriptomic relevance | Public reusable data at cutoff | Decision |
|---|---|---|---|
| Liu et al., *Multi-omics analysis reveals that sleep deprivation exacerbates intestinal barrier injury via the microbiota-dependent taurodeoxycholic acid/group 3 innate lymphoid/IL-22 axis* (2026), DOI 10.1186/s12967-026-08696-9 | Mouse CPW sleep-deprivation experiment (0–96 h) with intestinal RNA sequencing reported | No GEO, SRA, BioProject, ENA, or BioStudies accession and no reusable expression matrix located | **EXCLUDED — E6: data not publicly accessible at search cutoff** |

The report was published on 29 July 2026 and was therefore within the search cutoff. It changed the candidate and exclusion counts but did not change the eligible-dataset set, the frozen analyses, or the scientific conclusions.

Source: https://link.springer.com/article/10.1186/s12967-026-08696-9

## Related-article false positives checked

- The 2026 ApcMin/+ colorectal-cancer report (DOI 10.7150/ijbs.134241) included transcriptomic experiments in treated HCT116 cancer cells, not an intestinal-tissue sleep-deprivation-versus-control transcriptome. It was not promoted to the registry.
- The 2025 KynA/P4HA2/HIF-1alpha colon-cancer report (DOI 10.1016/j.molmet.2025.102109) used transcriptomics of KynA-treated LoVo cells rather than an intestinal-tissue sleep-deprivation-versus-control contrast. It was not promoted to the registry.
- Reports in ovary, heart, brain, adipose tissue, lacrimal gland, or microbiome-only assays were outside the predefined intestinal host-transcriptome eligibility frame; plausible records already represented in the registry retained their documented exclusions.

## Completeness conclusion

No fourth eligible public intestinal sleep-deprivation transcriptomic dataset was identified by the cutoff. The eligible studies remain GSE289089/PRJNA1221288, PRJNA1167170, and PRJNA862187. PRJNA1395688 remains the sole potentially eligible study. The newly identified Liu et al. report is retained as an E6 exclusion so that its absence from analysis is explicit rather than silent.

## Boundary

This was a targeted public-dataset identification and selection audit, not a PRISMA systematic review. Search-engine indexing can change after the cutoff; later release or correction of repository metadata would require a future registry update but does not alter this frozen cutoff assessment.
