# Public-dataset selection flow

Search cutoff: **2026-10-04**

This is a dataset-identification and selection audit, not a PRISMA systematic review.

```text
Database searches and accession/citation checks
  GEO (9 hits), BioProject (9), SRA focused searches
  (118 colon; 10 intestine; first 500/741 broad gut; 18 RNA-Seq-filtered),
  PubMed (27), Europe PMC (21), ENA (5 accession checks),
  and BioStudies/ArrayExpress (first 100 broad-query hits plus accession checks)
                              |
                              v
Plausible reports/accessions retained after title/metadata screening
and merged at the study level
                              |
                              v
23 unique candidate studies underwent full metadata assessment
                              |
               +--------------+--------------+
               |                             |
               v                             v
19 excluded                                  4 retained
  E1 not intestinal tissue: 7                  |
  E2 no SD/SR contrast: 2                      +--> 3 ELIGIBLE
  E4 non-transcriptomic only: 3                |      GSE289089 / PRJNA1221288
  E6 not publicly accessible: 5                |      PRJNA1167170
  E7 group labels unresolved: 1                |      PRJNA862187
  E10 no reusable expression data: 1           |
                                                +--> 1 POTENTIALLY_ELIGIBLE
                                                       PRJNA1395688
                                                       (library-metadata conflict)
```

Database hit counts overlap and therefore were not summed as independent records. GEO, SRA, BioProject, ENA, publication, and preprint records describing the same experiment were merged into one study-level record. No FASTQ files or candidate expression matrices were downloaded, and no candidate dataset was reanalysed.
