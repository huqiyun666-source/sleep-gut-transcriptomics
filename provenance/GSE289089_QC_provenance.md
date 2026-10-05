# GSE289089 QC provenance

Update date: 2026-10-04  
Scope: manuscript provenance correction based on the recovered contemporaneous analysis record. No QC filtering or manuscript result was recomputed.

## Provenance status

The QC parameters below were recovered from a record created during the original analysis period. They are therefore labelled **CONFIRMED_FROM_CONTEMPORANEOUS_RECORD**. The original preprocessing script is not present in the current project filesystem, and these parameters were not automatically recovered from a Seurat object.

## Confirmed parameters

| Method item | Value | Status | Manuscript-usable | Provenance note |
|---|---:|---|---|---|
| `CreateSeuratObject(min.cells)` | `3` | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Recovered from the contemporaneous analysis record; not from the current filesystem script collection or object history. |
| `CreateSeuratObject(min.features)` | `200` | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Recovered from the contemporaneous analysis record; not from object history. |
| Minimum `nFeature_RNA` boundary | `>= 500` | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Formal QC threshold. |
| Maximum `nFeature_RNA` boundary | `<= 6500` | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Formal QC threshold. |
| `percent.mt` cutoff | `< 25%` | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Formal QC threshold. |
| Hard `nCount_RNA` cutoff | None | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | No hard `nCount_RNA` threshold was used. |
| Thresholds across samples | Same thresholds for all four samples | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Applies to 0D_1, 0D_2, 2D_1, and 2D_2. |
| Post-QC cell set | 30,580 cells | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Formal post-QC cell count associated with the recovered filtering record. |
| Whether filtering was executed as separate per-sample subset calls | UNAVAILABLE | UNAVAILABLE | No | The recovered record confirms common thresholds, but does not preserve the exact execution structure. |
| `PercentageFeatureSet` pattern | UNAVAILABLE | UNAVAILABLE | No | The exact call/pattern is absent from the recovered contemporaneous record and current script collection. |

## Formal QC rule

The confirmed QC rule was:

```text
nFeature_RNA >= 500
nFeature_RNA <= 6500
percent.mt < 25
```

No hard `nCount_RNA` cutoff was applied, and the same thresholds were used for all four samples. This produced the formal 30,580-cell post-QC set.

## Object-based consistency check

The archived pre-QC and post-QC objects independently show 36,043 and 30,580 cells, respectively. Applying the confirmed three-variable rule to the pre-QC metadata reproduces the saved post-QC cell IDs overall and within each sample. This is corroborating evidence only; it is not the provenance source used to assign the confirmed status.

## Remaining unavailable information

The exact `PercentageFeatureSet` invocation and whether the QC expression was executed separately for each sample cannot be recovered from the archived record. Manuscript wording for either unavailable detail should state:

> The exact value or invocation could not be recovered from the archived analysis record.

