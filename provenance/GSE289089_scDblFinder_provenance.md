# GSE289089 scDblFinder provenance

Update date: 2026-10-04  
Scope: manuscript provenance correction based on the recovered contemporaneous analysis record. scDblFinder was not rerun.

## Recovered invocation

The contemporaneous analysis record preserves the following actual call:

```r
scDblFinder(sce, samples = "sample")
```

This call establishes that the input object was named `sce`, sample identity was supplied through the `samples` argument using the `sample` column, and all other arguments were left at package defaults. It does not preserve the package version or the numerical values selected internally by that installed package version.

## Parameter-level provenance

| Item | Value | Status | Manuscript-usable | Provenance note |
|---|---|---|---|---|
| Method / invocation | `scDblFinder(sce, samples = "sample")` | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Actual call recovered from the original analysis-period record; not reconstructed from cell counts. |
| Input object type | `SingleCellExperiment` (`sce`) | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | The `scDblFinder` call receives `sce`. |
| `samples` argument | metadata column `sample` | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | Explicitly specified as `samples = "sample"`. |
| Other arguments | Left at package defaults | CONFIRMED_FROM_CONTEMPORANEOUS_RECORD | Yes | No other arguments appear in the recovered call. |
| Package version | UNAVAILABLE | UNAVAILABLE | No | The exact version could not be recovered from the archived analysis record. |
| Numerical `dbr` / expected doublet rate | UNAVAILABLE | UNAVAILABLE | No | The argument was not manually supplied; the numerical default used by the unrecovered package version cannot be stated. |
| Numerical/explicit `clusters` value | UNAVAILABLE | UNAVAILABLE | No | The argument was not manually supplied; no resolved internal value was preserved. |
| Random seed | UNAVAILABLE | UNAVAILABLE | No | No seed-setting record associated with this call was recovered. |
| Output classification field | `scDblFinder.class` | CONFIRMED_FROM_OBJECT | Yes | Present in the archived clean and clustered Seurat-object metadata. |

## Cell-count consistency check

| Sample | Post-QC | Retained singlets | Doublets removed |
|---|---:|---:|---:|
| 0D_1 | 7,700 | 6,972 | 728 |
| 0D_2 | 7,274 | 6,672 | 602 |
| 2D_1 | 7,437 | 6,704 | 733 |
| 2D_2 | 8,169 | 7,323 | 846 |
| **Total** | **30,580** | **27,671** | **2,909** |

These values are an arithmetic and object-transition consistency check only. They were not used to infer any invocation parameter.

## Manuscript-ready wording

> Doublets were identified using scDblFinder with sample identity specified through the `samples` argument; other arguments were left at package defaults.

The package version, numerical default doublet rate, resolved `clusters` value, and random seed must not be stated because they remain unavailable.

