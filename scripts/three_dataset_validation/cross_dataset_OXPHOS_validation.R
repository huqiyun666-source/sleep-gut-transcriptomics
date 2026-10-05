#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: cross_dataset_OXPHOS_validation.R PROJECT_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(file.path(root, "reference", "R_library"), .libPaths()))

required <- c("Seurat", "SeuratObject", "edgeR", "fgsea", "msigdbr", "data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(data.table)
})

gse_dir <- file.path(root, "results", "GSE289089_celltype_OXPHOS")
bulk_dir <- file.path(root, "results", "PRJNA1167170_bulk")
final_dir <- file.path(root, "results", "cross_dataset_OXPHOS")
required_inputs <- c(
  file.path(gse_dir, "celltype_edgeR_objects.rds"),
  file.path(gse_dir, "celltype_OXPHOS_all_rankings.csv"),
  file.path(gse_dir, "celltype_sample_counts.csv"),
  file.path(bulk_dir, "tximport_object.rds"),
  file.path(bulk_dir, "edgeR_SD_vs_Control.csv"),
  file.path(bulk_dir, "primary_OXPHOS_validation.csv"),
  file.path(bulk_dir, "OXPHOS_leading_edge.csv")
)
if (!all(file.exists(required_inputs))) stop("Missing required input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
if (dir.exists(final_dir)) stop("Refusing to overwrite existing output directory: ", final_dir)
stage_dir <- tempfile(pattern = ".cross_dataset_OXPHOS_staging_", tmpdir = file.path(root, "results"))
dir.create(stage_dir, recursive = TRUE)
completed <- FALSE
on.exit(if (!completed) message("Incomplete staging output retained at: ", stage_dir), add = TRUE)

hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
reactome <- msigdbr(db_species = "MM", species = "Mus musculus",
                    collection = "M2", subcollection = "CP:REACTOME")
gocc <- msigdbr(db_species = "MM", species = "Mus musculus",
                collection = "M5", subcollection = "GO:CC")
msigdb_version <- unique(hallmark$db_version)
primary_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
primary_pathway <- split(hallmark$gene_symbol[hallmark$gs_name == primary_name], primary_name)
if (length(primary_pathway) != 1L) stop("Native-mouse Hallmark OXPHOS is unavailable")

make_rank <- function(de, method) {
  symbol_col <- if ("current_gene_symbol" %in% names(de)) "current_gene_symbol" else "gene_symbol"
  z <- copy(de[!is.na(get(symbol_col)) & get(symbol_col) != ""])
  score_col <- if (method == "logFC") "logFC" else "signed_stat"
  if (!score_col %in% names(z)) stop("Missing ranking column: ", score_col)
  z <- z[is.finite(get(score_col))]
  setorderv(z, c(symbol_col, "logCPM"), c(1L, -1L))
  z <- z[!duplicated(get(symbol_col))]
  ans <- z[[score_col]]
  names(ans) <- z[[symbol_col]]
  sort(ans, decreasing = TRUE)
}
run_one <- function(de, dataset, group, method, seed) {
  stats <- make_rank(de, method)
  set.seed(seed)
  fg <- as.data.table(fgseaMultilevel(primary_pathway, stats, minSize = 10,
                                      maxSize = 500, eps = 0))
  if (nrow(fg) != 1L) stop("OXPHOS fgsea did not return one result for ", dataset, " / ", group, " / ", method)
  data.table(
    dataset = dataset, cell_group = group, ranking_method = method,
    NES = fg$NES, pval = fg$pval, padj = fg$padj, size = fg$size,
    leadingEdge_count = length(fg$leadingEdge[[1]]),
    leadingEdge = paste(fg$leadingEdge[[1]], collapse = ";")
  )
}

# GSE289089: use the already fixed outputs from the cell-type script so the
# sensitivity table is numerically identical to the primary logFC analysis.
# The QL F statistic was retained directly from glmQLFTest in each edgeR object.
gse_objects <- readRDS(file.path(gse_dir, "celltype_edgeR_objects.rds"))
sensitivity_rows <- list(GSE289089 = fread(file.path(gse_dir, "celltype_OXPHOS_all_rankings.csv")))
for (i in seq_along(gse_objects)) {
  group_name <- names(gse_objects)[i]
  de <- as.data.table(gse_objects[[i]]$de)
  if (!"F" %in% names(de)) stop("GSE289089 QL F statistic is missing for ", group_name)
  expected_signed <- sign(de$logFC) * sqrt(de$F)
  if (!isTRUE(all.equal(de$signed_stat, expected_signed))) stop("Stored signed statistic is inconsistent for ", group_name)
}

# PRJNA1167170: refit the already imported 8-sample gene-count object to obtain
# the original edgeR QL F statistic; do not derive or invent it from p-values.
txi <- readRDS(file.path(bulk_dir, "tximport_object.rds"))
bulk_samples <- data.frame(
  sample = c("con1", "con2", "con3", "con4", "SD1", "SD2", "SD3", "SD4"),
  condition = factor(c(rep("Control", 4), rep("SD", 4)), levels = c("Control", "SD"))
)
rownames(bulk_samples) <- bulk_samples$sample
if (!identical(colnames(txi$counts), bulk_samples$sample)) stop("PRJNA1167170 count-column order differs from the prespecified design")
y_bulk <- DGEList(counts = txi$counts, samples = bulk_samples)
design_bulk <- model.matrix(~condition, data = bulk_samples)
keep_bulk <- filterByExpr(y_bulk, design = design_bulk)
y_bulk <- y_bulk[keep_bulk, , keep.lib.sizes = FALSE]
y_bulk <- calcNormFactors(y_bulk)
y_bulk <- estimateDisp(y_bulk, design_bulk)
fit_bulk <- glmQLFit(y_bulk, design_bulk)
qlf_bulk <- glmQLFTest(fit_bulk, coef = "conditionSD")
bulk_de <- as.data.table(topTags(qlf_bulk, n = Inf, sort.by = "none")$table, keep.rownames = "gene_id")
bulk_annot <- fread(file.path(bulk_dir, "edgeR_SD_vs_Control.csv"))[, .(gene_id, gene_symbol)]
bulk_de <- merge(bulk_annot, bulk_de, by = "gene_id", all.y = TRUE, sort = FALSE)
if (!"F" %in% names(bulk_de)) stop("PRJNA1167170 QL F statistic could not be extracted")
bulk_de[, signed_stat := sign(logFC) * sqrt(F)]
bulk_primary <- fread(file.path(bulk_dir, "primary_OXPHOS_validation.csv"))
if (nrow(bulk_primary) != 1L || bulk_primary$pathway != primary_name) stop("Unexpected PRJNA1167170 primary result format")
sensitivity_rows[["whole_colon logFC"]] <- data.table(
  dataset = "PRJNA1167170", cell_group = "whole_colon", ranking_method = "logFC",
  NES = bulk_primary$NES, pval = bulk_primary$pval, padj = bulk_primary$padj,
  size = bulk_primary$size,
  leadingEdge_count = length(strsplit(bulk_primary$leadingEdge, ";", fixed = TRUE)[[1]]),
  leadingEdge = bulk_primary$leadingEdge, status = "tested"
)
sensitivity_rows[["whole_colon signed_stat"]] <- run_one(bulk_de, "PRJNA1167170", "whole_colon", "signed_stat", 1167171L)

sensitivity <- rbindlist(sensitivity_rows, use.names = TRUE, fill = TRUE)
dataset_order <- c("GSE289089", "PRJNA1167170")
method_order <- c("logFC", "signed_stat")
sensitivity <- sensitivity[order(match(dataset, dataset_order), cell_group,
                                   match(ranking_method, method_order))]
fwrite(sensitivity, file.path(stage_dir, "OXPHOS_ranking_sensitivity.csv"))

# Leading edges are defined from the prespecified logFC ranking in each dataset.
gse_isc_result <- sensitivity[dataset == "GSE289089" & cell_group == "ISC" & ranking_method == "logFC"]
bulk_result <- sensitivity[dataset == "PRJNA1167170" & cell_group == "whole_colon" & ranking_method == "logFC"]
if (nrow(gse_isc_result) != 1L || nrow(bulk_result) != 1L) stop("Could not identify the two primary leading-edge results")
gse_lead <- strsplit(gse_isc_result$leadingEdge, ";", fixed = TRUE)[[1]]
bulk_lead <- strsplit(bulk_result$leadingEdge, ";", fixed = TRUE)[[1]]

gse_isc_de <- as.data.table(gse_objects[["ISC"]]$de)
gse_lead_dt <- merge(data.table(current_gene_symbol = gse_lead),
                     gse_isc_de[, .(source_gene_symbol = gene_symbol, current_gene_symbol,
                                    logFC, logCPM, F, signed_stat, PValue, FDR)],
                     by = "current_gene_symbol", all.x = TRUE, sort = FALSE)
setnames(gse_lead_dt, "current_gene_symbol", "gene_symbol")
bulk_lead_dt <- merge(data.table(gene_symbol = bulk_lead),
                      bulk_de[, .(gene_id, gene_symbol, logFC, logCPM, F, signed_stat, PValue, FDR)],
                      by = "gene_symbol", all.x = TRUE, sort = FALSE)
fwrite(gse_lead_dt, file.path(stage_dir, "GSE289089_ISC_OXPHOS_leading_edge.csv"))
fwrite(bulk_lead_dt, file.path(stage_dir, "PRJNA1167170_OXPHOS_leading_edge.csv"))

existing_bulk <- fread(file.path(bulk_dir, "OXPHOS_leading_edge.csv"))
existing_format_valid <- "gene_symbol" %in% names(existing_bulk)
existing_set_verified <- existing_format_valid && setequal(unique(existing_bulk$gene_symbol), unique(bulk_lead))

overlap <- intersect(unique(gse_lead), unique(bulk_lead))
union_genes <- union(unique(gse_lead), unique(bulk_lead))
overlap_summary <- data.table(
  n_GSE289089 = length(unique(gse_lead)),
  n_PRJNA1167170 = length(unique(bulk_lead)),
  n_overlap = length(overlap),
  n_union = length(union_genes),
  percentage_of_GSE289089_overlapped = 100 * length(overlap) / length(unique(gse_lead)),
  percentage_of_PRJNA1167170_overlapped = 100 * length(overlap) / length(unique(bulk_lead)),
  Jaccard = length(overlap) / length(union_genes),
  existing_PRJNA1167170_format_valid = existing_format_valid,
  existing_PRJNA1167170_leading_edge_set_verified = existing_set_verified
)
fwrite(overlap_summary, file.path(stage_dir, "OXPHOS_leading_edge_overlap.csv"))
fwrite(data.table(gene_symbol = union_genes,
                  in_GSE289089_ISC = union_genes %in% gse_lead,
                  in_PRJNA1167170_colon = union_genes %in% bulk_lead,
                  shared = union_genes %in% overlap),
       file.path(stage_dir, "OXPHOS_leading_edge_union.csv"))

# Functional categories are assigned from explicit native-mouse MSigDB GO/Reactome
# membership, not from gene-symbol prefixes alone. The four Complex II core
# subunits are an explicit Reactome-curated set (SDHA/B/C/D; R-HSA-70990,
# inferred to Mus musculus) because MSigDB 2026.1.Mm does not expose that
# Reactome complex as a standalone gene set.
get_set <- function(tbl, name) unique(tbl$gene_symbol[tbl$gs_name == name])
annotation_sets <- list(
  "Complex I / NDUF" = get_set(gocc, "GOCC_NADH_DEHYDROGENASE_COMPLEX"),
  "Complex III / UQCR" = get_set(gocc, "GOCC_RESPIRATORY_CHAIN_COMPLEX_III"),
  "Complex IV / COX" = get_set(gocc, "GOCC_RESPIRATORY_CHAIN_COMPLEX_IV"),
  "ATP synthase / ATP5" = get_set(gocc, "GOCC_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX"),
  "TCA" = get_set(reactome, "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"),
  "fatty acid oxidation" = get_set(reactome, "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION"),
  "mitochondrial translation" = get_set(reactome, "REACTOME_MITOCHONDRIAL_TRANSLATION")
)
complex_II_core <- c("Sdha", "Sdhb", "Sdhc", "Sdhd")
annotation_sets <- append(annotation_sets,
                          list("Complex II / SDH" = complex_II_core),
                          after = 1L)
hallmark_oxphos <- primary_pathway[[primary_name]]
priority <- c("Complex I / NDUF", "Complex II / SDH", "Complex III / UQCR",
              "Complex IV / COX", "ATP synthase / ATP5", "TCA",
              "fatty acid oxidation", "mitochondrial translation")
evidence_name <- c(
  "Complex I / NDUF" = "GOCC_NADH_DEHYDROGENASE_COMPLEX",
  "Complex II / SDH" = "Reactome SDH complex R-HSA-70990 (core SDHA/SDHB/SDHC/SDHD; inferred to Mus musculus)",
  "Complex III / UQCR" = "GOCC_RESPIRATORY_CHAIN_COMPLEX_III",
  "Complex IV / COX" = "GOCC_RESPIRATORY_CHAIN_COMPLEX_IV",
  "ATP synthase / ATP5" = "GOCC_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX",
  "TCA" = "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "fatty acid oxidation" = "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "mitochondrial translation" = "REACTOME_MITOCHONDRIAL_TRANSLATION"
)
functional <- rbindlist(lapply(overlap, function(gene) {
  matched <- priority[vapply(annotation_sets[priority], function(s) gene %in% s, logical(1))]
  primary_category <- if (length(matched)) matched[1] else if (gene %in% hallmark_oxphos) "other mitochondrial metabolism" else "unclassified"
  data.table(
    gene_symbol = gene,
    functional_category = primary_category,
    all_matching_categories = if (length(matched)) paste(matched, collapse = ";") else primary_category,
    annotation_evidence = if (length(matched)) paste(unname(evidence_name[matched]), collapse = ";") else if (gene %in% hallmark_oxphos) primary_name else NA_character_
  )
}))
functional <- merge(functional,
                    gse_lead_dt[, .(gene_symbol, GSE289089_ISC_logFC = logFC, GSE289089_ISC_FDR = FDR)],
                    by = "gene_symbol", all.x = TRUE, sort = FALSE)
functional <- merge(functional,
                    bulk_lead_dt[, .(gene_symbol, PRJNA1167170_logFC = logFC, PRJNA1167170_FDR = FDR)],
                    by = "gene_symbol", all.x = TRUE, sort = FALSE)
fwrite(functional, file.path(stage_dir, "OXPHOS_overlap_functional_annotation.csv"))

saveRDS(list(sensitivity = sensitivity, gse_leading_edge = gse_lead_dt,
             bulk_leading_edge = bulk_lead_dt, overlap = overlap_summary,
             functional_annotation = functional, bulk_edgeR_QLF = qlf_bulk),
        file.path(stage_dir, "cross_dataset_OXPHOS_results.rds"))

sink(file.path(stage_dir, "reproducibility_info.txt"))
cat("Analysis: cross-dataset OXPHOS ranking sensitivity and leading-edge overlap\n")
cat("GSE289089 Seurat object:", file.path(root, "processed", "GSE289089_clustered.rds"), "\n")
cat("GSE289089 edgeR objects:", file.path(gse_dir, "celltype_edgeR_objects.rds"), "\n")
cat("GSE289089 gene-symbol harmonization:", file.path(gse_dir, "gene_symbol_harmonization.csv"), "\n")
cat("Cluster-to-cell-group mapping:", file.path(gse_dir, "cluster_to_cell_group_mapping.csv"), "\n")
cat("PRJNA1167170 tximport object:", file.path(bulk_dir, "tximport_object.rds"), "\n")
cat("Biological replication: GSE289089 2 Control vs 2 SD; PRJNA1167170 4 Control vs 4 SD.\n")
cat("Ranking methods: edgeR logFC; sign(logFC)*sqrt(glmQLFTest F).\n")
cat("GSEA: fgseaMultilevel; HALLMARK_OXIDATIVE_PHOSPHORYLATION; minSize=10; maxSize=500; eps=0.\n")
cat("MSigDB release:", msigdb_version, "native mouse via msigdbr", as.character(packageVersion("msigdbr")), "\n")
cat("Functional annotation source: native-mouse MSigDB GO:CC and Reactome set membership; gene-name string matching was not used as the sole classifier.\n")
cat("Complex II core annotation: Reactome SDH complex R-HSA-70990, https://reactome.org/content/detail/R-HSA-70990 (SDHA/SDHB/SDHC/SDHD; inferred to Mus musculus).\n")
cat("Existing PRJNA1167170 leading-edge set verified against recomputation:", existing_set_verified, "\n")
cat("Procedural guidance: Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065\n")
cat("Package versions:\n")
for (pkg in required) cat("  ", pkg, ": ", as.character(packageVersion(pkg)), "\n", sep = "")
cat("\nSession information:\n")
print(sessionInfo())
sink()

if (!file.rename(stage_dir, final_dir)) stop("Could not atomically move staging output to final directory")
completed <- TRUE
message("Cross-dataset OXPHOS validation completed: ", final_dir)
