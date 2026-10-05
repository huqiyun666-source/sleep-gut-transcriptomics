#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: GSE289089_celltype_OXPHOS.R PROJECT_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(file.path(root, "reference", "R_library"), .libPaths()))

required <- c("Seurat", "SeuratObject", "Matrix", "edgeR", "fgsea",
              "msigdbr", "data.table", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(ggplot2)
})

object_path <- file.path(root, "processed", "GSE289089_clustered.rds")
features_path <- file.path(root, "data", "raw", "GSE289089", "GSM8783533_SD_0D_1_features.tsv.gz")
gtf_path <- file.path(root, "reference", "gencode_M39", "gencode.vM39.chr_patch_hapl_scaff.annotation.gtf.gz")
final_dir <- file.path(root, "results", "GSE289089_celltype_OXPHOS")
if (!file.exists(object_path)) stop("Seurat object is missing: ", object_path)
if (!file.exists(features_path)) stop("10x feature annotation is missing: ", features_path)
if (!file.exists(gtf_path)) stop("GENCODE M39 annotation is missing: ", gtf_path)
if (dir.exists(final_dir)) stop("Refusing to overwrite existing output directory: ", final_dir)
results_root <- file.path(root, "results")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
stage_dir <- tempfile(pattern = ".GSE289089_celltype_OXPHOS_staging_", tmpdir = results_root)
dir.create(stage_dir, recursive = TRUE)
dir.create(file.path(stage_dir, "DE"))
completed <- FALSE
on.exit(if (!completed) message("Incomplete staging output retained at: ", stage_dir), add = TRUE)

samples <- c("0D_1", "0D_2", "2D_1", "2D_2")
condition <- factor(c("Control", "Control", "SD", "SD"), levels = c("Control", "SD"))
names(condition) <- samples

# Prespecified mapping from the user request. Cluster 17 is recorded but not
# analyzed as an intestinal epithelial comparison.
mapping <- data.table(
  cell_group = c(
    "ISC", rep("TA_cycling_progenitor", 5), "Goblet_secretory",
    rep("Paneth_lineage", 2), rep("Enterocyte", 2), "Enteroendocrine",
    "Serotonergic_EC", "Tuft", "Stress_state_cluster0",
    "IFN_regenerative_cluster1", "EEC_progenitor_cluster14", "Immune_excluded"
  ),
  cluster = c(6, 2, 7, 10, 11, 13, 3, 5, 16, 4, 15, 9, 12, 8, 0, 1, 14, 17),
  analysis_role = c(
    rep("primary_epithelial_comparison", 14),
    rep("exploratory_state", 3), "excluded_non_epithelial"
  )
)
group_order <- unique(mapping[analysis_role != "excluded_non_epithelial", cell_group])
fwrite(mapping, file.path(stage_dir, "cluster_to_cell_group_mapping.csv"))

obj <- readRDS(object_path)
if (!inherits(obj, "Seurat")) stop("Object is not a Seurat object: ", object_path)
required_meta <- c("sample", "condition", "seurat_clusters")
if (!all(required_meta %in% colnames(obj[[]]))) {
  stop("Missing metadata field(s): ", paste(setdiff(required_meta, colnames(obj[[]])), collapse = ", "))
}
if (!setequal(unique(as.character(obj$sample)), samples)) stop("Unexpected sample labels in Seurat object")
observed_design <- unique(data.table(sample = as.character(obj$sample), condition = as.character(obj$condition)))
expected_design <- data.table(sample = samples, condition = as.character(condition))
setorder(observed_design, sample); setorder(expected_design, sample)
if (!identical(observed_design, expected_design)) stop("Sample-to-condition mapping differs from the prespecified design")
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA raw-count layer is unavailable")

# Harmonize legacy 10x feature symbols to the current GENCODE M39 symbols by
# stable Ensembl gene ID. Raw counts and the Seurat object are not modified.
feature_map <- fread(cmd = sprintf("gzip -dc %s", shQuote(features_path)), header = FALSE,
                     col.names = c("gene_id", "source_gene_symbol", "feature_type"))
# Seurat/Read10X uses make.unique() when source feature symbols are duplicated.
# Reproduce that deterministic naming so each Seurat row still maps to exactly
# one Ensembl gene ID.
feature_map[, seurat_gene_symbol := make.unique(source_gene_symbol)]
gtf_cmd <- sprintf("gzip -dc %s | awk -F '\\t' '$3 == \"gene\"'", shQuote(gtf_path))
gtf_gene <- fread(cmd = gtf_cmd, sep = "\t", header = FALSE, quote = "", fill = TRUE,
                  showProgress = FALSE)
current_map <- unique(data.table(
  gene_id = sub("\\..*$", "", sub('.*gene_id "([^"]+)".*', "\\1", gtf_gene$V9)),
  current_gene_symbol = sub('.*gene_name "([^"]+)".*', "\\1", gtf_gene$V9)
))
symbol_map <- merge(feature_map, current_map, by = "gene_id", all.x = TRUE, sort = FALSE)
symbol_map[is.na(current_gene_symbol) | current_gene_symbol == "", current_gene_symbol := source_gene_symbol]
if (symbol_map[, anyDuplicated(seurat_gene_symbol)] > 0L) stop("Seurat feature symbols remain ambiguous after make.unique")
symbol_lookup <- setNames(symbol_map$current_gene_symbol, symbol_map$seurat_gene_symbol)
unmapped_object_symbols <- setdiff(rownames(obj), names(symbol_lookup))
if (length(unmapped_object_symbols)) {
  symbol_lookup <- c(symbol_lookup, setNames(unmapped_object_symbols, unmapped_object_symbols))
}
symbol_map[, changed := source_gene_symbol != current_gene_symbol]
fwrite(symbol_map, file.path(stage_dir, "gene_symbol_harmonization.csv"))

obj$cell_group_OXPHOS <- NA_character_
for (i in seq_len(nrow(mapping))) {
  obj$cell_group_OXPHOS[as.character(obj$seurat_clusters) == as.character(mapping$cluster[i])] <- mapping$cell_group[i]
}

count_dt <- as.data.table(table(
  cell_group = factor(obj$cell_group_OXPHOS, levels = c(group_order, "Immune_excluded")),
  sample = factor(obj$sample, levels = samples)
))
setnames(count_dt, "N", "cell_count")
count_dt[, condition := ifelse(sample %in% c("0D_1", "0D_2"), "Control", "SD")]
count_dt <- merge(count_dt, unique(mapping[, .(cell_group, analysis_role)]), by = "cell_group", all.x = TRUE)
count_dt <- count_dt[order(match(cell_group, c(group_order, "Immune_excluded")),
                           match(sample, samples))]
fwrite(count_dt, file.path(stage_dir, "celltype_sample_counts.csv"))

hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
reactome <- msigdbr(db_species = "MM", species = "Mus musculus",
                    collection = "M2", subcollection = "CP:REACTOME")
msigdb_version <- unique(hallmark$db_version)
if (length(msigdb_version) != 1L) stop("Could not determine one MSigDB version")

primary_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
secondary_names <- c(
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION"
)
make_pathways <- function(tbl, wanted) {
  present <- intersect(wanted, unique(tbl$gs_name))
  split(tbl$gene_symbol[tbl$gs_name %in% present], tbl$gs_name[tbl$gs_name %in% present])
}
primary_pathway <- make_pathways(hallmark, primary_name)
secondary_pathways <- make_pathways(reactome, secondary_names)
if (length(primary_pathway) != 1L) stop("Primary OXPHOS set is unavailable")
if (!setequal(names(secondary_pathways), secondary_names)) stop("One or more secondary pathways are unavailable")
mitochondrial_pathways <- c(primary_pathway, secondary_pathways)

sanitize_aggregate_names <- function(x) gsub("-", "_", sub("^g", "", x))
make_rank <- function(de, method) {
  z <- copy(de[!is.na(current_gene_symbol) & current_gene_symbol != ""])
  score_col <- if (method == "logFC") "logFC" else "signed_stat"
  z <- z[is.finite(get(score_col))]
  setorder(z, current_gene_symbol, -logCPM)
  z <- z[!duplicated(current_gene_symbol)]
  ans <- z[[score_col]]
  names(ans) <- z$current_gene_symbol
  sort(ans, decreasing = TRUE)
}
run_fgsea <- function(pathways, stats, seed) {
  set.seed(seed)
  as.data.table(fgseaMultilevel(pathways, stats, minSize = 10, maxSize = 500, eps = 0))
}

pseudobulk_counts <- list()
edge_objects <- list()
library_rows <- list()
primary_rows <- list()
all_ranking_rows <- list()
secondary_rows <- list()
consistency_rows <- list()

for (g_idx in seq_along(group_order)) {
  group_name <- group_order[g_idx]
  group_cells <- colnames(obj)[obj$cell_group_OXPHOS == group_name]
  sub_obj <- subset(obj, cells = group_cells)
  pb <- AggregateExpression(sub_obj, assays = "RNA", group.by = "sample",
                            return.seurat = FALSE, verbose = FALSE)$RNA
  colnames(pb) <- sanitize_aggregate_names(colnames(pb))
  if (!setequal(colnames(pb), samples)) stop("Unexpected pseudobulk sample columns for ", group_name)
  pb <- pb[, samples, drop = FALSE]
  pseudobulk_counts[[group_name]] <- pb

  group_counts <- count_dt[cell_group == group_name]
  eligible <- nrow(group_counts) == 4L && all(group_counts$cell_count >= 20)
  status <- if (eligible) "tested" else "insufficient cell count"
  library_rows[[group_name]] <- data.table(
    cell_group = group_name,
    sample = samples,
    condition = as.character(condition),
    cell_count = group_counts$cell_count[match(samples, group_counts$sample)],
    pseudobulk_library_size = as.numeric(colSums(pb)),
    status = status
  )

  if (!eligible) {
    primary_rows[[group_name]] <- data.table(
      cell_group = group_name,
      n_control_cells = sum(group_counts[condition == "Control", cell_count]),
      n_SD_cells = sum(group_counts[condition == "SD", cell_count]),
      NES = NA_real_, pval = NA_real_, padj = NA_real_, size = NA_integer_,
      leadingEdge_count = NA_integer_, leadingEdge = NA_character_, status = status
    )
    next
  }

  y <- DGEList(counts = pb, samples = data.frame(sample = samples, condition = condition))
  design <- model.matrix(~condition)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "conditionSD")
  de <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene_symbol")
  if (!"F" %in% names(de)) stop("edgeR QL F statistic is absent for ", group_name)
  de[, current_gene_symbol := unname(symbol_lookup[gene_symbol])]
  de[is.na(current_gene_symbol) | current_gene_symbol == "", current_gene_symbol := gene_symbol]
  de[, signed_stat := sign(logFC) * sqrt(F)]
  setorder(de, PValue)
  fwrite(de[, .(gene_symbol, current_gene_symbol, logFC, logCPM, F, signed_stat, PValue, FDR)],
         file.path(stage_dir, "DE", paste0(group_name, "_edgeR_SD_vs_Control.csv")))
  edge_objects[[group_name]] <- list(y = y, design = design, fit = fit, qlf = qlf, de = de)

  ranks <- list(logFC = make_rank(de, "logFC"), signed_stat = make_rank(de, "signed_stat"))
  for (m_idx in seq_along(ranks)) {
    method <- names(ranks)[m_idx]
    fg <- run_fgsea(primary_pathway, ranks[[method]], 289089L + g_idx * 10L + m_idx)
    if (nrow(fg) != 1L) stop("OXPHOS fgsea did not return one row for ", group_name, " / ", method)
    lead <- fg$leadingEdge[[1]]
    all_ranking_rows[[paste(group_name, method)]] <- data.table(
      dataset = "GSE289089", cell_group = group_name, ranking_method = method,
      NES = fg$NES, pval = fg$pval, padj = fg$padj, size = fg$size,
      leadingEdge_count = length(lead), leadingEdge = paste(lead, collapse = ";"), status = status
    )
    if (method == "logFC") {
      primary_rows[[group_name]] <- data.table(
        cell_group = group_name,
        n_control_cells = sum(group_counts[condition == "Control", cell_count]),
        n_SD_cells = sum(group_counts[condition == "SD", cell_count]),
        NES = fg$NES, pval = fg$pval, padj = fg$padj, size = fg$size,
        leadingEdge_count = length(lead), leadingEdge = paste(lead, collapse = ";"), status = status
      )
      lead_de <- de[match(lead, current_gene_symbol)]
      leading_negative_fraction <- mean(lead_de$logFC < 0, na.rm = TRUE)
    }
  }

  sec <- run_fgsea(secondary_pathways, ranks$logFC, 389089L + g_idx)
  sec[, `:=`(
    dataset = "GSE289089", cell_group = group_name, ranking_method = "logFC",
    leadingEdge_count = lengths(leadingEdge),
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)),
    status = status
  )]
  primary_for_heat <- primary_rows[[group_name]][, .(
    pathway = primary_name, pval, padj, log2err = NA_real_, ES = NA_real_, NES,
    size, leadingEdge, dataset = "GSE289089", cell_group = group_name,
    ranking_method = "logFC", leadingEdge_count, status
  )]
  sec <- rbind(primary_for_heat, sec, fill = TRUE)
  missing_sets <- setdiff(names(mitochondrial_pathways), sec$pathway)
  if (length(missing_sets)) {
    missing_rows <- data.table(
      pathway = missing_sets,
      pval = NA_real_, padj = NA_real_, log2err = NA_real_, ES = NA_real_, NES = NA_real_,
      size = vapply(mitochondrial_pathways[missing_sets],
                    function(s) sum(unique(s) %in% names(ranks$logFC)), integer(1)),
      leadingEdge = NA_character_, dataset = "GSE289089", cell_group = group_name,
      ranking_method = "logFC", leadingEdge_count = NA_integer_,
      status = "not tested: matched size below prespecified minSize=10"
    )
    sec <- rbind(sec, missing_rows, fill = TRUE)
  }
  secondary_rows[[group_name]] <- sec[, .(dataset, cell_group, ranking_method, pathway,
                                           NES, pval, padj, size, leadingEdge_count,
                                           leadingEdge, status)]

  # Descriptive replicate-consistency check: mean normalized logCPM across the
  # fixed Hallmark OXPHOS genes, retaining biological samples as the units.
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  current_symbols <- unname(symbol_lookup[rownames(logcpm)])
  current_symbols[is.na(current_symbols) | current_symbols == ""] <- rownames(logcpm)[is.na(current_symbols) | current_symbols == ""]
  ox_rows <- rownames(logcpm)[current_symbols %in% primary_pathway[[primary_name]]]
  sample_score <- colMeans(logcpm[ox_rows, , drop = FALSE])
  pairwise_delta <- outer(sample_score[c("2D_1", "2D_2")],
                          sample_score[c("0D_1", "0D_2")], "-")
  consistency_rows[[group_name]] <- data.table(
    cell_group = group_name,
    matched_OXPHOS_genes = length(unique(current_symbols[current_symbols %in% primary_pathway[[primary_name]]])),
    mean_Control_OXPHOS_logCPM = mean(sample_score[c("0D_1", "0D_2")]),
    mean_SD_OXPHOS_logCPM = mean(sample_score[c("2D_1", "2D_2")]),
    SD_minus_Control_mean = mean(sample_score[c("2D_1", "2D_2")]) - mean(sample_score[c("0D_1", "0D_2")]),
    negative_pairwise_comparisons = sum(pairwise_delta < 0),
    total_pairwise_comparisons = length(pairwise_delta),
    pairwise_negative_fraction = mean(pairwise_delta < 0),
    leading_edge_negative_fraction = leading_negative_fraction,
    score_0D_1 = sample_score["0D_1"], score_0D_2 = sample_score["0D_2"],
    score_2D_1 = sample_score["2D_1"], score_2D_2 = sample_score["2D_2"]
  )
}

library_dt <- rbindlist(library_rows, use.names = TRUE, fill = TRUE)
primary_dt <- rbindlist(primary_rows, use.names = TRUE, fill = TRUE)
ranking_dt <- rbindlist(all_ranking_rows, use.names = TRUE, fill = TRUE)
secondary_dt <- rbindlist(secondary_rows, use.names = TRUE, fill = TRUE)
consistency_dt <- rbindlist(consistency_rows, use.names = TRUE, fill = TRUE)
library_dt <- library_dt[order(match(cell_group, group_order), match(sample, samples))]
primary_dt <- primary_dt[order(match(cell_group, group_order))]
ranking_dt <- ranking_dt[order(match(cell_group, group_order),
                               match(ranking_method, c("logFC", "signed_stat")))]
secondary_dt <- secondary_dt[order(match(cell_group, group_order),
                                   match(pathway, c(primary_name, secondary_names)))]
consistency_dt <- consistency_dt[order(match(cell_group, group_order))]

fwrite(library_dt, file.path(stage_dir, "pseudobulk_library_sizes.csv"))
fwrite(primary_dt, file.path(stage_dir, "celltype_OXPHOS_logFC_ranking.csv"))
fwrite(ranking_dt, file.path(stage_dir, "celltype_OXPHOS_all_rankings.csv"))
fwrite(secondary_dt, file.path(stage_dir, "celltype_mitochondrial_pathways.csv"))
fwrite(consistency_dt, file.path(stage_dir, "celltype_OXPHOS_replicate_consistency.csv"))
saveRDS(pseudobulk_counts, file.path(stage_dir, "pseudobulk_counts.rds"))
saveRDS(edge_objects, file.path(stage_dir, "celltype_edgeR_objects.rds"))
saveRDS(list(primary = primary_dt, all_rankings = ranking_dt,
             mitochondrial = secondary_dt, consistency = consistency_dt),
        file.path(stage_dir, "celltype_OXPHOS_results.rds"))

# General-purpose mechanism heatmap. Missing/insufficient values are explicitly
# grey and labeled; the diverging scale is centered at NES=0.
pathway_labels <- c(
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Hallmark OXPHOS",
  REACTOME_RESPIRATORY_ELECTRON_TRANSPORT = "Respiratory electron transport",
  REACTOME_MITOCHONDRIAL_TRANSLATION = "Mitochondrial translation",
  REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING = "ATP synthesis / chemiosmosis",
  REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE = "TCA cycle",
  REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION = "Mitochondrial FA beta-oxidation"
)
heat_grid <- CJ(cell_group = group_order,
                pathway = c(primary_name, secondary_names), unique = TRUE)
heat_dt <- merge(heat_grid, secondary_dt, by = c("cell_group", "pathway"), all.x = TRUE, sort = FALSE)
heat_dt[, cell_group := factor(cell_group, levels = rev(group_order))]
heat_dt[, pathway_label := factor(pathway_labels[pathway], levels = unname(pathway_labels))]
lim <- max(abs(heat_dt$NES), na.rm = TRUE)
p <- ggplot(heat_dt, aes(pathway_label, cell_group, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(
    is.na(NES),
    ifelse(grepl("matched size below", status), "gene set\n<10 genes", "insufficient\ncell count"),
    sprintf("%.2f", NES)
  )), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), na.value = "#BDBDBD",
                       name = "NES") +
  labs(x = NULL, y = NULL,
       title = "Cell-type pseudobulk mitochondrial pathway enrichment",
       subtitle = "GSE289089; SD vs Control; edgeR logFC ranking") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title.position = "plot")
ggsave(file.path(stage_dir, "celltype_mitochondrial_NES_heatmap.png"), p,
       width = 11, height = 7.5, dpi = 320, bg = "white")
ggsave(file.path(stage_dir, "celltype_mitochondrial_NES_heatmap.pdf"), p,
       width = 11, height = 7.5, device = "pdf", bg = "white")

sink(file.path(stage_dir, "GSE289089_celltype_reproducibility.txt"))
cat("Analysis: sample x cell-group pseudobulk OXPHOS specificity\n")
cat("Seurat object:", object_path, "\n")
cat("Gene-symbol harmonization: source 10x features", features_path, "to GENCODE M39", gtf_path, "by stable Ensembl gene ID.\n")
cat("Source symbols changed to current symbols:", sum(symbol_map$changed), "of", nrow(symbol_map), "features.\n")
cat("Biological samples: 0D_1, 0D_2 (Control); 2D_1, 2D_2 (SD)\n")
cat("No individual cell was treated as a biological replicate.\n")
cat("Minimum cell-count rule: every sample must have >=20 cells in a group for formal edgeR inference.\n")
cat("edgeR: DGEList -> filterByExpr -> calcNormFactors -> estimateDisp -> glmQLFit -> glmQLFTest\n")
cat("Design: ~ condition; Control reference; positive logFC means higher in SD.\n")
cat("GSEA ranking methods: logFC; sign(logFC)*sqrt(edgeR QL F).\n")
cat("GSEA engine: fgseaMultilevel; minSize=10; maxSize=500; eps=0.\n")
cat("Gene-set source: MSigDB", msigdb_version, "native mouse via msigdbr", as.character(packageVersion("msigdbr")), "\n")
cat("Cluster mapping file: cluster_to_cell_group_mapping.csv\n")
cat("Procedural guidance: Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065\n")
cat("Package versions:\n")
for (pkg in required) cat("  ", pkg, ": ", as.character(packageVersion(pkg)), "\n", sep = "")
cat("\nSession information:\n")
print(sessionInfo())
sink()

if (!file.rename(stage_dir, final_dir)) stop("Could not atomically move staging output to final directory")
completed <- TRUE
message("GSE289089 cell-type OXPHOS analysis completed: ", final_dir)
