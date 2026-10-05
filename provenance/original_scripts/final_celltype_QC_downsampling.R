#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: final_celltype_QC_downsampling.R PROJECT_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(file.path(root, "reference", "R_library"), .libPaths()))

required <- c("Seurat", "SeuratObject", "Matrix", "edgeR", "fgsea", "msigdbr", "data.table", "ggplot2")
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

seed <- 20261003L
iterations <- 20L
out_dir <- file.path(root, "results", "final_mito_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
outputs <- file.path(out_dir, c(
  "celltype_sample_QC_summary.csv", "celltype_QC_group_effects.csv",
  "downsampled_OXPHOS_sensitivity.csv", "downsampled_OXPHOS_iterations.csv",
  "Tuft_OXPHOS_robustness.csv", "celltype_QC_summary.png", "celltype_QC_summary.pdf",
  "downsampling_NES_stability.png", "downsampling_NES_stability.pdf"
))
if (any(file.exists(outputs))) stop("Refusing to overwrite existing output(s): ", paste(outputs[file.exists(outputs)], collapse = ", "))

object_path <- file.path(root, "processed", "GSE289089_clustered.rds")
mapping_path <- file.path(root, "results", "GSE289089_celltype_OXPHOS", "gene_symbol_harmonization.csv")
edge_path <- file.path(root, "results", "GSE289089_celltype_OXPHOS", "celltype_edgeR_objects.rds")
nuclear_path <- file.path(out_dir, "nuclear_only_OXPHOS_genes.csv")
inputs <- c(object_path, mapping_path, edge_path, nuclear_path)
if (!all(file.exists(inputs))) stop("Missing input(s): ", paste(inputs[!file.exists(inputs)], collapse = ", "))

samples <- c("0D_1", "0D_2", "2D_1", "2D_2")
sample_condition <- setNames(c("Control", "Control", "SD", "SD"), samples)
groups <- data.table(
  cell_group = c("ISC", "TA", "Tuft", "Serotonergic EC", "Enterocyte", "Goblet", "Paneth"),
  object_group = c("ISC", "TA_cycling_progenitor", "Tuft", "Serotonergic_EC", "Enterocyte", "Goblet_secretory", "Paneth_lineage"),
  clusters = c("6", "2;7;10;11;13", "8", "12", "4;15", "3", "5;16"),
  analysis_role = c(rep("primary", 6), "secondary")
)
downsample_groups <- groups[analysis_role == "primary"]

obj <- readRDS(object_path)
meta <- as.data.table(obj[[]], keep.rownames = "cell")
if (!all(c("sample", "condition", "seurat_clusters", "nFeature_RNA", "nCount_RNA", "percent.mt") %in% names(meta))) {
  stop("Required Seurat metadata fields are missing")
}
cluster_to_group <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
  data.table(cluster = strsplit(groups$clusters[i], ";", fixed = TRUE)[[1]], cell_group = groups$cell_group[i])
}))
meta[, cluster := as.character(seurat_clusters)]
meta <- merge(meta, cluster_to_group, by = "cluster", all.x = FALSE, all.y = FALSE, sort = FALSE)
meta[, sample := as.character(sample)]
meta[, condition := unname(sample_condition[sample])]

qc <- meta[, .(
  n_cells = .N,
  median_nFeature_RNA = as.numeric(median(nFeature_RNA, na.rm = TRUE)),
  median_nCount_RNA = as.numeric(median(nCount_RNA, na.rm = TRUE)),
  median_percent_mt = as.numeric(median(percent.mt, na.rm = TRUE)),
  median_S_score = if ("S.Score" %in% names(meta)) as.numeric(median(S.Score, na.rm = TRUE)) else NA_real_,
  median_G2M_score = if ("G2M.Score" %in% names(meta)) as.numeric(median(G2M.Score, na.rm = TRUE)) else NA_real_
), by = .(sample, condition, cell_group)]
qc <- merge(CJ(sample = samples, cell_group = groups$cell_group, unique = TRUE), qc,
            by = c("sample", "cell_group"), all.x = TRUE, sort = FALSE)
qc[is.na(condition), condition := unname(sample_condition[sample])]
qc <- qc[order(match(cell_group, groups$cell_group), match(sample, samples))]
fwrite(qc, outputs[1], na = "NA")

# Sample-level descriptive effect directions; no cell-level inferential tests are used.
effect_metrics <- c("n_cells", "median_nFeature_RNA", "median_nCount_RNA", "median_percent_mt")
effects <- rbindlist(lapply(groups$cell_group, function(g) {
  z <- qc[cell_group == g]
  rbindlist(lapply(effect_metrics, function(metric) {
    ctrl <- z[condition == "Control", get(metric)]
    sd <- z[condition == "SD", get(metric)]
    ctrl_med <- median(ctrl, na.rm = TRUE)
    sd_med <- median(sd, na.rm = TRUE)
    data.table(cell_group = g, metric = metric,
               Control_sample_median = ctrl_med, SD_sample_median = sd_med,
               SD_minus_Control = sd_med - ctrl_med,
               SD_to_Control_ratio = if (is.finite(ctrl_med) && ctrl_med != 0) sd_med / ctrl_med else NA_real_,
               direction = fifelse(sd_med > ctrl_med, "higher_in_SD", fifelse(sd_med < ctrl_med, "lower_in_SD", "no_change")))
  }))
}))
fwrite(effects, outputs[2], na = "NA")

qc_long <- melt(qc, id.vars = c("sample", "condition", "cell_group"),
                measure.vars = effect_metrics, variable.name = "metric", value.name = "value")
qc_long[, cell_group := factor(cell_group, levels = groups$cell_group)]
qc_long[, metric := factor(metric, levels = effect_metrics,
  labels = c("Cell count", "Median detected genes", "Median UMI count", "Median mitochondrial %"))]
p_qc <- ggplot(qc_long, aes(cell_group, value, color = condition, shape = condition)) +
  geom_point(position = position_dodge(width = 0.48), size = 2.2, alpha = 0.9) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(Control = "#0072B2", SD = "#D55E00")) +
  labs(x = NULL, y = NULL, color = "Condition", shape = "Condition",
       title = "Sample-level cell-group QC summaries",
       subtitle = "Each point is one biological sample; n=2 Control and n=2 SD") +
  theme_classic(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
ggsave(outputs[6], p_qc, width = 8.2, height = 6.2, dpi = 300, bg = "white")
ggsave(outputs[7], p_qc, width = 8.2, height = 6.2, device = "pdf", useDingbats = FALSE)

# Prepare immutable count and symbol inputs once, then equalize cell numbers independently 20 times.
counts <- LayerData(obj, assay = "RNA", layer = "counts")
symbol_map <- fread(mapping_path)
lookup <- setNames(symbol_map$current_gene_symbol, symbol_map$seurat_gene_symbol)
current_symbols <- unname(lookup[rownames(counts)])
current_symbols[is.na(current_symbols) | current_symbols == ""] <- rownames(counts)[is.na(current_symbols) | current_symbols == ""]
hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
setDT(hallmark)
set_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
oxphos <- unique(hallmark[gs_name == set_name, gene_symbol])
nuclear_genes <- fread(nuclear_path)[include_nuclear_only == TRUE, gene]

make_de <- function(pb) {
  y <- DGEList(counts = pb, samples = data.frame(sample = samples,
              condition = factor(unname(sample_condition[samples]), levels = c("Control", "SD"))))
  design <- model.matrix(~condition, data = y$samples)
  keep <- filterByExpr(y, design = design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "conditionSD")
  de <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "source_gene_symbol")
  de[, current_gene_symbol := current_symbols[match(source_gene_symbol, rownames(counts))]]
  de[is.na(current_gene_symbol) | current_gene_symbol == "", current_gene_symbol := source_gene_symbol]
  de[, signed_stat := sign(logFC) * sqrt(F)]
  de
}
make_rank <- function(de, method) {
  score <- if (method == "logFC") "logFC" else "signed_stat"
  z <- copy(de[is.finite(get(score)) & !is.na(current_gene_symbol) & current_gene_symbol != ""])
  setorder(z, current_gene_symbol, -logCPM)
  z <- z[!duplicated(current_gene_symbol)]
  ans <- z[[score]]
  names(ans) <- z$current_gene_symbol
  sort(ans, decreasing = TRUE)
}
run_one <- function(rank, genes, seed_value) {
  set.seed(seed_value)
  fg <- as.data.table(fgseaMultilevel(setNames(list(genes), set_name), rank,
                                     minSize = 10, maxSize = 500, eps = 0))
  if (nrow(fg) != 1L) stop("OXPHOS fgsea did not return one row")
  list(NES = fg$NES, pval = fg$pval, padj = fg$padj, size = fg$size,
       leadingEdge_count = length(fg$leadingEdge[[1]]))
}

iteration_rows <- vector("list", nrow(downsample_groups) * iterations * 2L)
k <- 0L
for (g_idx in seq_len(nrow(downsample_groups))) {
  g <- downsample_groups$cell_group[g_idx]
  group_cells <- meta[cell_group == g, .(cell, sample)]
  available <- group_cells[, .N, by = sample][match(samples, sample)]
  if (any(is.na(available$N))) stop("A sample is missing from cell group ", g)
  min_n <- min(available$N)
  for (iter in seq_len(iterations)) {
    sampling_seed <- seed + g_idx * 1000L + iter
    set.seed(sampling_seed)
    chosen <- unlist(lapply(samples, function(s) sample(group_cells[sample == s, cell], min_n, replace = FALSE)),
                     use.names = FALSE)
    chosen_sample <- rep(samples, each = min_n)
    pb <- do.call(cbind, lapply(samples, function(s) Matrix::rowSums(counts[, chosen[chosen_sample == s], drop = FALSE])))
    rownames(pb) <- rownames(counts)
    colnames(pb) <- samples
    de <- make_de(pb)
    for (method in c("logFC", "signed_stat")) {
      result <- run_one(make_rank(de, method), oxphos, sampling_seed + ifelse(method == "logFC", 10000L, 20000L))
      k <- k + 1L
      iteration_rows[[k]] <- data.table(cell_group = g, ranking_method = method,
        iteration = iter, seed = sampling_seed, min_n = min_n,
        NES = result$NES, pval = result$pval, padj = result$padj,
        size = result$size, leadingEdge_count = result$leadingEdge_count)
    }
  }
  message("Completed downsampling: ", g, " (min_n=", min_n, ")")
}
iterations_dt <- rbindlist(iteration_rows)
fwrite(iterations_dt, outputs[4])
summary_dt <- iterations_dt[, .(
  n_iterations = .N,
  min_n = unique(min_n),
  median_NES = median(NES),
  Q1_NES = quantile(NES, 0.25, names = FALSE),
  Q3_NES = quantile(NES, 0.75, names = FALSE),
  IQR_NES = IQR(NES),
  min_NES = min(NES),
  max_NES = max(NES),
  proportion_NES_negative = mean(NES < 0),
  proportion_FDR_lt_0.05 = mean(padj < 0.05)
), by = .(cell_group, ranking_method)]
summary_dt <- summary_dt[order(match(cell_group, downsample_groups$cell_group),
                               match(ranking_method, c("logFC", "signed_stat")))]
fwrite(summary_dt, outputs[3])

# Tuft sensitivity across full/nuclear-only sets and both ranking metrics, using all available cells.
edge_objects <- readRDS(edge_path)
tuft_de <- as.data.table(edge_objects[["Tuft"]]$de)
tuft_de[, signed_stat := sign(logFC) * sqrt(F)]
tuft_rows <- list()
k <- 0L
for (method in c("logFC", "signed_stat")) {
  method_seed <- seed + 50000L + match(method, c("logFC", "signed_stat"))
  rank <- make_rank(tuft_de, method)
  for (set_type in c("full", "nuclear_only")) {
    k <- k + 1L
    genes <- if (set_type == "full") oxphos else nuclear_genes
    z <- run_one(rank, genes, method_seed)
    tuft_rows[[k]] <- data.table(cell_group = "Tuft", ranking_method = method, set_type = set_type,
      NES = z$NES, pval = z$pval, padj = z$padj, size = z$size,
      leadingEdge_count = z$leadingEdge_count)
  }
}
fwrite(rbindlist(tuft_rows), outputs[5])

iterations_dt[, cell_group := factor(cell_group, levels = downsample_groups$cell_group)]
p_down <- ggplot(iterations_dt, aes(cell_group, NES, color = ranking_method, shape = ranking_method)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_boxplot(aes(group = interaction(cell_group, ranking_method)), position = position_dodge(width = 0.65),
               width = 0.52, outlier.shape = NA, alpha = 0.12, linewidth = 0.45) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65, seed = seed),
             size = 1.3, alpha = 0.65) +
  scale_color_manual(values = c(logFC = "#0072B2", signed_stat = "#D55E00"),
                     labels = c(logFC = "logFC", signed_stat = "signed statistic")) +
  scale_shape_manual(values = c(logFC = 16, signed_stat = 17),
                     labels = c(logFC = "logFC", signed_stat = "signed statistic")) +
  labs(x = NULL, y = "Hallmark OXPHOS NES", color = "Ranking", shape = "Ranking",
       title = "Equal-cell downsampling sensitivity",
       subtitle = "20 independent iterations per cell group; fixed master seed 20261003") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(outputs[8], p_down, width = 8.0, height = 4.8, dpi = 300, bg = "white")
ggsave(outputs[9], p_down, width = 8.0, height = 4.8, device = "pdf", useDingbats = FALSE)

message("QC and equal-cell downsampling completed: ", iterations, " iterations; seed=", seed)
