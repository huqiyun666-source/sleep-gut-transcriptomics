#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: final_OXPHOS_concordance.R PROJECT_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(file.path(root, "reference", "R_library"), .libPaths()))

required <- c("edgeR", "fgsea", "msigdbr", "data.table", "ggplot2", "ggrepel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

out_dir <- file.path(root, "results", "final_mito_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
outputs <- file.path(out_dir, c(
  "OXPHOS_gene_level_concordance.csv", "OXPHOS_concordance_statistics.csv",
  "nuclear_only_OXPHOS_genes.csv", "nuclear_only_OXPHOS_validation.csv",
  "OXPHOS_cross_dataset_logFC_scatter.png", "OXPHOS_cross_dataset_logFC_scatter.pdf",
  "nuclear_only_OXPHOS_comparison.png", "nuclear_only_OXPHOS_comparison.pdf"
))
if (any(file.exists(outputs))) stop("Refusing to overwrite existing output(s): ", paste(outputs[file.exists(outputs)], collapse = ", "))

gse_rds <- file.path(root, "results", "GSE289089_celltype_OXPHOS", "celltype_edgeR_objects.rds")
bulk_txi <- file.path(root, "results", "PRJNA1167170_bulk", "tximport_object.rds")
bulk_csv <- file.path(root, "results", "PRJNA1167170_bulk", "edgeR_SD_vs_Control.csv")
gtf_path <- file.path(root, "reference", "gencode_M39", "gencode.vM39.chr_patch_hapl_scaff.annotation.gtf.gz")
inputs <- c(gse_rds, bulk_txi, bulk_csv, gtf_path)
if (!all(file.exists(inputs))) stop("Missing input(s): ", paste(inputs[!file.exists(inputs)], collapse = ", "))

hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
setDT(hallmark)
msigdb_version <- unique(hallmark$db_version)
if (length(msigdb_version) != 1L) stop("Could not determine one MSigDB release")
set_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
oxphos <- unique(hallmark[gs_name == set_name, gene_symbol])
if (length(oxphos) < 10L) stop("Native-mouse Hallmark OXPHOS set is unavailable")

make_rank <- function(de, method, symbol_col) {
  z <- copy(de[!is.na(get(symbol_col)) & get(symbol_col) != ""])
  score_col <- if (method == "logFC") "logFC" else "signed_stat"
  z <- z[is.finite(get(score_col))]
  setorderv(z, c(symbol_col, "logCPM"), c(1L, -1L))
  z <- z[!duplicated(get(symbol_col))]
  ans <- z[[score_col]]
  names(ans) <- z[[symbol_col]]
  sort(ans, decreasing = TRUE)
}

# Reuse the fixed GSE289089 ISC edgeR result.
gse_objects <- readRDS(gse_rds)
gse_de <- as.data.table(gse_objects[["ISC"]]$de)
if (!all(c("current_gene_symbol", "logFC", "F") %in% names(gse_de))) stop("ISC edgeR object is incomplete")
gse_de[, signed_stat := sign(logFC) * sqrt(F)]

# Refit the stored PRJNA1167170 tximport counts to recover the original QL F statistic.
txi <- readRDS(bulk_txi)
bulk_samples <- data.frame(
  sample = c("con1", "con2", "con3", "con4", "SD1", "SD2", "SD3", "SD4"),
  condition = factor(c(rep("Control", 4), rep("SD", 4)), levels = c("Control", "SD"))
)
rownames(bulk_samples) <- bulk_samples$sample
if (!identical(colnames(txi$counts), bulk_samples$sample)) stop("Unexpected PRJNA1167170 sample order")
y <- DGEList(counts = txi$counts, samples = bulk_samples)
design <- model.matrix(~condition, data = bulk_samples)
keep <- filterByExpr(y, design = design)
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "conditionSD")
bulk_de <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene_id")
bulk_annot <- fread(bulk_csv)[, .(gene_id, gene_symbol)]
bulk_de <- merge(bulk_annot, bulk_de, by = "gene_id", all.y = TRUE, sort = FALSE)
bulk_de[, signed_stat := sign(logFC) * sqrt(F)]

# Complete Hallmark gene table; concordance statistics use only genes observed in both datasets.
gse_one <- copy(gse_de)
setorder(gse_one, current_gene_symbol, -logCPM)
gse_one <- gse_one[!duplicated(current_gene_symbol), .(gene = current_gene_symbol, GSE289089_ISC_logFC = logFC)]
bulk_one <- copy(bulk_de)
setorder(bulk_one, gene_symbol, -logCPM)
bulk_one <- bulk_one[!duplicated(gene_symbol), .(gene = gene_symbol, PRJNA1167170_logFC = logFC)]
concordance <- merge(data.table(gene = sort(oxphos)), gse_one, by = "gene", all.x = TRUE, sort = FALSE)
concordance <- merge(concordance, bulk_one, by = "gene", all.x = TRUE, sort = FALSE)
concordance[, `:=`(
  present_in_GSE289089 = !is.na(GSE289089_ISC_logFC),
  present_in_PRJNA1167170 = !is.na(PRJNA1167170_logFC)
)]
setorder(concordance, gene)
fwrite(concordance, outputs[1])

both <- concordance[present_in_GSE289089 & present_in_PRJNA1167170]
if (nrow(both) < 3L) stop("Too few shared OXPHOS genes for concordance analysis")
cor_rows <- rbindlist(lapply(c("spearman", "kendall", "pearson"), function(method) {
  z <- suppressWarnings(cor.test(both$GSE289089_ISC_logFC, both$PRJNA1167170_logFC,
                                 method = method, exact = FALSE))
  data.table(metric = paste0(method, "_correlation"), value = unname(z$estimate), p_value = z$p.value,
             numerator = NA_integer_, denominator = nrow(both))
}))
direction <- data.table(
  metric = c("double_negative_proportion", "double_positive_proportion", "opposite_direction_proportion", "zero_in_either_proportion"),
  numerator = c(
    sum(both$GSE289089_ISC_logFC < 0 & both$PRJNA1167170_logFC < 0),
    sum(both$GSE289089_ISC_logFC > 0 & both$PRJNA1167170_logFC > 0),
    sum(sign(both$GSE289089_ISC_logFC) * sign(both$PRJNA1167170_logFC) < 0),
    sum(both$GSE289089_ISC_logFC == 0 | both$PRJNA1167170_logFC == 0)
  ),
  denominator = nrow(both), p_value = NA_real_
)
direction[, value := numerator / denominator]
stats <- rbind(cor_rows, direction, fill = TRUE)
stats[, n_hallmark_total := length(oxphos)]
stats[, n_present_both := nrow(both)]
fwrite(stats, outputs[2])

# Labels are selected deterministically from five prespecified respiratory families.
label_pool <- both[grepl("^(Nduf|Uqcr|Cox|Atp5|Sdh)", gene)]
label_pool[, family := fifelse(grepl("^Nduf", gene), "NDUF",
                        fifelse(grepl("^Uqcr", gene), "UQCR",
                        fifelse(grepl("^Cox", gene), "COX",
                        fifelse(grepl("^Atp5", gene), "ATP5", "SDH"))))]
label_pool[, combined_magnitude := abs(scale(GSE289089_ISC_logFC)) + abs(scale(PRJNA1167170_logFC))]
setorder(label_pool, family, -combined_magnitude, gene)
labels <- label_pool[, head(.SD, 3L), by = family]

p_scatter <- ggplot(both, aes(GSE289089_ISC_logFC, PRJNA1167170_logFC)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_point(color = "#0072B2", alpha = 0.72, size = 1.8) +
  geom_text_repel(data = labels, aes(label = gene), size = 2.8, max.overlaps = Inf,
                  box.padding = 0.25, point.padding = 0.15, seed = 20261003) +
  labs(x = "GSE289089 ISC log2 fold-change (SD vs Control)",
       y = "PRJNA1167170 colon log2 fold-change (SD vs Control)",
       title = "Hallmark OXPHOS gene-level concordance",
       subtitle = sprintf("Shared genes n=%d; Spearman rho=%.3f", nrow(both), stats[metric == "spearman_correlation", value])) +
  theme_classic(base_size = 10)
ggsave(outputs[5], p_scatter, width = 6.4, height = 5.4, dpi = 300, bg = "white")
ggsave(outputs[6], p_scatter, width = 6.4, height = 5.4, device = "pdf", useDingbats = FALSE)

# Classify mitochondrial encoding from GENCODE M39 chromosome annotation, not symbol prefixes.
gtf_cmd <- sprintf("gzip -dc %s | awk -F '\\t' '$3 == \"gene\"'", shQuote(gtf_path))
gtf <- fread(cmd = gtf_cmd, sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE)
gtf_genes <- unique(data.table(
  gene = sub('.*gene_name "([^"]+)".*', "\\1", gtf$V9),
  chromosome = gtf$V1
))
annot <- gtf_genes[gene %in% oxphos, .(chromosomes = paste(sort(unique(chromosome)), collapse = ";"),
                                       is_mtDNA = any(chromosome == "chrM")), by = gene]
nuclear <- merge(data.table(gene = sort(oxphos)), annot, by = "gene", all.x = TRUE, sort = FALSE)
nuclear[, annotation_status := ifelse(is.na(chromosomes), "not_found_in_GENCODE_M39", "annotated")]
nuclear[is.na(is_mtDNA), is_mtDNA := NA]
nuclear[, encoding := fifelse(is.na(is_mtDNA), "unknown", fifelse(is_mtDNA, "mitochondrial", "nuclear"))]
nuclear[, include_nuclear_only := encoding == "nuclear"]
setorder(nuclear, gene)
fwrite(nuclear, outputs[3])
nuclear_genes <- nuclear[include_nuclear_only == TRUE, gene]
if (length(nuclear_genes) < 10L) stop("Nuclear-only OXPHOS set is unexpectedly small")

run_fgsea <- function(stats_vector, genes, dataset, method, set_type, seed) {
  set.seed(seed)
  fg <- as.data.table(fgseaMultilevel(setNames(list(genes), paste0(set_name, "_", set_type)),
                                     stats_vector, minSize = 10, maxSize = 500, eps = 0))
  if (nrow(fg) != 1L) stop("fgsea failed for ", dataset, " / ", method, " / ", set_type)
  data.table(dataset = dataset, ranking_method = method, set_type = set_type,
             NES = fg$NES, pval = fg$pval, padj = fg$padj, size = fg$size,
             leadingEdge_count = length(fg$leadingEdge[[1]]),
             leadingEdge = paste(fg$leadingEdge[[1]], collapse = ";"))
}
fg_rows <- list()
idx <- 0L
comparison_idx <- 0L
for (dataset in c("GSE289089_ISC", "PRJNA1167170_whole_colon")) {
  de <- if (dataset == "GSE289089_ISC") gse_de else bulk_de
  sym <- if (dataset == "GSE289089_ISC") "current_gene_symbol" else "gene_symbol"
  for (method in c("logFC", "signed_stat")) {
    comparison_idx <- comparison_idx + 1L
    rank <- make_rank(de, method, sym)
    for (set_type in c("full", "nuclear_only")) {
      idx <- idx + 1L
      genes <- if (set_type == "full") oxphos else nuclear_genes
      # The same seed within each dataset/ranking comparison prevents stochastic
      # differences when full and nuclear-only happen to contain identical genes.
      fg_rows[[idx]] <- run_fgsea(rank, genes, dataset, method, set_type, 20261003L + comparison_idx)
    }
  }
}
validation <- rbindlist(fg_rows)
fwrite(validation, outputs[4])

validation[, dataset_label := fifelse(dataset == "GSE289089_ISC", "GSE289089 ISC", "PRJNA1167170 colon")]
validation[, set_label := fifelse(set_type == "full", "Full OXPHOS", "Nuclear-only OXPHOS")]
p_nuclear <- ggplot(validation, aes(set_label, NES, color = ranking_method, shape = ranking_method, group = ranking_method)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_line(position = position_dodge(width = 0.25), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.25), size = 2.8) +
  facet_wrap(~dataset_label) +
  scale_color_manual(values = c(logFC = "#0072B2", signed_stat = "#D55E00"),
                     labels = c(logFC = "logFC", signed_stat = "signed statistic")) +
  scale_shape_manual(values = c(logFC = 16, signed_stat = 17),
                     labels = c(logFC = "logFC", signed_stat = "signed statistic")) +
  labs(x = NULL, y = "Normalized enrichment score (NES)", color = "Ranking", shape = "Ranking",
       title = "Full versus nuclear-only Hallmark OXPHOS",
       subtitle = "GENCODE M39: no mtDNA-encoded genes were present in the native-mouse Hallmark set") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "bottom")
ggsave(outputs[7], p_nuclear, width = 7.2, height = 4.4, dpi = 300, bg = "white")
ggsave(outputs[8], p_nuclear, width = 7.2, height = 4.4, device = "pdf", useDingbats = FALSE)

message("Concordance and nuclear-only analyses completed. MSigDB ", msigdb_version,
        "; Hallmark genes=", length(oxphos), "; nuclear-only=", length(nuclear_genes))
