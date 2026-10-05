#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: three_dataset_validation.R PROJECT_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(file.path(root, "reference", "R_library"), .libPaths()))
out_dir <- file.path(root, "results", "PRJNA862187_validation")
gse_dir <- file.path(root, "results", "GSE289089_celltype_OXPHOS")
bulk_dir <- file.path(root, "results", "PRJNA1167170_bulk")

required <- c("data.table", "msigdbr", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(msigdbr); library(ggplot2)})

targets <- file.path(out_dir, c("three_dataset_OXPHOS_concordance.csv",
  "three_dataset_OXPHOS_concordance_summary.csv", "three_dataset_pathway_summary.csv",
  "three_dataset_OXPHOS_summary.pdf", "three_dataset_OXPHOS_summary.png"))
if (any(file.exists(targets))) stop("Refusing overwrite: ", paste(targets[file.exists(targets)], collapse = ", "))

hallmark <- as.data.table(msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH"))
oxphos_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
oxphos_genes <- sort(unique(hallmark[gs_name == oxphos_name, gene_symbol]))
if (!length(oxphos_genes)) stop("Native-mouse Hallmark OXPHOS unavailable")

deduplicate <- function(dt, symbol_col) {
  z <- copy(dt[!is.na(get(symbol_col)) & get(symbol_col) != ""])
  if ("logCPM" %in% names(z)) setorderv(z, c(symbol_col, "logCPM"), c(1L, -1L)) else setorderv(z, symbol_col)
  z[!duplicated(get(symbol_col))]
}
gse <- deduplicate(fread(file.path(gse_dir, "DE", "ISC_edgeR_SD_vs_Control.csv")), "current_gene_symbol")
bulk <- deduplicate(fread(file.path(bulk_dir, "edgeR_SD_vs_Control.csv")), "gene_symbol")
third <- deduplicate(fread(file.path(out_dir, "edgeR_primary_SD_effect.csv")), "gene_symbol")

concordance <- data.table(gene = oxphos_genes)
concordance <- merge(concordance, gse[, .(gene = current_gene_symbol, GSE289089_logFC = logFC)], by = "gene", all.x = TRUE, sort = FALSE)
concordance <- merge(concordance, bulk[, .(gene = gene_symbol, PRJNA1167170_logFC = logFC)], by = "gene", all.x = TRUE, sort = FALSE)
concordance <- merge(concordance, third[, .(gene = gene_symbol, PRJNA862187_primary_SD_logFC = logFC)], by = "gene", all.x = TRUE, sort = FALSE)
concordance[, measured_in_all_three := complete.cases(GSE289089_logFC, PRJNA1167170_logFC, PRJNA862187_primary_SD_logFC)]
concordance[, negative_in_all_three := measured_in_all_three & GSE289089_logFC < 0 & PRJNA1167170_logFC < 0 & PRJNA862187_primary_SD_logFC < 0]
fwrite(concordance, file.path(out_dir, "three_dataset_OXPHOS_concordance.csv"))

pairs <- list(
  c("GSE289089_logFC", "PRJNA1167170_logFC"),
  c("GSE289089_logFC", "PRJNA862187_primary_SD_logFC"),
  c("PRJNA1167170_logFC", "PRJNA862187_primary_SD_logFC")
)
pair_summary <- rbindlist(lapply(pairs, function(cols) {
  ok <- complete.cases(concordance[, ..cols])
  test <- suppressWarnings(cor.test(concordance[[cols[1]]][ok], concordance[[cols[2]]][ok], method = "spearman", exact = FALSE))
  data.table(metric = "pairwise_spearman", dataset_1 = cols[1], dataset_2 = cols[2],
             n_genes = sum(ok), value = unname(test$estimate), p_value = test$p.value)
}))
n_complete <- sum(concordance$measured_in_all_three)
n_all_negative <- sum(concordance$negative_in_all_three)
overall <- data.table(
  metric = c("hallmark_gene_set_size", "measured_in_all_three", "negative_in_all_three",
             "negative_in_all_three_percent_of_complete", "negative_in_all_three_percent_of_hallmark"),
  dataset_1 = NA_character_, dataset_2 = NA_character_, n_genes = NA_integer_,
  value = c(length(oxphos_genes), n_complete, n_all_negative,
            100 * n_all_negative / n_complete, 100 * n_all_negative / length(oxphos_genes)),
  p_value = NA_real_
)
fwrite(rbind(pair_summary, overall, fill = TRUE), file.path(out_dir, "three_dataset_OXPHOS_concordance_summary.csv"))

gse_path <- fread(file.path(gse_dir, "celltype_mitochondrial_pathways.csv"))[
  cell_group == "ISC" & ranking_method == "logFC"]
gse_path <- gse_path[, .(dataset = "GSE289089", tissue = "small-intestinal ISC pseudobulk",
  design = "2 control vs 2 SD", n_control = 2L, n_SD = 2L, contrast = "SD-Control",
  pathway, OXPHOS_NES = NES, OXPHOS_FDR = padj)]

bulk_primary <- fread(file.path(bulk_dir, "primary_OXPHOS_validation.csv"))
bulk_secondary <- fread(file.path(bulk_dir, "secondary_mitochondrial_pathways.csv"))
bulk_path <- rbindlist(list(
  bulk_primary[, .(pathway, NES, padj)],
  bulk_secondary[, .(pathway, NES, padj)]
), use.names = TRUE)
bulk_path <- bulk_path[, .(dataset = "PRJNA1167170", tissue = "proximal/whole colon",
  design = "4 control vs 4 SD", n_control = 4L, n_SD = 4L, contrast = "SD-Control",
  pathway, OXPHOS_NES = NES, OXPHOS_FDR = padj)]

third_primary <- fread(file.path(out_dir, "primary_OXPHOS_validation.csv"))[, .(pathway, NES, padj)]
third_secondary <- fread(file.path(out_dir, "secondary_mitochondrial_pathways.csv"))[
  contrast == "average_SD_across_diets", .(pathway, NES, padj)]
third_path <- rbindlist(list(third_primary, third_secondary), use.names = TRUE)
third_path <- third_path[, .(dataset = "PRJNA862187", tissue = "whole colon / large intestine",
  design = "2x2 factorial; n=3 per cell", n_control = 6L, n_SD = 6L,
  contrast = "average SD effect across SCD and HFD", pathway,
  OXPHOS_NES = NES, OXPHOS_FDR = padj)]

path_summary <- rbindlist(list(gse_path, bulk_path, third_path), use.names = TRUE)
path_order <- c(oxphos_name,
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT", "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION")
if (!all(table(path_summary$dataset) == 6L) || !setequal(unique(path_summary$pathway), path_order)) {
  stop("Three-dataset pathway summary is incomplete")
}
fwrite(path_summary, file.path(out_dir, "three_dataset_pathway_summary.csv"))

plot_dt <- copy(path_summary)
plot_dt[, dataset_label := factor(dataset, levels = c("GSE289089", "PRJNA1167170", "PRJNA862187"),
  labels = c("GSE289089\nISC pseudobulk", "PRJNA1167170\ncolon bulk", "PRJNA862187\ncolon bulk, factorial"))]
short_names <- c(
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Hallmark OXPHOS",
  REACTOME_RESPIRATORY_ELECTRON_TRANSPORT = "Respiratory electron transport",
  REACTOME_MITOCHONDRIAL_TRANSLATION = "Mitochondrial translation",
  REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING = "ATP formation by chemiosmosis",
  REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE = "TCA cycle",
  REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION = "Mitochondrial FA beta-oxidation"
)
plot_dt[, pathway_label := factor(unname(short_names[pathway]), levels = rev(unname(short_names[path_order])))]
limit <- max(abs(plot_dt$OXPHOS_NES), na.rm = TRUE)
p <- ggplot(plot_dt, aes(dataset_label, pathway_label, fill = OXPHOS_NES)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f%s", OXPHOS_NES, ifelse(OXPHOS_FDR < 0.05, "*", ""))), size = 3.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-limit, limit), name = "NES") +
  labs(x = NULL, y = NULL, title = "Prespecified mitochondrial pathway validation across datasets",
       subtitle = "* FDR < 0.05; whole-colon datasets and ISC pseudobulk have different biological resolution") +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "right")
ggsave(file.path(out_dir, "three_dataset_OXPHOS_summary.pdf"), p, width = 8.2, height = 5.8)
ggsave(file.path(out_dir, "three_dataset_OXPHOS_summary.png"), p, width = 8.2, height = 5.8, dpi = 300)

message("Three-dataset OXPHOS validation outputs completed")
