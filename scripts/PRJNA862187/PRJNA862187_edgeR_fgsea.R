#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: PRJNA862187_edgeR_fgsea.R PROJECT_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- file.path(root, "results", "PRJNA862187_validation")
salmon_dir <- file.path(out_dir, "salmon")
metadata_file <- file.path(out_dir, "PRJNA862187_sample_metadata.csv")
gtf <- file.path(root, "reference", "gencode_M39", "gencode.vM39.chr_patch_hapl_scaff.annotation.gtf.gz")
.libPaths(c(file.path(root, "reference", "R_library"), .libPaths()))

required <- c("data.table", "tximport", "edgeR", "fgsea", "msigdbr", "ggplot2", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(tximport); library(edgeR); library(fgsea)
  library(msigdbr); library(ggplot2)
})

targets <- file.path(out_dir, c(
  "gene_counts.csv", "tximport_object.rds", "salmon_qc_summary.csv",
  "edgeR_primary_SD_effect.csv", "edgeR_SD_SCD.csv", "edgeR_SD_HFD.csv",
  "edgeR_sleep_diet_interaction.csv", "primary_OXPHOS_validation.csv",
  "secondary_mitochondrial_pathways.csv", "diet_specific_OXPHOS.csv",
  "interaction_pathways.csv", "PRJNA862187_PCA.pdf", "PRJNA862187_PCA.png",
  "PRJNA862187_primary_OXPHOS.pdf", "PRJNA862187_primary_OXPHOS.png",
  "PRJNA862187_diet_specific_OXPHOS.pdf", "PRJNA862187_diet_specific_OXPHOS.png"
))
if (any(file.exists(targets))) stop("Refusing to overwrite existing analysis output(s): ", paste(targets[file.exists(targets)], collapse = ", "))

meta <- fread(metadata_file)
expected_groups <- c("EC_SCD", "SD_SCD", "EC_HFD", "SD_HFD")
if (nrow(meta) != 12L || !all(table(meta$group) == 3L) || !setequal(unique(meta$group), expected_groups)) {
  stop("Metadata does not match frozen 4 groups x 3 biological replicates design")
}
meta[, group := factor(group, levels = expected_groups)]
meta[, group_order := match(group, expected_groups)]
setorder(meta, group_order, biological_replicate)
meta[, group_order := NULL]

quant_files <- file.path(salmon_dir, meta$sample, "quant.sf")
names(quant_files) <- meta$sample
if (!all(file.exists(quant_files))) stop("Missing quant.sf: ", paste(quant_files[!file.exists(quant_files)], collapse = ", "))
if (!file.exists(gtf)) stop("Missing GENCODE M39 GTF: ", gtf)

gtf_cmd <- sprintf("gzip -dc %s | awk -F '\\t' '$3 == \"transcript\"'", shQuote(gtf))
gtf_dt <- fread(cmd = gtf_cmd, sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE)
tx2gene <- unique(data.table(
  TXNAME = sub('.*transcript_id "([^"]+)".*', "\\1", gtf_dt$V9),
  GENEID = sub('.*gene_id "([^"]+)".*', "\\1", gtf_dt$V9),
  gene_symbol = sub('.*gene_name "([^"]+)".*', "\\1", gtf_dt$V9)
))
gene_annot <- unique(tx2gene[, .(gene_id = GENEID, gene_symbol)])

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene[, .(TXNAME, GENEID)], countsFromAbundance = "no")
saveRDS(txi, file.path(out_dir, "tximport_object.rds"))
counts_out <- data.table(gene_id = rownames(txi$counts), as.data.frame(txi$counts, check.names = FALSE))
counts_out <- merge(gene_annot, counts_out, by = "gene_id", all.y = TRUE, sort = FALSE)
setcolorder(counts_out, c("gene_id", "gene_symbol", meta$sample))
fwrite(counts_out, file.path(out_dir, "gene_counts.csv"))

salmon_qc <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
  sample <- meta$sample[i]
  info <- jsonlite::fromJSON(file.path(salmon_dir, sample, "aux_info", "meta_info.json"))
  fmt <- jsonlite::fromJSON(file.path(salmon_dir, sample, "lib_format_counts.json"))
  data.table(
    sample = sample, run = meta$SRR[i], group = as.character(meta$group[i]),
    diet = meta$diet[i], sleep = meta$sleep_condition[i],
    fragments = info$num_processed, mapped_fragments = info$num_mapped,
    mapping_rate = info$percent_mapped, library_type = info$detected_library_type,
    orphan_fragments = info$num_orphan,
    orphan_rate = 100 * info$num_orphan / info$num_processed,
    compatible_ratio = fmt$compatible_fragment_ratio,
    assigned_fragments = fmt$num_assigned_fragments,
    fragment_length_mean = info$frag_length_mean,
    fragment_length_sd = info$frag_length_sd
  )
}))
robust_low <- function(x) median(x) - 3 * mad(x, constant = 1.4826)
salmon_qc[, outlier_flag := mapping_rate < robust_low(mapping_rate) |
            compatible_ratio < robust_low(compatible_ratio) |
            fragments < robust_low(fragments)]
salmon_qc[, outlier_reason := fifelse(outlier_flag,
  paste0("technical metric below median-3*MAD; retained in analysis"), "none")]
fwrite(salmon_qc, file.path(out_dir, "salmon_qc_summary.csv"))

y <- DGEList(counts = txi$counts, samples = as.data.frame(meta))
design <- model.matrix(~0 + group, data = meta)
colnames(design) <- expected_groups
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

contrasts <- list(
  primary = makeContrasts(0.5 * SD_SCD + 0.5 * SD_HFD - 0.5 * EC_SCD - 0.5 * EC_HFD, levels = design),
  SD_SCD = makeContrasts(SD_SCD - EC_SCD, levels = design),
  SD_HFD = makeContrasts(SD_HFD - EC_HFD, levels = design),
  interaction = makeContrasts((SD_HFD - EC_HFD) - (SD_SCD - EC_SCD), levels = design)
)
output_names <- c(primary = "edgeR_primary_SD_effect.csv", SD_SCD = "edgeR_SD_SCD.csv",
                  SD_HFD = "edgeR_SD_HFD.csv", interaction = "edgeR_sleep_diet_interaction.csv")
de_results <- list()
for (nm in names(contrasts)) {
  test <- glmQLFTest(fit, contrast = contrasts[[nm]])
  dt <- as.data.table(topTags(test, n = Inf, sort.by = "none")$table, keep.rownames = "gene_id")
  dt <- merge(gene_annot, dt, by = "gene_id", all.y = TRUE, sort = FALSE)
  setorder(dt, PValue)
  fwrite(dt[, .(gene_id, gene_symbol, logFC, logCPM, F, PValue, FDR)], file.path(out_dir, output_names[[nm]]))
  de_results[[nm]] <- dt
}
saveRDS(list(y = y, design = design, fit = fit, contrasts = contrasts), file.path(out_dir, "edgeR_factorial_model.rds"))

logcpm <- cpm(y, log = TRUE, prior.count = 2)
pca <- prcomp(t(logcpm), scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_dt <- data.table(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_dt <- merge(pca_dt, meta[, .(sample, group, diet, sleep_condition)], by = "sample", sort = FALSE)
fwrite(pca_dt, file.path(out_dir, "PRJNA862187_PCA_coordinates.csv"))
palette <- c(EC_SCD = "#0072B2", SD_SCD = "#D55E00", EC_HFD = "#009E73", SD_HFD = "#CC79A7")
p_pca <- ggplot(pca_dt, aes(PC1, PC2, color = group, shape = group, label = sample)) +
  geom_point(size = 3.2) + geom_text(nudge_y = 0.25, size = 3, show.legend = FALSE) +
  scale_color_manual(values = palette) + scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]),
       title = "PRJNA862187: normalized gene-expression PCA") + theme_classic(base_size = 11)
ggsave(file.path(out_dir, "PRJNA862187_PCA.pdf"), p_pca, width = 7, height = 5.8)
ggsave(file.path(out_dir, "PRJNA862187_PCA.png"), p_pca, width = 7, height = 5.8, dpi = 300)

hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
reactome <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2", subcollection = "CP:REACTOME")
msigdb_version <- unique(hallmark$db_version)
if (length(msigdb_version) != 1L) stop("Could not resolve one MSigDB release")
primary_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
secondary_names <- c(
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION"
)
make_sets <- function(tbl, wanted) {
  present <- intersect(wanted, unique(tbl$gs_name))
  split(tbl$gene_symbol[tbl$gs_name %in% present], tbl$gs_name[tbl$gs_name %in% present])
}
primary_set <- make_sets(hallmark, primary_name)
secondary_sets <- make_sets(reactome, secondary_names)
if (length(primary_set) != 1L || !setequal(names(secondary_sets), secondary_names)) stop("A frozen pathway is unavailable in ", msigdb_version)

make_rank <- function(de) {
  z <- copy(de[!is.na(gene_symbol) & gene_symbol != "" & is.finite(logFC)])
  setorder(z, gene_symbol, -logCPM)
  z <- z[!duplicated(gene_symbol)]
  stats <- z$logFC; names(stats) <- z$gene_symbol
  sort(stats, decreasing = TRUE)
}
tidy_fgsea <- function(pathways, stats, contrast_name, analysis_class, seed) {
  set.seed(seed)
  ans <- as.data.table(fgseaMultilevel(pathways, stats, minSize = 10, maxSize = 500, eps = 0))
  ans[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  ans[, `:=`(contrast = contrast_name, analysis_class = analysis_class, msigdb_version = msigdb_version)]
  setcolorder(ans, c("analysis_class", "contrast", "pathway", "NES", "pval", "padj", "size", "leadingEdge", "msigdb_version"))
  ans
}
ranks <- lapply(de_results, make_rank)

primary_out <- tidy_fgsea(primary_set, ranks$primary, "average_SD_across_diets", "PRIMARY", 862187L)
primary_out[, interpretation := fifelse(NES < 0 & padj < 0.05, "supports",
  fifelse(NES < 0, "directionally supports", "does not support"))]
fwrite(primary_out, file.path(out_dir, "primary_OXPHOS_validation.csv"))

diet_out <- rbindlist(list(
  tidy_fgsea(primary_set, ranks$SD_SCD, "SD_SCD_minus_EC_SCD", "SECONDARY_DIET_SPECIFIC", 862188L),
  tidy_fgsea(primary_set, ranks$SD_HFD, "SD_HFD_minus_EC_HFD", "SECONDARY_DIET_SPECIFIC", 862189L)
))
fwrite(diet_out, file.path(out_dir, "diet_specific_OXPHOS.csv"))

secondary_out <- rbindlist(list(
  tidy_fgsea(secondary_sets, ranks$primary, "average_SD_across_diets", "SECONDARY", 862190L),
  tidy_fgsea(secondary_sets, ranks$SD_SCD, "SD_SCD_minus_EC_SCD", "SECONDARY", 862191L),
  tidy_fgsea(secondary_sets, ranks$SD_HFD, "SD_HFD_minus_EC_HFD", "SECONDARY", 862192L)
))
fwrite(secondary_out, file.path(out_dir, "secondary_mitochondrial_pathways.csv"))

interaction_sets <- c(primary_set, secondary_sets)
interaction_out <- tidy_fgsea(interaction_sets, ranks$interaction, "sleep_by_diet_interaction", "EXPLORATORY_INTERACTION", 862193L)
fwrite(interaction_out, file.path(out_dir, "interaction_pathways.csv"))

p_primary <- plotEnrichment(primary_set[[primary_name]], ranks$primary) +
  labs(title = "Primary external validation: Hallmark OXPHOS",
       subtitle = sprintf("Average SD effect across diets; NES %.3f, p %.3g, FDR %.3g", primary_out$NES, primary_out$pval, primary_out$padj),
       x = "Genes ranked by edgeR logFC", y = "Running enrichment score") + theme_classic(base_size = 11)
ggsave(file.path(out_dir, "PRJNA862187_primary_OXPHOS.pdf"), p_primary, width = 7, height = 5.5)
ggsave(file.path(out_dir, "PRJNA862187_primary_OXPHOS.png"), p_primary, width = 7, height = 5.5, dpi = 300)

p_diet <- ggplot(diet_out, aes(contrast, NES, color = NES < 0, shape = padj < 0.05)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.5) + geom_point(size = 4) +
  geom_text(aes(label = sprintf("FDR %.3g", padj)), nudge_y = 0.12, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00"), labels = c(`TRUE` = "NES < 0", `FALSE` = "NES >= 0")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05")) +
  labs(x = NULL, y = "Normalized enrichment score", color = "Direction", shape = "Evidence",
       title = "Diet-specific Hallmark OXPHOS validation") + theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(out_dir, "PRJNA862187_diet_specific_OXPHOS.pdf"), p_diet, width = 7, height = 5.5)
ggsave(file.path(out_dir, "PRJNA862187_diet_specific_OXPHOS.png"), p_diet, width = 7, height = 5.5, dpi = 300)

sink(file.path(out_dir, "reproducibility_info.txt"))
cat("Analysis: PRJNA862187 prespecified 2x2 factorial external validation\n")
cat("Frozen plan:", file.path(out_dir, "PRJNA862187_analysis_plan_frozen.md"), "\n")
cat("Reference: GRCm39 / GENCODE Mouse M39; Salmon index", file.path(root, "reference", "salmon_mouse_GRCm39"), "\n")
cat("Salmon:", system2("salmon", "--version", stdout = TRUE), "\n")
cat("tximport countsFromAbundance: no\n")
cat("edgeR design: ~0 + group; all 12 biological samples retained\n")
cat("Primary contrast: 0.5*SD_SCD + 0.5*SD_HFD - 0.5*EC_SCD - 0.5*EC_HFD\n")
cat("Secondary contrasts: SD_SCD-EC_SCD; SD_HFD-EC_HFD\n")
cat("Exploratory interaction: (SD_HFD-EC_HFD)-(SD_SCD-EC_SCD)\n")
cat("GSEA ranking: edgeR logFC, complete tested gene list; fgseaMultilevel minSize=10 maxSize=500 eps=0\n")
cat("MSigDB:", msigdb_version, "native mouse\n")
cat("Seeds: 862187-862193 as recorded in script\n")
cat("Technical outliers flagged, not removed:", sum(salmon_qc$outlier_flag), "\n")
cat("Package versions:\n")
for (pkg in required) cat("  ", pkg, ": ", as.character(packageVersion(pkg)), "\n", sep = "")
cat("\nSession information:\n"); print(sessionInfo())
sink()

message("PRJNA862187 edgeR and frozen pathway analysis completed")
