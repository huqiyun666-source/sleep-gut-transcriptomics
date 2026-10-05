#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: PRJNA1167170_tximport_edgeR_fgsea.R PROJECT_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- file.path(root, "results", "PRJNA1167170_bulk")
salmon_dir <- file.path(out_dir, "salmon")
gtf <- file.path(root, "reference", "gencode_M39", "gencode.vM39.chr_patch_hapl_scaff.annotation.gtf.gz")
.libPaths(c(file.path(root, "reference", "R_library"), .libPaths()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

protected_outputs <- file.path(out_dir, c(
  "gene_counts.csv", "tximport_object.rds",
  "salmon_qc_summary.csv", "edgeR_SD_vs_Control.csv",
  "primary_OXPHOS_validation.csv", "reproducibility_info.txt"
))
if (any(file.exists(protected_outputs))) {
  stop("Existing analysis output(s) would be overwritten: ",
       paste(protected_outputs[file.exists(protected_outputs)], collapse = ", "))
}

required <- c("data.table", "tximport", "edgeR", "fgsea", "msigdbr", "ggplot2", "pheatmap", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(tximport)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(pheatmap)
})

samples <- data.frame(
  sample = c("con1", "con2", "con3", "con4", "SD1", "SD2", "SD3", "SD4"),
  run = c("SRR30835575", "SRR30835574", "SRR30835573", "SRR30835572",
          "SRR30835571", "SRR30835570", "SRR30835569", "SRR30835568"),
  condition = factor(c(rep("Control", 4), rep("SD", 4)), levels = c("Control", "SD")),
  stringsAsFactors = FALSE
)
metadata_path <- file.path(out_dir, "sample_metadata.csv")
samples_out <- transform(samples, condition = as.character(condition))
if (file.exists(metadata_path)) {
  existing_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
  if (!identical(existing_metadata, samples_out)) {
    stop("Existing sample metadata differs from the prespecified design: ", metadata_path)
  }
} else {
  write.csv(samples_out, metadata_path, row.names = FALSE)
}

quant_files <- file.path(salmon_dir, samples$sample, "quant.sf")
names(quant_files) <- samples$sample
if (!all(file.exists(quant_files))) stop("Missing quant.sf: ", paste(quant_files[!file.exists(quant_files)], collapse = ", "))
if (!file.exists(gtf)) stop("GENCODE M39 GTF is missing: ", gtf)

# Parse versioned GENCODE transcript/gene identifiers and gene symbols from the
# same M39 release used to build the Salmon index.
gtf_cmd <- sprintf("gzip -dc %s | awk -F '\\t' '$3 == \"transcript\"'", shQuote(gtf))
gtf_dt <- fread(cmd = gtf_cmd, sep = "\t", header = FALSE,
                quote = "", fill = TRUE, showProgress = FALSE)
tx2gene <- unique(data.table(
  TXNAME = sub('.*transcript_id "([^"]+)".*', "\\1", gtf_dt$V9),
  GENEID = sub('.*gene_id "([^"]+)".*', "\\1", gtf_dt$V9),
  gene_symbol = sub('.*gene_name "([^"]+)".*', "\\1", gtf_dt$V9)
))
gene_annot <- unique(tx2gene[, .(gene_id = GENEID, gene_symbol)])

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene[, .(TXNAME, GENEID)],
                countsFromAbundance = "no")
saveRDS(txi, file.path(out_dir, "tximport_object.rds"))
count_dt <- data.table(gene_id = rownames(txi$counts), as.data.frame(txi$counts, check.names = FALSE))
count_dt <- merge(gene_annot, count_dt, by = "gene_id", all.y = TRUE, sort = FALSE)
setcolorder(count_dt, c("gene_id", "gene_symbol", samples$sample))
fwrite(count_dt, file.path(out_dir, "gene_counts.csv"))

# Summarize Salmon technical metrics before statistical testing.
salmon_qc <- rbindlist(lapply(seq_len(nrow(samples)), function(i) {
  meta <- jsonlite::fromJSON(file.path(salmon_dir, samples$sample[i], "aux_info", "meta_info.json"))
  fmt <- jsonlite::fromJSON(file.path(salmon_dir, samples$sample[i], "lib_format_counts.json"))
  data.table(
    sample = samples$sample[i], run = samples$run[i], condition = as.character(samples$condition[i]),
    processed_fragments = meta$num_processed, mapped_fragments = meta$num_mapped,
    mapping_rate_percent = meta$percent_mapped, inferred_library_type = meta$detected_library_type,
    orphan_fragments = meta$num_orphan,
    orphan_percent_of_processed = 100 * meta$num_orphan / meta$num_processed,
    compatible_fragment_ratio = fmt$compatible_fragment_ratio,
    assigned_fragments = fmt$num_assigned_fragments,
    fragment_length_mean = meta$frag_length_mean,
    fragment_length_sd = meta$frag_length_sd
  )
}))
fwrite(salmon_qc, file.path(out_dir, "salmon_qc_summary.csv"))

# edgeR quasi-likelihood analysis. No sample is excluded by this script.
y <- DGEList(counts = txi$counts, samples = samples)
design <- model.matrix(~condition, data = samples)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "conditionSD")
de <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene_id")
de <- merge(gene_annot, de, by = "gene_id", all.x = FALSE, all.y = TRUE, sort = FALSE)
setorder(de, PValue)
fwrite(de[, .(gene_id, gene_symbol, logFC, logCPM, PValue, FDR)],
       file.path(out_dir, "edgeR_SD_vs_Control.csv"))
saveRDS(y, file.path(out_dir, "edgeR_DGEList.rds"))

# Sample-level QC plots use normalized logCPM. Group is encoded by both color
# and plotting symbol; no PCA/MDS-based exclusion is performed.
group_cols <- c(Control = "#0072B2", SD = "#D55E00")
png(file.path(out_dir, "sample_MDS.png"), width = 1800, height = 1500, res = 220)
plotMDS(y, labels = samples$sample, col = group_cols[as.character(samples$condition)],
        pch = c(Control = 16, SD = 17)[as.character(samples$condition)], main = "MDS: PRJNA1167170")
legend("topright", legend = names(group_cols), col = group_cols, pch = c(16, 17), bty = "n")
dev.off()

logcpm <- cpm(y, log = TRUE, prior.count = 2)
pca <- prcomp(t(logcpm), scale. = FALSE)
pca_df <- data.frame(sample = rownames(pca$x), condition = samples$condition,
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2])
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
p <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = condition, label = sample)) +
  geom_point(size = 3) + geom_text(nudge_y = 0.15, show.legend = FALSE) +
  scale_color_manual(values = group_cols) + scale_shape_manual(values = c(Control = 16, SD = 17)) +
  labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]),
       title = "PCA of normalized logCPM") + theme_classic(base_size = 12)
ggsave(file.path(out_dir, "sample_PCA.png"), p, width = 7, height = 6, dpi = 300)
write.csv(pca_df, file.path(out_dir, "sample_PCA_coordinates.csv"), row.names = FALSE)

lib_df <- data.frame(sample = samples$sample, condition = samples$condition,
                     library_size = colSums(txi$counts))
p <- ggplot(lib_df, aes(sample, library_size / 1e6, fill = condition)) +
  geom_col() + scale_fill_manual(values = group_cols) +
  labs(x = NULL, y = "Estimated library size (million fragments)", title = "Library sizes") +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "sample_library_sizes.png"), p, width = 7, height = 5, dpi = 300)
write.csv(lib_df, file.path(out_dir, "sample_library_sizes.csv"), row.names = FALSE)

cor_mat <- cor(logcpm, method = "pearson")
ann <- data.frame(condition = samples$condition, row.names = samples$sample)
png(file.path(out_dir, "sample_correlation_heatmap.png"), width = 1800, height = 1600, res = 220)
pheatmap(cor_mat, annotation_col = ann, annotation_row = ann,
         annotation_colors = list(condition = group_cols), border_color = NA,
         main = "Sample correlation: normalized logCPM")
dev.off()
write.csv(cor_mat, file.path(out_dir, "sample_correlation_matrix.csv"))

# Use one symbol per tested gene for GSEA, retaining the most highly expressed
# gene when a symbol is duplicated. The requested ranking statistic is logFC.
rank_dt <- de[!is.na(gene_symbol) & gene_symbol != "" & is.finite(logFC)]
setorder(rank_dt, gene_symbol, -logCPM)
rank_dt <- rank_dt[!duplicated(gene_symbol)]
stats <- rank_dt$logFC
names(stats) <- rank_dt$gene_symbol
stats <- sort(stats, decreasing = TRUE)
set.seed(1167170)

hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
reactome <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2", subcollection = "CP:REACTOME")
gobp <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M5", subcollection = "GO:BP")
msigdb_version <- unique(hallmark$db_version)

make_pathways <- function(tbl, names_wanted) {
  present <- intersect(names_wanted, unique(tbl$gs_name))
  split(tbl$gene_symbol[tbl$gs_name %in% present], tbl$gs_name[tbl$gs_name %in% present])
}
tidy_fgsea <- function(res, analysis_class) {
  ans <- as.data.table(res)
  if (!nrow(ans)) return(data.table())
  ans[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  ans[, analysis_class := analysis_class]
  setcolorder(ans, c("analysis_class", "pathway", "NES", "pval", "padj", "size", "leadingEdge"))
  ans
}

# Primary endpoint: fixed before viewing this dataset.
primary_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
primary_pathway <- make_pathways(hallmark, primary_name)
if (!length(primary_pathway)) stop("Primary Hallmark OXPHOS set is unavailable in MSigDB ", msigdb_version)
primary_res <- fgseaMultilevel(primary_pathway, stats, minSize = 10, maxSize = 500, eps = 0)
primary_out <- tidy_fgsea(primary_res, "PRIMARY")
primary_out[, direction_vs_discovery := fifelse(NES < 0, "supports", fifelse(NES > 0, "opposite direction", "does not support"))]
fwrite(primary_out, file.path(out_dir, "primary_OXPHOS_validation.csv"))

png(file.path(out_dir, "primary_OXPHOS_enrichment.png"), width = 1900, height = 1400, res = 220)
print(plotEnrichment(primary_pathway[[primary_name]], stats) +
        labs(title = "Primary validation: HALLMARK_OXIDATIVE_PHOSPHORYLATION",
             subtitle = sprintf("SD vs Control; NES = %.3f, nominal p = %.3g, FDR = %.3g",
                                primary_out$NES[1], primary_out$pval[1], primary_out$padj[1]),
             x = "Genes ranked by edgeR logFC (SD vs Control)", y = "Running enrichment score") +
        theme_classic(base_size = 12))
dev.off()

secondary_names <- c(
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION"
)
secondary_pathways <- make_pathways(reactome, secondary_names)
if (!setequal(names(secondary_pathways), secondary_names)) {
  stop("A prespecified secondary pathway is absent from MSigDB ", msigdb_version)
}
secondary_out <- tidy_fgsea(fgseaMultilevel(secondary_pathways, stats, minSize = 10, maxSize = 500, eps = 0), "SECONDARY")
fwrite(secondary_out, file.path(out_dir, "secondary_mitochondrial_pathways.csv"))

exploratory_names <- list(
  hallmark = c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
  gobp = c("GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
           "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
           "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING",
           "GOBP_HYPOXIA_INDUCIBLE_FACTOR_1ALPHA_SIGNALING_PATHWAY"),
  reactome = c("REACTOME_MITOCHONDRIAL_UNFOLDED_PROTEIN_RESPONSE_UPRMT",
               "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
               "REACTOME_REGULATION_OF_GENE_EXPRESSION_BY_HYPOXIA_INDUCIBLE_FACTOR")
)
exploratory_pathways <- c(make_pathways(hallmark, exploratory_names$hallmark),
                         make_pathways(gobp, exploratory_names$gobp),
                         make_pathways(reactome, exploratory_names$reactome))
exploratory_out <- tidy_fgsea(fgseaMultilevel(exploratory_pathways, stats, minSize = 10, maxSize = 500, eps = 0), "EXPLORATORY")
excluded_small <- setdiff(names(exploratory_pathways), exploratory_out$pathway)
if (length(excluded_small)) {
  small_rows <- data.table(
    analysis_class = "EXPLORATORY",
    pathway = excluded_small,
    NES = NA_real_, pval = NA_real_, padj = NA_real_,
    size = vapply(exploratory_pathways[excluded_small], function(x) sum(unique(x) %in% names(stats)), integer(1)),
    leadingEdge = NA_character_,
    status = "not tested: matched size below prespecified minSize=10"
  )
} else {
  small_rows <- data.table()
}
# No ATF6-specific GO/Reactome set exists in the pinned native-mouse MSigDB
# release; record this explicitly instead of inventing a gene set.
atf6_unavailable <- data.table(analysis_class = "EXPLORATORY", pathway = "ATF6_MEDIATED_UPR",
                              NES = NA_real_, pval = NA_real_, padj = NA_real_, size = NA_integer_,
                              leadingEdge = NA_character_, status = paste0("not available as a specific native-mouse MSigDB set in ", msigdb_version))
exploratory_out[, status := "tested"]
exploratory_out <- rbind(exploratory_out, small_rows, atf6_unavailable, fill = TRUE)
fwrite(exploratory_out, file.path(out_dir, "exploratory_stress_pathways.csv"))

# Save primary leading edge with corresponding edgeR statistics.
leading <- strsplit(primary_out$leadingEdge[1], ";", fixed = TRUE)[[1]]
leading_dt <- merge(data.table(gene_symbol = leading),
                    de[, .(gene_id, gene_symbol, logFC, logCPM, PValue, FDR)],
                    by = "gene_symbol", all.x = TRUE, sort = FALSE)
fwrite(leading_dt, file.path(out_dir, "OXPHOS_leading_edge.csv"))

# Compare only when one unambiguous, explicitly named prior leading-edge table
# exists under results/ or processed/. Otherwise leave the comparison absent.
search_roots <- c(file.path(root, "results"), file.path(root, "processed"))
candidates <- unlist(lapply(search_roots[file.exists(search_roots)], function(d) {
  list.files(d, pattern = "GSE289089.*leading.*edge.*\\.(csv|tsv|txt)$", recursive = TRUE,
             full.names = TRUE, ignore.case = TRUE)
}), use.names = FALSE)
if (length(candidates) == 1L) {
  prior <- fread(candidates[[1]])
  symbol_col <- intersect(c("gene_symbol", "gene", "symbol"), names(prior))
  if (length(symbol_col) == 1L) {
    overlap <- data.table(gene_symbol = intersect(leading, unique(prior[[symbol_col]])))
    fwrite(overlap, file.path(out_dir, "OXPHOS_cross_dataset_overlap.csv"))
  }
}

key_pattern <- "^(Ndufa|Ndufb|Ndufs|Ndufv|Uqcr|Cox|Atp5)"
key_exact <- c("Sdha", "Sdhb", "Sdhc", "Sdhd", "Tfam", "Ppargc1a",
               "Xbp1", "Hspa5", "Atf6", "Atf4", "Ddit3", "Eif2ak3", "Ern1", "Ern2",
               "Hif1a", "Egln1", "Egln2", "Egln3")
key <- de[grepl(key_pattern, gene_symbol) | gene_symbol %in% key_exact,
          .(gene_id, gene_symbol, logFC, logCPM, PValue, FDR)]
key[, category := fifelse(grepl(key_pattern, gene_symbol) | gene_symbol %in% c("Sdha", "Sdhb", "Sdhc", "Sdhd", "Tfam", "Ppargc1a"),
                          "mitochondrial_OXPHOS",
                          fifelse(gene_symbol %in% c("Hif1a", "Egln1", "Egln2", "Egln3"), "HIF1", "stress"))]
setcolorder(key, c("category", "gene_id", "gene_symbol", "logFC", "logCPM", "PValue", "FDR"))
fwrite(key, file.path(out_dir, "key_genes.csv"))

# Reproducibility record, including the fixed interpretation boundaries.
sink(file.path(out_dir, "reproducibility_info.txt"))
cat("Analysis: independent colon-tissue transcriptomic validation\n")
cat("Comparison: SD vs Control; logFC > 0 means higher in SD\n")
cat("Primary endpoint: HALLMARK_OXIDATIVE_PHOSPHORYLATION; expected NES < 0\n")
cat("No samples were automatically excluded.\n")
cat("tximport countsFromAbundance: no\n")
cat("GSEA ranking statistic: edgeR logFC\n")
cat("GSEA seed: 1167170\n")
cat("MSigDB release:", msigdb_version, "(native mouse collections)\n")
cat("Reference: GENCODE Mouse M39, GRCm39, Ensembl 116\n")
cat("Salmon index:", file.path(root, "reference", "salmon_mouse_GRCm39"), "\n")
cat("Salmon version: 2.7.0\n")
cat("Procedural guidance: Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026). ",
    "Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. ",
    "arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065\n", sep = "")
cat("R package versions:\n")
for (pkg in required) cat("  ", pkg, ": ", as.character(packageVersion(pkg)), "\n", sep = "")
cat("\nSession information:\n")
print(sessionInfo())
sink()
