#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: ISC_TA_regeneration_analysis.R PROJECT_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(file.path(root, "reference", "R_library"), .libPaths()))

required <- c("Seurat", "edgeR", "fgsea", "msigdbr", "data.table", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(ggplot2)
})

seed <- 20261003L
out_dir <- file.path(root, "results", "final_mito_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
outputs <- file.path(out_dir, c(
  "ISC_TA_regeneration_pathways.csv", "ISC_TA_regeneration_gene_sets.csv",
  "ISC_TA_HIF_descriptive.csv", "ISC_TA_marker_gene_changes.csv",
  "ISC_TA_regeneration_pathway_heatmap.png", "ISC_TA_regeneration_pathway_heatmap.pdf",
  "reproducibility_info.txt"
))
if (any(file.exists(outputs))) stop("Refusing to overwrite existing output(s): ", paste(outputs[file.exists(outputs)], collapse = ", "))

edge_path <- file.path(root, "results", "GSE289089_celltype_OXPHOS", "celltype_edgeR_objects.rds")
object_path <- file.path(root, "processed", "GSE289089_clustered.rds")
gtf_path <- file.path(root, "reference", "gencode_M39", "gencode.vM39.chr_patch_hapl_scaff.annotation.gtf.gz")
if (!all(file.exists(c(edge_path, object_path, gtf_path)))) stop("Required input is missing")

hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
gobp <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M5", subcollection = "GO:BP")
reactome <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2", subcollection = "CP:REACTOME")
setDT(hallmark)
setDT(gobp)
setDT(reactome)
msigdb_version <- unique(hallmark$db_version)
if (length(msigdb_version) != 1L) stop("Could not determine one MSigDB release")

# Frozen before testing: one direct standard set for each prespecified biological family.
manifest <- data.table(
  family = c("WNT signaling", "intestinal stem-cell maintenance", "epithelial proliferation",
             "epithelial regeneration / wound repair", "MYC targets", "mTORC1 signaling",
             "ribosome biogenesis", "translation", "cell cycle", "DNA replication"),
  pathway = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "GOBP_INTESTINAL_STEM_CELL_HOMEOSTASIS",
              "GOBP_EPITHELIAL_CELL_PROLIFERATION", "GOBP_WOUND_HEALING",
              "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING",
              "GOBP_RIBOSOME_BIOGENESIS", "REACTOME_TRANSLATION",
              "REACTOME_CELL_CYCLE", "REACTOME_DNA_REPLICATION"),
  source = c("MSigDB Hallmark", "Gene Ontology BP", "Gene Ontology BP", "Gene Ontology BP",
             "MSigDB Hallmark", "MSigDB Hallmark", "Gene Ontology BP", "Reactome", "Reactome", "Reactome"),
  analysis_role = "primary regeneration coupling"
)
all_sets <- rbind(
  hallmark[, .(pathway = gs_name, gene_symbol)],
  gobp[, .(pathway = gs_name, gene_symbol)],
  reactome[, .(pathway = gs_name, gene_symbol)]
)
if (length(setdiff(manifest$pathway, all_sets$pathway))) stop("A prespecified pathway is unavailable: ", paste(setdiff(manifest$pathway, all_sets$pathway), collapse = ", "))
pathways <- split(all_sets[pathway %in% manifest$pathway, gene_symbol],
                  all_sets[pathway %in% manifest$pathway, pathway])
pathways <- lapply(pathways, unique)
manifest[, database_gene_set_size := lengths(pathways[pathway])]
manifest[, MSigDB_version := msigdb_version]
fwrite(manifest, outputs[2])

edge_objects <- readRDS(edge_path)
cell_groups <- c(ISC = "ISC", TA = "TA_cycling_progenitor")
make_rank <- function(de, method) {
  score <- if (method == "logFC") "logFC" else "signed_stat"
  z <- copy(de[is.finite(get(score)) & !is.na(current_gene_symbol) & current_gene_symbol != ""])
  setorder(z, current_gene_symbol, -logCPM)
  z <- z[!duplicated(current_gene_symbol)]
  ans <- z[[score]]
  names(ans) <- z$current_gene_symbol
  sort(ans, decreasing = TRUE)
}
run_collection <- function(rank, cell_group, method, seed_value) {
  set.seed(seed_value)
  fg <- as.data.table(fgseaMultilevel(pathways, rank, minSize = 5, maxSize = 1000, eps = 0))
  fg[, `:=`(cell_group = cell_group, ranking_method = method,
            leadingEdge_count = lengths(leadingEdge),
            leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)))]
  found <- fg[, .(pathway, cell_group, ranking_method, NES, pval, padj, size,
                  leadingEdge_count, leadingEdge, status = "tested")]
  missing_sets <- setdiff(manifest$pathway, found$pathway)
  if (length(missing_sets)) {
    matched <- vapply(pathways[missing_sets], function(x) sum(unique(x) %in% names(rank)), integer(1))
    found <- rbind(found, data.table(pathway = missing_sets, cell_group = cell_group,
      ranking_method = method, NES = NA_real_, pval = NA_real_, padj = NA_real_, size = matched,
      leadingEdge_count = NA_integer_, leadingEdge = NA_character_,
      status = "not tested: fewer than 5 matched genes"), fill = TRUE)
  }
  found
}

long_rows <- list()
markers <- c("Lgr5", "Olfm4", "Ascl2", "Axin2", "Myc", "Mki67", "Pcna")
marker_rows <- list()
k <- 0L
for (label in names(cell_groups)) {
  de <- as.data.table(edge_objects[[cell_groups[[label]]]]$de)
  de[, signed_stat := sign(logFC) * sqrt(F)]
  marker_rows[[label]] <- de[current_gene_symbol %in% markers,
    .(cell_group = label, gene = current_gene_symbol, logFC, logCPM, F, signed_stat, PValue, FDR)]
  for (method in c("logFC", "signed_stat")) {
    k <- k + 1L
    long_rows[[k]] <- run_collection(make_rank(de, method), label, method, seed + k)
  }
}
long <- rbindlist(long_rows)
long <- merge(long, manifest[, .(pathway, family, source)], by = "pathway", all.x = TRUE, sort = FALSE)
log_rows <- copy(long[ranking_method == "logFC"])
setnames(log_rows, c("NES", "pval", "padj", "size", "leadingEdge_count", "leadingEdge", "status"),
         c("NES_logFC", "p_logFC", "FDR_logFC", "size_logFC", "leadingEdge_count_logFC", "leadingEdge_logFC", "status_logFC"))
signed_rows <- copy(long[ranking_method == "signed_stat"])
setnames(signed_rows, c("NES", "pval", "padj", "size", "leadingEdge_count", "leadingEdge", "status"),
         c("NES_signedStat", "p_signedStat", "FDR_signedStat", "size_signedStat", "leadingEdge_count_signedStat", "leadingEdge_signedStat", "status_signedStat"))
wide <- merge(
  log_rows[, .(pathway, family, source, cell_group, NES_logFC, p_logFC, FDR_logFC, size_logFC,
               leadingEdge_count_logFC, leadingEdge_logFC, status_logFC)],
  signed_rows[, .(pathway, family, source, cell_group, NES_signedStat, p_signedStat, FDR_signedStat,
                  size_signedStat, leadingEdge_count_signedStat, leadingEdge_signedStat, status_signedStat)],
  by = c("pathway", "family", "source", "cell_group"), all = TRUE, sort = FALSE
)
wide <- merge(wide, manifest[, .(pathway, family_order = .I)], by = "pathway", all.x = TRUE)
wide <- wide[order(match(cell_group, c("ISC", "TA")), family_order)]
wide[, family_order := NULL]
fwrite(wide, outputs[1], na = "NA")

marker_dt <- rbindlist(marker_rows, fill = TRUE)
marker_dt <- merge(CJ(cell_group = c("ISC", "TA"), gene = markers, unique = TRUE), marker_dt,
                   by = c("cell_group", "gene"), all.x = TRUE, sort = FALSE)
marker_dt <- marker_dt[order(match(cell_group, c("ISC", "TA")), match(gene, markers))]
fwrite(marker_dt, outputs[4], na = "NA")

# Optional descriptive HIF/hypoxia check is isolated from the primary pathway family.
hif_name <- "HALLMARK_HYPOXIA"
hif_set <- unique(hallmark[gs_name == hif_name, gene_symbol])
hif_rows <- list()
k <- 0L
for (label in names(cell_groups)) {
  de <- as.data.table(edge_objects[[cell_groups[[label]]]]$de)
  de[, signed_stat := sign(logFC) * sqrt(F)]
  for (method in c("logFC", "signed_stat")) {
    k <- k + 1L
    rank <- make_rank(de, method)
    set.seed(seed + 100L + k)
    fg <- as.data.table(fgseaMultilevel(setNames(list(hif_set), hif_name), rank,
                                       minSize = 10, maxSize = 500, eps = 0))
    hif_rows[[k]] <- data.table(pathway = hif_name, cell_group = label, ranking_method = method,
      NES = fg$NES, pval = fg$pval, padj = fg$padj, size = fg$size,
      leadingEdge_count = length(fg$leadingEdge[[1]]),
      note = "descriptive only; Hif1a mRNA is not a measure of HIF1alpha protein activity")
  }
}
fwrite(rbindlist(hif_rows), outputs[3])

plot_dt <- copy(long)
plot_dt[, pathway_label := factor(family, levels = rev(manifest$family))]
plot_dt[, column := factor(paste(cell_group, fifelse(ranking_method == "logFC", "logFC", "signed stat"), sep = " | "),
                           levels = c("ISC | logFC", "ISC | signed stat", "TA | logFC", "TA | signed stat"))]
limit <- max(abs(plot_dt$NES), na.rm = TRUE)
p_heat <- ggplot(plot_dt, aes(column, pathway_label, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = ifelse(is.na(NES), "NA", sprintf("%.2f%s", NES, ifelse(padj < 0.05, "*", "")))), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-limit, limit), na.value = "#D9D9D9", name = "NES") +
  labs(x = NULL, y = NULL, title = "ISC and TA regeneration-associated programs",
       subtitle = "* FDR < 0.05 within the 10 prespecified pathway family") +
  theme_minimal(base_size = 9.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 28, hjust = 1),
        legend.position = "right")
ggsave(outputs[5], p_heat, width = 8.2, height = 5.8, dpi = 300, bg = "white")
ggsave(outputs[6], p_heat, width = 8.2, height = 5.8, device = "pdf", useDingbats = FALSE)

sink(outputs[7])
cat("Analysis: final OXPHOS concordance, nuclear-only sensitivity, sample-level QC, equal-cell downsampling, and ISC/TA regeneration coupling\n")
cat("Date/time:", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), "\n")
cat("Project root:", root, "\n")
cat("Random seed:", seed, "\n")
cat("Downsampling iterations: 20\n")
cat("Biological replication: GSE289089 2 Control vs 2 SD; PRJNA1167170 4 Control vs 4 SD.\n")
cat("Input Seurat object:", object_path, "\n")
cat("Input cell-type edgeR objects:", edge_path, "\n")
cat("Input PRJNA1167170 tximport object:", file.path(root, "results", "PRJNA1167170_bulk", "tximport_object.rds"), "\n")
cat("Input PRJNA1167170 DE annotation:", file.path(root, "results", "PRJNA1167170_bulk", "edgeR_SD_vs_Control.csv"), "\n")
cat("GENCODE release: M39; genome GRCm39; annotation:", gtf_path, "\n")
cat("MSigDB release:", msigdb_version, "native mouse; msigdbr", as.character(packageVersion("msigdbr")), "\n")
cat("GSEA: fgseaMultilevel, eps=0; OXPHOS minSize=10/maxSize=500; regeneration minSize=5/maxSize=1000.\n")
cat("Cell-cycle score fields were not present in the Seurat metadata and are reported as NA.\n")
cat("Procedural guidance: Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065\n")
cat("Package versions:\n")
for (pkg in required) cat("  ", pkg, ": ", as.character(packageVersion(pkg)), "\n", sep = "")
cat("R version:", R.version.string, "\n\nSession information:\n")
print(sessionInfo())
sink()

message("ISC/TA regeneration analysis completed. MSigDB ", msigdb_version)
