#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
object_file <- file.path(project_root, "processed", "GSE289089_clustered.rds")
output_dir <- file.path(project_root, "results", "manuscript_v1")

marker_file <- file.path(output_dir, "GSE289089_FindAllMarkers_reconstructed.csv")
annotation_file <- file.path(output_dir, "GSE289089_cluster_annotation_evidence_reconstructed.csv")
dotplot_png <- file.path(output_dir, "GSE289089_annotation_DotPlot.png")
dotplot_pdf <- file.path(output_dir, "GSE289089_annotation_DotPlot.pdf")

planned_outputs <- c(marker_file, annotation_file, dotplot_png, dotplot_pdf)
if (all(file.exists(planned_outputs))) {
  stop(
    "Refusing to overwrite the complete existing reconstruction output set."
  )
}

combined.clean <- readRDS(object_file)
if (!inherits(combined.clean, "Seurat")) {
  stop("Frozen input is not a Seurat object.")
}
if (ncol(combined.clean) != 27671L) {
  stop("Frozen input does not contain the expected 27,671 singlets.")
}
if (!"seurat_clusters" %in% colnames(combined.clean[[]])) {
  stop("Frozen input does not contain seurat_clusters.")
}

DefaultAssay(combined.clean) <- "RNA"
Idents(combined.clean) <- "seurat_clusters"
expected_clusters <- as.character(0:17)
observed_clusters <- as.character(sort(as.integer(unique(as.character(Idents(combined.clean))))))
if (!identical(observed_clusters, expected_clusters)) {
  stop("Frozen cluster identities are not exactly 0 through 17.")
}

# Exact marker invocation recovered from the contemporaneous analysis record.
# A completed marker file is reused when resuming after a graphics-device failure.
if (file.exists(marker_file)) {
  all_markers <- read.csv(marker_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  all_markers <- FindAllMarkers(
    combined.clean,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
}

if (!"gene" %in% colnames(all_markers)) {
  all_markers$gene <- rownames(all_markers)
}

marker_panels <- list(
  "Stem / ISC" = c("Lgr5", "Olfm4", "Ascl2", "Smoc2"),
  "Cycling" = c("Mki67", "Top2a", "Pcna"),
  "Paneth" = c("Lyz1", "Mmp7"),
  "Goblet" = c("Muc2", "Spdef", "Agr2", "Clca1"),
  "Enterocyte" = c("Alpi", "Apoa1", "Fabp1", "Krt20"),
  "EEC / EC" = c("Chga", "Chgb", "Neurog3", "Tph1"),
  "Tuft" = c("Dclk1", "Pou2f3", "Trpm5"),
  "Immune" = c("Ptprc", "Lst1", "Cd3d"),
  "Endothelial" = c("Pecam1", "Kdr")
)
canonical_markers <- unlist(marker_panels, use.names = FALSE)

working_annotations <- c(
  "0" = "stress-responsive epithelial",
  "1" = "IFN-responsive/regenerative epithelial",
  "2" = "cycling TA/progenitor",
  "3" = "secretory/Goblet-like",
  "4" = "mature enterocyte",
  "5" = "Paneth",
  "6" = "cycling ISC/early progenitor",
  "7" = "cycling TA/progenitor",
  "8" = "Tuft",
  "9" = "enteroendocrine",
  "10" = "cycling TA/progenitor",
  "11" = "cycling TA/progenitor",
  "12" = "Tph1+ serotonergic EC",
  "13" = "cycling immature epithelial progenitor",
  "14" = "EEC progenitor",
  "15" = "mature absorptive enterocyte",
  "16" = "cycling Paneth-lineage/Paneth progenitor",
  "17" = "immune"
)

effect_column <- if ("avg_log2FC" %in% colnames(all_markers)) {
  "avg_log2FC"
} else if ("avg_logFC" %in% colnames(all_markers)) {
  "avg_logFC"
} else {
  stop("Reconstructed marker output lacks an average log-fold-change column.")
}

summarize_cluster <- function(cluster_id) {
  cluster_markers <- all_markers[as.character(all_markers$cluster) == cluster_id, , drop = FALSE]
  marker_order <- order(
    -cluster_markers[[effect_column]],
    cluster_markers$p_val_adj,
    cluster_markers$gene
  )
  cluster_markers <- cluster_markers[marker_order, , drop = FALSE]
  representative <- head(unique(cluster_markers$gene), 10L)
  supporting <- canonical_markers[canonical_markers %in% cluster_markers$gene]

  data.frame(
    cluster = cluster_id,
    working_annotation = unname(working_annotations[[cluster_id]]),
    representative_top_markers = if (length(representative)) {
      paste(representative, collapse = "; ")
    } else {
      "none recovered"
    },
    canonical_supporting_markers = if (length(supporting)) {
      paste(supporting, collapse = "; ")
    } else {
      "none among reconstructed positive markers at the original thresholds"
    },
    annotation_note = paste(
      "Working annotation retained unchanged; evidence summarized from the",
      "reconstructed positive-marker table using the original parameters."
    ),
    stringsAsFactors = FALSE
  )
}

annotation_evidence <- do.call(rbind, lapply(expected_clusters, summarize_cluster))

write_atomic_csv <- function(x, destination) {
  temporary <- tempfile(pattern = paste0(basename(destination), "."), tmpdir = output_dir)
  on.exit(unlink(temporary), add = TRUE)
  write.csv(x, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, destination)) {
    stop("Could not move temporary CSV into place: ", destination)
  }
}

if (!file.exists(marker_file)) {
  write_atomic_csv(all_markers, marker_file)
}
if (!file.exists(annotation_file)) {
  write_atomic_csv(annotation_evidence, annotation_file)
}

present_panels <- lapply(marker_panels, function(features) {
  features[features %in% rownames(combined.clean)]
})
missing_features <- setdiff(canonical_markers, rownames(combined.clean))

annotation_dotplot <- DotPlot(
  combined.clean,
  features = present_panels,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  scale_color_viridis_c(option = "C", name = "Scaled average\nexpression") +
  labs(
    title = "GSE289089 canonical marker evidence by frozen cluster",
    x = NULL,
    y = "Frozen Seurat cluster",
    size = "Percent\nexpressed"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 9)
  )

png_temporary <- tempfile(pattern = "GSE289089_annotation_DotPlot_", tmpdir = output_dir, fileext = ".png")
pdf_temporary <- tempfile(pattern = "GSE289089_annotation_DotPlot_", tmpdir = output_dir, fileext = ".pdf")
on.exit(unlink(c(png_temporary, pdf_temporary)), add = TRUE)

if (!file.exists(dotplot_png)) {
  ggsave(
    filename = png_temporary,
    plot = annotation_dotplot,
    width = 15,
    height = 8.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

if (!file.exists(dotplot_pdf)) {
  ggsave(
    filename = pdf_temporary,
    plot = annotation_dotplot,
    device = grDevices::pdf,
    width = 15,
    height = 8.5,
    units = "in",
    bg = "white"
  )
}

if (!file.exists(dotplot_png) && !file.rename(png_temporary, dotplot_png)) {
  stop("Could not move temporary PNG into place.")
}
if (!file.exists(dotplot_pdf) && !file.rename(pdf_temporary, dotplot_pdf)) {
  stop("Could not move temporary PDF into place.")
}

message("Reconstructed marker rows: ", nrow(all_markers))
message("Clusters summarized: ", nrow(annotation_evidence))
message(
  "Canonical panel genes absent from the frozen RNA assay: ",
  if (length(missing_features)) paste(missing_features, collapse = ", ") else "none"
)
