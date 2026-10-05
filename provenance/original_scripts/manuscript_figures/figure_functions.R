suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(ggrepel)
  library(scales)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "processed", "GSE289089_clustered.rds"))) {
  stop("Run figure scripts from the sleep_gut_project root.")
}

main_dir <- file.path(project_root, "results", "manuscript_v1", "figures_main")
supp_dir <- file.path(project_root, "results", "manuscript_v1", "figures_supplement")

condition_colors <- c(Control = "#0072B2", SD = "#D55E00")
condition_shapes <- c(Control = 16, SD = 17)

cell_order <- c(
  "ISC", "TA_cycling_progenitor", "Goblet_secretory", "Paneth_lineage",
  "Enterocyte", "Enteroendocrine", "Serotonergic_EC", "Tuft",
  "Stress_state_cluster0", "IFN_regenerative_cluster1",
  "EEC_progenitor_cluster14", "Immune_excluded"
)
cell_labels <- c(
  ISC = "ISC", TA_cycling_progenitor = "TA", Goblet_secretory = "Goblet",
  Paneth_lineage = "Paneth", Enterocyte = "Enterocyte",
  Enteroendocrine = "EEC", Serotonergic_EC = "Serotonergic EC", Tuft = "Tuft",
  Stress_state_cluster0 = "Stress-responsive", IFN_regenerative_cluster1 = "IFN/regenerative",
  EEC_progenitor_cluster14 = "EEC progenitor", Immune_excluded = "Immune"
)
cell_colors <- c(
  ISC = "#0072B2", TA_cycling_progenitor = "#E69F00",
  Goblet_secretory = "#009E73", Paneth_lineage = "#CC79A7",
  Enterocyte = "#D55E00", Enteroendocrine = "#56B4E9",
  Serotonergic_EC = "#AA3377", Tuft = "#CCBB44",
  Stress_state_cluster0 = "#777777", IFN_regenerative_cluster1 = "#332288",
  EEC_progenitor_cluster14 = "#44AA99", Immune_excluded = "#000000"
)

pathway_order <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION"
)
pathway_labels <- c(
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = "OXPHOS",
  REACTOME_RESPIRATORY_ELECTRON_TRANSPORT = "Respiratory electron\ntransport",
  REACTOME_MITOCHONDRIAL_TRANSLATION = "Mitochondrial\ntranslation",
  REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING = "ATP synthesis",
  REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE = "TCA",
  REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION = "FA beta-oxidation"
)

theme_manuscript <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F2F2F2", color = "#777777", linewidth = 0.35),
      strip.text = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(size = rel(0.9), color = "#333333"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 14, family = "Helvetica"),
  plot.tag.position = c(0.01, 0.99)
)

save_pair <- function(plot, stem, directory, width, height) {
  png_file <- file.path(directory, paste0(stem, ".png"))
  pdf_file <- file.path(directory, paste0(stem, ".pdf"))
  if (file.exists(png_file) || file.exists(pdf_file)) {
    stop("Refusing to overwrite existing figure: ", stem)
  }
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  if (Sys.info()[["sysname"]] == "Darwin" && nzchar(Sys.which("sips"))) {
    status <- system2("sips", c("-s", "dpiWidth", "300", "-s", "dpiHeight", "300", png_file),
                      stdout = FALSE, stderr = FALSE)
    if (!identical(status, 0L)) stop("Could not set PNG DPI metadata: ", png_file)
  }
  ggsave(pdf_file, plot, width = width, height = height, units = "in", device = grDevices::pdf, bg = "white")
  invisible(c(pdf_file, png_file))
}

save_pair_replace <- function(plot, stem, directory, width, height) {
  png_file <- file.path(directory, paste0(stem, ".png"))
  pdf_file <- file.path(directory, paste0(stem, ".pdf"))
  png_tmp <- file.path(directory, paste0(".", stem, ".tmp.png"))
  pdf_tmp <- file.path(directory, paste0(".", stem, ".tmp.pdf"))
  on.exit(unlink(c(png_tmp, pdf_tmp)), add = TRUE)
  ggsave(png_tmp, plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  if (Sys.info()[["sysname"]] == "Darwin" && nzchar(Sys.which("sips"))) {
    status <- system2("sips", c("-s", "dpiWidth", "300", "-s", "dpiHeight", "300", png_tmp),
                      stdout = FALSE, stderr = FALSE)
    if (!identical(status, 0L)) stop("Could not set PNG DPI metadata: ", png_tmp)
  }
  ggsave(pdf_tmp, plot, width = width, height = height, units = "in", device = grDevices::pdf, bg = "white")
  if (!file.rename(png_tmp, png_file)) stop("Could not replace PNG: ", png_file)
  if (!file.rename(pdf_tmp, pdf_file)) stop("Could not replace PDF: ", pdf_file)
  invisible(c(pdf_file, png_file))
}

cache <- new.env(parent = emptyenv())
get_sc_object <- function() {
  if (!exists("sc", envir = cache, inherits = FALSE)) {
    assign("sc", readRDS(file.path(project_root, "processed", "GSE289089_clustered.rds")), envir = cache)
  }
  get("sc", envir = cache, inherits = FALSE)
}

get_oxphos_pathway <- function() {
  if (!exists("oxphos", envir = cache, inherits = FALSE)) {
    hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
    genes <- unique(hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"])
    assign("oxphos", genes, envir = cache)
  }
  get("oxphos", envir = cache, inherits = FALSE)
}

add_group_metadata <- function(obj) {
  mapping <- fread(file.path(project_root, "results", "GSE289089_celltype_OXPHOS", "cluster_to_cell_group_mapping.csv"))
  lookup <- setNames(mapping$cell_group, as.character(mapping$cluster))
  obj$manuscript_group <- unname(lookup[as.character(obj$seurat_clusters)])
  obj$manuscript_group_label <- unname(cell_labels[obj$manuscript_group])
  obj
}

composition_plot <- function(full = FALSE) {
  x <- fread(file.path(project_root, "results", "GSE289089_celltype_OXPHOS", "celltype_sample_counts.csv"))
  groups <- if (full) cell_order else cell_order[1:8]
  x <- x[cell_group %in% groups]
  x[, proportion := cell_count / sum(cell_count), by = sample]
  x[, cell_group := factor(cell_group, levels = groups)]
  x[, sample := factor(sample, levels = c("0D_1", "0D_2", "2D_1", "2D_2"))]
  ggplot(x, aes(sample, proportion, fill = cell_group)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = cell_colors, labels = cell_labels, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Cell proportion", fill = "Cell group",
         subtitle = "Descriptive cell proportions; no inferential test") +
    theme_manuscript() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.key.height = unit(3.5, "mm"))
}

workflow_plot <- function() {
  boxes <- data.frame(
    x = c(0.08, 0.25, 0.42, 0.60, 0.78, 0.93),
    y = rep(0.52, 6),
    label = c(
      "GSE289089 discovery\n4 samples\n27,671 singlets",
      "Cell-group pseudobulk\nand OXPHOS discovery",
      "PRJNA1167170\nproximal colon; 8 samples\nfirst validation",
      "Systematic search\n23 candidates\n3 eligible / 1 potential / 19 excluded",
      "PRJNA862187\ncolon; 12 samples\nfrozen-plan validation",
      "Cross-dataset\nconcordance and\nsensitivity"
    ),
    branch = c("single-cell", "single-cell", "bulk/validation", "search", "bulk/validation", "integration")
  )
  arrows <- data.frame(
    x = c(0.15, 0.32, 0.49, 0.67, 0.85),
    xend = c(0.18, 0.35, 0.53, 0.71, 0.88),
    y = rep(0.52, 5), yend = rep(0.52, 5)
  )
  ggplot() +
    geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(2.3, "mm")), linewidth = 0.45, color = "#555555") +
    geom_label(data = boxes, aes(x, y, label = label, fill = branch),
               size = 2.75, lineheight = 0.95, linewidth = 0.35,
               label.padding = unit(2.3, "mm"), color = "black") +
    scale_fill_manual(values = c("single-cell" = "#E8F3F8", "bulk/validation" = "#FBEBDD",
                                 "search" = "#EDF4F1", "integration" = "#EEE8F5"), guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() + theme(plot.margin = margin(8, 8, 8, 8))
}

annotation_dotplot <- function(base_size = 8) {
  obj <- get_sc_object()
  Idents(obj) <- "seurat_clusters"
  marker_panels <- list(
    "ISC" = c("Lgr5", "Olfm4", "Ascl2", "Smoc2"),
    "Cycling" = c("Mki67", "Top2a", "Pcna"),
    "Goblet" = c("Muc2", "Spdef", "Agr2", "Clca1"),
    "Paneth" = c("Lyz1", "Mmp7"),
    "Enterocyte" = c("Alpi", "Apoa1", "Fabp1", "Krt20"),
    "EEC" = c("Chga", "Chgb", "Neurog3"),
    "5-HT EC" = "Tph1",
    "Tuft" = c("Dclk1", "Pou2f3", "Trpm5"),
    "Immune" = c("Ptprc", "Lst1")
  )
  DotPlot(obj, features = marker_panels, group.by = "seurat_clusters", dot.scale = 5) +
    scale_color_gradient2(low = "#3B4CC0", mid = "#F2F2F2", high = "#B40426", midpoint = 0,
                          name = "Scaled average\nexpression") +
    labs(x = NULL, y = "Frozen cluster", size = "Percent\nexpressed") +
    theme_bw(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_line(color = "#E6E6E6", linewidth = 0.25),
      strip.background = element_rect(fill = "#F2F2F2", color = "#777777", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = base_size),
      legend.position = "right", plot.margin = margin(5, 5, 5, 5)
    )
}

nes_heatmap <- function(x, groups, show_values = TRUE) {
  x <- copy(x[cell_group %in% groups & pathway %in% pathway_order])
  x[, cell_group := factor(cell_group, levels = rev(groups))]
  x[, pathway := factor(pathway, levels = pathway_order)]
  lim <- max(abs(x$NES), na.rm = TRUE)
  x[, label := ifelse(is.na(NES), "NA", sprintf("%.2f%s", NES, ifelse(padj < 0.05, "  +", "")))]
  p <- ggplot(x, aes(pathway, cell_group, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lim, lim), na.value = "#BDBDBD", name = "NES") +
    scale_x_discrete(labels = pathway_labels) +
    scale_y_discrete(labels = cell_labels) +
    labs(x = NULL, y = NULL, caption = "+ FDR < 0.05") +
    theme_manuscript(8) +
    theme(axis.text.x = element_text(angle = 38, hjust = 1),
          panel.grid = element_blank(), plot.caption = element_text(hjust = 0, size = 7))
  if (show_values) p <- p + geom_text(aes(label = label), size = 2.25)
  p
}

make_celltype_rank <- function(group_name) {
  edges <- readRDS(file.path(project_root, "results", "GSE289089_celltype_OXPHOS", "celltype_edgeR_objects.rds"))
  de <- as.data.table(edges[[group_name]]$de)
  z <- copy(de[!is.na(current_gene_symbol) & current_gene_symbol != "" & is.finite(logFC)])
  setorder(z, current_gene_symbol, -logCPM)
  z <- z[!duplicated(current_gene_symbol)]
  ans <- z$logFC
  names(ans) <- z$current_gene_symbol
  sort(ans, decreasing = TRUE)
}

make_bulk_rank <- function() {
  de <- fread(file.path(project_root, "results", "PRJNA1167170_bulk", "edgeR_SD_vs_Control.csv"))
  de <- de[!is.na(gene_symbol) & gene_symbol != "" & is.finite(logFC)]
  setorder(de, gene_symbol, -logCPM)
  de <- de[!duplicated(gene_symbol)]
  ans <- de$logFC
  names(ans) <- de$gene_symbol
  sort(ans, decreasing = TRUE)
}

enrichment_plot <- function(rank, expected_nes, expected_fdr, title) {
  p <- fgsea::plotEnrichment(get_oxphos_pathway(), rank)
  p +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "#777777") +
    labs(title = title,
         subtitle = sprintf("NES = %.3f; FDR = %.3g", expected_nes, expected_fdr),
         x = "Genes ranked by edgeR logFC (SD vs Control)", y = "Running enrichment score") +
    theme_manuscript(8)
}

make_figure_1 <- function() {
  obj <- add_group_metadata(get_sc_object())
  Idents(obj) <- "seurat_clusters"
  p_b <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                 repel = TRUE, pt.size = 0.18, raster = FALSE) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster") + theme_manuscript(8)
  p_c <- DimPlot(obj, reduction = "umap", group.by = "manuscript_group", label = FALSE,
                 pt.size = 0.18, raster = FALSE, cols = cell_colors) +
    scale_color_manual(values = cell_colors, labels = cell_labels, drop = FALSE) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Manuscript cell group") + theme_manuscript(8)
  p <- workflow_plot() /
    (p_b | p_c) /
    annotation_dotplot(8) /
    composition_plot(FALSE) +
    plot_layout(heights = c(0.65, 1, 1.05, 0.8), guides = "keep") +
    plot_annotation(tag_levels = "A", theme = tag_theme)
  save_pair_replace(p, "Figure_1", main_dir, 16, 15)
}

make_figure_2 <- function() {
  primary_groups <- cell_order[1:8]
  x <- fread(file.path(project_root, "results", "GSE289089_celltype_OXPHOS", "celltype_OXPHOS_logFC_ranking.csv"))
  x <- x[cell_group %in% primary_groups]
  x[, cell_group := factor(cell_group, levels = rev(primary_groups))]
  x[, fdr_label := sprintf("FDR %.2g", padj)]
  p_a <- ggplot(x, aes(NES, cell_group, fill = cell_group)) +
    geom_col(width = 0.62, color = "black", linewidth = 0.25) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    geom_point(data = x[padj < 0.05], aes(x = NES - 0.08, y = cell_group), inherit.aes = FALSE,
               shape = 18, size = 2.1, color = "black") +
    geom_text(aes(x = 0.04, label = sprintf("NES %.2f; %s", NES, fdr_label)),
              hjust = 0, size = 2.4) +
    scale_fill_manual(values = cell_colors, guide = "none") +
    scale_y_discrete(labels = cell_labels) +
    scale_x_continuous(limits = c(min(x$NES) - 0.25, 0.9)) +
    labs(x = "Hallmark OXPHOS NES", y = NULL, caption = "Diamond: FDR < 0.05") +
    theme_manuscript(8)

  mito <- fread(file.path(project_root, "results", "GSE289089_celltype_OXPHOS", "celltype_mitochondrial_pathways.csv"))
  p_b <- nes_heatmap(mito, primary_groups, TRUE)

  curve_groups <- c("ISC", "TA_cycling_progenitor", "Tuft", "Serotonergic_EC")
  formal <- fread(file.path(project_root, "results", "GSE289089_celltype_OXPHOS", "celltype_OXPHOS_logFC_ranking.csv"))
  curves <- lapply(curve_groups, function(g) {
    row <- formal[cell_group == g]
    enrichment_plot(make_celltype_rank(g), row$NES, row$padj, unname(cell_labels[g]))
  })
  p <- (p_a | p_b) / wrap_plots(curves, nrow = 1) +
    plot_layout(heights = c(1.12, 0.88), widths = c(0.82, 1.18)) +
    plot_annotation(tag_levels = "A", theme = tag_theme)
  save_pair(p, "Figure_2", main_dir, 16, 10.5)
}

select_shared_genes <- function() {
  ann <- fread(file.path(project_root, "results", "cross_dataset_OXPHOS", "OXPHOS_overlap_functional_annotation.csv"))
  quota <- c(
    "Complex I / NDUF" = 3, "Complex II / SDH" = 1, "Complex III / UQCR" = 3,
    "Complex IV / COX" = 3, "ATP synthase / ATP5" = 3,
    "mitochondrial translation" = 3, "TCA" = 3, "fatty acid oxidation" = 2,
    "other mitochondrial metabolism" = 2
  )
  selected <- rbindlist(lapply(names(quota), function(category) {
    z <- ann[functional_category == category][order(gene_symbol)]
    head(z, quota[[category]])
  }))
  selected[, display_category := factor(functional_category, levels = names(quota))]
  selected
}

bulk_leading_edge_heatmap <- function() {
  selected <- select_shared_genes()
  dge <- readRDS(file.path(project_root, "results", "PRJNA1167170_bulk", "edgeR_DGEList.rds"))
  lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  de <- fread(file.path(project_root, "results", "PRJNA1167170_bulk", "edgeR_SD_vs_Control.csv"))
  setorder(de, gene_symbol, -logCPM)
  de <- de[!duplicated(gene_symbol)]
  ids <- de$gene_id[match(selected$gene_symbol, de$gene_symbol)]
  mat <- lcpm[match(ids, rownames(lcpm)), , drop = FALSE]
  rownames(mat) <- selected$gene_symbol
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  long <- as.data.table(as.table(z))
  setnames(long, c("gene", "sample", "z"))
  long[, gene := factor(gene, levels = rev(selected$gene_symbol))]
  conditions <- setNames(as.character(dge$samples$condition), rownames(dge$samples))
  long[, sample_label := paste0(sample, "\n", conditions[as.character(sample)])]
  long[, sample_label := factor(sample_label,
                                levels = paste0(colnames(mat), "\n", conditions[colnames(mat)]))]
  lim <- max(abs(long$z), na.rm = TRUE)
  ggplot(long, aes(sample_label, gene, fill = z)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lim, lim), name = "Row z-score\nof logCPM") +
    labs(x = NULL, y = NULL, subtitle = "Module-balanced alphabetical selection from 90 shared leading-edge genes") +
    theme_manuscript(7.5) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
}

make_figure_3 <- function() {
  pca <- fread(file.path(project_root, "results", "PRJNA1167170_bulk", "sample_PCA_coordinates.csv"))
  p_a <- ggplot(pca, aes(PC1, PC2, color = condition, shape = condition, label = sample)) +
    geom_hline(yintercept = 0, color = "#DDDDDD", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "#DDDDDD", linewidth = 0.3) +
    geom_point(size = 3) + geom_text_repel(size = 2.6, show.legend = FALSE) +
    scale_color_manual(values = condition_colors) + scale_shape_manual(values = condition_shapes) +
    labs(x = "PC1", y = "PC2", color = NULL, shape = NULL,
         subtitle = "Four biological replicates per condition") + theme_manuscript()

  primary <- fread(file.path(project_root, "results", "PRJNA1167170_bulk", "primary_OXPHOS_validation.csv"))
  p_b <- enrichment_plot(make_bulk_rank(), primary$NES[1], primary$padj[1], "Hallmark OXPHOS")

  sec <- fread(file.path(project_root, "results", "PRJNA1167170_bulk", "secondary_mitochondrial_pathways.csv"))
  sec[, pathway := factor(pathway, levels = rev(pathway_order[-1]))]
  p_c <- ggplot(sec, aes(NES, pathway)) +
    geom_segment(aes(x = 0, xend = NES, yend = pathway), color = "#777777", linewidth = 0.6) +
    geom_point(aes(fill = padj < 0.05), shape = 21, size = 3, color = "black") +
    geom_text(aes(x = 0.03, label = sprintf("NES %.2f; FDR %.2g", NES, padj)),
              hjust = 0, size = 2.4) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    scale_y_discrete(labels = pathway_labels) +
    scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "white"), name = "FDR < 0.05") +
    scale_x_continuous(limits = c(min(sec$NES) - 0.2, 0.85)) +
    labs(x = "NES", y = NULL) + theme_manuscript(8)

  p <- (p_a | p_b) / (p_c | bulk_leading_edge_heatmap()) +
    plot_annotation(tag_levels = "A", theme = tag_theme)
  save_pair(p, "Figure_3", main_dir, 14, 10.5)
}

representative_shared_15 <- function() {
  wanted <- c(
    "Ndufa3", "Ndufa1", "Ndufs4", "Uqcrb", "Uqcr11", "Cox17", "Cox7c", "Cox6c",
    "Atp5me", "Atp5mf", "Atp5f1e", "Sdha", "Mrpl34", "Cpt1a", "Cs"
  )
  ann <- fread(file.path(project_root, "results", "cross_dataset_OXPHOS", "OXPHOS_overlap_functional_annotation.csv"))
  found <- ann[match(wanted, gene_symbol)]
  if (anyNA(found$gene_symbol)) stop("A prespecified representative shared gene is absent.")
  found
}

make_figure_4 <- function() {
  x <- fread(file.path(project_root, "results", "final_mito_validation", "OXPHOS_gene_level_concordance.csv"))
  x <- x[present_in_GSE289089 & present_in_PRJNA1167170]
  labels <- representative_shared_15()
  stats <- fread(file.path(project_root, "results", "final_mito_validation", "OXPHOS_concordance_statistics.csv"))
  rho <- stats[metric == "spearman_correlation"]
  p_a <- ggplot(x, aes(GSE289089_ISC_logFC, PRJNA1167170_logFC)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "#555555") +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "#555555") +
    geom_point(size = 1.4, alpha = 0.65, color = "#666666") +
    geom_point(data = labels, aes(GSE289089_ISC_logFC, PRJNA1167170_logFC),
               inherit.aes = FALSE, size = 2, color = "#0072B2") +
    geom_text_repel(data = labels,
                    aes(GSE289089_ISC_logFC, PRJNA1167170_logFC, label = gene_symbol),
                    inherit.aes = FALSE, size = 2.2, max.overlaps = Inf,
                    box.padding = 0.25, min.segment.length = 0) +
    annotate("text", x = Inf, y = -Inf,
             label = sprintf("Spearman rho = %.3f\nP = %.3g", rho$value, rho$p_value),
             hjust = 1.05, vjust = -0.4, size = 3) +
    labs(x = "GSE289089 ISC logFC", y = "PRJNA1167170 colon logFC") + theme_manuscript()

  direction_source <- fread(file.path(project_root, "results", "PRJNA862187_validation",
                                      "three_dataset_OXPHOS_concordance.csv"))
  complete <- direction_source[measured_in_all_three == TRUE]
  if (nrow(complete) != 191L || sum(complete$negative_in_all_three) != 82L) {
    stop("Frozen three-dataset direction counts do not match the manuscript.")
  }
  direction <- data.table(
    category = factor(c("Negative in all three", "Other directional patterns"),
                      levels = c("Negative in all three", "Other directional patterns")),
    n = c(82L, 109L)
  )
  direction[, proportion := 100 * n / 191]
  p_b <- ggplot(direction, aes(category, proportion, fill = category)) +
    geom_col(width = 0.62, color = "black", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%d / 191\n%.2f%%", n, proportion)),
              vjust = -0.35, size = 3.2) +
    scale_fill_manual(values = c("Negative in all three" = "#0072B2",
                                 "Other directional patterns" = "#BDBDBD"), guide = "none") +
    scale_y_continuous(limits = c(0, 66), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Hallmark OXPHOS genes (%)",
         subtitle = "Directional consistency across three datasets; not effect-size equivalence") +
    theme_manuscript() +
    theme(axis.text.x = element_text(angle = 12, hjust = 1))

  p <- (p_a | p_b) + plot_annotation(tag_levels = "A", theme = tag_theme)
  save_pair_replace(p, "Figure_4", main_dir, 13, 5.3)
}

make_figure_5 <- function() {
  primary_groups <- cell_order[1:8]
  ranks <- fread(file.path(project_root, "results", "GSE289089_celltype_OXPHOS", "celltype_OXPHOS_all_rankings.csv"))
  ranks <- ranks[cell_group %in% primary_groups]
  ranks[, cell_group := factor(cell_group, levels = rev(primary_groups))]
  ranks[, ranking_method := factor(ranking_method, levels = c("logFC", "signed_stat"),
                                   labels = c("logFC ranking", "Signed-stat ranking"))]
  p_a <- ggplot(ranks, aes(NES, cell_group, group = cell_group, color = cell_group)) +
    geom_line(color = "#999999", linewidth = 0.6) +
    geom_point(aes(shape = ranking_method), size = 2.8, stroke = 0.8) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    scale_color_manual(values = cell_colors, guide = "none") +
    scale_shape_manual(values = c(16, 1), name = NULL) +
    scale_y_discrete(labels = cell_labels) +
    labs(x = "Hallmark OXPHOS NES", y = NULL) + theme_manuscript(8)

  ds <- fread(file.path(project_root, "results", "final_mito_validation", "downsampled_OXPHOS_iterations.csv"))
  ds <- ds[ranking_method == "logFC"]
  ds[, cell_group := factor(cell_group, levels = c("ISC", "TA", "Goblet", "Enterocyte", "Serotonergic EC", "Tuft"))]
  short_colors <- c(ISC = unname(cell_colors["ISC"]), TA = unname(cell_colors["TA_cycling_progenitor"]),
                    Goblet = unname(cell_colors["Goblet_secretory"]), Enterocyte = unname(cell_colors["Enterocyte"]),
                    `Serotonergic EC` = unname(cell_colors["Serotonergic_EC"]), Tuft = unname(cell_colors["Tuft"]))
  p_b <- ggplot(ds, aes(cell_group, NES, color = cell_group)) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    geom_boxplot(width = 0.58, outlier.shape = NA, fill = "white", color = "#555555", linewidth = 0.5) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = 20261004), size = 1.2, alpha = 0.75) +
    scale_color_manual(values = short_colors, guide = "none") +
    labs(x = NULL, y = "Hallmark OXPHOS NES", subtitle = "20 equal-cell downsampling iterations; logFC ranking") +
    theme_manuscript(8) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  qc <- fread(file.path(project_root, "results", "final_mito_validation", "celltype_sample_QC_summary.csv"))
  qc <- qc[cell_group %in% c("ISC", "TA", "Tuft", "Serotonergic EC")]
  qlong <- melt(qc,
                id.vars = c("sample", "cell_group", "condition", "n_cells"),
                measure.vars = c("median_nFeature_RNA", "median_nCount_RNA", "median_percent_mt"),
                variable.name = "metric", value.name = "median")
  qlong[, metric := factor(metric,
                           levels = c("median_nFeature_RNA", "median_nCount_RNA", "median_percent_mt"),
                           labels = c("Median nFeature_RNA", "Median nCount_RNA", "Median percent.mt"))]
  p_c <- ggplot(qlong, aes(cell_group, median, color = condition, shape = condition, group = sample)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.35)) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    scale_color_manual(values = condition_colors) + scale_shape_manual(values = condition_shapes) +
    labs(x = NULL, y = "Sample-level median", color = NULL, shape = NULL,
         subtitle = "Descriptive sample-level QC summaries; each point is one biological sample") +
    theme_manuscript(8) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  p <- (p_a | p_b) / p_c + plot_layout(heights = c(1, 0.88)) +
    plot_annotation(tag_levels = "A", theme = tag_theme)
  save_pair(p, "Figure_5", main_dir, 14, 10)
}

make_supplementary_figure_1 <- function() {
  pre <- readRDS(file.path(project_root, "processed", "GSE289089_combined_preQC.rds"))[[]]
  post <- readRDS(file.path(project_root, "processed", "GSE289089_combined_postQC.rds"))[[]]
  pre$stage <- "Before filtering"; post$stage <- "After filtering"
  x <- rbind(pre[, c("nFeature_RNA", "nCount_RNA", "percent.mt", "stage")],
             post[, c("nFeature_RNA", "nCount_RNA", "percent.mt", "stage")])
  x$cell <- seq_len(nrow(x))
  long <- melt(as.data.table(x), id.vars = c("cell", "stage"),
               measure.vars = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               variable.name = "metric", value.name = "value")
  long[, stage := factor(stage, levels = c("Before filtering", "After filtering"))]
  long[, metric := factor(metric, levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"))]
  thresholds <- data.table(
    metric = factor(c("nFeature_RNA", "nFeature_RNA", "percent.mt"),
                    levels = c("nFeature_RNA", "nCount_RNA", "percent.mt")),
    value = c(500, 6500, 25)
  )
  p <- ggplot(long, aes(stage, value, fill = stage)) +
    geom_violin(scale = "width", trim = TRUE, color = "#555555", linewidth = 0.35) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 2, fill = "white") +
    geom_hline(data = thresholds, aes(yintercept = value), inherit.aes = FALSE,
               linetype = "dashed", linewidth = 0.45, color = "#D55E00") +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("Before filtering" = "#BBBBBB", "After filtering" = "#0072B2"), guide = "none") +
    labs(x = NULL, y = "Cell-level value",
         subtitle = "Archived pre-QC and post-QC objects; thresholds were not re-estimated",
         caption = "Dashed lines: confirmed nFeature_RNA and percent.mt thresholds; no hard nCount_RNA cutoff") +
    theme_manuscript() + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    plot_annotation(tag_levels = "A", theme = tag_theme)
  save_pair(p, "Supplementary_Figure_1", supp_dir, 12, 4.5)
}

make_supplementary_figure_2 <- function() {
  p <- annotation_dotplot(9) + labs(title = "GSE289089 full reconstructed canonical-marker DotPlot")
  save_pair(p, "Supplementary_Figure_2", supp_dir, 15, 8.5)
}

make_supplementary_figure_4 <- function() {
  p <- composition_plot(TRUE) + labs(title = "Full sample-level cell-group composition")
  save_pair(p, "Supplementary_Figure_4", supp_dir, 10, 7)
}

make_supplementary_figure_5 <- function() {
  x <- fread(file.path(project_root, "results", "PRJNA862187_validation", "salmon_qc_summary.csv"))
  sample_order <- c("EC_SCD_1", "EC_SCD_2", "EC_SCD_3", "SD_SCD_1", "SD_SCD_2", "SD_SCD_3",
                    "EC_HFD_1", "EC_HFD_2", "EC_HFD_3", "SD_HFD_1", "SD_HFD_2", "SD_HFD_3")
  x[, sample := factor(sample, levels = sample_order)]
  x[, sleep := factor(sleep, levels = c("EC", "SD"))]
  p1 <- ggplot(x, aes(sample, mapping_rate, color = sleep, shape = sleep)) +
    geom_point(size = 2.7) +
    scale_color_manual(values = c(EC = condition_colors[["Control"]], SD = condition_colors[["SD"]])) +
    scale_shape_manual(values = c(EC = condition_shapes[["Control"]], SD = condition_shapes[["SD"]])) +
    labs(x = NULL, y = "Mapping rate (%)", color = "Sleep", shape = "Sleep") +
    theme_manuscript(8) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p2 <- ggplot(x, aes(sample, assigned_fragments / 1e6, color = sleep, shape = sleep)) +
    geom_point(size = 2.7) +
    scale_color_manual(values = c(EC = condition_colors[["Control"]], SD = condition_colors[["SD"]])) +
    scale_shape_manual(values = c(EC = condition_shapes[["Control"]], SD = condition_shapes[["SD"]])) +
    labs(x = NULL, y = "Assigned fragments (millions)", color = "Sleep", shape = "Sleep") +
    theme_manuscript(8) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  meta <- melt(x[, .(sample, Diet = diet, Sleep = as.character(sleep), `Library type` = library_type)],
               id.vars = "sample", variable.name = "metadata", value.name = "value")
  meta[, metadata := factor(metadata, levels = c("Diet", "Sleep", "Library type"))]
  meta[, fill_key := paste(metadata, value, sep = ":")]
  p3 <- ggplot(meta, aes(sample, metadata, fill = fill_key)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = value), size = 2.8) +
    scale_fill_manual(values = c("Diet:SCD" = "#E8E8E8", "Diet:HFD" = "#CC79A7",
                                 "Sleep:EC" = condition_colors[["Control"]],
                                 "Sleep:SD" = condition_colors[["SD"]],
                                 "Library type:ISR" = "#F0E442"), guide = "none") +
    labs(x = NULL, y = NULL) + theme_manuscript(8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks.y = element_blank())
  p <- (p1 | p2) / p3 + plot_layout(heights = c(1, 0.62), guides = "collect") +
    plot_annotation(tag_levels = "A", theme = tag_theme)
  save_pair_replace(p, "Supplementary_Figure_5", supp_dir, 13, 8)
}

make_supplementary_figure_6 <- function() {
  mito <- fread(file.path(project_root, "results", "GSE289089_celltype_OXPHOS", "celltype_mitochondrial_pathways.csv"))
  groups <- cell_order[1:11]
  p <- nes_heatmap(mito, groups, TRUE) +
    labs(title = "Full cell-group by mitochondrial-pathway NES matrix")
  save_pair(p, "Supplementary_Figure_6", supp_dir, 12, 8)
}

make_supplementary_figure_7 <- function() {
  pathways <- fread(file.path(project_root, "results", "PRJNA862187_validation", "interaction_pathways.csv"))
  pathways[, pathway := factor(pathway, levels = rev(pathway_order))]
  pathways[, result_label := sprintf("NES %.2f; FDR %.3g", NES, padj)]
  p1 <- ggplot(pathways, aes(NES, pathway)) +
    geom_segment(aes(x = 0, xend = NES, yend = pathway), color = "#777777", linewidth = 0.6) +
    geom_point(aes(fill = padj < 0.05), shape = 21, size = 3.1, color = "black") +
    geom_text(aes(label = result_label), hjust = ifelse(pathways$NES < 0, 1.05, -0.05), size = 2.5) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    scale_y_discrete(labels = pathway_labels) +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white"), name = "FDR < 0.05") +
    scale_x_continuous(limits = c(-3.05, 2.45)) +
    labs(x = "Interaction NES", y = NULL,
         subtitle = "Exploratory sleep-by-diet pathway interaction") + theme_manuscript(8)

  genes <- fread(file.path(project_root, "results", "PRJNA862187_validation", "edgeR_sleep_diet_interaction.csv"))
  if (nrow(genes) != 17900L || sum(genes$FDR < 0.05, na.rm = TRUE) != 0L) {
    stop("Frozen gene-level interaction results do not match the manuscript.")
  }
  p2 <- ggplot(genes, aes(logFC, -log10(FDR))) +
    geom_point(size = 0.75, alpha = 0.4, color = "#777777") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.45, color = "#D55E00") +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    annotate("text", x = Inf, y = Inf, label = "0 / 17,900 genes at FDR < 0.05",
             hjust = 1.05, vjust = 1.5, size = 3) +
    labs(x = "Interaction logFC", y = expression(-log[10](FDR)),
         subtitle = "Gene-level exploratory interaction") + theme_manuscript(8)
  p <- (p1 | p2) + plot_annotation(tag_levels = "A", theme = tag_theme)
  save_pair_replace(p, "Supplementary_Figure_7", supp_dir, 13, 5.8)
}

make_supplementary_figure_8 <- function() {
  x <- fread(file.path(project_root, "results", "final_mito_validation", "ISC_TA_regeneration_pathways.csv"))
  keep <- c("WNT signaling", "epithelial proliferation", "mTORC1 signaling", "translation",
            "ribosome biogenesis", "cell cycle", "DNA replication", "epithelial regeneration / wound repair")
  x <- x[family %in% keep]
  long1 <- x[, .(family, cell_group, ranking = "logFC ranking", NES = NES_logFC, FDR = FDR_logFC, status = status_logFC)]
  long2 <- x[, .(family, cell_group, ranking = "Signed-stat ranking", NES = NES_signedStat, FDR = FDR_signedStat, status = status_signedStat)]
  long <- rbind(long1, long2)
  long[, family := factor(family, levels = rev(keep))]
  long[, cell_group := factor(cell_group, levels = c("ISC", "TA"))]
  lim <- max(abs(long$NES), na.rm = TRUE)
  long[, label := ifelse(is.na(NES), "Not tested", sprintf("%.2f", NES))]
  p <- ggplot(long, aes(cell_group, family, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.7) +
    geom_point(data = long[!is.na(FDR) & FDR < 0.05], shape = 21, size = 1.8, fill = "black", color = "white",
               position = position_nudge(x = 0.32, y = 0.32)) +
    facet_wrap(~ranking, nrow = 1) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lim, lim), na.value = "#BDBDBD", name = "NES") +
    labs(x = NULL, y = NULL, caption = "Corner dot: FDR < 0.05; gray: not tested") +
    theme_manuscript() + theme(panel.grid = element_blank())
  save_pair(p, "Supplementary_Figure_8", supp_dir, 11, 6.5)
}

make_all_figures <- function() {
  make_figure_1()
  make_figure_2()
  make_figure_3()
  make_figure_4()
  make_figure_5()
  make_supplementary_figure_1()
  make_supplementary_figure_2()
  make_supplementary_figure_4()
  make_supplementary_figure_5()
  make_supplementary_figure_6()
  make_supplementary_figure_7()
  make_supplementary_figure_8()
}
