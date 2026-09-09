# ===================================================================================
# 01. allelic divergence of regulator classes in mesenchyme, endoderm, PNS (Fig.4a)
# ===================================================================================
library(ggplot2)
dir <- "/fraction_comparison"
tissues <- c("endo", "mesen", "PNS")
class_levels <- c("upstream", "intermediate", "downstream", "target")
tissue_cols <- c(endo = "#70c5c7", mesen = "#f9cf29", PNS = "#317ec3")
class_shapes <- c(upstream = 21, intermediate = 22, downstream = 24, target = 25)
outline_col <- "#4D4D4D"

format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("p<0.001")
  paste0("p=", signif(p, 2))
}

star_label <- function(p) {
  if (is.na(p) || p >= 0.05) return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  "*"
}

read_one <- function(tissue) {
  path <- file.path(dir, paste0(tissue, "_regulator_class_allelic_divergence_fraction_tissue_level.csv"))
  if (!file.exists(path)) {
    stop("Missing input CSV: ", path)
  }
  x <- read.csv(path, check.names = FALSE)
  required_cols <- c("OG", "regulator_class", "allelic_divergence_fraction")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }
  x$tissue <- tissue
  x
}

make_pairwise_stats <- function(plot_data) {
  pairs <- combn(class_levels, 2, simplify = FALSE)
  rows <- list()
  for (tissue in tissues) {
    sub <- plot_data[plot_data$tissue == tissue, , drop = FALSE]
    for (pair in pairs) {
      x <- sub$allelic_divergence_fraction[sub$regulator_class == pair[[1]]]
      y <- sub$allelic_divergence_fraction[sub$regulator_class == pair[[2]]]
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      p_value <- NA_real_
      statistic <- NA_real_
      if (length(x) > 0 && length(y) > 0) {
        wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
        p_value <- wt$p.value
        statistic <- unname(wt$statistic)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        tissue = tissue,
        group1 = pair[[1]],
        group2 = pair[[2]],
        comparison = paste(pair[[1]], pair[[2]], sep = "_vs_"),
        n_group1 = length(x),
        n_group2 = length(y),
        median_group1 = ifelse(length(x) > 0, median(x), NA_real_),
        median_group2 = ifelse(length(y) > 0, median(y), NA_real_),
        wilcox_w = statistic,
        p_value = p_value,
        p_label = format_p(p_value),
        significance = star_label(p_value),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

make_annotations <- function(stats) {
  sig <- stats[!is.na(stats$p_value) & stats$p_value < 0.05, , drop = FALSE]
  if (nrow(sig) == 0) return(sig)
  dodge_width <- 0.72
  offsets <- setNames(seq(-dodge_width / 3, dodge_width / 3, length.out = length(tissues)), tissues)
  sig$x <- match(sig$group1, class_levels) + offsets[sig$tissue]
  sig$xend <- match(sig$group2, class_levels) + offsets[sig$tissue]
  sig <- sig[order(sig$tissue, sig$x, sig$xend), , drop = FALSE]
  sig$rank <- seq_len(nrow(sig))
  sig$y <- 1.035 + (sig$rank - 1) * 0.06
  sig$y_tick <- sig$y - 0.035
  sig
}

all_data <- do.call(rbind, lapply(tissues, read_one))
plot_data <- all_data[
  all_data$regulator_class %in% class_levels &
    !is.na(all_data$allelic_divergence_fraction),
  ,
  drop = FALSE
]
plot_data$tissue <- factor(plot_data$tissue, levels = tissues)
plot_data$regulator_class <- factor(plot_data$regulator_class, levels = class_levels)

stats <- make_pairwise_stats(plot_data)
ann <- make_annotations(stats)
y_top <- if (nrow(ann) > 0) max(1.08, max(ann$y, na.rm = TRUE) + 0.08) else 1.08

count_df <- as.data.frame(table(plot_data$regulator_class, plot_data$tissue), stringsAsFactors = FALSE)
names(count_df) <- c("regulator_class", "tissue", "n")
count_line <- paste(
  apply(count_df, 1, function(row) paste0(row[["regulator_class"]], "-", row[["tissue"]], " n=", row[["n"]])),
  collapse = "   "
)

p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = regulator_class, y = allelic_divergence_fraction, fill = tissue)
) +
  ggplot2::geom_boxplot(
    width = 0.62,
    position = ggplot2::position_dodge(width = 0.72),
    outlier.shape = NA,
    linewidth = 0.35,
    alpha = 0.78,
    color = outline_col
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = regulator_class),
    position = ggplot2::position_jitterdodge(jitter.width = 0.08, dodge.width = 0.72),
    size = 1.45,
    stroke = 0.24,
    alpha = 0.72,
    color = outline_col
  ) +
  ggplot2::scale_fill_manual(values = tissue_cols, drop = FALSE) +
  ggplot2::scale_shape_manual(values = class_shapes, drop = FALSE) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(order = 1, nrow = 1),
    shape = ggplot2::guide_legend(order = 2, nrow = 1)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    expand = ggplot2::expansion(mult = c(0.03, 0.08))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, y_top), clip = "off") +
  ggplot2::labs(
    title = "Allelic divergence fraction by regulator class and tissue",
    subtitle = "Only significant within-tissue class comparisons are annotated (Wilcoxon p<0.05)",
    x = "regulator_class",
    y = "Allelic divergence fraction",
    fill = "tissue",
    shape = "regulator_class"
  ) +
  ggplot2::theme_classic(base_size = 9.6) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 10.8),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8.2),
    axis.text.x = ggplot2::element_text(angle = 24, hjust = 1, size = 8.4),
    axis.text.y = ggplot2::element_text(size = 8.4),
    axis.title = ggplot2::element_text(size = 9.4),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(size = 8.6),
    legend.text = ggplot2::element_text(size = 8.4),
    plot.margin = ggplot2::margin(8, 9, 8, 8)
  )

if (nrow(ann) > 0) {
  p <- p +
    ggplot2::geom_segment(
      data = ann,
      ggplot2::aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.28
    ) +
    ggplot2::geom_segment(
      data = ann,
      ggplot2::aes(x = x, xend = x, y = y_tick, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.28
    ) +
    ggplot2::geom_segment(
      data = ann,
      ggplot2::aes(x = xend, xend = xend, y = y_tick, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.28
    ) +
    ggplot2::geom_text(
      data = ann,
      ggplot2::aes(x = (x + xend) / 2, y = y + 0.03, label = significance),
      inherit.aes = FALSE,
      size = 3.0,
      fontface = "bold"
    )
}
out_pdf <- file.path(dir, "regulator_class_tissue_grouped_boxplot_tissue.pdf")
ggplot2::ggsave(out_pdf, p, width = 8, height = 5)



# ===================================================================================
# 02. Cis- and trans- regulatory decomposition of cardiopharyngeal GRN genes (Fig.4d)
# ===================================================================================
library(ggplot2)
library(ggrepel)
plot_data <- read.csv(
  "heart_7gene_list_transition_categories_with_cis_trans.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
plot_data$transition_categories <- factor(
  plot_data$transition_categories,
  levels = c(
    "switched",
    "propagated",
    "buffered",
    "attenuated",
    "exposed"
  )
)
transition_shapes <- c(
  "propagated" = 22,  
  "buffered" = 24 
)
class_colors <- c(
  "trans_only" = "#1391eb",
  "cis_dominant" = "#eeac4f",
  "trans_dominant" = "#b745ef"
)
p <- ggplot(
  plot_data,
  aes(
    x = log2_cis,
    y = log2_trans,
    fill = class_detailed,
    shape = transition_categories
  )
) +
  # y = x
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.6
  ) +
  # y = -x
  geom_abline(
    slope = -1,
    intercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.3
  ) +
  geom_point(
    aes(fill = class_detailed),
    size = 4,
    colour = "black",
    stroke = 0.8
  ) +
  geom_text_repel(
    aes(label = gene),
    color = "black",
    size = 6,
    fontface = "italic",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.linewidth = 0.4,
    max.overlaps = Inf
  ) +
  scale_shape_manual(
    values = transition_shapes,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = class_colors,
    drop = FALSE
  ) +
  labs(
    x = expression(log[2]~cis),
    y = expression(log[2]~trans),
    fill = "class_detailed",
    shape = "transition_categories"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 22,
        colour = "black",
        size = 5
      )
    )
  )+
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.3
    ),
    legend.position = "right"
  ) +
  coord_equal()
ggsave(
  filename = "heart_7gene_cis_trans_four_quadrant_plot.pdf",
  plot = p,
  width = 9,
  height = 7
)



# ========================================================================================
# 03. chromatin accessibility and expression level of Foxf in epidermis and heart (Fig.4e)
# ========================================================================================
### chromatin accessibility (coverage plot) ###
mtb.atac <- readRDS("CsCi_MTB_ATAC_anno.rds")
DefaultAssay(mtb.atac) <- "peaks"
# OG05402 - foxf - Cint.KY.KY2019:KY.Chr3.1005 Csav.dovetail-Aniseed-MSTRG.831
Idents(mtb.atac) <- "node_anno.predicted"
idents_to_plot <- c("heart_B_heart")
# idents_to_plot <- c("epidermis")

ann <- Annotation(mtb.atac)
Ci_gene <- "Cint.KY.KY2019:KY.Chr3.1005"
gene.gr.Ci <- ann[ann$gene_name == Ci_gene]
chrom.Ci <- as.character(seqnames(gene.gr.Ci)[1])    
start_pos.Ci <- min(start(gene.gr.Ci))                
end_pos.Ci   <- max(end(gene.gr.Ci))   
region.Ci <- paste0(chrom.Ci, "-", start_pos.Ci, "-", end_pos.Ci)
pdf("CsCi_MTB_foxf_Ci_heart_peak.pdf", width = 7, height = 6)
CoveragePlot(
  object = mtb.atac,
  region = region.Ci,
  idents = idents_to_plot,
  extend.upstream = 2000,
  extend.downstream = 2000,
  ymax = 30
)
dev.off()

Cs_gene <- "Csav.dovetail-Aniseed_MSTRG.831"
gene.gr.Cs <- ann[ann$gene_name == Cs_gene]
chrom.Cs <- as.character(seqnames(gene.gr.Cs)[1])    
start_pos.Cs <- min(start(gene.gr.Cs))                
end_pos.Cs   <- max(end(gene.gr.Cs))   
region.Cs <- paste0(chrom.Cs, "-", start_pos.Cs, "-", end_pos.Cs)
pdf("/CsCi_MTB_foxf_Cs_heart_peak.pdf", width = 7, height = 6)
CoveragePlot(
  object = mtb.atac,
  region = region.Cs,
  idents = idents_to_plot,
  extend.upstream = 2000,
  extend.downstream = 2000,
  ymax = 30
)
dev.off()



### expression level (violin plot) ###
library(Seurat)
library(SeuratObject)
library(ggplot2)
work_dir <- "/foxf_exp"
rds_path <- file.path(work_dir, "CsCi.MTB.RNA.allele.exp.anno.rds")
genes <- c(
  "Cint.KY.KY2019:KY.Chr3.1005",
  "Csav.dovetail-Aniseed-MSTRG.831"
)
gene_colors <- c(
  "Cint.KY.KY2019:KY.Chr3.1005" = "#e34f5b",
  "Csav.dovetail-Aniseed-MSTRG.831" = "#639cd5"
)
tissues <- c("heart", "epidermis")
out_pdf <- file.path(work_dir, "CsCi_MTB_RNA_allele_exp_violin.pdf")
out_summary <- file.path(work_dir, "CsCi_MTB_RNA_allele_exp_summary.csv")

get_metadata <- function(obj) {
  if (inherits(obj, "Seurat")) {
    return(obj[[]])
  }
  if (!is.null(obj@meta.data)) {
    return(obj@meta.data)
  }
  stop("Could not find Seurat-style metadata on the object.")
}

get_assay_names <- function(obj) {
  if (inherits(obj, "Seurat")) {
    return(names(obj@assays))
  }
  character(0)
}

get_layer <- function(obj, assay, layer) {
  assay_obj <- obj[[assay]]
  if ("layers" %in% slotNames(assay_obj) && layer %in% names(assay_obj@layers)) {
    if (exists("LayerData", where = asNamespace("SeuratObject"), mode = "function")) {
      return(SeuratObject::LayerData(assay_obj, layer = layer))
    }
  }
  if (layer %in% slotNames(assay_obj)) {
    return(slot(assay_obj, layer))
  }
  if (exists("GetAssayData", where = asNamespace("SeuratObject"), mode = "function")) {
    return(SeuratObject::GetAssayData(obj, assay = assay, layer = layer))
  }
  stop(sprintf("Could not retrieve %s layer for assay %s.", layer, assay))
}

get_data_matrix <- function(obj, assay) {
  get_layer(obj, assay, "data")
}

obj <- readRDS(rds_path)
if (!inherits(obj, "Seurat")) {
  stop(sprintf("Unexpected object class: %s", paste(class(obj), collapse = ", ")))
}
meta <- get_metadata(obj)
if (!"tissue.plot.big" %in% colnames(meta)) {
  stop("metadata column tissue.plot.big was not found.")
}

assays <- get_assay_names(obj)
candidate_assays <- unique(c(SeuratObject::DefaultAssay(obj), "RNA", assays))
candidate_assays <- candidate_assays[candidate_assays %in% assays]

selected_assay <- NULL
selected_mat <- NULL
selected_counts <- NULL
presence_by_assay <- data.frame()

for (assay in candidate_assays) {
  mat <- tryCatch(get_data_matrix(obj, assay), error = function(e) e)
  if (inherits(mat, "error")) {
    message("Skipping assay ", assay, ": ", conditionMessage(mat))
    next
  }
  present <- genes %in% rownames(mat)
  presence_by_assay <- rbind(
    presence_by_assay,
    data.frame(assay = assay, gene = genes, present = present)
  )
  if (all(present) && is.null(selected_assay)) {
    selected_assay <- assay
    selected_mat <- mat
    selected_counts <- get_layer(obj, assay, "counts")
  }
}

if (nrow(presence_by_assay) == 0) {
  stop("No readable data layer was found in any assay.")
}

if (is.null(selected_assay)) {
  missing_summary <- subset(presence_by_assay, !present)
  stop(sprintf(
    "The requested genes were not both found in the same data layer. Missing checks:\n%s",
    paste(capture.output(print(missing_summary, row.names = FALSE)), collapse = "\n")
  ))
}

cells <- rownames(meta)[meta$tissue.plot.big %in% tissues]
if (length(cells) == 0) {
  stop("No cells found for requested tissues in tissue.plot.big.")
}
cells <- intersect(cells, colnames(selected_mat))
if (length(cells) == 0) {
  stop("Requested tissue cells were not present in the expression matrix.")
}

expr <- as.matrix(selected_mat[genes, cells, drop = FALSE])
counts_expr <- selected_counts[genes, cells, drop = FALSE]
plot_df <- data.frame(
  expression_level = as.vector(expr),
  tissue = rep(meta[cells, "tissue.plot.big"], each = length(genes)),
  gene = rep(genes, times = length(cells)),
  stringsAsFactors = FALSE
)
plot_df$tissue <- factor(plot_df$tissue, levels = tissues)
plot_df$gene <- factor(plot_df$gene, levels = genes)

summary_df <- do.call(rbind, lapply(tissues, function(tissue_name) {
  tissue_cells <- intersect(
    rownames(meta)[meta$tissue.plot.big == tissue_name],
    colnames(selected_mat)
  )
  do.call(rbind, lapply(genes, function(gene_name) {
    values <- as.numeric(selected_mat[gene_name, tissue_cells, drop = TRUE])
    counts_values <- as.numeric(selected_counts[gene_name, tissue_cells, drop = TRUE])
    data.frame(
      tissue = tissue_name,
      gene = gene_name,
      cells = length(tissue_cells),
      raw_count_sum = sum(counts_values),
      data_layer_sum = sum(values),
      data_layer_mean = mean(values),
      data_layer_median = median(values),
      data_layer_q1 = unname(quantile(values, 0.25)),
      data_layer_q3 = unname(quantile(values, 0.75)),
      data_layer_max = max(values),
      expressing_cells_counts_gt0 = sum(counts_values > 0),
      expressing_cells_data_gt0 = sum(values > 0),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(summary_df, out_summary, row.names = FALSE, quote = TRUE)

p <- ggplot(plot_df, aes(x = expression_level, y = tissue, fill = gene)) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    width = 0.82,
    position = position_dodge(width = 0.95)
  ) +
  scale_fill_manual(values = gene_colors, drop = FALSE) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(x = "expression level", y = "tissue", fill = "gene") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black")
  )
ggsave(out_pdf, p, width = 8, height = 5)
