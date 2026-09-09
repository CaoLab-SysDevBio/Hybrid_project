# ==============================================================================
# 01. transition categories definition and lineage tree plot (Fig.2a)
# We used nervous system lineages as examples.
# ==============================================================================
base_dir <- "/CNS/allelic_rates"
default_lineages <- c("CNS_A8.7_8_15_16", "CNS_a8.17_20_25_26", "CNS_b8.17_20")
args <- commandArgs(trailingOnly = TRUE)
lineages <- if (length(args) > 0) args else default_lineages

bias_levels <- c("maternal", "balanced", "paternal")
transition_levels <- c("bi_bi_switch", "bi_bi_stable", "bi_ba", "ba_ba", "ba_bi") # switched, propagated, attenuated, buffered, exposed
transition_alpha <- c(
  bi_bi_switch = 0.80,
  bi_bi_stable = 0.10,
  bi_ba = 0.60,
  ba_ba = 0.30,
  ba_bi = 1.00
)
stage_colors <- c(
  "a.iniG" = "#952b47",
  "b.midG" = "#c64550",
  "c.earN" = "#e16c45",
  "d.latN" = "#eca764",
  "e.ETB" = "#a6cb9d",
  "f.MTB" = "#68b59c",
  "g.LTB1" = "#3981b0",
  "h.LTB2" = "#b07caa",
  "i.larva" = "#5f5a9b"
)
pie_alpha <- c(maternal = 1.00, balanced = 0.75, paternal = 0.50)
boundary_col <- "#58595b"
alpha_colors <- function(base_col, alphas) {
  cols <- vapply(alphas, function(alpha) {
    grDevices::adjustcolor(base_col, alpha.f = alpha)
  }, character(1))
  names(cols) <- names(alphas)
  cols
}

normalize_node_key <- function(x) {
  x <- as.character(x)
  sub("^i\\.lv_", "i.larva_", x)
}

stage_key <- function(node) {
  node <- as.character(node)
  hits <- names(stage_colors)[vapply(names(stage_colors), function(prefix) {
    startsWith(node, prefix)
  }, logical(1))]
  if (length(hits) == 0) return(names(stage_colors)[1])
  hits[which.max(nchar(hits))]
}

node_color <- function(node) {
  stage_colors[[stage_key(node)]]
}

clean_gene_vector <- function(x) {
  x <- unique(as.character(unlist(x, use.names = FALSE)))
  x[!is.na(x) & x != ""]
}

classify_bias <- function(rate) {
  out <- rep(NA_character_, length(rate))
  out[!is.na(rate) & rate > 0.55] <- "maternal"
  out[!is.na(rate) & rate < 0.45] <- "paternal"
  out[!is.na(rate) & rate >= 0.45 & rate <= 0.55] <- "balanced"
  factor(out, levels = bias_levels)
}

classify_transition <- function(parent_bias, child_bias) {
  out <- rep(NA_character_, length(parent_bias))
  parent_biased <- parent_bias %in% c("maternal", "paternal")
  child_biased <- child_bias %in% c("maternal", "paternal")
  parent_balanced <- parent_bias == "balanced"
  child_balanced <- child_bias == "balanced"

  out[parent_biased & child_balanced] <- "bi_ba" #attenuated
  out[parent_balanced & child_balanced] <- "ba_ba" #buffered
  out[parent_balanced & child_biased] <- "ba_bi" #exposed
  out[parent_biased & child_biased & parent_bias == child_bias] <- "bi_bi_stable" #propagated
  out[parent_biased & child_biased & parent_bias != child_bias] <- "bi_bi_switch" #switched
  factor(out, levels = transition_levels)
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

sector_points <- function(start, end, r_outer, r_inner = 0, n = 120) {
  if (end < start) {
    tmp <- start
    start <- end
    end <- tmp
  }
  theta <- seq(start, end, length.out = max(8, ceiling(n * abs(end - start) / (2 * pi))))
  if (r_inner <= 0) {
    list(x = c(0, r_outer * cos(theta), 0), y = c(0, r_outer * sin(theta), 0))
  } else {
    theta_inner <- rev(theta)
    list(
      x = c(r_outer * cos(theta), r_inner * cos(theta_inner)),
      y = c(r_outer * sin(theta), r_inner * sin(theta_inner))
    )
  }
}

draw_pie <- function(counts, colors, radius = 0.58, center = c(0, 0), start = pi / 2) {
  counts <- counts[counts > 0]
  if (length(counts) == 0) return(invisible(NULL))
  colors <- colors[names(counts)]
  props <- counts / sum(counts)
  angles <- c(start, start - cumsum(props) * 2 * pi)
  for (i in seq_along(counts)) {
    pts <- sector_points(angles[i + 1], angles[i], radius, 0)
    polygon(center[1] + pts$x, center[2] + pts$y, col = colors[i], border = NA)
  }
  invisible(NULL)
}

draw_donut <- function(counts, colors, r_inner = 0.66, r_outer = 0.98, center = c(0, 0), start = pi / 2, add_radial_boundaries = TRUE) {
  counts <- counts[counts > 0]
  if (length(counts) == 0) return(invisible(NULL))
  colors <- colors[names(counts)]
  props <- counts / sum(counts)
  angles <- c(start, start - cumsum(props) * 2 * pi)
  for (i in seq_along(counts)) {
    pts <- sector_points(angles[i + 1], angles[i], r_outer, r_inner)
    polygon(center[1] + pts$x, center[2] + pts$y, col = colors[i], border = NA)
  }
  if (add_radial_boundaries && length(counts) > 1) {
    internal <- angles[-c(1, length(angles))]
    for (theta in internal) {
      lines(
        center[1] + c(r_inner, r_outer) * cos(theta),
        center[2] + c(r_inner, r_outer) * sin(theta),
        col = boundary_col,
        lwd = 0.6,
        lend = "butt"
      )
    }
  }
  invisible(NULL)
}

plot_node_chart <- function(out_pdf, node, bias_counts, transition_counts = NULL) {
  base_col <- node_color(node)
  bias_cols <- alpha_colors(base_col, pie_alpha[bias_levels])
  transition_cols <- alpha_colors(base_col, transition_alpha[transition_levels])

  grDevices::pdf(out_pdf, width = 2.2, height = 2.2, useDingbats = FALSE)
  op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), asp = 1)
  if (is.null(transition_counts) || sum(transition_counts) == 0) {
    draw_pie(bias_counts, bias_cols, radius = 0.92)
  } else {
    draw_donut(transition_counts[transition_levels], transition_cols, r_inner = 0.66, r_outer = 0.98)
    draw_pie(bias_counts[bias_levels], bias_cols, radius = 0.56)
  }
  par(op)
  grDevices::dev.off()
}

compute_tree_layout <- function(nodes, parents) {
  parent_map <- parents
  names(parent_map) <- nodes
  root <- nodes[is.na(parent_map) | parent_map == ""]
  if (length(root) != 1) stop("Expected exactly one root, found: ", paste(root, collapse = ", "))
  children <- split(nodes[parent_map != ""], parent_map[parent_map != ""])
  leaf_order <- character(0)
  depth <- setNames(rep(NA_integer_, length(nodes)), nodes)

  walk_depth <- function(node, d) {
    depth[[node]] <<- d
    kids <- children[[node]]
    if (is.null(kids) || length(kids) == 0) {
      leaf_order <<- c(leaf_order, node)
    } else {
      for (kid in kids) walk_depth(kid, d + 1)
    }
  }
  walk_depth(root, 0)

  x <- setNames(rep(NA_real_, length(nodes)), nodes)
  x[leaf_order] <- seq_along(leaf_order)
  assign_x <- function(node) {
    kids <- children[[node]]
    if (is.null(kids) || length(kids) == 0) return(x[[node]])
    kid_x <- vapply(kids, assign_x, numeric(1))
    x[[node]] <<- mean(kid_x)
    x[[node]]
  }
  assign_x(root)
  data.frame(node = nodes, parent_node = parent_map, x = x[nodes], y = -depth[nodes], stringsAsFactors = FALSE)
}

plot_lineage_tree <- function(out_pdf, rel, node_bias_counts, node_transition_counts) {
  layout <- compute_tree_layout(rel$node, rel$parent_node)
  max_depth <- max(-layout$y)
  leaf_n <- sum(!(layout$node %in% rel$parent_node[rel$parent_node != ""]))
  width <- max(7, min(18, leaf_n * 0.55 + 1.2))
  height <- max(7, min(22, (max_depth + 1) * 0.78 + 1.0))
  node_radius <- min(0.24, max(0.12, 0.18 + 0.02 * (12 / max(leaf_n, 12))))

  grDevices::pdf(out_pdf, width = width, height = height, useDingbats = FALSE)
  op <- par(mar = c(0.15, 0.15, 0.15, 0.15), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(
    xlim = range(layout$x) + c(-0.7, 0.7),
    ylim = range(layout$y) + c(-0.5, 0.5),
    asp = 1
  )

  for (parent in unique(rel$parent_node[rel$parent_node != ""])) {
    kids <- rel$node[rel$parent_node == parent]
    parent_row <- layout[layout$node == parent, , drop = FALSE]
    kid_rows <- layout[match(kids, layout$node), , drop = FALSE]
    if (nrow(kid_rows) == 1) {
      lines(c(parent_row$x, kid_rows$x), c(parent_row$y - node_radius, kid_rows$y + node_radius), col = "#6d6e71", lwd = 0.6)
    } else {
      mid_y <- parent_row$y - 0.48
      lines(c(parent_row$x, parent_row$x), c(parent_row$y - node_radius, mid_y), col = "#6d6e71", lwd = 0.6)
      lines(range(kid_rows$x), c(mid_y, mid_y), col = "#6d6e71", lwd = 0.6)
      for (i in seq_len(nrow(kid_rows))) {
        lines(c(kid_rows$x[i], kid_rows$x[i]), c(mid_y, kid_rows$y[i] + node_radius), col = "#6d6e71", lwd = 0.6)
      }
    }
  }

  for (i in seq_len(nrow(layout))) {
    node <- layout$node[i]
    base_col <- node_color(node)
    bias_cols <- alpha_colors(base_col, pie_alpha[bias_levels])
    transition_cols <- alpha_colors(base_col, transition_alpha[transition_levels])
    center <- c(layout$x[i], layout$y[i])
    transition_counts <- node_transition_counts[[node]]
    if (!is.null(transition_counts) && sum(transition_counts) > 0) {
      draw_donut(transition_counts[transition_levels], transition_cols, r_inner = node_radius * 0.66, r_outer = node_radius, center = center)
      draw_pie(node_bias_counts[[node]][bias_levels], bias_cols, radius = node_radius * 0.56, center = center)
    } else {
      draw_pie(node_bias_counts[[node]][bias_levels], bias_cols, radius = node_radius * 0.92, center = center)
    }
  }
  par(op)
  grDevices::dev.off()
}

read_inputs <- function(lineage) {
  avg_path <- file.path(base_dir, paste0(lineage, "_cutoff30_joint_avg_psmean_10stage_nodes.csv"))
  rds_path <- file.path(base_dir, paste0(lineage, "_gene_totalCount_cutoff30_byNode.rds"))
  rel_path <- file.path(base_dir, paste0("mother_daughter_relationship_", lineage, ".csv"))
  og_path <- file.path(base_dir, "OG_genelist_without_mlt_data2017_blast.csv")

  avg <- read.csv(avg_path, check.names = FALSE, stringsAsFactors = FALSE)
  names(avg)[1] <- "OG"
  avg$OG <- as.character(avg$OG)
  rel <- read.csv(rel_path, check.names = FALSE, stringsAsFactors = FALSE)
  names(rel)[1:2] <- c("node", "parent_node")
  rel$node <- as.character(rel$node)
  rel$parent_node <- as.character(rel$parent_node)
  rel$parent_node[is.na(rel$parent_node)] <- ""
  rds_obj <- readRDS(rds_path)
  names(rds_obj) <- normalize_node_key(names(rds_obj))
  og <- read.csv(og_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Ci_OG" %in% names(og)) stop("Ci_OG column not found in ", og_path)
  og_genes <- clean_gene_vector(og$Ci_OG)

  list(avg = avg, rel = rel, rds_obj = rds_obj, og_genes = og_genes)
}

run_lineage <- function(lineage) {
  message("Processing ", lineage)
  inputs <- read_inputs(lineage)
  avg <- inputs$avg
  rel <- inputs$rel
  rds_obj <- inputs$rds_obj
  og_genes <- inputs$og_genes
  plot_dir <- file.path(output_base_dir, paste0(lineage, "_plot"))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  old_node_pdfs <- list.files(
    plot_dir,
    pattern = paste0("^", regex_escape(lineage), "_.*_same_color_custom_alpha_border_donut[.]pdf$"),
    full.names = TRUE
  )
  if (length(old_node_pdfs) > 0) unlink(old_node_pdfs)
  old_summaries <- file.path(
    plot_dir,
    c(
      paste0(lineage, "_cutoff30_same_color_custom_alpha_border_donut_bias_summary.csv"),
      paste0(lineage, "_cutoff30_same_color_custom_alpha_border_donut_transition_summary.csv")
    )
  )
  unlink(old_summaries[file.exists(old_summaries)])

  missing_avg <- setdiff(rel$node, names(avg))
  missing_rds <- setdiff(rel$node, names(rds_obj))
  if (length(missing_avg) > 0) stop(lineage, " nodes missing in avg table: ", paste(missing_avg, collapse = ", "))
  if (length(missing_rds) > 0) stop(lineage, " nodes missing in RDS: ", paste(missing_rds, collapse = ", "))

  avg_gene_set <- clean_gene_vector(avg$OG)
  avg_index <- match(avg$OG, avg$OG)
  names(avg_index) <- avg$OG

  node_genes <- list()
  node_rates <- list()
  node_bias <- list()
  node_bias_counts <- list()
  node_gene_counts <- list()

  for (node in rel$node) {
    genes <- Reduce(intersect, list(og_genes, avg_gene_set, clean_gene_vector(rds_obj[[node]])))
    idx <- match(genes, avg$OG)
    rates <- suppressWarnings(as.numeric(avg[idx, node]))
    valid <- !is.na(rates)
    genes <- genes[valid]
    rates <- rates[valid]
    bias <- as.character(classify_bias(rates))
    names(rates) <- genes
    names(bias) <- genes
    counts <- table(factor(bias, levels = bias_levels))

    node_genes[[node]] <- genes
    node_rates[[node]] <- rates
    node_bias[[node]] <- bias
    node_bias_counts[[node]] <- counts
    node_gene_counts[[node]] <- data.frame(
      lineage = lineage,
      node = node,
      intersect_gene_n = length(genes),
      maternal_n = unname(counts[["maternal"]]),
      balanced_n = unname(counts[["balanced"]]),
      paternal_n = unname(counts[["paternal"]]),
      stringsAsFactors = FALSE
    )
  }

  node_transition_counts <- list()
  node_transition_gene_lists <- list()
  transition_summary <- list()

  for (node in rel$node) {
    parent <- rel$parent_node[rel$node == node][1]
    empty_lists <- setNames(vector("list", length(transition_levels)), transition_levels)
    if (is.na(parent) || parent == "") {
      node_transition_counts[[node]] <- table(factor(character(0), levels = transition_levels))
      node_transition_gene_lists[[node]] <- list(parent_node = parent, shared_gene_n = 0, categories = empty_lists)
      next
    }

    shared <- intersect(node_genes[[parent]], node_genes[[node]])
    parent_bias <- node_bias[[parent]][shared]
    child_bias <- node_bias[[node]][shared]
    transition <- as.character(classify_transition(parent_bias, child_bias))
    names(transition) <- shared
    valid <- !is.na(transition)
    transition <- transition[valid]
    counts <- table(factor(transition, levels = transition_levels))
    gene_lists <- lapply(transition_levels, function(cat) names(transition)[transition == cat])
    names(gene_lists) <- transition_levels

    node_transition_counts[[node]] <- counts
    node_transition_gene_lists[[node]] <- list(parent_node = parent, shared_gene_n = length(transition), categories = gene_lists)
    for (cat in transition_levels) {
      transition_summary[[length(transition_summary) + 1]] <- data.frame(
        lineage = lineage,
        node = node,
        parent_node = parent,
        category = cat,
        n = length(gene_lists[[cat]]),
        stringsAsFactors = FALSE
      )
    }
  }

  for (node in rel$node) {
    out_pdf <- file.path(plot_dir, paste0(lineage, "_", safe_name(node), "_same_color_custom_alpha_border_donut.pdf"))
    plot_node_chart(out_pdf, node, node_bias_counts[[node]], node_transition_counts[[node]])
  }

  tree_pdf <- file.path(output_base_dir, paste0("lineage_tree_plot_", lineage, "_same_color_custom_alpha_border_donut.pdf"))
  plot_lineage_tree(tree_pdf, rel, node_bias_counts, node_transition_counts)

  transition_summary_df <- if (length(transition_summary) > 0) {
    do.call(rbind, transition_summary)
  } else {
    data.frame(lineage = character(0), node = character(0), parent_node = character(0), category = character(0), n = integer(0))
  }
  bias_summary_df <- do.call(rbind, node_gene_counts)

  write.csv(
    bias_summary_df,
    file.path(plot_dir, paste0(lineage, "_cutoff30_same_color_custom_alpha_border_donut_bias_summary.csv")),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    transition_summary_df,
    file.path(plot_dir, paste0(lineage, "_cutoff30_same_color_custom_alpha_border_donut_transition_summary.csv")),
    row.names = FALSE,
    quote = TRUE
  )

  transition_rds <- list(
    metadata = list(
      lineage = lineage,
      base_dir = base_dir,
      generated_at = as.character(Sys.time()),
      bias_thresholds = c(paternal_lt = 0.45, balanced_min = 0.45, balanced_max = 0.55, maternal_gt = 0.55),
      transition_order = transition_levels,
      transition_alpha = transition_alpha,
      note = "Each node stores parent_node, shared_gene_n, and transition category gene lists calculated from parent-child shared intersected OG genes."
    ),
    node_transition_gene_lists = node_transition_gene_lists,
    count_summary = transition_summary_df,
    node_bias_summary = bias_summary_df
  )
  rds_out <- file.path(output_base_dir, paste0(lineage, "_cutoff30_same_color_custom_alpha_border_donut_transition_gene_lists.rds"))
  saveRDS(transition_rds, rds_out)

  data.frame(
    lineage = lineage,
    nodes = nrow(rel),
    node_pdf_n = length(list.files(plot_dir, pattern = "_same_color_custom_alpha_border_donut[.]pdf$", full.names = FALSE)),
    tree_pdf = tree_pdf,
    transition_rds = rds_out,
    plot_dir = plot_dir,
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(lineages, run_lineage))
summary_path <- if (setequal(lineages, default_lineages) && length(lineages) == length(default_lineages)) {
  file.path(base_dir, "CNS_cutoff30_same_color_custom_alpha_border_donut_run_summary.csv")
} else {
  file.path(output_base_dir, paste0(paste(lineages, collapse = "_"), "_cutoff30_same_color_custom_alpha_border_donut_rerun_summary.csv"))
}
write.csv(results, summary_path, row.names = FALSE, quote = TRUE)



# ============================================================================
# 02. transition categories statistic and bias persistence (Fig.2b, 2c)
# ============================================================================
work_dir <- "/all_node_transition_mode_statistic"
input_csv <- file.path(work_dir, "transition_category_node_percentages_from_transition_gene_lists.csv")
transition_pdf <- file.path(work_dir, "transition_category_boxplot_all_categories.pdf")
bias_pdf <- file.path(work_dir, "bias_persistence_by_lineage.pdf")
library(ggplot2)
lineage_colors <- c(
  "CNS_a8.17_20_25_26" = "#317ec3",
  "CNS_A8.7_8_15_16" = "#0472f8",
  "CNS_b8.17_20" = "#6695c4",
  "endoderm_A7.1_2_5" = "#17dedb",
  "endoderm_B7.1_2" = "#70c5c7",
  "epidermis_a8.21_24_27_32" = "#6658ff",
  "epidermis_b8.21_32" = "#605a9a",
  "germ_B7.6" = "#386236",
  "heart_B7.5" = "#ec5624",
  "mesenchyme_A7.6" = "#e7cf21",
  "mesenchyme_B7.7_8.5" = "#a08b11",
  "muscle_B8.7_8_15_16" = "#ff9500",
  "notochord_A8.5_6_13_14" = "#4bc051",
  "notochord_B8.6" = "#029115"
)
lineage_order <- names(lineage_colors)
exclude_nodes <- c(
  "i.larva_33",
  "i.larva_21",
  "i.larva_24",
  "i.larva_27",
  "i.larva_35",
  "i.larva_10",
  "a.iniG_b_line",
  "i.larva_Tll1+ mesenchyme",
  "i.larva_noto2"
) # less than 6 cells

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
plot_data <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("lineage", "node", "category_label", "n", "percent")
missing_cols <- setdiff(required_cols, names(plot_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

original_rows <- nrow(plot_data)
matched_excluded_nodes <- sort(unique(plot_data$node[plot_data$node %in% exclude_nodes]))
missing_excluded_nodes <- setdiff(exclude_nodes, matched_excluded_nodes)
plot_data <- plot_data[!(plot_data$node %in% exclude_nodes), ]
removed_rows <- original_rows - nrow(plot_data)

plot_data <- plot_data[plot_data$lineage %in% lineage_order, ]
plot_data$n <- as.numeric(plot_data$n)
plot_data$percent <- as.numeric(plot_data$percent)

### Transition category boxplot ###
category_top_to_bottom <- c("propagated", "switched", "attenuated", "buffered", "exposed")
low_max <- 10
high_min <- 80
gap_width <- 4
segment_width <- 10
x_break_transform <- function(x) {
  ifelse(
    x <= low_max,
    x,
    ifelse(
      x >= high_min,
      low_max + gap_width + (x - high_min) * segment_width / (100 - high_min),
      low_max + gap_width * (x - low_max) / (high_min - low_max)
    )
  )
}

transition_data <- plot_data[plot_data$category_label %in% category_top_to_bottom & !is.na(plot_data$percent), ]
box_rows <- list()
i <- 0L
for (category in category_top_to_bottom) {
  for (lineage in lineage_order) {
    values <- transition_data$percent[transition_data$category_label == category & transition_data$lineage == lineage]
    values <- values[!is.na(values)]
    if (length(values) == 0) next
    stats <- boxplot.stats(values, coef = 1.5, do.conf = FALSE, do.out = FALSE)$stats
    i <- i + 1L
    box_rows[[i]] <- data.frame(
      category_label = category,
      lineage = lineage,
      whisker_low = stats[1],
      q1 = stats[2],
      median = stats[3],
      q3 = stats[4],
      whisker_high = stats[5],
      n = length(values),
      stringsAsFactors = FALSE
    )
  }
}
box_df <- do.call(rbind, box_rows)
box_df$category_index <- match(box_df$category_label, rev(category_top_to_bottom))
box_df$lineage_index <- match(box_df$lineage, lineage_order)
lineage_n <- length(lineage_order)
band_height <- 0.82
slot_height <- band_height / lineage_n
box_df$y <- box_df$category_index - band_height / 2 + (box_df$lineage_index - 0.5) * slot_height
box_df$box_half_height <- slot_height * 0.36
for (col in c("whisker_low", "q1", "median", "q3", "whisker_high")) {
  box_df[[paste0("x_", col)]] <- x_break_transform(box_df[[col]])
}
axis_breaks <- c(0, 2.5, 5, 7.5, 10, 80, 85, 90, 95, 100)
axis_breaks_t <- x_break_transform(axis_breaks)
x_max_t <- x_break_transform(100)

p_transition <- ggplot(box_df) +
  geom_segment(aes(x = x_whisker_low, xend = x_q1, y = y, yend = y, color = lineage), linewidth = 0.34) +
  geom_segment(aes(x = x_q3, xend = x_whisker_high, y = y, yend = y, color = lineage), linewidth = 0.34) +
  geom_rect(
    aes(xmin = x_q1, xmax = x_q3, ymin = y - box_half_height, ymax = y + box_half_height, fill = lineage),
    color = "grey20",
    linewidth = 0.24,
    alpha = 0.96
  ) +
  geom_segment(aes(x = x_median, xend = x_median, y = y - box_half_height, yend = y + box_half_height), color = "grey10", linewidth = 0.28) +
  annotate("text", x = low_max + gap_width / 2, y = 0.31, label = "//", size = 3.2, color = "grey20") +
  scale_fill_manual(values = lineage_colors, drop = FALSE) +
  scale_color_manual(values = lineage_colors, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, x_max_t),
    breaks = axis_breaks_t,
    labels = paste0(axis_breaks, "%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0.25, 5.55),
    breaks = seq_along(rev(category_top_to_bottom)),
    labels = rev(category_top_to_bottom),
    expand = c(0, 0)
  ) +
  labs(x = "Percentage", y = "Category", fill = "Lineage") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 6.2),
    legend.title = element_text(size = 7),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    plot.margin = margin(8, 10, 8, 8)
  ) +
  guides(
    color = "none",
    fill = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(alpha = 0.96))
  )

ggsave(transition_pdf, p_transition, width = 8.6, height = 7.1, units = "in", device = cairo_pdf)

### Bias persistence boxplot ###
target_categories <- c("propagated", "attenuated", "switched")
bias_sub <- plot_data[plot_data$category_label %in% target_categories, c("lineage", "node", "category_label", "n")]
wide <- reshape(
  bias_sub,
  idvar = c("lineage", "node"),
  timevar = "category_label",
  direction = "wide"
)
count_cols <- paste0("n.", target_categories)
for (col in count_cols) {
  if (!col %in% names(wide)) wide[[col]] <- 0
  wide[[col]][is.na(wide[[col]])] <- 0
}
wide$denominator <- wide[["n.propagated"]] + wide[["n.attenuated"]] + wide[["n.switched"]]
wide$bias_persistence <- ifelse(
  wide$denominator > 0,
  100 * wide[["n.propagated"]] / wide$denominator,
  NA_real_
)
bias_df <- wide[!is.na(wide$bias_persistence), c("lineage", "node", "bias_persistence")]
lineage_medians <- aggregate(bias_persistence ~ lineage, data = bias_df, FUN = median)
lineage_medians <- lineage_medians[order(lineage_medians$bias_persistence, decreasing = FALSE), ]
bias_df$lineage <- factor(bias_df$lineage, levels = rev(lineage_medians$lineage))
pdf_height <- max(4.8, min(10, 0.31 * length(unique(bias_df$lineage)) + 1.8))

p_bias <- ggplot(bias_df, aes(x = bias_persistence, y = lineage, color = lineage)) +
  geom_boxplot(
    fill = "white",
    width = 0.58,
    outlier.shape = NA,
    linewidth = 0.42
  ) +
  scale_color_manual(values = lineage_colors, drop = FALSE) +
  scale_y_discrete(labels = capitalize_first) +
  scale_x_continuous(
    breaks = seq(80, 100, 5),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(xlim = c(80, 100)) +
  labs(x = "Bias persistence (%)", y = "Lineage") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 8),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(bias_pdf, p_bias, width = 6.2, height = pdf_height, units = "in", device = cairo_pdf)




# ==============================================================================
# 03. Cp26+ endoderm lineage (Fig.2d)
# ==============================================================================
library(dplyr)
library(tidyverse) 
library(igraph)
library(ggraph)
library(ggplot2)
library(scales)
pattern_df <- read.csv("/endo_cp26_avg_psmean_10stage_MP_nmlt.csv")
df <- pattern_df[, c("d.latN","f.MTB","h.LTB2","j.larva24")]
colnames(df) <- c("latN","MTB","LTB2","larva24")
df$latN <- paste0("latN_", df$latN)
df$MTB <- paste0("MTB_", df$MTB)
df$LTB2 <- paste0("LTB2_", df$LTB2)
df$larva24 <- paste0("larva24_", df$larva24)
# Build lineage nodes
df$node1 <- df$latN
df$node2 <- paste(df$latN, df$MTB, sep = "|")
df$node3 <- paste(df$latN, df$MTB, df$LTB2, sep = "|")
df$node4 <- paste(df$latN, df$MTB, df$LTB2, df$larva24, sep = "|")
# Create edges
edge1 <- df %>%
  count(node1, node2) %>%
  group_by(node1) %>%
  mutate(freq = n, weight = n / sum(n)) %>%
  ungroup() %>%
  rename(from = node1, to = node2) %>%
  select(from, to, weight, freq)

edge2 <- df %>%
  count(node2, node3) %>%
  group_by(node2) %>%
  mutate(freq = n, weight = n / sum(n)) %>%
  ungroup() %>%
  rename(from = node2, to = node3) %>%
  select(from, to, weight, freq)

edge3 <- df %>%
  count(node3, node4) %>%
  group_by(node3) %>%
  mutate(freq = n, weight = n / sum(n)) %>%
  ungroup() %>%
  rename(from = node3, to = node4) %>%
  select(from, to, weight, freq)

# Add root node
root_edges <- df %>%
  count(node1) %>%
  mutate(freq = n, weight = n / sum(n)) %>%
  mutate(from = "ROOT") %>%
  rename(to = node1) %>%
  select(from, to, weight, freq)
edges <- bind_rows(root_edges, edge1, edge2, edge3)
# Add color by endpoint (M/P)
edges$to_state <- ifelse(
  grepl("_M$", edges$to),
  "M",
  ifelse(grepl("_P$", edges$to), "P", "OTHER")
)
# node freq
node_freq <- edges %>%
  group_by(to) %>%
  summarise(freq = sum(freq), .groups = "drop")
# Build graph
g <- graph_from_data_frame(edges)
# Simplify node labels for plotting
V(g)$label <- sapply(V(g)$name, function(x) {
  if (x == "ROOT") return("ROOT")
  tail(strsplit(x, "\\|")[[1]], 1)
})
# # Node frequency
V(g)$freq <- node_freq$freq[match(V(g)$name, node_freq$to)]
# Assign node colors
V(g)$color <- ifelse(
  grepl("_M$", V(g)$label),
  "#ed2424",   
  ifelse(
    grepl("_P$", V(g)$label),
    "#0b61ec", 
    "grey50"
  )
)
# Plot lineage tree
pdf("endo_cp26_4stage_treeplot.nmlt.pdf", width=8, height=6)
p <- ggraph(g, layout = "tree") +
  geom_edge_diagonal(
    aes(width = weight, color = to_state, label = paste0(round(weight*100,1), "%")),
    angle_calc = "along", # along with edge
    label_dodge = unit(0.3, "cm"),
    # label_push = unit(-1, "cm"),
    label_size = 3,
    alpha = 0.8
  ) +
  scale_edge_width(range = c(0.3, 4)) +
  scale_edge_color_manual(
    values = c(
      "M" = "#f98785",   # red
      "P" = "#82b1fc",   # blue
      "OTHER" = "grey70"
    )
  ) +
  geom_node_point(
    aes(color = I(color)),
    size = 4) +
  geom_node_text(
    # aes(label = label),
    aes(label = freq),
    repel = TRUE,
    size = 4
  ) +
  theme_void() +
  ggtitle("Swith of allelic bias across stages in endo_cp26 lineage") +  
     theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
print(p)
dev.off()




# ==============================================================================
# 04. dynamic change across development and lineages (Fig.2e,2f,2g and Ext.Data.Fig.10)
# We used Myt1l as an example.
# ==============================================================================
# OG01536 - Myt1l
mesen.fringe7 <- read.csv("allelic.rate.new/mesen_A_Fringe7_posterior_mean.csv", row.names = 1)
mesen.hox5 <- read.csv("allelic.rate.new/mesen_A_Hox5_posterior_mean.csv", row.names = 1)
MYT1L.fringe7 <- mesen.fringe7[,"OG01536", drop = FALSE]
MYT1L.fringe7$lineage <- "mesen_A_Fringe7"
MYT1L.hox5 <- mesen.hox5[,"OG01536", drop = FALSE]
MYT1L.hox5$lineage <- "mesen_A_Hox5"

MYT1L.df <- rbind(MYT1L.fringe7, MYT1L.hox5)
MYT1L.df$cell <- rownames(MYT1L.df)
rownames(MYT1L.df) <- NULL
MYT1L.df <- merge(MYT1L.df, cell.info[, c("cell", "stage")], by = "cell", sort = FALSE)
MYT1L.df$OG01536 <- ifelse(MYT1L.df$OG01536 > 1, 1,
                           ifelse(MYT1L.df$OG01536 < 0, 0,
                                  MYT1L.df$OG01536))
write.csv(MYT1L.df, "MYT1L.mesen.lineage.allelic.rate.csv", row.names = FALSE)
pdf("MYT1L.lineage.dodge.boxplot.smooth.pdf", width=8, height=5)
ggplot(MYT1L.df, aes(stage, OG01536)) +
  geom_boxplot(aes(fill = lineage),
               position = position_dodge(0.8),
               alpha = 0.7,
               outlier.shape = NA) +
  geom_point(aes(shape = lineage, color = lineage),
             position = position_dodge(0.8),
             alpha = 0.8, size = 0.6) +
  geom_smooth(aes(group = lineage, color = lineage),
              method = "loess",
              se = FALSE,
              linewidth = 0.7) +
  scale_color_manual(values = c(
    "mesen_A_Fringe7" = "#f7666f",  # red
    "mesen_A_Hox5"  = "#458cef"   # blue
  )) +
  scale_fill_manual(values = c(
    "mesen_A_Fringe7" = "#f7666f",
    "mesen_A_Hox5"  = "#458cef"
  )) +
    theme_classic() +
  labs(
    x = "Stage",
    y = "Estimated allelic rate",
    title = "Allelic rate of MYT1L across stages"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
dev.off()