# ===================================================================================
# 01. Percentage of embryos with normal morphology (Ext.Data.Fig.2)
# ===================================================================================
stages <- c("iniG", "midG", "earN", "latN", "ETB", "MTB", "LatTI", "LatTII", "Larva18", "Larva24")
crosses <- list(
  list(
    key = "C_intestinalis_female_x_C_intestinalis_male",
    label = "C. intestinalis female x C. intestinalis male",
    color = "#009B48",
    ratios = c("85/94", "83/91", "63/75", "82/91", "92/104", "110/121", "128/151", "137/153", "89/108", "73/104")
  ),
  list(
    key = "C_savignyi_female_x_C_savignyi_male",
    label = "C. savignyi female x C. savignyi male",
    color = "#E85B1F",
    ratios = c("99/102", "118/121", "128/135", "116/127", "114/126", "129/141", "138/144", "104/109", "97/116", "101/110")
  ),
  list(
    key = "C_intestinalis_female_x_C_savignyi_male",
    label = "C. intestinalis female x C. savignyi male",
    color = "#3658B5",
    ratios = c("125/128", "173/190", "114/129", "131/143", "122/136", "161/173", "186/204", "140/165", "110/126", "92/104")
  ),
  list(
    key = "C_savignyi_female_x_C_intestinalis_male",
    label = "C. savignyi female x C. intestinalis male",
    color = "#FF1493",
    ratios = c("121/144", "85/133", "46/114", "37/169", "15/146", "26/188", "23/184", "6/69", "5/57", "4/60")
  )
)

ratio_percent <- function(ratio) {
  parts <- strsplit(ratio, "/", fixed = TRUE)[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2]) * 100
}

survival <- data.frame(stage = stages, stringsAsFactors = FALSE)
for (cross in crosses) {
  survival[[cross$key]] <- cross$ratios
  survival[[paste0(cross$key, "_percent")]] <- round(vapply(cross$ratios, ratio_percent, numeric(1)), 2)
}

csv_path <- file.path("cross_stage_survival_rates.csv")
pdf_path <- file.path( "cross_stage_survival_rates_lineplot.pdf")
write.csv(survival, csv_path, row.names = FALSE, quote = TRUE)

x <- seq_along(stages)
y_values <- lapply(crosses, function(cross) vapply(cross$ratios, ratio_percent, numeric(1)))

pdf(pdf_path, width = 8, height = 5)
par(mar = c(5, 4, 4, 6), xpd = NA)
plot(
  x,
  y_values[[1]],
  type = "n",
  xaxt = "n",
  xlim = c(0.75, length(stages) + 1.15),
  ylim = c(0, 105),
  xlab = "Stage",
  ylab = "Survival rate (%)",
  main = "Survival rate by cross and stage"
)
axis(1, at = x, labels = stages, las = 2)
grid(nx = NA, ny = NULL, col = "#D9D9D9", lty = "solid")
box(bty = "l")

for (i in seq_along(crosses)) {
  lines(x, y_values[[i]], type = "o", pch = 16, lwd = 2.2, col = crosses[[i]]$color)
  text(
    x = length(stages) + 0.12,
    y = y_values[[i]][length(stages)],
    labels = sprintf("%.1f%%", y_values[[i]][length(stages)]),
    col = crosses[[i]]$color,
    adj = 0,
    cex = 0.85,
    font = 2
  )
}

legend(
  "bottomleft",
  inset = c(0, 1.02),
  legend = vapply(crosses, function(cross) cross$label, character(1)),
  col = vapply(crosses, function(cross) cross$color, character(1)),
  lty = 1,
  lwd = 2.2,
  pch = 16,
  bty = "n",
  ncol = 2,
  cex = 0.78
)
dev.off()



# ===================================================================================
# 02. Tissue composition of hybrid and pure C. intestinalis embryos  (Ext.Data.Fig.3)
# ===================================================================================
library(ggplot2)
library(patchwork)
library(scales)
input_dir <- "/pure_hybrid_cell_fraction"
output_pdf <- file.path(input_dir, "pure_hybrid_tissue_composition.pdf")
hybrid_file <- file.path(input_dir, "CsCi.10s.cell.anno.region.csv")
pure_file <- file.path(input_dir, "pureCi.10stage.integrated.metadata.csv")
hybrid <- read.csv(hybrid_file, check.names = FALSE, stringsAsFactors = FALSE)
pure <- read.csv(pure_file, check.names = FALSE, stringsAsFactors = FALSE)

required_hybrid <- "tissue.beau.adj"
required_pure <- "tissue.big"
stopifnot(required_hybrid %in% names(hybrid), required_pure %in% names(pure))

clean_tissue <- function(x) {
  x <- trimws(x)
  x[x == "muscle&heart"] <- "muscle & heart"
  x
}

summarize_tissue <- function(x, condition) {
  x <- clean_tissue(x)
  x <- x[!is.na(x) & nzchar(x)]
  counts <- as.data.frame(table(x), stringsAsFactors = FALSE)
  names(counts) <- c("tissue", "cells")
  counts$condition <- condition
  counts$proportion <- counts$cells / sum(counts$cells)
  counts
}

hybrid_sum <- summarize_tissue(hybrid[[required_hybrid]], "Hybrid (CsCi)")
pure_sum <- summarize_tissue(pure[[required_pure]], "Pure (Ci)")

tissue_order <- c(
  "epidermis",
  "nervous system",
  "mesenchyme",
  "endoderm",
  "notochord",
  "muscle & heart",
  "germ"
)

all_tissues <- union(hybrid_sum$tissue, pure_sum$tissue)
unexpected <- setdiff(all_tissues, tissue_order)
if (length(unexpected) > 0) {
  tissue_order <- c(tissue_order, sort(unexpected))
}

tissue_colors <- c(
  "epidermis" = "#605a9a",
  "nervous system" = "#317ec3",
  "mesenchyme" = "#e7cf21",
  "endoderm" = "#70c5c7",
  "notochord" = "#4bc051",
  "muscle & heart" = "#ff9500",
  "germ" = "#386236"
)
if (length(unexpected) > 0) {
  extra_cols <- hue_pal()(length(unexpected))
  names(extra_cols) <- unexpected
  tissue_colors <- c(tissue_colors, extra_cols)
}

hybrid_sum$tissue <- factor(hybrid_sum$tissue, levels = rev(tissue_order))
pure_sum$tissue <- factor(pure_sum$tissue, levels = rev(tissue_order))

x_max <- max(c(hybrid_sum$proportion, pure_sum$proportion))
x_limit <- max(0.45, ceiling((x_max + 0.07) * 20) / 20)

n_label <- function(n) paste0("n = ", comma(n), " cells")

make_panel <- function(df, title_text) {
  ggplot(df, aes(x = proportion, y = tissue, fill = tissue)) +
    geom_col(width = 0.68, color = NA) +
    geom_text(
      aes(label = paste0(
        "N = ", comma(cells),
        "\n",
        percent(proportion, accuracy = 0.1)
      )),
      hjust = -0.14,
      lineheight = 0.92,
      size = 3.0,
      family = "Arial",
      color = "#202020"
    ) +
    scale_fill_manual(values = tissue_colors, guide = "none") +
    scale_x_continuous(
      limits = c(0, x_limit),
      breaks = seq(0, x_limit, by = 0.1),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = title_text,
      subtitle = n_label(sum(df$cells)),
      x = "Proportion of annotated cells",
      y = NULL
    ) +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = "#202020", size = 9.5),
      axis.text.x = element_text(color = "#202020", size = 8.5),
      axis.title.x = element_text(size = 9.5, margin = margin(t = 7)),
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 8.5, color = "#555555", margin = margin(b = 8)),
      plot.margin = margin(8, 18, 8, 7)
    )
}

p_hybrid <- make_panel(hybrid_sum, "Hybrid (CsCi)")
p_pure <- make_panel(pure_sum, "Pure (Ci)")
combined <- p_hybrid + p_pure +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 12, family = "Arial"),
    plot.tag.position = c(0, 1)
  )
ggsave(filename = output_pdf, plot = combined, width = 10, height = 5)



# ====================================================================================================
# 03. Reproducibility of tissue composition and allelic bias across LatTII replicates (Ext.Data.Fig.4)
# ====================================================================================================
### tissue composition (Ext.Data.Fig.4a) ###
input_dir <- "/LTB2_replicates"
files <- list.files(input_dir, pattern = "_cell_anno\\.csv$", full.names = TRUE)
required_col <- "tissue.beau.adj"
read_one <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!required_col %in% colnames(x)) {
    stop("Column '", required_col, "' not found in: ", basename(path))
  }

  replicate_id <- sub("_cell_anno\\.csv$", "", basename(path))
  replicate_label <- sub("^LTB2_", "", replicate_id)
  tissue <- trimws(as.character(x[[required_col]]))
  tissue[tissue %in% c("", "NA", "NaN")] <- NA_character_
  tissue <- gsub("^muscle\\s*&\\s*heart$", "muscle & heart", tissue)

  data.frame(
    replicate = replicate_label,
    replicate_file = basename(path),
    tissue = tissue,
    stringsAsFactors = FALSE
  )
}

dat <- do.call(rbind, lapply(files, read_one))
dat <- dat[!is.na(dat$tissue), , drop = FALSE]

if (nrow(dat) == 0) {
  stop("All tissue labels are empty or NA.")
}

count_df <- as.data.frame(xtabs(~ replicate + tissue, data = dat), stringsAsFactors = FALSE)
names(count_df) <- c("replicate", "tissue", "n")

preferred_tissue_order <- c(
  "epidermis",
  "nervous system",
  "mesenchyme",
  "endoderm",
  "notochord",
  "muscle & heart",
  "germ"
)
observed_tissues <- unique(count_df$tissue)
extra_tissues <- sort(setdiff(observed_tissues, preferred_tissue_order))
tissue_order <- c(preferred_tissue_order[preferred_tissue_order %in% observed_tissues], extra_tissues)

replicate_order <- sub("^LTB2_", "", sub("_cell_anno\\.csv$", "", basename(files)))

replicate_total <- aggregate(n ~ replicate, data = count_df, sum)
names(replicate_total)[2] <- "replicate_total_n"
count_df <- merge(count_df, replicate_total, by = "replicate", all.x = TRUE)
count_df$proportion <- count_df$n / count_df$replicate_total_n
count_df$percent_label <- paste0(sprintf("%.1f", 100 * count_df$proportion), "%")

count_df$replicate <- factor(count_df$replicate, levels = replicate_order)
count_df$tissue <- factor(count_df$tissue, levels = tissue_order)
count_df <- count_df[order(count_df$replicate, count_df$tissue), , drop = FALSE]

summary_out <- count_df[, c("replicate", "tissue", "n", "replicate_total_n", "proportion", "percent_label")]
write.csv(
  summary_out,
  file.path(input_dir, "LTB2_replicates_tissue_beau_adj_circular_barplot_summary.csv"),
  row.names = FALSE
)

empty_bar <- 2
plot_df <- data.frame()
for (rep_name in replicate_order) {
  one <- count_df[count_df$replicate == rep_name, , drop = FALSE]
  empty <- one[rep(1, empty_bar), , drop = FALSE]
  empty[, c("tissue", "n", "replicate_total_n", "proportion", "percent_label")] <- NA
  empty$replicate <- rep_name
  plot_df <- rbind(plot_df, one, empty)
}
plot_df$id <- seq_len(nrow(plot_df))

label_df <- plot_df[!is.na(plot_df$proportion), , drop = FALSE]
number_of_bar <- nrow(plot_df)
label_df$angle <- 90 - 360 * (label_df$id - 0.5) / number_of_bar
label_df$hjust <- ifelse(label_df$angle < -90, 1, 0)
label_df$angle <- ifelse(label_df$angle < -90, label_df$angle + 180, label_df$angle)
label_df$label_y <- pmax(label_df$proportion + 0.035, 0.075)

base_start <- aggregate(id ~ replicate, data = plot_df, FUN = min)
base_end <- aggregate(id ~ replicate, data = plot_df, FUN = max)
names(base_start)[2] <- "start"
names(base_end)[2] <- "end"
base_df <- merge(base_start, base_end, by = "replicate", all = TRUE)
base_df$end <- base_df$end - empty_bar
base_df$replicate_mid <- (base_df$start + base_df$end) / 2
base_df <- merge(base_df, replicate_total, by = "replicate", all.x = TRUE)
base_df$replicate_label <- paste0(base_df$replicate, "\nN = ", scales::comma(base_df$replicate_total_n))
base_df$angle <- 90 - 360 * (base_df$replicate_mid - 0.5) / number_of_bar
base_df$hjust <- ifelse(base_df$angle < -90, 1, 0)
base_df$angle <- ifelse(base_df$angle < -90, base_df$angle + 180, base_df$angle)

grid_df <- data.frame(y = c(0.10, 0.20, 0.30), label = c("10%", "20%", "30%"))

library(ggplot2)
library(scales)
tissue_palette <- c(
  "epidermis" = "#C8B6D9",
  "nervous system" = "#94B8D7",
  "mesenchyme" = "#edd93e",
  "endoderm" = "#70c5c7",
  "notochord" = "#4bc051",
  "muscle & heart" = "#ff9500",
  "germ" = "#386236"
)
present_tissues <- levels(count_df$tissue)
missing_cols <- setdiff(present_tissues, names(tissue_palette))
if (length(missing_cols) > 0) {
  fallback_cols <- hue_pal(l = 75, c = 45)(length(missing_cols))
  names(fallback_cols) <- missing_cols
  tissue_palette <- c(tissue_palette, fallback_cols)
}
tissue_palette <- tissue_palette[present_tissues]

p <- ggplot(plot_df, aes(x = factor(id), y = proportion, fill = tissue)) +
  geom_hline(
    data = grid_df,
    aes(yintercept = y),
    inherit.aes = FALSE,
    color = "#E7E2DB",
    linewidth = 0.35
  ) +
  geom_col(width = 0.82, color = "white", linewidth = 0.28, na.rm = TRUE) +
  geom_text(
    data = label_df,
    aes(x = factor(id), y = label_y, label = percent_label, angle = angle, hjust = hjust),
    inherit.aes = FALSE,
    size = 2.55,
    color = "#333333",
    family = "sans"
  ) +
  geom_text(
    data = base_df,
    aes(x = replicate_mid, y = 0.425, label = replicate_label, angle = angle, hjust = hjust),
    inherit.aes = FALSE,
    size = 3.4,
    lineheight = 0.92,
    fontface = "bold",
    color = "#333333",
    family = "sans"
  ) +
  annotate(
    "text",
    x = 1,
    y = grid_df$y,
    label = grid_df$label,
    hjust = 0,
    size = 2.6,
    color = "#8A8178",
    family = "sans"
  ) +
  scale_fill_manual(values = tissue_palette, name = "Cell type", na.translate = FALSE) +
  scale_y_continuous(limits = c(-0.04, 0.45), expand = c(0, 0)) +
  coord_polar(start = 0, clip = "off") +
  labs(
    title = "Tissue composition across LTB2 replicates",
    subtitle = "Circular barplot; each bar is the within-replicate proportion of one tissue.beau.adj class",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 10.5, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 13.5, color = "#222222", hjust = 0.5, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 9.3, color = "#555555", hjust = 0.5, margin = margin(b = 9)),
    legend.position = "bottom",
    legend.title = element_text(size = 9.5, color = "#333333"),
    legend.text = element_text(size = 9.0, color = "#333333"),
    legend.key.size = unit(0.42, "cm"),
    legend.margin = margin(t = 0),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 12, r = 20, b = 10, l = 20)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

pdf_out <- file.path(input_dir, "LTB2_replicates_tissue_beau_adj_circular_barplot.pdf")
ggsave(pdf_out, p, width = 7.2, height = 7.2, units = "in", device = cairo_pdf, bg = "white")


### within replicates comparison (Ext.Data.Fig.4b) ###
input_dir <- "/LTB2_replicates"
summary_path <- file.path(input_dir, "LTB2_replicates_tissue_beau_adj_stacked_summary.csv")
df <- read.csv(summary_path, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("tissue", "replicate", "within_replicate_prop")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df$within_replicate_prop <- as.numeric(df$within_replicate_prop)
df <- df[
  !is.na(df$within_replicate_prop) & !is.na(df$replicate) & !is.na(df$tissue),
  ,
  drop = FALSE
]

replicate_order <- unique(df$replicate)
tissue_order <- unique(df$tissue)
df$replicate <- factor(df$replicate, levels = replicate_order)
df$tissue <- factor(df$tissue, levels = tissue_order)

if (length(replicate_order) < 2) {
  stop("At least two replicate groups are required.")
}

kw <- kruskal.test(within_replicate_prop ~ replicate, data = df)

complete_tissues <- names(which(table(df$tissue) == length(replicate_order)))
friedman_df <- df[df$tissue %in% complete_tissues, , drop = FALSE]
fried <- friedman.test(within_replicate_prop ~ replicate | tissue, data = friedman_df)

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  as.character(signif(p, 3))
}

kw_label <- paste0(
  "Kruskal-Wallis: H = ",
  sprintf("%.3f", unname(kw$statistic)),
  ", df = ",
  unname(kw$parameter),
  ", p = ",
  format_p(kw$p.value)
)

friedman_label <- paste0(
  "Friedman paired by tissue: chi-squared = ",
  sprintf("%.3f", unname(fried$statistic)),
  ", df = ",
  unname(fried$parameter),
  ", p = ",
  format_p(fried$p.value)
)

group_stats <- do.call(
  rbind,
  lapply(split(df$within_replicate_prop, df$replicate), function(x) {
    data.frame(
      n_tissues = length(x),
      median = median(x),
      q1 = unname(quantile(x, 0.25)),
      q3 = unname(quantile(x, 0.75)),
      min = min(x),
      max = max(x),
      mean = mean(x),
      stringsAsFactors = FALSE
    )
  })
)
group_stats$replicate <- rownames(group_stats)
rownames(group_stats) <- NULL
group_stats <- group_stats[, c("replicate", "n_tissues", "median", "q1", "q3", "min", "max", "mean")]

test_out <- rbind(
  data.frame(
    test = "Kruskal-Wallis rank sum test",
    statistic = unname(kw$statistic),
    df = unname(kw$parameter),
    p_value = kw$p.value,
    paired_by_tissue = FALSE,
    n_groups = length(replicate_order),
    total_observations = nrow(df),
    stringsAsFactors = FALSE
  ),
  data.frame(
    test = "Friedman rank sum test",
    statistic = unname(fried$statistic),
    df = unname(fried$parameter),
    p_value = fried$p.value,
    paired_by_tissue = TRUE,
    n_groups = length(replicate_order),
    total_observations = nrow(friedman_df),
    stringsAsFactors = FALSE
  )
)

write.csv(
  group_stats,
  file.path(input_dir, "LTB2_replicates_within_replicate_prop_boxplot_group_stats.csv"),
  row.names = FALSE
)
write.csv(
  test_out,
  file.path(input_dir, "LTB2_replicates_within_replicate_prop_boxplot_test_results.csv"),
  row.names = FALSE
)

library(ggplot2)
library(scales)

replicate_palette <- c(
  "Cs1Ci1" = "#9DB6A3",
  "Cs1Ci2" = "#D6A6A1",
  "Cs2Ci1" = "#A9B7D1"
)

missing_cols <- setdiff(replicate_order, names(replicate_palette))
if (length(missing_cols) > 0) {
  fallback_cols <- hue_pal(l = 70, c = 55)(length(missing_cols))
  names(fallback_cols) <- missing_cols
  replicate_palette <- c(replicate_palette, fallback_cols)
}
replicate_palette <- replicate_palette[replicate_order]

y_top <- max(df$within_replicate_prop, na.rm = TRUE) * 1.34
y_top <- max(y_top, 0.42)

p <- ggplot(df, aes(x = replicate, y = within_replicate_prop, fill = replicate)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "#333333",
    linewidth = 0.55,
    alpha = 0.92
  ) +
  annotate(
    "text",
    x = mean(seq_along(replicate_order)),
    y = y_top * 0.965,
    label = paste(kw_label, friedman_label, sep = "\n"),
    size = 3.25,
    lineheight = 1.05,
    color = "#333333",
    family = "sans"
  ) +
  scale_fill_manual(values = replicate_palette, guide = "none") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, y_top),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Distribution of within-replicate tissue proportions",
    subtitle = "Each box summarizes the seven tissue.beau.adj classes within one replicate",
    x = "Replicate",
    y = "Within-replicate proportion"
  ) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#222222"),
    plot.subtitle = element_text(size = 9.8, color = "#555555", margin = margin(b = 8)),
    axis.title.x = element_text(color = "#333333", margin = margin(t = 8)),
    axis.title.y = element_text(color = "#333333", margin = margin(r = 8)),
    axis.text = element_text(color = "#333333"),
    axis.line = element_line(linewidth = 0.45, color = "#333333"),
    axis.ticks = element_line(linewidth = 0.35, color = "#333333"),
    panel.grid.major.y = element_line(color = "#E8E2DA", linewidth = 0.35),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 12, r = 16, b = 12, l = 16)
  )

pdf_out <- file.path(input_dir, "LTB2_replicates_within_replicate_prop_boxplot_tests.pdf")
ggsave(pdf_out, p, width = 5.8, height = 4.6, units = "in", device = cairo_pdf, bg = "white")



### allelic rates correlation between replicates -- density plot (Ext.Data.Fig.4c) ###
input_dir <- "/LTB2_replicates"
rate_pattern <- "_cutoff30_joint_avg_psmean_allcells\\.csv$"
og_file <- file.path(input_dir, "OG_genelist_without_mlt.csv")

files <- sort(list.files(input_dir, pattern = rate_pattern, full.names = TRUE))
og_table <- read.csv(og_file, stringsAsFactors = FALSE, check.names = FALSE)
ci_og <- unique(og_table$Ci_OG[!is.na(og_table$Ci_OG) & nzchar(og_table$Ci_OG)])

read_rate_table <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("gene", "mean_rate")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("Missing column(s) in ", basename(path), ": ", paste(missing, collapse = ", "))
  }
  x <- x[, required]
  x$mean_rate <- suppressWarnings(as.numeric(x$mean_rate))
  x <- x[!is.na(x$gene) & nzchar(x$gene) & is.finite(x$mean_rate), ]
  x <- x[!duplicated(x$gene), ]
  names(x)[2] <- tools::file_path_sans_ext(basename(path))
  x
}

short_name <- function(path) {
  sub(rate_pattern, "", basename(path))
}

format_p <- function(p) {
  if (!is.finite(p)) {
    return("NA")
  }
  format.pval(p, digits = 3, eps = 1e-300)
}

plot_pair <- function(file_x, file_y) {
  name_x <- short_name(file_x)
  name_y <- short_name(file_y)
  x <- read_rate_table(file_x)
  y <- read_rate_table(file_y)
  merged <- merge(x, y, by = "gene", all = FALSE)
  names(merged) <- c("gene", "rate_x", "rate_y")
  shared_genes_before_nmlt <- nrow(merged)
  merged <- merged[merged$gene %in% ci_og, ]
  merged <- merged[is.finite(merged$rate_x) & is.finite(merged$rate_y), ]

  if (nrow(merged) < 3) {
    stop("Too few shared finite nmlt genes for ", name_x, " vs ", name_y, ": ", nrow(merged))
  }

  cor_res <- suppressWarnings(cor.test(merged$rate_x, merged$rate_y, method = "spearman", exact = FALSE))
  rho <- unname(cor_res$estimate)
  annotation <- sprintf("Spearman rho = %.3f    P = %s    shared nmlt genes = %s",
                        rho, format_p(cor_res$p.value), format(nrow(merged), big.mark = ","))

  point_cols <- densCols(merged$rate_x, merged$rate_y,
                         colramp = colorRampPalette(c("#25408F", "#2C7FB8", "#41B6C4", "#A1DAB4", "#FFFFCC")))
  order_idx <- order(col2rgb(point_cols)[1, ])

  out_file <- file.path(input_dir, paste0(name_x, "_vs_", name_y, "_allelic_rate_density_scatter_nmlt.pdf"))
  pdf(out_file, width = 5.4, height = 5.7, useDingbats = FALSE)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mar = c(4.7, 4.9, 3.7, 1.2), mgp = c(2.8, 0.8, 0))
  plot(merged$rate_x[order_idx], merged$rate_y[order_idx],
       pch = 16, cex = 0.55, col = point_cols[order_idx],
       xlim = c(0, 1), ylim = c(0, 1), asp = 1,
       xlab = paste0(name_x, " allelic rate"),
       ylab = paste0(name_y, " allelic rate"),
       main = paste0(name_x, " vs ", name_y, " nmlt"))
  abline(0, 1, lwd = 1.3, lty = 2, col = "gray35")
  grid(col = adjustcolor("gray80", alpha.f = 0.55), lty = 1)
  box(lwd = 1.1)
  mtext(annotation, side = 3, line = 0.35, cex = 0.82)

  invisible(data.frame(
    comparison = paste(name_x, name_y, sep = "_vs_"),
    shared_genes_before_nmlt = shared_genes_before_nmlt,
    shared_nmlt_genes = nrow(merged),
    spearman_rho = rho,
    p_value = cor_res$p.value,
    output_pdf = out_file,
    stringsAsFactors = FALSE
  ))
}
pairs <- combn(files, 2, simplify = FALSE)
summary_df <- do.call(rbind, lapply(pairs, function(pair) plot_pair(pair[1], pair[2])))
write.csv(summary_df,
          file.path(input_dir, "LTB2_replicates_pairwise_allelic_rate_spearman_summary_nmlt.csv"),
          row.names = FALSE)



### allelic rates correlation between replicates within each tissue -- heatmap plot (Ext.Data.Fig.4d) ###
input_dir <- "/LTB2_replicates"
rate_pattern <- "_cutoff30_joint_avg_psmean_7tissue\\.csv$"
og_file <- file.path(input_dir, "OG_genelist_without_mlt.csv")

files <- sort(list.files(input_dir, pattern = rate_pattern, full.names = TRUE))
og_table <- read.csv(og_file, stringsAsFactors = FALSE, check.names = FALSE)
ci_og <- unique(og_table$Ci_OG[!is.na(og_table$Ci_OG) & nzchar(og_table$Ci_OG)])

short_name <- function(path) {
  sub(rate_pattern, "", basename(path))
}

read_tissue_table <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(x) < 2) {
    stop("Expected a gene column plus tissue columns in ", basename(path))
  }

  tissue_cols <- names(x)[-1]
  out <- x[, seq_len(ncol(x)), drop = FALSE]
  names(out)[1] <- "gene"
  out <- out[!is.na(out$gene) & nzchar(out$gene), , drop = FALSE]
  out <- out[!duplicated(out$gene), , drop = FALSE]

  for (col in tissue_cols) {
    out[[col]] <- suppressWarnings(as.numeric(out[[col]]))
  }
  out
}

star_label <- function(p) {
  if (!is.finite(p)) {
    return("")
  }
  if (p < 0.001) {
    return("***")
  }
  if (p < 0.01) {
    return("**")
  }
  if (p < 0.05) {
    return("*")
  }
  ""
}

tables <- lapply(files, read_tissue_table)
names(tables) <- vapply(files, short_name, character(1))
tissue_cols <- names(tables[[1]])[-1]

for (nm in names(tables)) {
  this_cols <- names(tables[[nm]])[-1]
  if (!identical(this_cols, tissue_cols)) {
    stop("Tissue columns differ for ", nm, ". Expected: ",
         paste(tissue_cols, collapse = ", "), "; found: ", paste(this_cols, collapse = ", "))
  }
}

pairs <- combn(names(tables), 2, simplify = FALSE)
pair_labels <- vapply(pairs, function(pair) paste(pair[1], pair[2], sep = " vs "), character(1))

results <- list()
for (i in seq_along(pairs)) {
  pair <- pairs[[i]]
  x <- tables[[pair[1]]]
  y <- tables[[pair[2]]]
  merged <- merge(x, y, by = "gene", suffixes = c(".x", ".y"), all = FALSE)
  shared_before_nmlt <- nrow(merged)
  merged <- merged[merged$gene %in% ci_og, , drop = FALSE]
  shared_nmlt <- nrow(merged)

  for (tissue in tissue_cols) {
    x_col <- paste0(tissue, ".x")
    y_col <- paste0(tissue, ".y")
    keep <- is.finite(merged[[x_col]]) & is.finite(merged[[y_col]])
    dat <- merged[keep, , drop = FALSE]

    if (nrow(dat) >= 3) {
      cor_res <- suppressWarnings(cor.test(dat[[x_col]], dat[[y_col]], method = "spearman", exact = FALSE))
      rho <- unname(cor_res$estimate)
      p_value <- cor_res$p.value
    } else {
      rho <- NA_real_
      p_value <- NA_real_
    }

    results[[length(results) + 1]] <- data.frame(
      comparison = pair_labels[i],
      tissue = tissue,
      shared_genes_before_nmlt = shared_before_nmlt,
      shared_nmlt_genes = shared_nmlt,
      finite_genes_for_tissue = nrow(dat),
      spearman_rho = rho,
      p_value = p_value,
      stars = star_label(p_value),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, results)
summary_file <- file.path(input_dir, "LTB2_7tissue_pairwise_spearman_heatmap_nmlt_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

rho_mat <- matrix(NA_real_, nrow = length(pair_labels), ncol = length(tissue_cols),
                  dimnames = list(pair_labels, tissue_cols))
stars_mat <- matrix("", nrow = length(pair_labels), ncol = length(tissue_cols),
                    dimnames = list(pair_labels, tissue_cols))

for (row_i in seq_len(nrow(summary_df))) {
  rho_mat[summary_df$comparison[row_i], summary_df$tissue[row_i]] <- summary_df$spearman_rho[row_i]
  stars_mat[summary_df$comparison[row_i], summary_df$tissue[row_i]] <- summary_df$stars[row_i]
}

out_pdf <- file.path(input_dir, "LTB2_7tissue_pairwise_spearman_heatmap_nmlt.pdf")
pdf(out_pdf, width = 10.2, height = 4.6, useDingbats = FALSE)
old_par <- par(no.readonly = TRUE)
on.exit({
  par(old_par)
  dev.off()
}, add = TRUE)

palette_fun <- colorRampPalette(c("#FFFFFF", "#B2182B"))
palette_cols <- palette_fun(101)
value_to_col <- function(v) {
  if (!is.finite(v)) {
    return("#E5E5E5")
  }
  scaled <- (v - 0.5) / 0.5
  palette_cols[pmin(101, pmax(1, floor(scaled * 100) + 1))]
}

n_rows <- nrow(rho_mat)
n_cols <- ncol(rho_mat)
par(mar = c(5.2, 10.8, 3.4, 5.8), xpd = NA)
plot(NA, xlim = c(0, n_cols + 1.35), ylim = c(0, n_rows),
     xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "",
     main = "LTB2 7-tissue pairwise Spearman correlation (nmlt genes)")

for (i in seq_len(n_rows)) {
  y_bottom <- n_rows - i
  for (j in seq_len(n_cols)) {
    v <- rho_mat[i, j]
    rect(j - 1, y_bottom, j, y_bottom + 1,
         col = value_to_col(v), border = "white", lwd = 1.4)
    value_label <- if (is.finite(v)) sprintf("%.3f", v) else "NA"
    text(j - 0.5, y_bottom + 0.60, value_label, cex = 0.86, font = 2)
    text(j - 0.5, y_bottom + 0.38, stars_mat[i, j], cex = 0.82, font = 2)
  }
}

axis(2, at = n_rows - seq_len(n_rows) + 0.5, labels = rownames(rho_mat), las = 2, tick = FALSE, cex.axis = 0.86)
text(seq_len(n_cols) - 0.5, -0.18, labels = tissue_cols, srt = 30, adj = 1, cex = 0.86)

legend_x0 <- n_cols + 0.35
legend_y0 <- 0.15
legend_y1 <- n_rows - 0.15
legend_steps <- 101
legend_ys <- seq(legend_y0, legend_y1, length.out = legend_steps + 1)
for (k in seq_len(legend_steps)) {
  rect(legend_x0, legend_ys[k], legend_x0 + 0.22, legend_ys[k + 1],
       col = palette_cols[k], border = NA)
}
axis(4, at = seq(legend_y0, legend_y1, length.out = 6),
     labels = sprintf("%.1f", seq(0.5, 1, length.out = 6)),
     las = 1, tick = FALSE, cex.axis = 0.75)
text(legend_x0 + 0.11, legend_y1 + 0.22, "spearman correlation", cex = 0.78, font = 2)



# =============================================================================================================
# 04. Developmental dynamics of allelic bias composition (Ext.Data.Fig.6)
# We used panel a and b as examples. (The plotting code of panel c and Ext.Data.Fig.8 are similar to panel b.)
# =============================================================================================================
input_csv <- "/lineage_bias_statistic/stage_bias_all_node_bias_with_stage.csv"
output_dir <- "/lineage_bias_statistic"

stage_levels <- c("iniG", "midG", "earN", "latN", "ETB", "MTB", "LTB1", "LTB2", "larva")
bias_levels <- c("maternal", "balanced", "paternal")
bias_colors <- c(
  maternal = "#e34f5b",
  balanced = "#8A8F98",
  paternal = "#0072B2"
)

remove_nodes <- c(
  "i.larva_33",
  "i.larva_21",
  "i.larva_24",
  "i.larva_27",
  "i.larva_35",
  "i.larva_10",
  "a.iniG_b_line",
  "i.larva_Tll1+ mesenchyme",
  "i.larva_noto2"
)

nodes <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("source_file", "lineage", "node", "stage", "maternal", "balanced", "paternal")
missing <- setdiff(needed, names(nodes))
if (length(missing) > 0) {
  stop("Missing columns in input CSV: ", paste(missing, collapse = ", "))
}

removed_rows <- nodes[nodes$node %in% remove_nodes, ]
filtered <- nodes[!(nodes$node %in% remove_nodes), ]

numeric_cols <- c("maternal", "balanced", "paternal")
filtered[numeric_cols] <- lapply(filtered[numeric_cols], as.numeric)
filtered$stage <- factor(filtered$stage, levels = stage_levels)
filtered <- filtered[order(filtered$stage, filtered$lineage, filtered$node), ]
filtered$stage <- as.character(filtered$stage)

write.csv(
  removed_rows,
  file.path(output_dir, "stage_bias_removed_nodes.csv"),
  row.names = FALSE
)
write.csv(
  filtered,
  file.path(output_dir, "stage_bias_all_node_bias_with_stage_filtered.csv"),
  row.names = FALSE
)

long_counts <- rbind(
  data.frame(stage = filtered$stage, node = filtered$node, lineage = filtered$lineage, bias = "maternal", count = filtered$maternal),
  data.frame(stage = filtered$stage, node = filtered$node, lineage = filtered$lineage, bias = "balanced", count = filtered$balanced),
  data.frame(stage = filtered$stage, node = filtered$node, lineage = filtered$lineage, bias = "paternal", count = filtered$paternal)
)
long_counts$stage <- factor(long_counts$stage, levels = stage_levels)
long_counts$bias <- factor(long_counts$bias, levels = bias_levels)

stage_bias <- aggregate(count ~ stage + bias, data = long_counts, FUN = sum, na.rm = TRUE)
stage_totals <- aggregate(count ~ stage, data = stage_bias, FUN = sum, na.rm = TRUE)
names(stage_totals)[2] <- "stage_total"
stage_bias <- merge(stage_bias, stage_totals, by = "stage", all.x = TRUE)
stage_bias$percent <- ifelse(stage_bias$stage_total > 0, stage_bias$count / stage_bias$stage_total * 100, NA_real_)
stage_bias$label <- ifelse(is.na(stage_bias$percent), "", sprintf("%.1f%%", stage_bias$percent))
stage_bias$stage <- factor(stage_bias$stage, levels = stage_levels)
stage_bias$bias <- factor(stage_bias$bias, levels = bias_levels)
stage_bias <- stage_bias[order(stage_bias$stage, stage_bias$bias), ]

write.csv(
  stage_bias,
  file.path(output_dir, "stage_bias_summary_counts_percent_filtered_long.csv"),
  row.names = FALSE
)

node_totals <- filtered$maternal + filtered$balanced + filtered$paternal
node_percent <- rbind(
  data.frame(stage = filtered$stage, node = filtered$node, lineage = filtered$lineage, bias = "maternal", percent = filtered$maternal / node_totals * 100),
  data.frame(stage = filtered$stage, node = filtered$node, lineage = filtered$lineage, bias = "balanced", percent = filtered$balanced / node_totals * 100),
  data.frame(stage = filtered$stage, node = filtered$node, lineage = filtered$lineage, bias = "paternal", percent = filtered$paternal / node_totals * 100)
)
node_percent <- node_percent[is.finite(node_percent$percent), ]
node_percent$stage <- factor(node_percent$stage, levels = stage_levels)
node_percent$bias <- factor(node_percent$bias, levels = bias_levels)
node_percent <- node_percent[order(node_percent$stage, node_percent$bias, node_percent$lineage, node_percent$node), ]

write.csv(
  node_percent,
  file.path(output_dir, "stage_bias_node_level_percent_filtered_long.csv"),
  row.names = FALSE
)

plot_stacked <- function(path, stage_bias) {
  plot_data <- xtabs(percent ~ bias + stage, data = stage_bias)
  pdf(path, width = 10, height = 6)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mar = c(7.5, 5.5, 4, 8), mgp = c(3.4, 0.8, 0), xpd = FALSE)
  mids <- barplot(
    plot_data,
    col = bias_colors[rownames(plot_data)],
    border = "white",
    ylim = c(0, 100),
    xlim = c(0, ncol(plot_data) * 1.35 + 1.8),
    ylab = "Percentage",
    xlab = "Stage",
    main = "Allelic bias composition by stage after node filtering",
    names.arg = rep("", ncol(plot_data)),
    yaxt = "n"
  )
  axis(2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
  abline(h = seq(20, 80, 20), col = "#E5E7EB", lwd = 0.8)
  axis(1, at = mids, labels = colnames(plot_data), las = 2, tick = FALSE)
  box()

  par(xpd = NA)
  legend(
    x = max(mids) + 1.0,
    y = 96,
    legend = bias_levels,
    fill = bias_colors[bias_levels],
    border = NA,
    bty = "n",
    xjust = 0
  )

  for (j in seq_len(ncol(plot_data))) {
    y0 <- 0
    label_positions <- numeric(nrow(plot_data))
    for (i in seq_len(nrow(plot_data))) {
      val <- plot_data[i, j]
      label_positions[i] <- y0 + val / 2
      y0 <- y0 + ifelse(is.na(val), 0, val)
    }
    for (i in 2:length(label_positions)) {
      if (label_positions[i] - label_positions[i - 1] < 3.0) {
        label_positions[i] <- label_positions[i - 1] + 3.0
      }
    }
    label_positions <- pmin(label_positions, 99)
    for (i in seq_len(nrow(plot_data))) {
      val <- plot_data[i, j]
      if (!is.na(val) && val > 0) {
        text(mids[j], label_positions[i], sprintf("%.1f%%", val), cex = 0.72)
      }
    }
  }
}

plot_boxplot <- function(path, node_percent) {
  stages <- stage_levels[stage_levels %in% as.character(unique(node_percent$stage))]
  stage_centers <- seq_along(stages) * 4
  bias_offsets <- c(maternal = -0.85, balanced = 0, paternal = 0.85)

  box_values <- list()
  positions <- numeric()
  box_colors <- character()
  for (stage in stages) {
    for (bias in bias_levels) {
      values <- node_percent$percent[node_percent$stage == stage & node_percent$bias == bias]
      box_values[[paste(stage, bias, sep = "_")]] <- values
      positions <- c(positions, stage_centers[match(stage, stages)] + bias_offsets[bias])
      box_colors <- c(box_colors, bias_colors[bias])
    }
  }

  pdf(path, width = 11, height = 6.2)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mar = c(7.5, 5.5, 4, 8), mgp = c(3.4, 0.8, 0), xpd = FALSE)
  boxplot(
    box_values,
    at = positions,
    col = box_colors,
    border = "#333333",
    outline = TRUE,
    pch = 16,
    cex = 0.45,
    xaxt = "n",
    yaxt = "n",
    names = rep("", length(box_values)),
    staplewex = 0,
    ylim = c(0, 100),
    xlim = c(min(stage_centers) - 2.1, max(stage_centers) + 5.0),
    ylab = "Node-level percentage",
    xlab = "Stage",
    main = "Node-level allelic bias percentage by stage after node filtering"
  )
  axis(2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
  abline(h = seq(20, 80, 20), col = "#E5E7EB", lwd = 0.8)
  axis(1, at = stage_centers, labels = stages, las = 2, tick = FALSE)

  for (bias in bias_levels) {
    subset <- node_percent[node_percent$bias == bias, ]
    if (length(unique(subset$stage_index)) < 2) {
        next
    }
    fit <- lm(percent ~ stage_index, data = subset)
    pred <- predict(fit, newdata = data.frame(stage_index = seq_along(stages)))
    line_x <- stage_centers + bias_offsets[bias]
    lines(line_x, pred, col = bias_colors[bias], lwd = 2.4)
    points(line_x, pred, col = bias_colors[bias], pch = 19, cex = 0.65)
  }
  par(xpd = NA)
  legend(
    x = max(stage_centers) + 1.8,
    y = 96,
    legend = bias_levels,
    fill = bias_colors[bias_levels],
    border = NA,
    bty = "n",
    xjust = 0
  )
}

stacked_pdf <- file.path(output_dir, "stage_bias_filtered_percent_stacked.pdf")
boxplot_pdf <- file.path(output_dir, "stage_bias_filtered_node_percent_boxplot.pdf")

plot_stacked(stacked_pdf, stage_bias)
plot_boxplot(boxplot_pdf, node_percent)


# =======================================================================================================
# 05. Allelic bias composition across related larval lineages within individual tissues (Ext.Data.Fig.7)
# =======================================================================================================
input_csv <- "/stage_bias_all_node_bias_with_stage_filtered.csv"
output_dir <- "/lineage_bias_statistic"

bias_levels <- c("maternal", "balanced", "paternal")
bias_colors <- c(
  maternal = "#e34f5b",
  balanced = "#8A8F98",
  paternal = "#0072B2"
)
tissues <- c("endoderm", "mesenchyme", "CNS")

nodes <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("lineage", "node", "stage", "maternal", "balanced", "paternal")
missing <- setdiff(needed, names(nodes))
if (length(missing) > 0) {
  stop("Missing columns in input CSV: ", paste(missing, collapse = ", "))
}

nodes$maternal <- as.numeric(nodes$maternal)
nodes$balanced <- as.numeric(nodes$balanced)
nodes$paternal <- as.numeric(nodes$paternal)
nodes$total <- nodes$maternal + nodes$balanced + nodes$paternal
nodes$tissue <- sub("_.*$", "", nodes$lineage)

larva <- nodes[nodes$stage == "larva" & nodes$tissue %in% tissues & nodes$total > 0, ]
if (nrow(larva) == 0) {
  stop("No larva-stage rows found for requested tissues.")
}

larva$maternal_percent <- larva$maternal / larva$total * 100
larva$balanced_percent <- larva$balanced / larva$total * 100
larva$paternal_percent <- larva$paternal / larva$total * 100

larva <- larva[order(factor(larva$tissue, levels = tissues), larva$lineage, larva$node), ]
write.csv(
  larva,
  file.path(output_dir, "larva_tissue_node_bias_percent_filtered.csv"),
  row.names = FALSE
)

plot_tissue <- function(tissue_name) {
  x <- larva[larva$tissue == tissue_name, ]
  if (nrow(x) == 0) {
    warning("No rows for tissue: ", tissue_name)
    return(NULL)
  }

  x <- x[order(x$lineage, x$node), ]
  percent_mat <- rbind(
    maternal = x$maternal_percent,
    balanced = x$balanced_percent,
    paternal = x$paternal_percent
  )
  colnames(percent_mat) <- paste(x$lineage, x$node, sep = "\n")

  pdf_file <- file.path(
    output_dir,
    paste0("larva_", tissue_name, "_node_bias_percent_horizontal_stacked.pdf")
  )

  height <- max(5.2, 2.8 + nrow(x) * 0.34)
  width <- if (tissue_name == "CNS") 12.5 else 13.2

  pdf(pdf_file, width = width, height = height)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mar = c(5, 16, 4, 4), mgp = c(3, 0.8, 0), xpd = FALSE)
  mids <- barplot(
    percent_mat,
    horiz = TRUE,
    col = bias_colors[rownames(percent_mat)],
    border = "white",
    xlim = c(0, 132),
    xaxt = "n",
    las = 1,
    cex.names = 0.62,
    xlab = "Percentage",
    ylab = "",
    main = paste0("Larva ", tissue_name, " node allelic bias composition")
  )
  axis(1, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"))
  abline(v = seq(20, 80, 20), col = "#E5E7EB", lwd = 0.8)
  box()

  for (j in seq_len(ncol(percent_mat))) {
    x0 <- 0
    for (i in seq_len(nrow(percent_mat))) {
      value <- percent_mat[i, j]
      if (!is.na(value) && value >= 5) {
        text(
          x = x0 + value / 2,
          y = mids[j],
          labels = sprintf("%.1f%%", value),
          cex = 0.62,
          col = "black"
        )
      }
      x0 <- x0 + ifelse(is.na(value), 0, value)
    }
  }

  legend(
    x = 104,
    y = max(mids),
    legend = bias_levels,
    fill = bias_colors[bias_levels],
    border = NA,
    bty = "n",
    cex = 0.9
  )

  pdf_file
}
outputs <- unlist(lapply(tissues, plot_tissue), use.names = FALSE)



# =======================================================================================================
# 06. Representative allelic state transitions of Ank3  (Ext.Data.Fig.9)
# =======================================================================================================
rate_dir <- "/lineage_node_allelic_rate"
output_dir <- "/lineage_bias_statistic"
gene_id <- "OG00689" # Ank3

transitions <- data.frame(
  category = c("bi_bi_switch", "bi_bi_stable", "bi_ba", "ba_ba", "ba_bi"), # switched, propagated, attenuated, buffered, exposed
  lineage = c(
    "CNS_a8.17_20_25_26",
    "CNS_A8.7_8_15_16",
    "endoderm_A7.1_2_5",
    "CNS_b8.17_20",
    "mesenchyme_B7.7_8.5"
  ),
  node = c(
    "e.ETB_22",
    "h.LTB2_5",
    "g.LTB1_cl0.3",
    "e.ETB_21",
    "d.latN_cl1"
  ),
  parent_node = c(
    "d.latN_11",
    "g.LTB1_2",
    "f.MTB_cl0-2.3",
    "d.latN_cl1-2",
    "c.earN_cl4.5"
  ),
  stringsAsFactors = FALSE
)

transitions$rate_csv <- file.path(
  rate_dir,
  paste0(transitions$lineage, "_cutoff30_joint_avg_psmean_10stage_nodes.csv")
)

read_gene_rates <- function(path, gene_id, parent_node, node) {
  if (!file.exists(path)) {
    stop("Missing allelic-rate CSV: ", path)
  }

  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(x)[1] <- "gene"
  missing_nodes <- setdiff(c(parent_node, node), names(x))
  if (length(missing_nodes) > 0) {
    stop("Missing node columns in ", basename(path), ": ", paste(missing_nodes, collapse = ", "))
  }

  gene_row <- x[x$gene == gene_id, , drop = FALSE]
  if (nrow(gene_row) != 1) {
    stop("Expected exactly one row for ", gene_id, " in ", basename(path), "; found ", nrow(gene_row))
  }

  c(
    parent_rate = as.numeric(gene_row[[parent_node]]),
    node_rate = as.numeric(gene_row[[node]])
  )
}

rates <- do.call(rbind, lapply(seq_len(nrow(transitions)), function(i) {
  row <- transitions[i, ]
  values <- read_gene_rates(row$rate_csv, gene_id, row$parent_node, row$node)
  data.frame(
    gene = gene_id,
    category = row$category,
    lineage = row$lineage,
    parent_node = row$parent_node,
    node = row$node,
    parent_rate = values[["parent_rate"]],
    node_rate = values[["node_rate"]],
    delta = values[["node_rate"]] - values[["parent_rate"]],
    rate_csv = basename(row$rate_csv),
    stringsAsFactors = FALSE
  )
}))

rate_table <- file.path(output_dir, "OG00689_selected_transition_allelic_rates.csv")
write.csv(rates, rate_table, row.names = FALSE)

pdf_file <- file.path(output_dir, "OG00689_five_transition_categories_allelic_rate_changes.pdf")

category_colors <- c(
  bi_bi_switch = "#7B4FA3",
  bi_bi_stable = "#4D4D4D",
  bi_ba = "#009E73",
  ba_ba = "#CC79A7",
  ba_bi = "#E69F00"
)

category_labels <- c(
  bi_bi_switch = "bi_bi_switch",
  bi_bi_stable = "bi_bi_stable",
  bi_ba = "bi_ba",
  ba_ba = "ba_ba",
  ba_bi = "ba_bi"
)

pdf(pdf_file, width = 10.8, height = 5.9)
old_par <- par(no.readonly = TRUE)
on.exit({
  par(old_par)
  dev.off()
}, add = TRUE)

par(mar = c(5.8, 5.3, 4.2, 2), mgp = c(3.2, 0.8, 0), xpd = FALSE)

x_parent <- seq_len(nrow(rates)) - 0.18
x_node <- seq_len(nrow(rates)) + 0.18
y_values <- c(rates$parent_rate, rates$node_rate)
y_min <- max(0, min(y_values, na.rm = TRUE) - 0.08)
y_max <- min(1, max(y_values, na.rm = TRUE) + 0.08)
if (y_max - y_min < 0.25) {
  mid <- mean(c(y_min, y_max))
  y_min <- max(0, mid - 0.125)
  y_max <- min(1, mid + 0.125)
}

plot(
  NA,
  xlim = c(0.45, nrow(rates) + 0.55),
  ylim = c(y_min, y_max),
  xaxt = "n",
  xlab = "",
  ylab = "Estimated allelic rate",
  main = paste0(gene_id, " allelic rate changes across transition categories")
)
abline(h = seq(0, 1, 0.1), col = "#E5E7EB", lwd = 0.7)
abline(h = c(0.45, 0.55), col = "#BDBDBD", lty = 2, lwd = 0.9)

for (i in seq_len(nrow(rates))) {
  category <- rates$category[i]
  col <- category_colors[[category]]
  arrows(
    x0 = x_parent[i],
    y0 = rates$parent_rate[i],
    x1 = x_node[i],
    y1 = rates$node_rate[i],
    col = col,
    lwd = 2.2,
    length = 0.08,
    angle = 20
  )
  points(
    c(x_parent[i], x_node[i]),
    c(rates$parent_rate[i], rates$node_rate[i]),
    pch = c(21, 19),
    bg = c("white", col),
    col = col,
    cex = 1.15,
    lwd = 1.4
  )
  text(
    c(x_parent[i], x_node[i]),
    c(rates$parent_rate[i], rates$node_rate[i]),
    labels = sprintf("%.3f", c(rates$parent_rate[i], rates$node_rate[i])),
    pos = c(2, 4),
    cex = 0.75,
    col = "#222222"
  )
}

axis(
  1,
  at = seq_len(nrow(rates)),
  labels = category_labels[rates$category],
  las = 2,
  tick = FALSE
)
mtext("Transition category", side = 1, line = 4.4)

legend(
  "topright",
  legend = c("parent node", "child node", "balanced thresholds (0.45, 0.55)"),
  pch = c(21, 19, NA),
  pt.bg = c("white", "#666666", NA),
  col = c("#666666", "#666666", "#BDBDBD"),
  lty = c(NA, NA, 2),
  lwd = c(NA, NA, 1),
  bty = "n",
  cex = 0.85
)
box()



# =======================================================================================================
# 07. Concordance between RNA and ATAC allelic bias across tissues (Ext.Data.Fig.12)
# =======================================================================================================
library(data.table)
library(ggplot2)
work_dir <- "/MTB_ATAC_RNA/nmlt_gene_by_tissue"
tissue_order <- c(
  "endoderm", "epidermis", "germ", "heart",
  "mesenchyme", "muscle", "nervous_system", "notochord"
)
input_files <- file.path(
  work_dir,
  paste0("CsCi_MTB_RNA_ATAC_", tissue_order, "_common_expressed_genes_allelic_rate.csv")
)
required_cols <- c("gene", "rna_bias", "compare")

all_dt <- rbindlist(
  lapply(seq_along(input_files), function(i) {
    path <- input_files[[i]]
    x <- fread(path, header = TRUE, check.names = FALSE)
    missing_cols <- setdiff(required_cols, names(x))
    if (length(missing_cols) > 0L) {
      stop(
        "Missing required column(s) in ", basename(path), ": ",
        paste(missing_cols, collapse = ", ")
      )
    }
    x[, .(
      tissue = tissue_order[[i]],
      gene,
      rna_bias,
      compare
    )]
  }),
  use.names = TRUE
)

all_dt <- all_dt[
  tissue %in% tissue_order &
    rna_bias %in% c("maternal", "paternal") &
    compare %in% c("same_class", "different_class")
]

summary_dt <- all_dt[, .N, by = .(tissue, rna_bias, compare)]
complete_grid <- CJ(
  tissue = tissue_order,
  rna_bias = c("maternal", "paternal"),
  compare = c("same_class", "different_class"),
  unique = TRUE
)
summary_dt <- merge(
  complete_grid,
  summary_dt,
  by = c("tissue", "rna_bias", "compare"),
  all.x = TRUE,
  sort = FALSE
)
summary_dt[is.na(N), N := 0L]
summary_dt[, total_N := sum(N), by = .(tissue, rna_bias)]
summary_dt[, fraction := fifelse(total_N > 0, N / total_N, NA_real_)]
summary_dt[, percent := fraction * 100]
summary_dt[, compare_label := fifelse(compare == "same_class", "Concordance", "Discordance")]
summary_dt[, compare_label := factor(compare_label, levels = c("Concordance", "Discordance"))]
summary_dt[, tissue_plot := factor(tissue, levels = rev(tissue_order))]
summary_dt[, label := fifelse(
  !is.na(fraction) & fraction >= 0.08,
  sprintf("%.1f%%\n(n=%d)", percent, N),
  ""
)]

plot_one_bias <- function(bias_value, title, out_path) {
  pdt <- copy(summary_dt[rna_bias == bias_value])
  pdt[, total_label := sprintf("%s (N=%d)", tissue, total_N)]
  total_labels <- unique(pdt[, .(tissue_plot, total_label)])
  label_map <- setNames(total_labels$total_label, total_labels$tissue_plot)
  p <- ggplot(pdt, aes(x = tissue_plot, y = fraction, fill = compare_label)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.35, na.rm = TRUE) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3.4,
      lineheight = 0.9,
      color = "white",
      fontface = "bold",
      na.rm = TRUE
    ) +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = label_map) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = function(x) paste0(round(x * 100), "%"),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_fill_manual(
      values = c("Concordance" = "#faa638", "Discordance" = "#4476af"),
      drop = FALSE
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Percentage of genes",
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.margin = margin(10, 20, 10, 10)
    )

  ggsave(out_path, p, width = 8.5, height = 5)
}

maternal_pdf <- file.path(
  work_dir,
  "CsCi_MTB_RNA_ATAC_rna_maternal_bias_compare_fraction_by_tissue_total_count_gt30.pdf"
)
paternal_pdf <- file.path(
  work_dir,
  "CsCi_MTB_RNA_ATAC_rna_paternal_bias_compare_fraction_by_tissue_total_count_gt30.pdf"
)

plot_one_bias(
  "maternal",
  "RNA maternal-biased common expressed genes: Concordance vs Discordance by tissue",
  maternal_pdf
)
plot_one_bias(
  "paternal",
  "RNA paternal-biased common expressed genes: Concordance vs Discordance by tissue",
  paternal_pdf
)