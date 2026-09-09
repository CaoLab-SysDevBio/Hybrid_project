# ==============================================================================
# 01. species and parent-of-origin associated biases (Fig.3a)
# ==============================================================================
CiCs_ai <- read.csv("CiCs_pc20_posterior_mean.csv", row.names = 1)
CiCs_ai_f_clean <- CiCs_ai[ , colSums(is.na(CiCs_ai)) == 0]
df_CiCs_ai_f_clean <- CiCs_ai_f_clean
df_CiCs_ai_f_clean <- as.data.frame(lapply(df_CiCs_ai_f_clean, function(x) {
    pmin(pmax(x, 0), 1)
}))
any(is.na(df_CiCs_ai_f_clean))
rownames(df_CiCs_ai_f_clean) <- rownames(CiCs_ai_f_clean)
write.csv(df_CiCs_ai_f_clean, "CiCs_pc20_posterior_mean_processed.csv")

cols_all_gt_05 <- sapply(df_CiCs_ai_f_clean, function(x) all(x > 0.5))
df_gt_05 <- df_CiCs_ai_f_clean[, cols_all_gt_05]
df_gt_05_m <- colMeans(df_gt_05)
df_gt_05_m <- as.data.frame(df_gt_05_m)
colnames(df_gt_05_m) <- "avg_allelic_rate"
df_gt_05_m$class <- "maternal_bias"

cols_all_lt_05 <- sapply(df_CiCs_ai_f_clean, function(x) all(x < 0.5))
df_lt_05 <- df_CiCs_ai_f_clean[, cols_all_lt_05]
df_lt_05_m <- colMeans(df_lt_05)
df_lt_05_m <- as.data.frame(df_lt_05_m)
colnames(df_lt_05_m) <- "avg_allelic_rate"
df_lt_05_m$class <- "paternal_bias"
df_CiCs_all <- rbind(df_gt_05_m, df_lt_05_m)
df_CiCs_all$OG <- rownames(df_CiCs_all)
rownames(df_CiCs_all) <- NULL
colnames(df_CiCs_all)[1:2] <- c("CiCs_avg_ai_rate", "CiCs_bias")

## CsCi
CsCi_ai <- read.csv("CsCi_pc20_posterior_mean.csv", row.names = 1)
CsCi_ai_clean <- CsCi_ai[ , colSums(is.na(CsCi_ai)) == 0]
df_CsCi_ai_clean <- CsCi_ai_clean
df_CsCi_ai_clean <- as.data.frame(lapply(df_CsCi_ai_clean, function(x) {
    pmin(pmax(x, 0), 1)
}))
any(is.na(df_CsCi_ai_clean))
df_CsCi_ai_clean[1:5, 1:5]
rownames(df_CsCi_ai_clean) <- rownames(CsCi_ai_clean)
write.csv(df_CsCi_ai_clean, "CsCi_pc20_posterior_mean_processed.csv")

cols_all_gt_05_r <- sapply(df_CsCi_ai_clean, function(x) all(x > 0.5))
df_gt_05_r <- df_CsCi_ai_clean[, cols_all_gt_05_r]
df_gt_05_r_m <- colMeans(df_gt_05_r)
df_gt_05_r_m <- as.data.frame(df_gt_05_r_m)
colnames(df_gt_05_r_m) <- "avg_allelic_rate"
df_gt_05_r_m$class <- "maternal_bias"

cols_all_lt_05_r <- sapply(df_CsCi_ai_clean, function(x) all(x < 0.5))
df_lt_05_r <- df_CsCi_ai_clean[, cols_all_lt_05_r]
df_lt_05_r_m <- colMeans(df_lt_05_r)
df_lt_05_r_m <- as.data.frame(df_lt_05_r_m)
colnames(df_lt_05_r_m) <- "avg_allelic_rate"
df_lt_05_r_m$class <- "paternal_bias"

df_CsCi_all <- rbind(df_gt_05_r_m, df_lt_05_r_m)
df_CsCi_all$OG <- rownames(df_CsCi_all)
rownames(df_CsCi_all) <- NULL
colnames(df_CsCi_all)[1:2] <- c("CsCi_avg_ai_rate", "CsCi_bias")

df_2cross <- merge(df_CsCi_all, df_CiCs_all, by = "OG")
library(dplyr)
df <- df_2cross
df <- df %>%
  mutate(class = case_when(
    CsCi_bias == "maternal_bias" & CiCs_bias == "maternal_bias" ~ "maternal_specific",
    CsCi_bias == "maternal_bias" & CiCs_bias == "paternal_bias" ~ "Ci_specific",
    CsCi_bias == "paternal_bias" & CiCs_bias == "paternal_bias" ~ "paternal_specific",
    CsCi_bias == "paternal_bias" & CiCs_bias == "maternal_bias" ~ "Cs_specific",
    TRUE ~ NA_character_
  ))
write.csv(df, "CsCi_CiCs_specific_class.csv", row.names = FALSE)

### plot ###
library(ggplot2)
quad_colors <- c(
"Cs_specific" = "#D55E00",
"Ci_specific" = "#0072B2",
"paternal_specific"= "#009E73",
"maternal_specific" = "#CC79A7"
)
df_plot <- df
df_pie <- df_plot %>%
count(class, name = "n") %>%
    mutate(
    freq = n / sum(n),
    pct  = freq * 100,
    label = paste0(round(pct, 1), "%")
    ) %>%
    arrange(desc(class)) %>%
    mutate(
    ypos = cumsum(freq) - 0.5 * freq
    )
pdf("CsCi_CiCs_specific_class_nmlt_counts.cutoff.30_pieplot.pdf", width = 6, height = 5)
    ggplot(df_pie, aes(x = "", y = freq, fill = class)) +
        geom_bar(
        stat = "identity",
        width = 1,
        color = "white"
        ) +
        coord_polar(theta = "y") +
        scale_fill_manual(
        values = quad_colors,
        drop = FALSE
        ) +
        geom_text(
        aes(y = ypos, label = label),
        color = "white",
        size = 4
        ) +
        theme_void() +
        labs(fill = "specific_effects")
dev.off()



# ==============================================================================
# 02. variance partition (Fig.3b)
# ==============================================================================
############################# varPar###################################
CsCi.rate <- read.csv("CsCi_allelic_rate_count>30_nmlt.csv", row.names = 1, check.names = FALSE)
gene.anno.CsCi <- read.csv("CsCi_gene_anno.csv")

library(variancePartition)
library(ggplot2)
library(dplyr)
epsilon <- 1e-3  # avoid 0 or 1
logit_rate_mat <- log((CsCi.rate + epsilon) / (1 - CsCi.rate + epsilon))
rownames(gene.anno.CsCi) <- gene.anno.CsCi$gene
gene.df.CsCi <- gene.anno.CsCi
gene.df.CsCi[, c(2:5,7:25)] <- lapply(gene.df.CsCi[, c(2:5,7:25)], factor) 
form <- ~ avg_exp + (1|tissue.marker) + (1|TF) + (1|species_effect_all) + (1|parental_effect_all)
all(colnames(logit_rate_mat) == rownames(gene.df.CsCi))
run_varPart_pipeline <- function(logit_rate_mat, form, gene.df,
                                 out_csv,
                                 out_pdf,
                                 pdf_width = 12,
                                 pdf_height = 6,
                                 y_max = 100) {
  varPart <- fitExtractVarPartModel(logit_rate_mat, form, gene.df)
  write.csv(varPart, out_csv, row.names = TRUE)
  p <- plotVarPart(sortCols(varPart)) +
    ylim(0, y_max) +
    ggtitle("Variance Partition") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13),
    legend.position = "none")
  pdf(out_pdf, width = pdf_width, height = pdf_height)
  print(p)
  dev.off()
  return(varPart)
}
varPart <- run_varPart_pipeline(
  logit_rate_mat = logit_rate_mat,
  form = form,
  gene.df = gene.df.CsCi,
  out_csv = "CsCi_rate_varPart.csv",
  out_pdf = "CsCi_rate_varPart.pdf",
  y_max = 70
)



# ==============================================================================
# 03. cis- and trans-associated regulatory effects (Fig.3e and Ext.Data.Fig.11)
# ==============================================================================
library(dplyr)
library(tidyr)
library(Matrix)
library(Seurat)
library(ggplot2)
combined.10hpf.all <- readRDS("combined.10hpf.all.anno.rds")
counts <- GetAssayData(combined.10hpf.all, layer = "counts")
meta <- combined.10hpf.all@meta.data
datasets <- c("H_Ci", "H_Cs", "P_Ci", "P_Cs")
count_sum_list <- lapply(datasets, function(ds) {
  cells <- rownames(meta)[meta$dataset == ds]
  rowSums(counts[, cells, drop = FALSE])
})

count_sum_mat <- do.call(cbind, count_sum_list)
colnames(count_sum_mat) <- datasets
write.csv(count_sum_mat, "combined.10hpf.gene.sum.counts.csv")

# CPM normalization
cpm_mat <- apply(count_sum_mat, 2, function(x) {
  x / sum(x) * 1e6
})
cpm_mat <- as.data.frame(cpm_mat)
colSums(cpm_mat)

# filter low-count gene
keep_genes <- rownames(count_sum_mat)[
  (count_sum_mat[, "H_Ci"] + count_sum_mat[, "H_Cs"] > 30) &
  (count_sum_mat[, "P_Ci"] + count_sum_mat[, "P_Cs"] > 30)
]

cpm_filtered <- cpm_mat[keep_genes, ]
write.csv(cpm_filtered, "combined.10hpf.gene.totalCPM.exp.csv")

# calculate cis, trans, total
pseudo <- 1 
cis <- (cpm_filtered$H_Ci + pseudo) / (cpm_filtered$H_Cs + pseudo)
total <- (cpm_filtered$P_Ci + pseudo) / (cpm_filtered$P_Cs + pseudo)
trans <- total / cis
log2_cis <- log2(cis)
log2_total <- log2(total)
log2_trans <- log2(trans)

result <- data.frame(
  gene = rownames(cpm_filtered),
  cis = cis,
  total = total,
  trans = trans,
  log2_cis = log2_cis,
  log2_total = log2_total,
  log2_trans = log2_trans
)
write.csv(result, "10hpf.gene.cis.trans.log2.CPM.csv", row.names = FALSE)

# classification
result$class <- "uncategorized"
cis_mid   <- result$cis   >= 0.8 & result$cis   <= 1.2
trans_mid <- result$trans >= 0.8 & result$trans <= 1.2
total_mid <- result$total >= 0.8 & result$total <= 1.2
cis_out   <- result$cis   < 0.8 | result$cis   > 1.2
trans_out <- result$trans < 0.8 | result$trans > 1.2
total_out <- result$total < 0.8 | result$total > 1.2
# conserved → cis_only → trans_only → interaction
result$class[cis_mid & total_mid] <- "conserved"
result$class[cis_out & trans_mid & result$class == "uncategorized"] <- "cis_only"
result$class[cis_mid & total_out & result$class == "uncategorized"] <- "trans_only"
result$class[cis_out & trans_out & result$class == "uncategorized"] <- "cis_trans_interaction"
write.csv(result, "10hpf.gene.cis.trans.log2.CPM.class.0.8_1.2.csv", row.names = FALSE)

### after excluding maternal loaded transcript ###
result.nmlt <- subset(result, result$mlt == "No")
class_colors <- c(
  "conserved" = "grey70",
  "cis_only" = "#f24242",
  "trans_only" = "#1391eb",
  "cis_dominant" = "#eeac4f",
  "trans_dominant" = "#b745ef"
)
pdf("cis_trans_dominant_classification_plots.nmlt.gene.pdf", width = 7, height = 7)
    # Pie chart
    pie_data <- bar_data %>%
    mutate(percentage = count / sum(count) * 100,
            label = paste0(class_detailed, " (", round(percentage,1), "%)"))
    print(
    ggplot(pie_data, aes(x = "", y = percentage, fill = class_detailed)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y") +
        theme_void() +
        ggtitle("Gene proportion per cis_trans category") +
        geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values = class_colors)
    )

    # Scatter plot log2cis vs log2trans
    print(ggplot(result.nmlt, aes(x = log2_cis, y = log2_trans, color = class_detailed)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = class_colors) +
    theme_minimal() +
    ggtitle("log2cis vs log2trans") +
    xlab("log2(cis)") +
    ylab("log2(trans)")+
    theme(legend.title = element_blank()))
dev.off()



# ==============================================================================
# 04. correlation of scRNA and scATAC allelic rates (Fig.3h)
# ==============================================================================
library(data.table)
library(ggplot2)
work_dir <- "/MTB_ATAC_RNA/"
input_path <- file.path(work_dir, "CsCi_MTB_RNA_ATAC_common_genes_mean_allelic_rate.csv")
out_pdf <- file.path(work_dir, "CsCi_MTB_RNA_ATAC_common_genes_mean_allelic_rate_contour.pdf")
stats_csv <- file.path(work_dir, "CsCi_MTB_RNA_ATAC_common_genes_mean_allelic_rate_spearman_statistics.csv")
dt <- fread(input_path)
setnames(dt, 1, "gene")
needed <- c("gene", "rna.mean", "atac.mean")
missing <- setdiff(needed, names(dt))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))
dt[, rna.mean := as.numeric(rna.mean)]
dt[, atac.mean := as.numeric(atac.mean)]
dt <- dt[is.finite(rna.mean) & is.finite(atac.mean)]

ct <- suppressWarnings(cor.test(dt$rna.mean, dt$atac.mean, method = "spearman", exact = FALSE))
rho <- unname(ct$estimate)
pval <- ct$p.value
sig <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns")))
p_text <- if (pval < 2.2e-16) "P < 2.2e-16" else paste0("P = ", formatC(pval, format = "e", digits = 2))
stat_text <- sprintf("Spearman rho = %.3f; %s (%s); N = %d", rho, p_text, sig, nrow(dt))

fwrite(
  data.table(
    comparison = "RNA_ATAC_mean_allelic_rate",
    N = nrow(dt),
    spearman_rho = rho,
    p_value = pval,
    significance = sig
  ),
  stats_csv
)

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

p_contour <- ggplot(dt, aes(x = rna.mean, y = atac.mean)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    bins = 7,
    alpha = 0.78,
    color = NA
  ) +
  geom_point(color = "black", alpha = 0.35, size = 0.45) +
  geom_smooth(method = "lm", se = TRUE, color = "#E66101", fill = "grey75", linewidth = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  scale_fill_gradientn(colors = c("#ffff33", "#fdae61", "#d53e4f", "#3f007d"), name = "Density") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    title = "Density contour",
    subtitle = stat_text,
    x = "Mean RNA allelic rate",
    y = "Mean ATAC allelic rate"
  ) +
  base_theme +
  theme(legend.position = "right")
ggsave(out_pdf, p_contour, width = 6, height = 5)



# ========================================================================================================
# 05. the odds-ratio analysis (Fig.3c,3d,3f,3g,3i,3j)
# We used the relationships of species/parental biases and transition categories (Fig.3c-d) as an example.
# ========================================================================================================
### Odds-ratio bubble plot ###
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
base_dir <- "/midG_all_lineage/"
stats_file <- file.path(
  base_dir,
  "midG_all_lineage_transition_vs_specific_class_stats.csv"
)
out_dir <- file.path(base_dir, "species_vs_parental")

transition_levels <- c("switched", "propagated", "attenuated", "buffered", "exposed")
specific_levels <- c(
  "Ci_specific",
  "Cs_specific",
  "maternal_specific",
  "paternal_specific",
  "unclassified"
)

score_defs <- list(
  "species specific propagation" = list(
    associated_classes = c("Ci_specific", "Cs_specific"),
    target_categories = c("propagated")
  ),
  "species specific change" = list(
    associated_classes = c("Ci_specific", "Cs_specific"),
    target_categories = c("switched", "exposed", "attenuated")
  ),
  "parental specific propagation" = list(
    associated_classes = c("maternal_specific", "paternal_specific"),
    target_categories = c("propagated")
  ),
  "parental specific change" = list(
    associated_classes = c("maternal_specific", "paternal_specific"),
    target_categories = c("switched", "exposed", "attenuated")
  )
)
score_levels <- names(score_defs)
pstar <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

calc_or <- function(mat, node_meta, score_name, def) {
  assoc_rows <- rownames(mat) %in% def$associated_classes
  target_cols <- colnames(mat) %in% def$target_categories

  a <- sum(mat[assoc_rows, target_cols, drop = FALSE], na.rm = TRUE)
  b <- sum(mat[assoc_rows, !target_cols, drop = FALSE], na.rm = TRUE)
  c <- sum(mat[!assoc_rows, target_cols, drop = FALSE], na.rm = TRUE)
  d <- sum(mat[!assoc_rows, !target_cols, drop = FALSE], na.rm = TRUE)

  odds_ratio <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))

  data.frame(
    node = node_meta$node,
    tissue = node_meta$tissue,
    lineage = node_meta$lineage,
    parent_node = node_meta$parent_node,
    score = score_name,
    associated_classes = paste(def$associated_classes, collapse = ";"),
    target_transition_category = paste(def$target_categories, collapse = ";"),
    intersection_n = a,
    odds_ratio = odds_ratio,
    log2_odds_ratio = log2(odds_ratio),
    fisher_p_value = ft$p.value,
    neg_log10_fisher_p = -log10(max(ft$p.value, .Machine$double.xmin)),
    stringsAsFactors = FALSE
  )
}

stats <- read.csv(stats_file, check.names = FALSE) %>%
  mutate(
    node = as.character(node),
    tissue = as.character(tissue),
    lineage = as.character(lineage),
    parent_node = as.character(parent_node),
    specific_class = factor(specific_class, levels = specific_levels),
    transition_category = factor(transition_category, levels = transition_levels),
    intersection_n = as.numeric(intersection_n)
  )

node_meta <- stats %>%
  distinct(node, tissue, lineage, parent_node) %>%
  arrange(tissue, lineage, node)
node_levels <- node_meta$node
score_rows <- list()
for (i in seq_len(nrow(node_meta))) {
  meta <- node_meta[i, , drop = FALSE]
  mat_df <- stats %>%
    filter(node == meta$node) %>%
    select(specific_class, transition_category, intersection_n) %>%
    pivot_wider(
      names_from = transition_category,
      values_from = intersection_n,
      values_fill = 0
    )
  mat <- as.matrix(mat_df[, transition_levels, drop = FALSE])
  rownames(mat) <- as.character(mat_df$specific_class)
  mat <- mat[specific_levels, transition_levels, drop = FALSE]
  storage.mode(mat) <- "numeric"
  score_rows[[i]] <- bind_rows(lapply(score_levels, function(score_name) {
    calc_or(mat, meta, score_name, score_defs[[score_name]])
  }))
}
score_df <- bind_rows(score_rows) %>%
  mutate(
    score = factor(score, levels = score_levels),
    node = factor(node, levels = node_levels),
    sig_label = pstar(fisher_p_value)
  )
write.csv(
  score_df,
  file.path(out_dir, "midG_all_lineage_species_vs_parental_OR_bubble_data.csv"),
  row.names = FALSE
)

bubble <- ggplot(score_df, aes(x = node, y = score)) +
  geom_point(
    aes(color = log2_odds_ratio, size = neg_log10_fisher_p),
    alpha = 0.9
  ) +
  geom_text(
    data = score_df %>% filter(sig_label != ""),
    aes(label = sig_label),
    color = "black",
    size = 4,
    fontface = "bold",
    vjust = 0.5
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-4, 4),
    oob = squish,
    name = "log2(OR)"
  ) +
  scale_size_continuous(
    range = c(2.2, 8.8),
    name = "-log10(Fisher p)"
  ) +
  scale_y_discrete(limits = rev(score_levels), drop = FALSE) +
  labs(
    x = NULL,
    y = NULL,
    caption = "* p < 0.05    ** p < 0.01    *** p < 0.001"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.35),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black", size = 9),
    plot.caption = element_text(hjust = 0, size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )
ggsave(
  file.path(out_dir, "midG_all_lineage_species_vs_parental_OR_bubble_pstar_marked.pdf"),
  bubble,
  width = 10,
  height = 5,
  useDingbats = FALSE
)



### Odds-ratio violin plot ###
library(ggplot2)
library(dplyr)
library(tidyr)
in_file <- "/midG_all_lineage_species_vs_parental_score_OR.csv"
out_stage <- file.path("/midG_species_vs_parental/")

score_levels <- c(
  "species specific propagation",
  "species specific change",
  "parental specific propagation",
  "parental specific change"
)

score_labels <- c(
  "species\npropagation",
  "species\nchange",
  "parental\npropagation",
  "parental\nchange"
)
names(score_labels) <- score_levels

pair_colors <- c(
  "species specific propagation" = "#C8586B",
  "species specific change" = "#F2A07D",
  "parental specific propagation" = "#5A9BD4",
  "parental specific change" = "#9AD4E8"
)

pstar <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", signif(p, 3))
}

paired_p <- function(dat, group_a, group_b) {
  wide <- dat %>%
    filter(score %in% c(group_a, group_b)) %>%
    select(node, score, log2_odds_ratio) %>%
    pivot_wider(names_from = score, values_from = log2_odds_ratio) %>%
    filter(!is.na(.data[[group_a]]), !is.na(.data[[group_b]]))
  if (nrow(wide) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(wide[[group_a]], wide[[group_b]], paired = TRUE, exact = FALSE)$p.value)
}

df <- read.csv(in_file, check.names = FALSE) %>%
  mutate(
    node = as.character(node),
    score = factor(score, levels = score_levels),
    log2_odds_ratio = ifelse(is.finite(log2_odds_ratio), log2_odds_ratio, NA_real_)
  ) %>%
  filter(!is.na(log2_odds_ratio))

species_p <- paired_p(df, "species specific propagation", "species specific change")
parental_p <- paired_p(df, "parental specific propagation", "parental specific change")

y_max <- max(df$log2_odds_ratio, na.rm = TRUE)
y_min <- min(df$log2_odds_ratio, na.rm = TRUE)
y_span <- max(y_max - y_min, 1)
anno_df <- data.frame(
  x = c(1.5, 3.5),
  x_start = c(1, 3),
  x_end = c(2, 4),
  y = c(y_max + 0.12 * y_span, y_max + 0.29 * y_span),
  label = c(
    paste0(format_p(species_p), "  ", pstar(species_p)),
    paste0(format_p(parental_p), "  ", pstar(parental_p))
  )
)

p <- ggplot(df, aes(x = score, y = log2_odds_ratio, fill = score)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35, linetype = "dashed") +
  geom_violin(
    trim = FALSE,
    width = 0.68,
    color = "grey20",
    linewidth = 0.32,
    alpha = 0.94
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    color = "grey20",
    fill = "white",
    linewidth = 0.28
  ) +
  geom_segment(
    data = anno_df,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_segment(
    data = anno_df,
    aes(x = x_start, xend = x_start, y = y, yend = y - 0.035 * y_span),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_segment(
    data = anno_df,
    aes(x = x_end, xend = x_end, y = y, yend = y - 0.035 * y_span),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_text(
    data = anno_df,
    aes(x = x, y = y + 0.045 * y_span, label = label),
    inherit.aes = FALSE,
    size = 3.1
  ) +
  scale_x_discrete(labels = score_labels, drop = FALSE, expand = expansion(add = 0.56)) +
  scale_fill_manual(values = pair_colors, guide = "none") +
  coord_cartesian(ylim = c(y_min - 0.08 * y_span, y_max + 0.43 * y_span), clip = "off") +
  labs(x = NULL, y = "log2(OR)") +
  theme_classic(base_size = 10) +
  theme(
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.32),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey20", linewidth = 0.35),
    axis.ticks = element_line(color = "grey20", linewidth = 0.3),
    axis.text.x = element_text(color = "black", size = 9.5),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(size = 10.5),
    plot.margin = margin(9, 6, 4, 6)
  )

ggsave(
  file.path(out_stage, "midG_all_lineage_species_vs_parental_log2OR_violin.pdf"),
  p,
  width = 4.7,
  height = 3.6,
  useDingbats = FALSE
)
write.csv(
  data.frame(
    comparison = c("species propagation vs species change", "parental propagation vs parental change"),
    test = "paired Wilcoxon by node",
    p_value = c(species_p, parental_p),
    p_label = pstar(c(species_p, parental_p))
  ),
  file.path(out_stage, "midG_all_lineage_species_vs_parental_log2OR_violin_pvalues.csv"),
  row.names = FALSE
)