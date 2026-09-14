# ============================================================
# 01. allelic bias score (Fig.1d, 1e)
# ============================================================
library(ggplot2)
library(ggforce)
library(dplyr)
# all stage cells
df <- read.csv("/CsCi.10s.allelic.rates.1e-3_1e-300.processed.nmlt.anno.tissue.beau.adj.csv", row.names = 1) # cell x gene
mat <- as.matrix(df[, -ncol(df)])
group <- df$cell_type
sum_mat <- rowsum(mat, group)
avg_mat <- sum_mat / as.vector(table(group))
# binary
binary_mat <- avg_mat
binary_mat[binary_mat > 0.5] <- 1
binary_mat[binary_mat < 0.5] <- -1
gene_mean <- colMeans(binary_mat, na.rm = TRUE)
df_plot <- data.frame(gene_mean = gene_mean)
df_plot <- df_plot %>%
  mutate(class = case_when(
    gene_mean == 1  ~ "Stable_maternal_bias",
    gene_mean == -1 ~ "Stable_paternal_bias",
    TRUE ~ "Switching_allelic_bias"
  ))
write.csv(df_plot, "10s.allelic.score.df.plot.csv")
# pie chart
pie_df <- df_plot %>%
  count(class) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0("n = ", n, "\n", sprintf("%.1f%%", percent))
  )
plot_colors <- c(
  Stable_maternal_bias = "#fa7569",
  Stable_paternal_bias = "#82a9fc",
  Switching_allelic_bias = "#fcb95c"
)
pdf("10s.allelic.celltype.class.pie.plot.pdf", width = 6, height = 5)
ggplot(pie_df, aes(x = "", y = n, fill = class)) +
  geom_bar(stat = "identity", width = 1, color = "white", alpha = 0.8) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = plot_colors) +
  theme_void() +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4) +
  labs(title = "Distribution of allelic bias category") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# only plot switching class
df_plot_change <- subset(df_plot, df_plot$class == "Switching_allelic_bias")
xmin <- min(df_plot_change$gene_mean)
xmax <- max(df_plot_change$gene_mean)
dens <- density(df_plot_change$gene_mean)
ymax <- max(dens$y)
pdf("10s.allelic.hetero.density.plot.pdf", width = 6, height = 4)
ggplot(df_plot_change, aes(x = gene_mean)) +
  geom_density(fill = "#fcb95c", alpha = 0.8) +
  theme_classic() +
  xlab("Allelic bias score") +
  ylab("Density") +
  labs(title = "Allelic bias switch across cell types") +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, 1.2)) +
  scale_x_continuous(
    breaks = c(
      xmin,
      pretty(c(xmin, xmax), n = 5),
      xmax
    )
  )
dev.off()



# ============================================================
# 02. allelic marker heatmap (Fig.1f)
# ============================================================
### find top 10 allelic markers ###
df <- read.csv("CsCi.10s.allelic.rates.1e-3_1e-300.processed.nmlt.anno.tissue.beau.adj.csv", row.names = 1)
cell_type <- df[, ncol(df)]
names(cell_type) <- rownames(df)
mat <- df[, -ncol(df)]

# create function for marker identification
get_markers_one_ct <- function(mat, cell_type, ct, top_n = 10,
                              min_diff = 0.1, padj_cutoff = 0.05) {
  cells_in  <- names(cell_type)[cell_type == ct]
  cells_out <- names(cell_type)[cell_type != ct]
  genes <- colnames(mat)
  pvals <- numeric(length(genes))
  diff  <- numeric(length(genes))
  direction <- character(length(genes))
  for (i in seq_along(genes)) {
    g <- genes[i]
    x <- mat[cells_in, g]
    y <- mat[cells_out, g]
    # Wilcoxon test
    wt <- wilcox.test(x, y)
    pvals[i] <- wt$p.value
    diff[i] <- median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
    direction[i] <- ifelse(diff[i] > 0, "greater", "less")
  }
  res <- data.frame(
    gene = genes,
    pval = pvals,
    padj = p.adjust(pvals, method = "BH"),
    diff = diff,
    direction = direction
  )
  res <- res[
    res$padj < padj_cutoff & abs(res$diff) > min_diff,
  ]
  if (nrow(res) == 0) return(NULL) # avoid error
  res <- res[order(-abs(res$diff), res$padj), ] # first order as padj, second order as diff
  res_top <- head(res, top_n)
  res_top$cell_type <- ct
  return(res_top)
}

# each cell type
celltypes <- unique(cell_type)
marker_list <- lapply(celltypes, function(ct) {
  get_markers_one_ct(mat, cell_type, ct,
                     top_n = 10,
                     min_diff = 0.25,
                     padj_cutoff = 0.05)
})
marker_list <- Filter(Negate(is.null), marker_list)
marker_df <- do.call(rbind, marker_list)
write.csv(marker_df, "CsCi.10s.celltype_specific_allelic_top10_markers.csv", row.names = FALSE)

### plot ###
df <- read.csv("CsCi.10s.allelic.rates.1e-3_1e-300.processed.nmlt.anno.tissue.beau.adj.csv", row.names = 1)
cell_type <- df[, ncol(df)]
names(cell_type) <- rownames(df)
mat <- df[, -ncol(df)]

marker_df <- read.csv("CsCi.10s.celltype_specific_allelic_top10_markers.csv")
length(unique(marker_df$gene))
marker_genes <- unique(marker_df$gene)
mat_marker <- mat[, marker_genes, drop = FALSE]
mat_marker <- t(mat_marker)

celltype_order <- c( "muscle&heart", "notochord", "endoderm", "germ", "mesenchyme", "nervous_system", "epidermis")
cell_order <- order(factor(cell_type, levels = celltype_order))
mat_marker <- mat_marker[, cell_order, drop = FALSE]

cell_type_df <- as.data.frame(cell_type)
head(cell_type_df)
cell_type_df <- cell_type_df[cell_order, , drop = FALSE]
all(rownames(cell_type_df) == colnames(mat_marker))
annotation_col <- data.frame(
  Cell.type = cell_type_df[,1]
)
rownames(annotation_col) <- rownames(cell_type_df)

ann_colors <- list(
  Cell.type = setNames(
    colorRampPalette(c(
      "#f49852", "#4bc051", "#66bfc1", "#386236", "#e9db43", "#41b3f1", "#7966a2"))
      (length(celltype_order)), 
    celltype_order
  ))
gene.homolog <- read.csv("OG.pairs_human.homolog.csv")
gene_map <- gene.homolog[, c("Ci_OG", "gene_homolog")]
rownames(gene_map) <- gene_map$Ci_OG
labels_row <- gene_map[rownames(mat_marker), "gene_homolog"]
labels_row[is.na(labels_row)] <- "KH.C1.391"

library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)
cell_type_factor <- factor(
  annotation_col$Cell.type,
  levels = celltype_order
)
celltype_colors <- ann_colors$Cell.type[celltype_order]
allelic_colors <- colorRamp2(
  c(0,   0.1,       0.2,       0.3,       0.4,       0.5,
    0.6, 0.7,       0.8,       0.9,       1),
  c("#4575B4", "#6699C5", "#91BFDB", "#BFDDED", "#E5F3EE", "#FFFFBF",
    "#FEE8A5", "#FDC980", "#F88D5A", "#E6533A", "#C92525")
)
top_anno <- HeatmapAnnotation(
  `Cell type` = anno_block(
    gp = gpar(
      fill = celltype_colors,
      col = NA
    ),
    labels = celltype_order,
    labels_rot = 10,
    labels_gp = gpar(
      fontsize = 10,
      fontface = "bold"
    )
  ),
  annotation_name = NULL,
  show_legend = FALSE,
  height = unit(5, "mm")
)

pdf(
  "allelic.top10.gene.padj.diff.heatmap.clusteing.label.pdf",
  width = 15,
  height = 10,
  useDingbats = FALSE
)

ht <- Heatmap(
  mat_marker,
  name = "Estimated allelic rate",
  col = allelic_colors,
  na_col = "#808080",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = cell_type_factor,
  cluster_column_slices = FALSE,
  column_gap = unit(0, "mm"),
  column_title = NULL,
  top_annotation = top_anno,
  row_labels = labels_row,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(
    fontsize = 12,
    fontface = "italic"
  ),
  rect_gp = gpar(
    col = NA,
    lwd = 0
  ),
  # allelic rate legend
  heatmap_legend_param = list(
  title = "Estimated allelic rate",
  title_position = "topcenter",
  legend_direction = "vertical",
  at = seq(0, 1, by = 0.2),
  labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")
  ),
  use_raster = TRUE,
  raster_device = "png",
  raster_quality = 4
)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()


# ============================================================
# 03. GO enrichment (Fig.1g, 1h)
# ============================================================
### find markers
library(data.table)
base_dir <- "allelic_rates_GO"
marker_dir <- file.path(base_dir, "greater_less_markers")
input_file <- file.path(base_dir, "CsCi.10s.allelic.rates.1e-3_1e-300.processed.nmlt.anno.tissue.beau.adj.csv")
output_dir <- file.path(marker_dir, "median_rate_filtered")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("^_+|_+$", "", x)
}
parse_marker_file <- function(path) {
  nm <- basename(path)
  m <- regexec("^(.*)_(greater|less)_markers_minDiff0\\.25_padj0\\.05\\.csv$", nm)
  hit <- regmatches(nm, m)[[1]]
  if (length(hit) != 3L) {
    stop("Cannot parse marker file name: ", nm)
  }
  list(cell_type_safe = hit[2], direction = hit[3])
}
message("Finding marker files")
marker_files <- list.files(
  marker_dir,
  pattern = "_(greater|less)_markers_minDiff0\\.25_padj0\\.05\\.csv$",
  full.names = TRUE
)
marker_files <- marker_files[!grepl("^all_", basename(marker_files))]
marker_files <- marker_files[order(basename(marker_files))]
if (length(marker_files) == 0L) {
  stop("No marker files found in: ", marker_dir)
}
message("Marker files: ", length(marker_files))

marker_tables <- lapply(marker_files, fread)
names(marker_tables) <- vapply(marker_files, function(f) sub("\\.csv$", "", basename(f)), character(1))
all_marker_genes <- unique(unlist(lapply(marker_tables, function(x) x$gene)))
message("Unique marker genes: ", length(all_marker_genes))

message("Reading header from original rate matrix")
header <- names(fread(input_file, nrows = 0))
gene_cols <- intersect(all_marker_genes, header)
missing_genes <- setdiff(all_marker_genes, gene_cols)
if (length(missing_genes) > 0L) {
  fwrite(data.table(gene = missing_genes), file.path(output_dir, "marker_genes_missing_from_rate_matrix.csv"))
  warning("Marker genes missing from rate matrix: ", length(missing_genes))
}

message("Reading original rate matrix columns: cell_type + ", length(gene_cols), " marker genes")
dt <- fread(input_file, select = c("cell_type", gene_cols), showProgress = TRUE)
celltypes <- sort(unique(dt$cell_type))
message("Cell types in matrix: ", paste(celltypes, collapse = ", "))

message("Computing per-cell-type median rates")
median_list <- vector("list", length(celltypes))
for (i in seq_along(celltypes)) {
  ct <- celltypes[i]
  message("  median: ", ct)
  idx <- dt$cell_type == ct
  med <- vapply(gene_cols, function(g) median(dt[[g]][idx], na.rm = TRUE), numeric(1))
  median_list[[i]] <- data.table(cell_type = ct, gene = gene_cols, median_rate_in_cell_type = med)
}
median_dt <- rbindlist(median_list)
median_file <- file.path(output_dir, "marker_gene_median_rates_by_cell_type.csv")
fwrite(median_dt, median_file)
message("Wrote medians: ", median_file)

summary_rows <- list()
filtered_tables <- list()

for (marker_file in marker_files) {
  meta <- parse_marker_file(marker_file)
  markers <- fread(marker_file)
  if (!all(c("gene", "direction", "cell_type") %in% names(markers))) {
    stop("Required columns gene/direction/cell_type not found in: ", marker_file)
  }

  ct <- unique(markers$cell_type)
  dir_label <- unique(markers$direction)
  if (length(ct) != 1L || length(dir_label) != 1L) {
    stop("Expected one cell_type and one direction in: ", marker_file)
  }
  if (dir_label != meta$direction) {
    warning("Direction in file name and table differ for: ", marker_file)
  }

  med_sub <- median_dt[cell_type == ct, .(gene, median_rate_in_cell_type)]
  out <- merge(markers, med_sub, by = "gene", all.x = TRUE, sort = FALSE)

  n_before <- nrow(out)
  if (dir_label == "greater") {
    out_filtered <- out[!is.na(median_rate_in_cell_type) & median_rate_in_cell_type > 0.5]
    median_filter <- "median_rate_in_cell_type > 0.5"
    suffix <- "medianRateGT0.5"
  } else if (dir_label == "less") {
    out_filtered <- out[!is.na(median_rate_in_cell_type) & median_rate_in_cell_type < 0.5]
    median_filter <- "median_rate_in_cell_type < 0.5"
    suffix <- "medianRateLT0.5"
  } else {
    stop("Unknown direction: ", dir_label)
  }

  setcolorder(out_filtered, c("gene", setdiff(names(out_filtered), "gene")))
  out_filtered <- out_filtered[order(-abs(diff), padj)]

  out_file <- file.path(
    output_dir,
    sprintf("%s_%s_markers_minDiff0.25_padj0.05_%s.csv", safe_name(ct), dir_label, suffix)
  )
  fwrite(out_filtered, out_file)
  message("Wrote: ", out_file, " (", nrow(out_filtered), "/", n_before, ")")

  key <- sprintf("%s_%s", safe_name(ct), dir_label)
  filtered_tables[[key]] <- out_filtered
  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    cell_type = ct,
    direction = dir_label,
    n_before = n_before,
    n_after = nrow(out_filtered),
    n_removed = n_before - nrow(out_filtered),
    kept_fraction = if (n_before > 0L) nrow(out_filtered) / n_before else NA_real_,
    median_filter = median_filter,
    output_file = out_file
  )
}

summary_dt <- rbindlist(summary_rows)
summary_file <- file.path(output_dir, "median_rate_filtered_marker_summary.csv")
fwrite(summary_dt, summary_file)

all_filtered <- rbindlist(filtered_tables, use.names = TRUE, fill = TRUE, idcol = "marker_set")
all_file <- file.path(output_dir, "all_greater_less_markers_minDiff0.25_padj0.05_medianRateFiltered.csv")
fwrite(all_filtered, all_file)

### GO big category definition
go_categories <- 
list(transcription = c("transcription", "gene expression", "rna polymerase ii", 
"transcription factor", "chromatin", "epigenetic", "regulation of transcription"
), metabolism = c("metabolic", "biosynthetic", "catabolic", "oxidation-reduction", 
"oxidative", "lipid", "fatty acid", "cholesterol", "eicosanoid", 
"lipoprotein", "plasma lipoprotein particle", "carbohydrate", 
"glucose", "amino acid", "energy", "small molecule", "cytochrome", 
"p450", "hydroxylase", "omega-hydroxylase", "homeostasis", "ion homeostasis", 
"sodium", "potassium", "calcium", "chemical homeostasis", "temperature homeostasis", 
"body fluid"), signaling = c("signal transduction", "signaling pathway", 
"cell communication", "receptor signaling", "ligand", "g-protein coupled receptor", 
"wnt", "notch", "hedgehog", "tgf", "bmp", "fgf", "mapk", "erk", 
"jak-stat", "calcium signaling", "calcium ion", "cytosolic calcium", 
"hormone", "hormonal"), phosphorylation = c("phosphorylation", 
"dephosphorylation", "kinase"), development = c("development", 
"developmental process", "differentiation", "cell fate", "pattern specification", 
"morphogenesis", "organogenesis", "growth", "organ growth", "tissue growth", 
"maturation", "remodeling", "circulatory system", "blood circulation", 
"blood vessel", "vascular", "angiogenesis", "vasodilation", "vasoconstriction", 
"blood pressure", "heart", "cardiac", "heart contraction", "heart rate", 
"cardiac muscle", "muscle development", "muscle system", "muscle contraction", 
"striated muscle", "striated muscle contraction", "muscle hypertrophy", 
"hypertrophy", "renal", "kidney", "glomerular", "filtration", 
"urine", "neurogenesis", "bone", "mineralization", "ossification", 
"system process", "multicellular organismal", "multi-multicellular organism", 
"reproductive", "pregnancy", "spermatogenesis", "hair cycle", 
"molting cycle", "tube", "tube size", "tube diameter"), neural_process = c("synaptic", 
"trans-synaptic", "chemical synaptic", "postsynaptic", "presynaptic", 
"neurotransmitter", "glutamate", "gaba", "gabaergic", "cholinergic", 
"serotonin", "neuronal", "neuron", "nervous system", "nerve impulse", 
"action potential", "membrane potential", "sensory perception", 
"temperature stimulus", "pain", "behavior", "circadian", "sleep", 
"wake", "signal release"), translation_rna_processing = c("translation", 
"ribosome", "ribosomal", "mrna processing", "rna splicing", "rrna processing", 
"trna"), cell_cycle = c("cell cycle", "cell division", "mitosis", 
"mitotic", "cell proliferation", "cell population proliferation", 
"dna replication", "chromosome segregation"), dna_repair = c("dna repair", 
"dna damage", "double-strand break", "recombination"), transport = c("transport", 
"transmembrane transport", "vesicle", "vesicle-mediated", "trafficking", 
"import", "export", "uptake", "reuptake", "secretion", "secretion by cell", 
"exocytosis", "endocytosis", "localization"), protein_turnover = c("proteolysis", 
"protein catabolic", "ubiquitin", "proteasome", "protein processing", 
"protein maturation", "protein modification", "protein folding", 
"protein complex", "protein assembly", "protein trimerization", 
"protein homotrimerization", "protein oligomerization", "glycosylation", 
"fucosylation", "amyloid", "protein polymerization"), cell_adhesion = c("cell adhesion", 
"cell-substrate adhesion", "substrate adhesion", "extracellular matrix", 
"extracellular structure", "external encapsulating structure", 
"collagen", "collagen fibril", "matrix organization", "cell junction", 
"cell migration", "cell motility", "locomotion", "smooth muscle cell migration", 
"vascular associated smooth muscle cell migration", "leukocyte migration", 
"cell-matrix adhesion"), stress_response = c("response to", "stress", 
"immune", "inflammatory", "cytokine", "acute-phase", "acute phase", 
"wound", "healing", "detoxification", "oxidative stress", "platelet", 
"coagulation", "blood coagulation", "hemostasis"), apoptosis = c("apoptotic", 
"programmed cell death", "cell death"), viral_process = c("viral", 
"virus"))


### GO biological process (BP) enrichemnt 
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(scales)

base_dir <- "/allelic_rates_GO"
marker_dir <- file.path(base_dir, "greater_less_markers", "median_rate_filtered")
anno_file <- file.path(base_dir, "all.OG.KY.KH.homolog.anno.duplicatedOG.new.csv")
out_dir <- file.path(base_dir, "greater_less_GO_BP_gene_filtered")
category_file <- file.path(out_dir, "go_categories_definition.R")

padj_cutoff <- 0.05
min_gs_size <- 5
max_gs_size <- 500

subdirs <- c(
  "enrichGO_full",
  "enrichGO_padj_lt_0.05",
  "marker_entrez_mapping",
  "big_category/GOterm_with_big_category",
  "big_category/big_category_proportion",
  "bubble_plots",
  "stacked_bar_plots"
)
for (d in subdirs) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("^_+|_+$", "", x)
}

split_symbols <- function(x) {
  x <- x[!is.na(x)]
  x <- unlist(strsplit(x, "[;,|]+"))
  x <- trimws(x)
  unique(x[nzchar(x) & toupper(x) != "NA"])
}

parse_marker_name <- function(path) {
  nm <- basename(path)
  m <- regexec("^(.*)_(greater|less)_markers_minDiff0\\.25_padj0\\.05_medianRate(GT|LT)0\\.5\\.csv$", nm)
  hit <- regmatches(nm, m)[[1]]
  if (length(hit) != 4L) stop("Cannot parse marker file name: ", nm)
  list(cell_type = hit[2], direction = hit[3])
}

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2L) return(NA_real_)
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

wrap_label <- function(x, width = 34) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

write_or_backup <- function(x, path) {
  tryCatch({
    fwrite(x, path)
    path
  }, error = function(e) {
    backup <- sub("\\.csv$", "_new.csv", path)
    fwrite(x, backup)
    warning("Could not overwrite ", path, "; wrote ", backup, " instead. Error: ", conditionMessage(e))
    backup
  })
}

message("Loading big category definition: ", category_file)
category_env <- new.env(parent = baseenv())
source(category_file, local = category_env)
go_categories <- category_env$go_categories

classify_description <- function(description) {
  desc <- tolower(description)
  matched_categories <- character()
  matched_keywords <- character()
  for (category in names(go_categories)) {
    keys <- go_categories[[category]]
    hits <- keys[vapply(keys, function(k) grepl(k, desc, fixed = TRUE), logical(1))]
    if (length(hits) > 0L) {
      matched_categories <- c(matched_categories, category)
      matched_keywords <- c(matched_keywords, paste0(category, ":", paste(hits, collapse = "|")))
    }
  }
  data.table(
    big_category = if (length(matched_categories) > 0L) matched_categories[1] else "",
    all_matched_categories = paste(matched_categories, collapse = ";"),
    matched_keywords = paste(matched_keywords, collapse = ";")
  )
}

message("Reading homolog annotation: ", anno_file)
anno <- fread(anno_file, select = c("OG", "gene_homolog"))
anno <- anno[!is.na(OG) & OG != "" & !is.na(gene_homolog) & gene_homolog != ""]
anno_long <- anno[, .(human_symbol_raw = split_symbols(gene_homolog)), by = OG]
anno_long <- unique(anno_long[nzchar(human_symbol_raw)])
anno_long[, human_symbol_upper := toupper(human_symbol_raw)]

symbol_keys <- unique(c(anno_long$human_symbol_raw, anno_long$human_symbol_upper))
symbol_keys <- symbol_keys[nzchar(symbol_keys)]

message("Mapping human homolog symbols to Entrez IDs")
symbol_map <- as.data.table(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = symbol_keys,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
))
symbol_map <- unique(symbol_map[!is.na(ENTREZID) & ENTREZID != ""])

map_raw <- merge(
  anno_long[, .(OG, human_symbol_raw, SYMBOL_QUERY = human_symbol_raw)],
  symbol_map,
  by.x = "SYMBOL_QUERY",
  by.y = "SYMBOL",
  allow.cartesian = TRUE
)
map_upper <- merge(
  anno_long[, .(OG, human_symbol_raw, SYMBOL_QUERY = human_symbol_upper)],
  symbol_map,
  by.x = "SYMBOL_QUERY",
  by.y = "SYMBOL",
  allow.cartesian = TRUE
)
og_to_entrez <- unique(rbindlist(list(map_raw, map_upper), fill = TRUE)[
  !is.na(ENTREZID) & ENTREZID != "",
  .(OG, human_symbol = human_symbol_raw, human_symbol_mapped = SYMBOL_QUERY, ENTREZID)
])
universe_entrez <- unique(og_to_entrez$ENTREZID)
fwrite(og_to_entrez, file.path(out_dir, "OG_to_human_homolog_ENTREZ_mapping.csv"))
message("Annotated OG with Entrez mapping: ", uniqueN(og_to_entrez$OG))
message("Background Entrez IDs: ", length(universe_entrez))

marker_files <- list.files(
  marker_dir,
  pattern = "_(greater|less)_markers_minDiff0\\.25_padj0\\.05_medianRate(GT|LT)0\\.5\\.csv$",
  full.names = TRUE
)
marker_files <- marker_files[!grepl("^all_", basename(marker_files))]
marker_files <- marker_files[order(basename(marker_files))]
if (length(marker_files) == 0L) stop("No median-rate-filtered marker files found in: ", marker_dir)
message("Marker files: ", length(marker_files))

enrich_summary <- list()
all_full <- list()
all_sig <- list()
all_annotated <- list()
all_stats <- list()
all_class_summary <- list()

for (marker_file in marker_files) {
  meta <- parse_marker_name(marker_file)
  markers <- fread(marker_file)
  marker_genes <- unique(markers$gene)
  mapped <- unique(og_to_entrez[OG %in% marker_genes])
  entrez <- unique(mapped$ENTREZID)
  prefix <- sprintf("%s_%s", safe_name(meta$cell_type), meta$direction)

  fwrite(mapped, file.path(out_dir, "marker_entrez_mapping", sprintf("%s_marker_OG_to_human_homolog_ENTREZ.csv", prefix)))
  message("Running enrichGO BP: ", prefix, " (OG=", length(marker_genes), ", Entrez=", length(entrez), ")")

  if (length(entrez) == 0L) {
    full_df <- data.table()
  } else {
    ego <- enrichGO(
      gene = entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      readable = TRUE
    )
    full_df <- as.data.table(as.data.frame(ego))
    if (nrow(full_df) > 0L) {
      full_df[, `:=`(ontology = "BP", cell_type = meta$cell_type, direction = meta$direction)]
      setcolorder(full_df, c("ontology", "cell_type", "direction", setdiff(names(full_df), c("ontology", "cell_type", "direction"))))
    }
  }

  sig_df <- full_df[p.adjust < padj_cutoff]
  full_file <- file.path(out_dir, "enrichGO_full", sprintf("%s_enrichGO_BP_full.csv", prefix))
  sig_file <- file.path(out_dir, "enrichGO_padj_lt_0.05", sprintf("%s_enrichGO_BP_padj_lt_0.05.csv", prefix))
  fwrite(full_df, full_file)
  fwrite(sig_df, sig_file)

  if (nrow(sig_df) > 0L) {
    cls <- rbindlist(lapply(sig_df$Description, classify_description))
    annotated <- cbind(sig_df, cls)
  } else {
    annotated <- copy(sig_df)
    annotated[, `:=`(big_category = character(), all_matched_categories = character(), matched_keywords = character())]
  }

  annotated_file <- file.path(out_dir, "big_category", "GOterm_with_big_category", sprintf("%s_GOterm_big_category.csv", prefix))
  fwrite(annotated, annotated_file)

  total_terms <- nrow(annotated)
  classified_terms <- if (total_terms > 0L) sum(annotated$big_category != "") else 0L
  unclassified_terms <- total_terms - classified_terms

  if (total_terms > 0L) {
    stats <- annotated[, .N, by = .(ontology, cell_type, direction, big_category)][order(-N, big_category)]
    stats[, total_terms := total_terms]
    stats[, classified_terms := classified_terms]
    stats[, unclassified_terms := unclassified_terms]
    stats[, prop_of_all_terms := N / total_terms]
    stats[, prop_of_classified_terms := fifelse(big_category == "" | classified_terms == 0L, NA_real_, N / classified_terms)]
    setcolorder(stats, c("ontology", "cell_type", "direction", "big_category", "N", "total_terms", "classified_terms", "unclassified_terms", "prop_of_all_terms", "prop_of_classified_terms"))
  } else {
    stats <- data.table(
      ontology = "BP", cell_type = meta$cell_type, direction = meta$direction,
      big_category = character(), N = integer(), total_terms = integer(),
      classified_terms = integer(), unclassified_terms = integer(),
      prop_of_all_terms = numeric(), prop_of_classified_terms = numeric()
    )
  }
  stats_file <- file.path(out_dir, "big_category", "big_category_proportion", sprintf("%s_big_category_proportion.csv", prefix))
  fwrite(stats, stats_file)

  enrich_summary[[prefix]] <- data.table(
    ontology = "BP",
    cell_type = meta$cell_type,
    direction = meta$direction,
    marker_OG = length(marker_genes),
    mapped_Entrez = length(entrez),
    enriched_terms = nrow(full_df),
    significant_terms = nrow(sig_df),
    full_file = full_file,
    significant_file = sig_file
  )
  all_full[[prefix]] <- full_df
  all_sig[[prefix]] <- sig_df
  all_annotated[[prefix]] <- annotated
  all_stats[[prefix]] <- stats
  all_class_summary[[prefix]] <- data.table(
    ontology = "BP",
    cell_type = meta$cell_type,
    direction = meta$direction,
    total_terms = total_terms,
    classified_terms = classified_terms,
    unclassified_terms = unclassified_terms,
    classified_fraction = if (total_terms > 0L) classified_terms / total_terms else NA_real_,
    unclassified_fraction = if (total_terms > 0L) unclassified_terms / total_terms else NA_real_,
    annotated_file = annotated_file,
    stats_file = stats_file
  )
}

summary_dt <- rbindlist(enrich_summary, use.names = TRUE, fill = TRUE)
write_or_backup(summary_dt, file.path(out_dir, "GO_BP_enrichment_summary.csv"))
write_or_backup(rbindlist(all_full, use.names = TRUE, fill = TRUE), file.path(out_dir, "all_enrichGO_BP_full.csv"))
write_or_backup(rbindlist(all_sig, use.names = TRUE, fill = TRUE), file.path(out_dir, "all_enrichGO_BP_padj_lt_0.05.csv"))
write_or_backup(rbindlist(all_annotated, use.names = TRUE, fill = TRUE), file.path(out_dir, "big_category", "all_GOterms_with_big_category.csv"))
write_or_backup(rbindlist(all_stats, use.names = TRUE, fill = TRUE), file.path(out_dir, "big_category", "all_big_category_proportions.csv"))
write_or_backup(rbindlist(all_class_summary, use.names = TRUE, fill = TRUE), file.path(out_dir, "big_category", "big_category_classification_summary.csv"))

### plots 
library(data.table)
library(ggplot2)
library(scales)
library(grid)

out_dir <- "/greater_less_GO_BP_gene_filtered"
sig_dir <- file.path(out_dir, "enrichGO_padj_lt_0.05")
category_file <- file.path(out_dir, "go_categories_definition.R")
big_dir <- file.path(out_dir, "big_category")
bubble_dir <- file.path(out_dir, "bubble_plots")
stack_dir <- file.path(out_dir, "stacked_bar_plots")

dir.create(file.path(big_dir, "GOterm_with_big_category"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(big_dir, "big_category_proportion"), recursive = TRUE, showWarnings = FALSE)
dir.create(bubble_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stack_dir, recursive = TRUE, showWarnings = FALSE)

cell_type_levels <- c("endoderm", "epidermis", "germ", "mesenchyme", "muscle_heart", "nervous_system", "notochord")

message("Loading categories: ", category_file)
category_env <- new.env(parent = baseenv())
source(category_file, local = category_env)
go_categories <- category_env$go_categories
category_levels <- c(names(go_categories), "unclassified")

category_colors <- c(
  transcription = "#0072B2",
  metabolism = "#E69F00",
  signaling = "#009E73",
  phosphorylation = "#D55E00",
  development = "#CC79A7",
  neural_process = "#56B4E9",
  translation_rna_processing = "#F0E442",
  cell_cycle = "#332288",
  dna_repair = "#117733",
  transport = "#88CCEE",
  protein_turnover = "#AA4499",
  cell_adhesion = "#44AA99",
  stress_response = "#DDCC77",
  apoptosis = "#882255",
  viral_process = "#999933",
  unclassified = "#BDBDBD"
)

parse_file_name <- function(path) {
  nm <- basename(path)
  m <- regexec("^(.*)_(greater|less)_enrichGO_BP_padj_lt_0\\.05\\.csv$", nm)
  hit <- regmatches(nm, m)[[1]]
  if (length(hit) != 3L) stop("Cannot parse file name: ", nm)
  list(cell_type = hit[2], direction = hit[3])
}

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2L) return(NA_real_)
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

wrap_label <- function(x, width = 34) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

classify_description <- function(description) {
  desc <- tolower(description)
  matched_categories <- character()
  matched_keywords <- character()
  for (category in names(go_categories)) {
    keys <- go_categories[[category]]
    hits <- keys[vapply(keys, function(k) grepl(k, desc, fixed = TRUE), logical(1))]
    if (length(hits) > 0L) {
      matched_categories <- c(matched_categories, category)
      matched_keywords <- c(matched_keywords, paste0(category, ":", paste(hits, collapse = "|")))
    }
  }
  data.table(
    big_category = if (length(matched_categories) > 0L) matched_categories[1] else "",
    all_matched_categories = paste(matched_categories, collapse = ";"),
    matched_keywords = paste(matched_keywords, collapse = ";")
  )
}

write_or_backup <- function(x, path) {
  tryCatch({
    fwrite(x, path)
    path
  }, error = function(e) {
    backup <- sub("\\.csv$", "_new.csv", path)
    fwrite(x, backup)
    warning("Could not overwrite ", path, "; wrote ", backup, " instead. Error: ", conditionMessage(e))
    backup
  })
}

message("Classifying significant GO terms from: ", sig_dir)
files <- list.files(sig_dir, pattern = "_(greater|less)_enrichGO_BP_padj_lt_0\\.05\\.csv$", full.names = TRUE)
files <- files[order(basename(files))]
if (length(files) == 0L) stop("No significant enrichGO CSV files found in: ", sig_dir)

all_annotated <- list()
all_stats <- list()
all_summary <- list()

for (f in files) {
  meta <- parse_file_name(f)
  prefix <- sprintf("%s_%s", meta$cell_type, meta$direction)
  message("  ", prefix)
  dt <- fread(f)
  if (nrow(dt) > 0L) {
    cls <- rbindlist(lapply(dt$Description, classify_description))
    dt <- cbind(dt, cls)
  } else {
    dt[, `:=`(big_category = character(), all_matched_categories = character(), matched_keywords = character())]
  }

  annotated_file <- file.path(big_dir, "GOterm_with_big_category", sprintf("%s_GOterm_big_category.csv", prefix))
  fwrite(dt, annotated_file)

  total_terms <- nrow(dt)
  classified_terms <- if (total_terms > 0L) sum(dt$big_category != "") else 0L
  unclassified_terms <- total_terms - classified_terms

  if (total_terms > 0L) {
    stats <- dt[, .N, by = .(ontology, cell_type, direction, big_category)][order(-N, big_category)]
    stats[, total_terms := total_terms]
    stats[, classified_terms := classified_terms]
    stats[, unclassified_terms := unclassified_terms]
    stats[, prop_of_all_terms := N / total_terms]
    stats[, prop_of_classified_terms := fifelse(big_category == "" | classified_terms == 0L, NA_real_, N / classified_terms)]
    setcolorder(stats, c("ontology", "cell_type", "direction", "big_category", "N", "total_terms", "classified_terms", "unclassified_terms", "prop_of_all_terms", "prop_of_classified_terms"))
  } else {
    stats <- data.table(
      ontology = "BP", cell_type = meta$cell_type, direction = meta$direction,
      big_category = character(), N = integer(), total_terms = integer(),
      classified_terms = integer(), unclassified_terms = integer(),
      prop_of_all_terms = numeric(), prop_of_classified_terms = numeric()
    )
  }
  stats_file <- file.path(big_dir, "big_category_proportion", sprintf("%s_big_category_proportion.csv", prefix))
  fwrite(stats, stats_file)

  all_annotated[[prefix]] <- dt
  all_stats[[prefix]] <- stats
  all_summary[[prefix]] <- data.table(
    ontology = "BP", cell_type = meta$cell_type, direction = meta$direction,
    total_terms = total_terms, classified_terms = classified_terms, unclassified_terms = unclassified_terms,
    classified_fraction = if (total_terms > 0L) classified_terms / total_terms else NA_real_,
    unclassified_fraction = if (total_terms > 0L) unclassified_terms / total_terms else NA_real_,
    annotated_file = annotated_file, stats_file = stats_file
  )
}

go <- rbindlist(all_annotated, use.names = TRUE, fill = TRUE)
props <- rbindlist(all_stats, use.names = TRUE, fill = TRUE)
class_summary <- rbindlist(all_summary, use.names = TRUE, fill = TRUE)
write_or_backup(go, file.path(big_dir, "all_GOterms_with_big_category.csv"))
write_or_backup(props, file.path(big_dir, "all_big_category_proportions.csv"))
write_or_backup(class_summary, file.path(big_dir, "big_category_classification_summary.csv"))

message("Drawing bubble plots")
go <- go[!is.na(Description) & Description != ""]
go[, GeneRatio_num := parse_ratio(GeneRatio)]
go[, p.adjust := as.numeric(p.adjust)]
props2 <- props[!is.na(big_category) & big_category != ""]
props2[, N := as.integer(N)]
props2[, prop_of_all_terms := as.numeric(prop_of_all_terms)]

top_categories <- props2[order(cell_type, direction, -prop_of_all_terms, -N, big_category), head(.SD, 2), by = .(cell_type, direction)]
if (nrow(top_categories) > 0L) {
  setnames(top_categories, c("N", "prop_of_all_terms"), c("category_term_count", "category_prop_of_all_terms"))
}

representatives <- if (nrow(top_categories) > 0L) {
  merge(
    go,
    top_categories[, .(cell_type, direction, big_category, category_term_count, category_prop_of_all_terms)],
    by = c("cell_type", "direction", "big_category"),
    allow.cartesian = TRUE
  )
} else data.table()

if (nrow(representatives) > 0L) {
  representatives <- representatives[order(cell_type, direction, -category_prop_of_all_terms, -category_term_count, big_category, p.adjust, ID), head(.SD, 1), by = .(cell_type, direction, big_category)]
}
fwrite(representatives, file.path(bubble_dir, "selected_top2_category_representative_GOterms.csv"))

make_bubble_plot <- function(direction_label) {
  selected_ids <- unique(representatives[direction == direction_label, ID])
  selected_meta <- representatives[direction == direction_label, .(
    representative_description = Description[1],
    selected_big_category = paste(unique(big_category), collapse = ";")
  ), by = .(ID)]
  plot_dt <- go[direction == direction_label & ID %in% selected_ids]
  if (nrow(plot_dt) == 0L) {
    warning("No bubble plot data for direction: ", direction_label)
    return(invisible(NULL))
  }
  plot_dt <- merge(plot_dt, selected_meta, by = "ID", all.x = TRUE, allow.cartesian = TRUE)
  plot_dt[, selected_big_category := fifelse(is.na(selected_big_category) | selected_big_category == "", "unclassified", selected_big_category)]
  plot_dt[, term_label := wrap_label(paste0(selected_big_category, " | ", Description), width = 42)]
  term_order <- plot_dt[, .(min_padj = min(p.adjust, na.rm = TRUE), max_ratio = max(GeneRatio_num, na.rm = TRUE), selected_rank = min(match(ID, selected_ids), na.rm = TRUE)), by = .(ID, term_label)][order(selected_rank, min_padj, -max_ratio)]
  plot_dt[, term_label := factor(term_label, levels = rev(term_order$term_label))]
  present_cells <- intersect(cell_type_levels, unique(plot_dt$cell_type))
  if (length(present_cells) == 0L) present_cells <- unique(plot_dt$cell_type)
  plot_dt[, cell_type := factor(cell_type, levels = present_cells)]

  plot_data_file <- file.path(bubble_dir, sprintf("%s_big_category_bubble_plot_data.csv", direction_label))
  fwrite(plot_dt, plot_data_file)
  n_terms <- uniqueN(plot_dt$term_label)
  plot_height <- max(6, min(18, 0.36 * n_terms + 2.5))
  plot_width <- max(6.8, 0.55 * length(present_cells) + 4.4)

  p <- ggplot(plot_dt, aes(x = cell_type, y = term_label)) +
    geom_point(aes(size = GeneRatio_num, fill = p.adjust), shape = 21, color = "grey25", stroke = 0.25, alpha = 0.95) +
    scale_size_continuous(name = "GeneRatio", range = c(2.2, 8.5), labels = percent_format(accuracy = 1)) +
    scale_fill_gradient(name = "p.adjust", low = "#d73027", high = "#4575b4", trans = "reverse", labels = scientific_format(digits = 2)) +
    scale_x_discrete(expand = expansion(add = 0.25)) +
    labs(title = sprintf("GO BP enrichment bubble plot (%s)", direction_label), x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
  pdf_file <- file.path(out_dir, sprintf("%s_GO_BP_big_category_bubble_plot.pdf", direction_label))
  sub_pdf_file <- file.path(bubble_dir, sprintf("%s_GO_BP_big_category_bubble_plot.pdf", direction_label))
  ggsave(pdf_file, p, width = plot_width, height = plot_height, units = "in", limitsize = FALSE)
  ggsave(sub_pdf_file, p, width = plot_width, height = plot_height, units = "in", limitsize = FALSE)
  message("  bubble: ", direction_label)
}
make_bubble_plot("greater")
make_bubble_plot("less")

message("Drawing stacked bar plots")
stack_dt <- copy(props)
stack_dt[, big_category := fifelse(is.na(big_category) | big_category == "", "unclassified", big_category)]
stack_dt <- stack_dt[N > 0]
stack_dt <- stack_dt[, .(N = sum(N)), by = .(ontology, cell_type, direction, big_category)]
stack_dt[, total_N := sum(N), by = .(cell_type, direction)]
stack_dt[, percent := N / total_N]

present_categories <- c(intersect(category_levels, unique(stack_dt$big_category)), setdiff(unique(stack_dt$big_category), category_levels))
stack_dt[, big_category := factor(big_category, levels = present_categories)]
present_cells <- intersect(cell_type_levels, unique(stack_dt$cell_type))
stack_dt[, cell_type := factor(cell_type, levels = rev(present_cells))]
fwrite(stack_dt, file.path(stack_dir, "BP_big_category_stacked_bar_plot_data.csv"))

category_colors <- category_colors[names(category_colors) %in% present_categories]
make_stack_plot <- function(direction_label) {
  pdt <- stack_dt[direction == direction_label]
  if (nrow(pdt) == 0L) {
    warning("No stacked bar data for direction: ", direction_label)
    return(invisible(NULL))
  }
  p <- ggplot(pdt, aes(x = percent, y = cell_type, fill = big_category)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.18) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.01))) +
    scale_fill_manual(values = category_colors, breaks = present_categories, drop = FALSE, name = "Big category") +
    labs(title = sprintf("BP big category composition (%s)", direction_label), x = "GO term proportion", y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.32, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
  pdf_file <- file.path(out_dir, sprintf("%s_BP_big_category_stacked_bar.pdf", direction_label))
  sub_pdf_file <- file.path(stack_dir, sprintf("%s_BP_big_category_stacked_bar.pdf", direction_label))
  ggsave(pdf_file, p, width = 7.1, height = 3.8, units = "in")
  ggsave(sub_pdf_file, p, width = 7.1, height = 3.8, units = "in")
  message("  stacked: ", direction_label)
}
make_stack_plot("greater")
make_stack_plot("less")


