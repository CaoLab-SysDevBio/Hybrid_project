# ===========================================================================================
# Construct peaks count at chromatin level. We used data at mid-tailbud stage as an example.
# ===========================================================================================
library(Signac)
library(Seurat)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(data.table)
library(patchwork)
library(dplyr)
# MTB peak counts
mtb.obj <- readRDS("CsCi_MTB_ATAC.rds")
DefaultAssay(mtb.obj) <- "peaks"
anno <- Annotation(mtb.obj)
exon <- anno[anno$type == "exon"]
exon.list <- split(exon, exon$gene_id)
tss <- lapply(exon.list, function(x) {
  if (as.character(strand(x)[1]) == "+") {
    min(start(x))
  } else {
    max(end(x))
  }
})
gene.ids <- names(exon.list)
tss <- unlist(tss)
seqs <- sapply(exon.list, function(x) as.character(seqnames(x)[1]))
strs <- sapply(exon.list, function(x) as.character(strand(x)[1]))
tss.gr <- GRanges(
  seqnames = seqs,
  ranges = IRanges(start = tss, width = 1),
  strand = strs,
  gene_id = gene.ids
)
prom <- GRanges(
  seqnames = seqnames(tss.gr),
  ranges = IRanges(
    start = ifelse(
      as.character(strand(tss.gr)) == "+",
      start(tss.gr) - 2000,
      start(tss.gr)
    ),
    end = ifelse(
      as.character(strand(tss.gr)) == "+",
      start(tss.gr),
      start(tss.gr) + 2000
    )
  ),
  strand = strand(tss.gr),
  gene_id = tss.gr$gene_id
)
gene_body <- lapply(exon.list, function(x) {
  GRanges(
    seqnames = unique(seqnames(x)),
    ranges = IRanges(
      start = min(start(x)),
      end = max(end(x))
    ),
    strand = unique(strand(x)),
    gene_id = x$gene_id[1]
  )
})
gene_body.gr <- GenomicRanges::GRangesList(gene_body)
gene_body.gr <- unlist(gene_body.gr)
names(gene_body.gr) <- NULL
gene_region <- c(prom, gene_body.gr)

# by gene group and then reduce
gene_region_list <- split(gene_region, gene_region$gene_id)
gene_region_final <- lapply(gene_region_list, function(gr) {
  gr_reduced <- reduce(gr)
  gr_reduced$gene_id <- unique(gr$gene_id)
  gr_reduced
})
gene_region_final <- unlist(GRangesList(gene_region_final))
class(gene_region_final)
names(gene_region_final) <- NULL

# overlap
mtb.peaks <- granges(mtb.obj) #get peaks
hits <- findOverlaps(mtb.peaks, gene_region_final)

peak_gene_df <- data.frame(
  peak = mtb.peaks[queryHits(hits)],
  gene_id = gene_region_final$gene_id[subjectHits(hits)],
  gene_strand = as.character(strand(gene_region_final)[subjectHits(hits)])
)
peak_gene_df$peak_name <- paste(
  peak_gene_df$peak.seqnames,
  peak_gene_df$peak.start,
  peak_gene_df$peak.end,
  sep = "-"
)
# match with OG
OG.pairs <- read.csv("OG.Ci.Cs.pairs.csv")
OG.pairs$Cint_gene <- gsub("Cint\\.KY-{15}", "Cint.KY.", OG.pairs$Cint_gene)
OG.pairs$Csav_gene <- gsub("^(Csav[^-]*-[^-]+)-", "\\1_", OG.pairs$Csav_gene)
OG.long <- rbind(
  data.frame(
    OG = OG.pairs$Ci_OG,
    gene_id = OG.pairs$Cint_gene
  ),
  data.frame(
    OG = OG.pairs$Cs_OG,
    gene_id = OG.pairs$Csav_gene
  )
)
peak_gene_df_OG <- peak_gene_df %>%
  left_join(OG.long, by = "gene_id")

peak_gene_df_OG <- peak_gene_df_OG %>%
  group_by(peak_name) %>%
  mutate(overlap_peak = ifelse(n() > 1, "overlap", "unique")) %>%
  ungroup()
peak_gene_df_OG <- as.data.frame(peak_gene_df_OG)
write.csv(peak_gene_df_OG, "CsCi_MTB_all_peak_allele_OG.csv", row.names = FALSE)

peak_gene_df_OG_f <- peak_gene_df_OG[!is.na(peak_gene_df_OG$OG), ]
peak_gene_df_OG_f$allele <- ifelse(grepl("^Cint", peak_gene_df_OG_f$peak_name),
                             "maternal", "paternal")
write.csv(peak_gene_df_OG_f, "CsCi_MTB_all_peak_allele_OG_noNA.csv", row.names = FALSE)

DefaultAssay(mtb.obj) <- "peaks"
frag_total <- GetAssayData(mtb.obj, layer = "counts")
library(Matrix)
split_peaks_by_OG <- split(peak_gene_df_OG_f$peak_name, peak_gene_df_OG_f$OG)
head(split_peaks_by_OG)
OG_counts_total <- sapply(split_peaks_by_OG, function(peaks) {
  Matrix::colSums(frag_total[peaks, , drop = FALSE])
})
write.csv(OG_counts_total, "CsCi_MTB_atac_all_peak_counts_OG_total.csv")

# filter
gene_total <- colSums(OG_counts_total)
keep_genes <- names(gene_total)[gene_total >= 30]
OG_counts_filtered <- OG_counts_total[, keep_genes]
write.csv(OG_counts_filtered, "CsCi_MTB_atac_all_peak_counts_OG_total_filter.csv")

OG_counts_maternal <- sapply(split_peaks_by_OG, function(peaks) {
  maternal_peaks <- peaks[grepl("^Cint", peaks)]
  if (length(maternal_peaks) == 0) return(rep(0, ncol(frag_total)))
  Matrix::colSums(frag_total[maternal_peaks, , drop = FALSE])
})
all(rownames(OG_counts_maternal) == rownames(OG_counts_total))
write.csv(OG_counts_maternal, "CsCi_MTB_atac_all_peak_counts_OG_maternal.csv")
OG_counts_maternal_filtered <- OG_counts_maternal[, keep_genes]
write.csv(OG_counts_maternal_filtered, "CsCi_MTB_atac_all_peak_counts_OG_maternal_filter.csv")
