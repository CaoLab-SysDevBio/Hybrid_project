# ============================================================
# We used multiomic data at LatTII stage as an example.
# ============================================================
library(Signac)
library(Seurat)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(data.table)
# load the RNA and ATAC data
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
fragpath <- "atac_fragments.tsv.gz"

# get gene annotations using gtf file
gtf_df <- fread("hybrid.genes.gtf", header = FALSE, sep = "\t")
gtf_df[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', V9)]
gtf_df[, transcript_id := sub('.*transcript_id "([^"]+)".*', '\\1', V9)]
# Cint gene id and chr id
gtf_df[, gene_id := gsub("Cint\\.KY_{15}", "Cint.KY.", gene_id)]
gtf_df[, transcript_id := gsub("Cint\\.KY_{15}", "Cint.KY.", transcript_id)]
gtf_df$V1 <- gsub("Cint\\.KY_{15}", "Cint.KY.", gtf_df$V1)
# Csav 
gtf_df$V1 <- sub("^Csav\\.dovetail[-_]", "Csav.dovetail.", gtf_df$V1)
gtf_df$V1 <- gsub("_", ".", gtf_df$V1)
# construct GRanges
annotation <- GRanges(
  seqnames = gtf_df$V1,
  ranges   = IRanges(start = gtf_df$V4, end = gtf_df$V5),
  strand   = gtf_df$V7
)
mcols(annotation)$type <- gtf_df$V3
mcols(annotation)$gene_id <- gtf_df$gene_id
mcols(annotation)$tx_id <- gtf_df$transcript_id
mcols(annotation)$gene_name <- gtf_df$gene_id
mcols(annotation)$gene_biotype <- ifelse(
  gtf_df$V3 == "gene",
  "gene",
  "transcript"
)
colnames(mcols(annotation))
# [1] "type"         "gene_id"      "tx_id"        "gene_name"    "gene_biotype"
saveRDS(annotation, file = "data/CsCi_signac_gene_annotation.rds")

annotation <- readRDS("data/CsCi_signac_gene_annotation.rds")
frags <- read.table(
  gzfile("atac_fragments.tsv.gz"),
  comment.char = "#",
  header = FALSE,
  stringsAsFactors = FALSE
)
anno.atac <- annotation[annotation$type %in% c("exon", "CDS")]
frag_chrs <- unique(frags$V1)
anno.atac <- anno.atac[seqnames(anno.atac) %in% frag_chrs]
seqlevels(anno.atac) <- frag_chrs
seqinfo(anno.atac) <- Seqinfo(seqnames = frag_chrs)

anno.atac$gene_biotype <- "protein_coding"
anno.atac$type <- factor(anno.atac$type, levels = c("exon", "CDS"))
saveRDS(anno.atac, file = "CsCi_signac_gene_annotation_atac.rds")

# create a Seurat object containing the RNA adata
CsCi.obj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
# create ATAC assay and add it to the object
CsCi.obj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = anno.atac
)

# Quality control
# We can compute per-cell quality control metrics using the DNA accessibility data and remove cells that are outliers for these metrics, as well as cells with low or unusually high counts for either the RNA or ATAC assay.
DefaultAssay(CsCi.obj) <- "ATAC"
CsCi.obj <- NucleosomeSignal(CsCi.obj)
CsCi.obj <- TSSEnrichment(CsCi.obj)
CsCi.obj$orig.ident <- gsub("SeuratProject", "CsCi_LTB2", CsCi.obj$orig.ident)

# The relationship between variables stored in the object metadata can be visualized using the DensityScatter() function. This can also be used to quickly find suitable cutoff values for different QC metrics by setting quantiles=TRUE:
# density means the density of cell number
pdf("CsCi_LTB2_QC_plots.pdf", width = 8, height = 6)
DensityScatter(
  CsCi.obj, 
  x = 'nCount_ATAC', 
  y = 'TSS.enrichment', 
  log_x = TRUE, 
  quantiles = TRUE
)
DensityScatter(
  CsCi.obj,
  x = "nCount_ATAC",
  y = "nucleosome_signal",
  log_x = TRUE
)
VlnPlot(
  object = CsCi.obj,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
dev.off()

# filter out low quality cells
CsCi.obj <- subset(
  x = CsCi.obj,
  subset = nCount_ATAC < 10000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 200 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 0.5
)

# Gene expression data processing
# We can normalize the gene expression data using SCTransform, and reduce the dimensionality using PCA.
DefaultAssay(CsCi.obj) <- "RNA"
CsCi.obj <- SCTransform(CsCi.obj)
CsCi.obj <- RunPCA(CsCi.obj)

# DNA accessibility data processing
# Here we process the DNA accessibility assay the same way we would process a scATAC-seq dataset, by performing latent semantic indexing (LSI).
DefaultAssay(CsCi.obj) <- "ATAC"
CsCi.obj <- FindTopFeatures(CsCi.obj, min.cutoff = 'q0') #Cutoff for feature to be included in the VariableFeatures for the object.
CsCi.obj <- RunTFIDF(CsCi.obj) #TF-IDF normalization
CsCi.obj <- RunSVD(CsCi.obj)
saveRDS(CsCi.obj, file = "CsCi_LTB2_multiome.rds")