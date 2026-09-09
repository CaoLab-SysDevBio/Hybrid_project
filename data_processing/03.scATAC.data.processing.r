# ============================================================
# The scATAC-seq data is from mid-tailbud stage.
# ============================================================
library(Signac)
library(Seurat)
library(ggplot2)
library(ggforce)
library(GenomicRanges)
library(data.table)
library(patchwork)
counts <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

mtb <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

### QC ###
mtb <- NucleosomeSignal(mtb)
mtb <- TSSEnrichment(mtb)
mtb$pct_reads_in_peaks <- mtb$peak_region_fragments / mtb$passed_filters * 100
pdf("CsCi_MTB_QC_plots.pdf", width = 8, height = 6)
DensityScatter(
  mtb, 
  x = 'nCount_peaks', 
  y = 'TSS.enrichment', 
  log_x = TRUE, 
  quantiles = TRUE
)
DensityScatter(
  mtb,
  x = "nCount_peaks",
  y = "nucleosome_signal",
  log_x = TRUE
)
VlnPlot(
  object = mtb,
  features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"),
  ncol = 4,
  pt.size = 0.1
)
dev.off()

mtb <- subset(
  x = mtb,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 30 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)
mtb <- RunTFIDF(mtb)
mtb <- FindTopFeatures(mtb, min.cutoff = 'q0')
mtb <- RunSVD(mtb)
pdf("CsCi_MTB_cor_LSI.pdf", width = 8, height = 6)
DepthCor(mtb)
dev.off()
mtb <- RunUMAP(mtb, reduction = 'lsi', dims = 2:30)
mtb <- FindNeighbors(mtb, reduction = 'lsi', dims = 2:30)
mtb <- FindClusters(mtb, verbose = FALSE, algorithm = 3)
pdf("CsCi_MTB_cluster.pdf", width = 8, height = 6)
DimPlot(mtb, label = TRUE) + NoLegend()
dev.off()

### Create a gene activity matrix ###
gene.activities <- GeneActivity(mtb)
gene.activities <- gene.activities[!duplicated(rownames(gene.activities)), ]
mtb[['Activity']] <- CreateAssayObject(counts = gene.activities)
mtb[["RNA"]] <- NULL
mtb$nCount_RNA <- NULL
mtb$nFeature_RNA <- NULL
mtb <- NormalizeData(
  object = mtb,
  assay = 'Activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(mtb$nCount_Activity)
)
saveRDS(mtb, "CsCi_MTB_ATAC.rds")



# =====================================================================================
# integrate RNA and ATAC objects for the same cell states of scDALI input
# =====================================================================================
# MTB RNA object
mtb.rna <- readRDS("data/CsCi.MTB.RNA.OG.rds")
mtb.gene.activity <- readRDS("data/CsCi.MTB.gene.activity.OG.rds")
# integration
features <- SelectIntegrationFeatures(
  object.list = list(mtb.rna, mtb.gene.activity),
  nfeatures = 2000
)
mtb.rna.merge <- mtb.rna
mtb.rna.merge <- ScaleData(mtb.rna.merge, features = features)
mtb.rna.merge <- RunPCA(mtb.rna.merge, features = features)
mtb.gene.activity.merge <- mtb.gene.activity
mtb.gene.activity.merge <- ScaleData(mtb.gene.activity.merge, features = features)
mtb.gene.activity.merge <- RunPCA(mtb.gene.activity.merge, features = features)

anchors.integrate <- FindIntegrationAnchors(
  object.list = list(mtb.rna.merge, mtb.gene.activity.merge),
  anchor.features = features,
  reduction = "cca"
)
integrated  <- IntegrateData(
  anchorset = anchors.integrate
)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

pc30 <- Embeddings(integrated, reduction = "pca")[, 1:30]
write.csv(pc30, file = "CsCi_mtb_integrate_RNA_gene.activity_OG_PCA30.csv")
cell.atac <- colnames(integrated)[integrated$orig.ident == "CsCi_MTB_OG"]
pc30.ATAC <- pc30[cell.atac, , drop = FALSE]
cell.rna <- colnames(integrated)[integrated$orig.ident == "CsCi_MTB"]
pc30.RNA <- pc30[cell.rna, , drop = FALSE]
write.csv(pc30.ATAC, file = "CsCi_MTB_integrate_ATAC_gene.activity_OG_PCA30.csv")
write.csv(pc30.RNA, file = "CsCi_MTB_integrate_RNA_OG_PCA30.csv")
saveRDS(integrated, "data/CsCi.MTB.integrated.RNA.ATAC.gene.activity.OG.rds")