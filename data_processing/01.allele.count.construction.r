# ============================================================
# We used three replicates at LatTII stage as an example.
# ============================================================
### allele matrix construction ###
library(Seurat)
library(Matrix)
input_root <- paste0("/allelic_expression_matrices/")
output_dir <- file.path("/count_matrices/")

samples <- c("LTB2_Cs2Ci1", "LTB2_Cs1Ci1", "LTB2_Cs1Ci2")
for (sample_name in samples) {
  message("Processing: ", sample_name)
  matrix_dir <- file.path(input_root, sample_name, "filtered_feature_bc_matrix")
  # read 10x raw UMI count matrix
  count_matrix <- Read10X(data.dir = matrix_dir, gene.column = 1, unique.features = TRUE)
  count_matrix <- as(count_matrix, "dgCMatrix")
  output_file <- file.path(output_dir, paste0(sample_name, ".raw_allelic_count_matrix.rds"))
  saveRDS(count_matrix, file = output_file, compress = TRUE)
}

### mapping OG ###
library(dplyr)
file_dir <- paste0("/count_matrices/")
OG_gene <- read.table("/data/Ci_KY-Cs_dovetail-Aniseed.matchID.OG.all.txt", sep='\t')
colnames(OG_gene) <- c("ortho_grp","Gene")
samples <- c("LTB2_Cs2Ci1", "LTB2_Cs1Ci1", "LTB2_Cs1Ci2")
for(sample_name in samples){
    message("Processing: ", sample_name)
    count_df <- readRDS(paste0(file_dir, sample_name, ".raw_allelic_count_matrix.rds"))
    counts.matrix <- as.data.frame(as.matrix(count_df))
    counts.matrix$Gene<-rownames(counts.matrix)
    #### transfer gene name '_' to '-'
    counts.matrix$Gene <- gsub("_", "-", counts.matrix$Gene)
    counts.matrix.OG<-merge(counts.matrix,OG_gene,by="Gene",sort = FALSE)
    message("After OG mapping: ",nrow(counts.matrix.OG), " rows × ", ncol(counts.matrix.OG), " columns")
    counts.matrix.OG.sum <- counts.matrix.OG[,-1] %>% group_by(ortho_grp) %>% summarise_all(funs(sum))
    counts.matrix.OG.sum <- as.data.frame(counts.matrix.OG.sum)
    rownames(counts.matrix.OG.sum) <- counts.matrix.OG.sum[,1]
    counts.matrix.OG.sum <- counts.matrix.OG.sum[,-1]
    message("dim of counts.matrix.OG.sum: ",nrow(counts.matrix.OG.sum)," rows × ",  ncol(counts.matrix.OG.sum)," columns")
    colnames(counts.matrix.OG.sum) <- gsub("\\.", "-", colnames(counts.matrix.OG.sum)) 
    # . is Regular Expression or .\ is Escape Character
    write.table(counts.matrix.OG.sum, file=paste0(file_dir, sample_name, ".counts.matrix.OG.total.txt"), sep='\t')
    counts.matrix.OG <- as.data.frame(counts.matrix.OG)
    message("dim of counts.matrix.OG: ",nrow(counts.matrix.OG), " rows × ", ncol(counts.matrix.OG)," columns") #include gene symbol and OG columns
    write.table(counts.matrix.OG,file=paste0(file_dir, sample_name, ".raw.counts.matrix.OG.txt"), sep='\t', row.names = FALSE)
    ### generate maternal matrix ###
    counts.matrix.OG.maternal <- counts.matrix.OG[grepl("^Cint\\.KY", counts.matrix.OG$Gene), , drop = FALSE]
    message("dim of counts.matrix.OG.maternal: ",nrow(counts.matrix.OG.maternal), "x", ncol(counts.matrix.OG.maternal)) 
    counts.matrix.OG.maternal <- counts.matrix.OG.maternal[,-1]
    counts.matrix.OG.maternal <- as.data.frame(counts.matrix.OG.maternal)
    rownames(counts.matrix.OG.maternal) <- counts.matrix.OG.maternal$ortho_grp
    counts.matrix.OG.maternal <- counts.matrix.OG.maternal[, -ncol(counts.matrix.OG.maternal)]
    message("dim of counts.matrix.OG.maternal after deleting: ",nrow(counts.matrix.OG.maternal), "x", ncol(counts.matrix.OG.maternal)) 
    colnames(counts.matrix.OG.maternal) <- gsub("\\.", "-", colnames(counts.matrix.OG.maternal)) 
    missing_genes <- setdiff(rownames(counts.matrix.OG.sum),rownames(counts.matrix.OG.maternal))
    message("Number of genes missing from maternal matrix: ", length(missing_genes))

    if (length(missing_genes) > 0) {
    missing_gene_matrix <- matrix(
        0,
        nrow = length(missing_genes),
        ncol = ncol(counts.matrix.OG.maternal),
        dimnames = list(
        missing_genes,
        colnames(counts.matrix.OG.maternal)
        )
    )
    ### add missing genes to maternal matrix ###
    counts.matrix.OG.maternal <- rbind(
        counts.matrix.OG.maternal,
        missing_gene_matrix
    )
    }
    counts.matrix.OG.maternal <- counts.matrix.OG.maternal[rownames(counts.matrix.OG.sum), colnames(counts.matrix.OG.sum), drop = FALSE]
    counts.matrix.OG.maternal <- as.data.frame( counts.matrix.OG.maternal)
    stopifnot(
    identical(
        row.names(counts.matrix.OG.maternal),
        row.names(counts.matrix.OG.sum)
    )
    )
    stopifnot(
    identical(
        colnames(counts.matrix.OG.maternal),
        colnames(counts.matrix.OG.sum)
    )
    )
    message(
    "dim of aligned maternal matrix: ",
    nrow(counts.matrix.OG.maternal),
    " × ",
    ncol(counts.matrix.OG.maternal)
    )
    message(
    "Gene order identical: ",
    identical(
        rownames(counts.matrix.OG.maternal),
        rownames(counts.matrix.OG.sum)
    )
    )
    message(
    "Cell order identical: ",
    identical(
        colnames(counts.matrix.OG.maternal),
        colnames(counts.matrix.OG.sum)
    )
    )
    write.table(counts.matrix.OG.maternal,file=paste0(file_dir, sample_name, ".counts.matrix.OG.maternal.txt"), sep = "\t")

    ### generate paternal matrix ###
    counts.matrix.OG.paternal <- counts.matrix.OG[grepl("^Csav\\.dovetail", counts.matrix.OG$Gene), , drop = FALSE]
    message("dim of counts.matrix.OG.paternal: ",nrow(counts.matrix.OG.paternal), "x", ncol(counts.matrix.OG.paternal)) 
    counts.matrix.OG.paternal <- counts.matrix.OG.paternal[,-1]
    counts.matrix.OG.paternal <- as.data.frame(counts.matrix.OG.paternal)
    rownames(counts.matrix.OG.paternal) <- counts.matrix.OG.paternal$ortho_grp
    counts.matrix.OG.paternal <- counts.matrix.OG.paternal[, -ncol(counts.matrix.OG.paternal)]
    message("dim of counts.matrix.OG.paternal after deleting: ",nrow(counts.matrix.OG.paternal), "x", ncol(counts.matrix.OG.paternal)) 
    colnames(counts.matrix.OG.paternal) <- gsub("\\.", "-", colnames(counts.matrix.OG.paternal)) 
    missing_genes <- setdiff(rownames(counts.matrix.OG.sum),rownames(counts.matrix.OG.paternal))
    message("Number of genes missing from maternal matrix: ", length(missing_genes))

    if (length(missing_genes) > 0) {
    missing_gene_matrix <- matrix(
        0,
        nrow = length(missing_genes),
        ncol = ncol(counts.matrix.OG.paternal),
        dimnames = list(
        missing_genes,
        colnames(counts.matrix.OG.paternal)
        )
    )
    ### add missing genes to paternal matrix ###
    counts.matrix.OG.paternal <- rbind(
        counts.matrix.OG.paternal,
        missing_gene_matrix
    )
    }
    counts.matrix.OG.paternal <- counts.matrix.OG.paternal[rownames(counts.matrix.OG.sum), colnames(counts.matrix.OG.sum), drop = FALSE]
    counts.matrix.OG.paternal <- as.data.frame( counts.matrix.OG.paternal)
    stopifnot(
    identical(
        row.names(counts.matrix.OG.paternal),
        row.names(counts.matrix.OG.sum)
    )
    )
    stopifnot(
    identical(
        colnames(counts.matrix.OG.paternal),
        colnames(counts.matrix.OG.sum)
    )
    )
    message(
    "dim of aligned paternal matrix: ",
    nrow(counts.matrix.OG.paternal),
    " × ",
    ncol(counts.matrix.OG.paternal)
    )
    message(
    "Gene order identical: ",
    identical(
        rownames(counts.matrix.OG.paternal),
        rownames(counts.matrix.OG.sum)
    )
    )
    message(
    "Cell order identical: ",
    identical(
        colnames(counts.matrix.OG.paternal),
        colnames(counts.matrix.OG.sum)
    )
    )
    write.table(counts.matrix.OG.paternal,file=paste0(file_dir, sample_name, ".counts.matrix.OG.paternal.txt"), sep = "\t")
}

# ============================================================
# filter gene and cells
# ============================================================
cell.info <- read.csv("/data/CsCi.10s.cell.anno.region.csv")
colnames(cell.info)[1] <- "cell"
replicates <- c("LTB2_Cs1Ci1", "LTB2_Cs1Ci2", "LTB2_Cs2Ci1") 
file_dir <- paste0("/count_matrices/")
for(rep in replicates){
  message("Processing: ", rep)
  total.count <- read.table(paste0(file_dir, rep, ".counts.matrix.OG.total.txt"), sep = "\t", row.names = 1, check.names = FALSE)
  maternal.count <- read.table(paste0(file_dir, rep, ".counts.matrix.OG.maternal.txt"), sep = "\t", row.names = 1, check.names = FALSE)
  cell_info_use <- cell.info[cell.info$orig.ident == rep, , drop = FALSE]
  cell_info_use$cell <- sub("_[0-9]+$", "", cell_info_use$cell)
  write.csv(cell_info_use, file = paste0(file_dir, rep, "_cell_info.csv"), row.names = FALSE, quote = FALSE)
  cells_use <- cell_info_use$cell
  sub_mat <- total.count[, cells_use, drop = FALSE]
  sub_mat_maternal <- maternal.count[, cells_use, drop = FALSE]
  message("mapping cells: ", "total: ", nrow(sub_mat), " x ", ncol(sub_mat))
  message("mapping cells: ", "maternal: ", nrow(sub_mat_maternal), " x ", ncol(sub_mat_maternal))
  # filter
  gene_sum <- rowSums(sub_mat)
  sub_mat <- sub_mat[gene_sum > 30, , drop = FALSE]
  sub_mat_maternal <- sub_mat_maternal[rownames(sub_mat), , drop = FALSE]
  message("filtering genes: ", "total: ", nrow(sub_mat), " x ", ncol(sub_mat))
  message("filtering genes: ", "maternal: ", nrow(sub_mat_maternal), " x ", ncol(sub_mat_maternal))
  # output
  write.table(sub_mat,
              file = paste0(file_dir, rep, "_cutoff30.counts.f.total.txt"),
              sep = "\t", quote = FALSE)
  write.table(sub_mat_maternal,
              file = paste0(file_dir, rep, "_cutoff30.counts.f.maternal.txt"),
              sep = "\t", quote = FALSE)
}


# ============================================================
# create seurat object
# ============================================================
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
replicates <- c("LTB2_Cs1Ci1", "LTB2_Cs1Ci2", "LTB2_Cs2Ci1") 
file_dir <- paste0("/count_matrices/")
obj_output_dir <- paste0("/exp_seurat_obj/")

for(rep in replicates){
    # seurat obj
    message("Processing: ", rep)
    mat.count <- read.table(paste0(file_dir, rep,  "_cutoff30.counts.f.total.txt"), sep = "\t", row.names = 1, check.names = FALSE)
    seu <- CreateSeuratObject(counts = mat.count, project = rep, min.cells = 3, min.features = 200)
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu)
    pca.coord <- Embeddings(seu, reduction = "pca") #50PCs
    write.csv(pca.coord, file = paste0(obj_output_dir, rep, "_cutoff30_PCA_dim50_coordinates.csv"), quote = FALSE)
    seu <- RunUMAP(seu, dims = 1:30)
    # check pca
    cell.anno <- read.csv(paste0(file_dir, rep, "_cell_info.csv"))
    rownames(cell.anno) <- cell.anno$cell
    cell.anno <- cell.anno[colnames(seu), ]
    seu$tissue.beau.adj <- cell.anno$tissue.beau.adj
    pdf(paste0(obj_output_dir, rep, "_cutoff30_pca_scatter.pdf"), width = 10, height = 8)
    print(DimPlot(object = seu, reduction = "pca", group.by = "tissue.beau.adj", pt.size = 0.5, label = TRUE, repel = TRUE) + ggtitle(rep))
    dev.off()
    pdf(paste0(obj_output_dir, rep, "_cutoff30_tissue_umap.pdf"), width = 10, height = 8)
    print(DimPlot(object = seu, reduction = "umap", group.by = "tissue.beau.adj", pt.size = 0.5, label = TRUE, repel = TRUE) + ggtitle(rep))
    dev.off()
    saveRDS(seu, file = paste0(obj_output_dir, rep, "_cutoff30_obj.rds"))
    rm(seu)
    gc()
}
 



