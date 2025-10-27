library(Seurat)
library(ggplot2)
library(writexl)
library(tidyverse)
library(qs)
library(Matrix)

# Example folder with data
# Need to mkdir QC and RDS (actually QS)

setwd ("~/Desktop/sc_TBI/new_data/GSE262317_RAW")

file_list <- list.files(pattern = ".barcodes.tsv.gz$")
file_list <- sub(".barcodes.tsv.gz$", "", file_list)


mad_outlier = function(X, metric, nmads){
  M = X@meta.data[[metric]]
  median_M = median(M, na.rm = TRUE)
  mad_M = mad(M, na.rm = TRUE)
  outlier = (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))
  return(outlier)}

Seurat_Object_QC_list <- function(Samples) {
  for (Sample in Samples) {
    message("Processing sample: ", Sample)
    
    # --- Step 1: Load data ---
    X <- ReadMtx(
      mtx = paste0(Sample, "_matrix.mtx.gz"),
      cells = paste0(Sample, "_barcodes.tsv.gz"),
      features = paste0(Sample, "_features.tsv.gz")
    )
    seurat_obj <- CreateSeuratObject(counts = X, project = Sample)
    rm(X); gc()
    
    # --- Step 2: Deduplicate genes (sum across .1, .2, etc.) ---
    counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
    dup_genes <- grep("\\.", rownames(counts), value = TRUE)
    
    if (length(dup_genes) > 0) {
      base_names <- sub("\\..*", "", dup_genes)
      
      for (gene in unique(base_names)) {
        matching_rows <- grep(paste0("^", gene, "(\\.|$)"), rownames(counts), value = TRUE)
        
        if (length(matching_rows) > 1) {
          summed_counts <- Matrix::colSums(counts[matching_rows, , drop = FALSE])
          counts <- counts[!(rownames(counts) %in% matching_rows), , drop = FALSE]
          counts <- rbind(counts, Matrix::Matrix(summed_counts, nrow = 1, sparse = TRUE, dimnames = list(gene, NULL)))
        }
      }
      
      # Update Seurat object with deduplicated counts
      seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "counts", new.data = counts)
      seurat_obj <- subset(seurat_obj, features = rownames(counts))
    }
    
    message("Deduplication done. Number of genes: ", nrow(seurat_obj))
    
    # --- Step 3: QC metrics ---
    Cells_Before <- length(Cells(seurat_obj))
    seurat_obj$log1p_total_counts <- log1p(seurat_obj@meta.data$nCount_RNA)
    seurat_obj$log1p_n_genes_by_counts <- log1p(seurat_obj@meta.data$nFeature_RNA)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
    
    # --- Step 4: MAD filtering and mitochondrial threshold ---
    bool_vector <- !mad_outlier(seurat_obj, 'log1p_total_counts', 5) &
      !mad_outlier(seurat_obj, 'log1p_n_genes_by_counts', 5)
    seurat_obj <- subset(seurat_obj, cells = which(bool_vector))
    seurat_obj <- subset(seurat_obj, subset = percent.mt < 2)
    
    # --- Step 5: Save QC summary and filtered object ---
    Count <- data.frame(
      Base = Cells_Before,
      QC = length(Cells(seurat_obj)),
      ID = Sample
    )
    
    dir.create("QC", showWarnings = FALSE)
    dir.create("RDS", showWarnings = FALSE)
    
    write.csv(Count, paste0("QC/", Sample, "_QC_Cells.csv"), row.names = FALSE)
    qsave(seurat_obj, paste0("RDS/", Sample, "_AfterQC.qs"))
    
    message("Finished sample: ", Sample, " ✓")
  }
}

Seurat_Object_QC_list (file_list)

