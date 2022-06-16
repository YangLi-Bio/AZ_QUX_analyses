#####################################################################
#                                                                   #
#                   Extract matrix and metadata                     #
#                                                                   #
#####################################################################


work.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/20220524_DITQUX_145_INSERM_Raw/"
ifnb.integrated <- qs::qread(paste0(work.dir, "PCA_UMAP_rPCA_integrated_object.qsave"))


#####################################################################
#                                                                   #
#                 Get expression matrix and metadata                #
#                                                                   #
#####################################################################


library(Seurat)
names(ifnb.integrated@reductions)
# rna.mat <- GetAssayData(object = ifnb.integrated, assay = "RNA", slot = "data")


# Save matrix
integrated.mat <- GetAssayData(object = ifnb.integrated, assay = "integrated", slot = "data")
class(integrated.mat)
dim(integrated.mat)
integrated.mat[1:5, 1:5]
Seurat3.dir <- paste0(work.dir, "/Seurat3_files/")
dir.create(Seurat3.dir)
qs::qsave(integrated.mat, paste0(Seurat3.dir, "Integrated_matrix.qsave"))


# Save embeddings
umap.mat <- ifnb.integrated@reductions$umap@cell.embeddings
qs::qsave(umap.mat, paste0(Seurat3.dir, "Integrated_UMAP.qsave"))


# Load projects
dir.list <- list.dirs(work.dir)
dir.list <- dir.list[c(3, 5, seq(7, 49, by = 3))]
head(dir.list)
library(dplyr)
project.list <- strsplit(dir.list, split = "\\/") %>% sapply(., "[", 10)


# Build correspondence between barcodes and projects
obj.list <- qs::qread(paste0(work.dir, "Object_list.qsave"))
barcode.proj <- Reduce("c", lapply(seq_along(obj.list), function(i) {
  projects <- rep(project.list[[i]], ncol(obj.list[[i]]))
  names(projects) <- colnames((obj.list[[i]]))
  projects
}))
qs::qsave(barcode.proj, paste0(Seurat3.dir, "Barcodes_project.qsave"))
ncol(ifnb.integrated) == length(barcode.proj)
rm(obj.list)


# Write csv files
# dir.create(paste0(work.dir, "csv_files/"))
# csv.dir <- paste0(work.dir, "csv_files/")
# write.csv(integrated.mat, paste0(csv.dir, "Integrated_matrix.csv"), quote = F)
# write.csv(umap.mat, paste0(csv.dir, "UMAP_matrix.csv"), quote = F)
# 
# barcode.df <- data.frame(Barocdes = names(barcode.proj), Project = barcode.proj)
# write.csv(barcode.df, paste0(csv.dir, "Barcodes_matrix.csv"), quote = F)


# Prepare annData for scanpy
library(anndata)
