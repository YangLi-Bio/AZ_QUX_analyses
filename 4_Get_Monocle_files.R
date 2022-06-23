#####################################################################
#                                                                   #
#                  Get files in Monocle3 format                     #
#                                                                   #
#####################################################################


library(Seurat)
work.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/20220524_DITQUX_145_INSERM_Raw/"
Seurat3.dir <- paste0(work.dir, "/Seurat3_files/")
ifnb.integrated <- qs::qread(paste0(work.dir, "PCA_UMAP_rPCA_integrated_object.qsave"))


# Add project data
barcode.proj <- qs::qread(paste0(Seurat3.dir, "Barcodes_project.qsave"))
barcode.proj <- factor(barcode.proj, levels = unique(barcode.proj))
ifnb.integrated <- AddMetaData(ifnb.integrated, metadata = barcode.proj, col.name = "Project")


# Get Monocle data
library(monocle3)
library(SeuratWrappers)
DefaultAssay(ifnb.integrated) <- "RNA" # integrated
cds <- as.cell_data_set(ifnb.integrated)
dim(cds)
monocle.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/Monocle_RNA_dir/"
dir.create(monocle.dir)
library(Matrix)
qs::qsave(cds, paste0(monocle.dir, "cds_RNA.qsave"))

writeMM(exprs(cds), paste0(monocle.dir, "matrix.mtx"))

write.table(as.data.frame(cbind(rownames(exprs(cds)), 
                                rownames(exprs(cds)))), 
            file = paste0(monocle.dir, "features.tsv"), sep = "\t", 
            row.names = F, col.names = F, quote = F)

write(colnames(exprs(cds)), file = paste0(monocle.dir, "barcodes.tsv"))

write.table(data.frame(cell = names(ifnb.integrated$Project), 
                       project = ifnb.integrated$Project), 
            file=paste0(monocle.dir, "meta.tsv"), 
            quote=FALSE, sep='\t', col.names = NA)

write.table(ifnb.integrated@reductions$umap@cell.embeddings, 
            file=paste0(monocle.dir, "coords.tsv"), 
            quote=FALSE, sep='\t', col.names = NA)
