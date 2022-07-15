###############################################################
#                                                             #
# This script fulfills the following tasks:                   #
# 1. Find cell-cluster-specific markers and save them to      #
#    csv files                                                #
#                                                             #
###############################################################


# Scripts, libraries, and parameters
library(Seurat)
library(dplyr)
library(future)
R.dir <- "/fs/ess/scratch/PCON0022/liyang/AstraZeneca/QUX/Rfiles/"
integrated.file <- "Obj_Inserm_Jose_2batches.qsave"


# # Load the Seurat object
# integrated <- qs::qread(paste0(R.dir, integrated.file))
# integrated$seurat_clusters
# length(unique(integrated$seurat_clusters))
# unique(integrated@meta.data$orig.ident)
# unique(integrated$Project)
# head(integrated$Project)
# 
# 
# # Set the sample IDs for Jose's data
# Jose.list <- integrated@meta.data$orig.ident
# length(Jose.list)
# names(Jose.list) <- rownames(integrated@meta.data)
# head(Jose.list)
# Jose.list <- Jose.list[which(Jose.list != "QUX")]
# length(Jose.list)
# head(Jose.list)
# 
# 
# # Set the sample IDs for Inserm data
# Inserm.list <- integrated$Project
# length(Inserm.list)
# head(Inserm.list)
# unique(Inserm.list)
# Inserm.list <- Inserm.list[which(!is.na(Inserm.list))]
# length(Inserm.list)
# unique(Inserm.list)
# head(Inserm.list)
# 
# 
# # Check whether cell sets are consistent
# integrated.meta <- c(Jose.list, Inserm.list)
# length(integrated.meta) == ncol(integrated)
# setdiff(names(integrated.meta), colnames(integrated))
# setdiff(colnames(integrated), names(integrated.meta))
# 
# 
# # Reset the sample IDs
# integrated <- AddMetaData(integrated, metadata = integrated.meta, 
#                           col.name = "Project")
# length(unique(integrated$Project))
# unique(integrated$Project)
# qs::qsave(integrated, paste0(R.dir, "Inserm_Jose_integrated_add_sample_IDs.qsave"))
  

# DEG analysis
integrated <- qs::qread(paste0(R.dir, integrated.file))
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- integrated$seurat_clusters
table(Idents(integrated))
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)
integrated.marmers <- FindAllMarkers(integrated, only.pos = F)
dim(integrated.marmers)
qs::qsave(integrated.marmers, paste0(R.dir, "Inserm_Jose_23samples_markers.qsave"))


# Session information
sessionInfo()