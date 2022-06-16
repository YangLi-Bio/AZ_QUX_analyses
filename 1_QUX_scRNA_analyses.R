####################################################################
#                                                                  #
#           Perform scRNA-Seq analyses including:                  #
#           1. Integration of multiple samples                     #
#           2. Dimension reduction                                 #
#           3. Cell clustering                                     #
#                                                                  #
####################################################################


# Parameters
work.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/20220524_DITQUX_145_INSERM_Raw/"
library(Seurat)


# # # Load data into matrices
# # dir.list <- list.dirs(work.dir)
# # dir.list
# # length(dir.list)
# # dir.list <- dir.list[c(3, 5, seq(7, 49, by = 3))]
# # length(dir.list)
# # ifnb.list <- c()
# # for (i in seq_along(dir.list)) {
# #   message ("Loading the ", i, "-th matrix ...\n")
# #   ifnb.list[[i]] <- CreateSeuratObject(counts = Seurat::Read10X(dir.list[i]),
# #                                     project = "QUX", min.cells = 3, min.features = 200)
# # }
# # library(qs)
# # qs::qsave(ifnb.list, paste0(work.dir, "Object_list.qsave"))
# ifnb.list <- qs::qread(paste0(work.dir, "Object_list.qsave"))
# 
# 
# # # Size of all matrices
# # sapply(ifnb.list, dim)
# # sum(sapply(ifnb.list, ncol))
# 
# 
# # # Performing integration on datasets normalized with SCTransform
# # ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
# 
# 
# # After acquiring the data, we first perform standard normalization and variable feature selection.
# # library(future)
# # library(future.apply)
# # plan("multicore", workers = 4)
# # options(future.globals.maxSize = 20000 * 1024^2)
# library(pbapply)
# ifnb.list <- pblapply(X = ifnb.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = FALSE)
#   x <- FindVariableFeatures(x, verbose = FALSE)
# })
# qs::qsave(ifnb.list, paste0(work.dir, "Normalized_HVG_object_list.qsave"))
# 
# 
# # features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
# # qs::qsave(features, paste0(work.dir, "Integration_features.qsave"))
# # ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
# # qs::qsave(ifnb.list, paste0(work.dir, "Object_list.qsave"))
# # # ifnb.list <- qs::qread(paste0(work.dir, "Object_list.qsave"))
# 
# 
# # Next, select features for downstream integration, and run PCA on each object 
# # in the list, which is required for running the reciprocal PCA workflow.
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# qs::qsave(features, paste0(work.dir, "rPCA_integration_features.qsave"))
# ifnb.list <- pblapply(X = ifnb.list, FUN = function(x) {
#   x <- ScaleData(x, features = features, verbose = FALSE)
#   x <- RunPCA(x, features = features, verbose = FALSE)
# })
# qs::qsave(ifnb.list, paste0(work.dir, "PCA_object_list.qsave"))
# 
# 
# # rPCA integration
# anchors <- FindIntegrationAnchors(object.list = ifnb.list, reference = c(1, 2), reduction = "rpca", 
#                                   dims = 1:30)
# ifnb.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
# qs::qsave(ifnb.integrated, paste0(work.dir, "rPCA_integrated_object.qsave"))
# 
# 
# # # CCA integration
# # immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
# #                                          anchor.features = features)
# # qs::qsave(immune.anchors, paste0(work.dir, "Integration_anchors.qsave"))
# # immune.anchors <- qs::qread(paste0(work.dir, "Integration_anchors.qsave"))
# # immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
# # qs::qsave(immune.combined.sct, paste0(work.dir, "Integrated_object.qsave"))


# Downstream analyses
ifnb.integrated <- qs::qread(paste0(work.dir, "rPCA_integrated_object.qsave"))
ifnb.integrated <- ScaleData(ifnb.integrated, verbose = FALSE)
ifnb.integrated <- RunPCA(ifnb.integrated, verbose = FALSE)
ifnb.integrated <- RunUMAP(ifnb.integrated, dims = 1:30)
qs::qsave(ifnb.integrated, paste0(work.dir, "PCA_UMAP_rPCA_integrated_object.qsave"))
work.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/20220524_DITQUX_145_INSERM_Raw/"
ifnb.integrated <- qs::qread(paste0(work.dir, "PCA_UMAP_rPCA_integrated_object.qsave"))


# p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
# p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
#               repel = TRUE)
# p1 + p2


# Session information
sessionInfo()
