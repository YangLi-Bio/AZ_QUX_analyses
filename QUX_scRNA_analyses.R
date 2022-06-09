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


# Load data into matrices
dir.list <- list.dirs(work.dir)
dir.list
length(dir.list)
dir.list <- dir.list[c(3, 5, seq(7, 49, by = 3))]
length(dir.list)
ifnb.list <- c()
for (i in seq_along(dir.list)) {
  message ("Loading the ", i, "-th matrix ...\n")
  ifnb.list[[i]] <- CreateSeuratObject(counts = Seurat::Read10X(dir.list[i]),
                                    project = "QUX", min.cells = 3, min.features = 200)
}
library(qs)
qs::qsave(ifnb.list, paste0(work.dir, "Object_list.qsave"))


# Size of all matrices
sapply(ifnb.list, dim)
sum(sapply(ifnb.list, ncol))


# Performing integration on datasets normalized with SCTransform
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
qs::qsave(features, paste0(work.dir, "Integration_features.qsave"))
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
qs::qsave(ifnb.list, paste0(work.dir, "Object_list.qsave"))
# ifnb.list <- qs::qread(paste0(work.dir, "Object_list.qsave"))



immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
qs::qsave(immune.anchors, paste0(work.dir, "Integration_anchors.qsave"))
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
qs::qsave(immune.combined.sct, paste0(work.dir, "Integrated_object.qsave"))


immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
qs::qsave(immune.combined.sct, paste0(work.dir, "Integrated_object.qsave"))


# p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
# p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
#               repel = TRUE)
# p1 + p2


# Session information
sessionInfo()
