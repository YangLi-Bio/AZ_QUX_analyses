###############################################################
#                                                             #
# This script fulfills the following tasks:                   #
# 1. Integrate Inserm and Jose's data                         #
# 2. Perform dimension reduction and UMAP embedding           #
# 3. Perform cell clustering                                  #
# 4. DEG analysis                                             #
# 5. Prepare the data files for Cell Browser                  #
#                                                             #
###############################################################


# Parameters
library(Seurat)
library(pbapply)
library(dplyr)
code.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/scratch/PCON0022/liyang/AstraZeneca/QUX/Rfiles/"
Inserm.file <- "Cell_clusters_res_0.5.qsave"
Jose.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/Filtered_data/Jose_COPD/"

source(paste0(code.dir, "transcriptome_tools.R"))


# # Load Jose's data
# dir.list <- list.dirs(path = Jose.dir)
# dir.list
# dir.list <- dir.list[c(3, 5, 7)]
# dir.list
# names <- strsplit(dir.list, split = "\\/") %>% sapply(., "[[", (11)) %>% paste0("Jose_", .)
# names
# obj.list <- pblapply(seq_along(dir.list), function(i) {
#   CreateSeuratObject(counts = Read10X(data.dir = dir.list[[i]]), 
#                      project = names[[i]], min.cells = 3, min.features = 200)
# })
# qs::qsave(obj.list, paste0(R.dir, "Jose_object_list.qsave"))
# length(obj.list)
# sapply(obj.list, ncol)
# sapply(obj.list, nrow)
# 
# 
# # Integrate Jose's data
# Jose.integrated <- integrate_scRNA(ifnb.list = obj.list, name.list = names)
# qs::qsave(Jose.integrated, paste0(R.dir, "Jose_integrated.qsave"))


# Integrate Inserm and Jose's data
Jose.integrated <- qs::qread(paste0(R.dir, "Jose_integrated.qsave"))
Inserm.integrated <- qs::qread(paste0(R.dir, Inserm.file))
dim(Inserm.integrated)
DefaultAssay(Inserm.integrated) <- "RNA"
length(unique(Inserm.integrated$Project))
DefaultAssay(Jose.integrated) <- "RNA"
Jose.integrated@project.name
integrated <- integrate_scRNA(ifnb.list = list(Inserm = Inserm.integrated, 
                                                    Jose = Jose.integrated), 
                              clustering = T, umap = T)
qs::qsave(integrated, paste0(R.dir, "Inserm_Jose_integrated_cluster_0.5.qsave"))
