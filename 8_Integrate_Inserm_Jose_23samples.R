###############################################################
#                                                             #
# Integrate the 17 Inserm samples, 3 COPD samples and 3       #
# healthy samples (23 samples)                                #
#                                                             #
###############################################################


# Libraries
library(Seurat)
library(pbapply)


# Gloabl parameters
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/R_files/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/Tables/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"


# Local parameters
integrated.file <- "Inserm_Jose_integrated_add_sample_IDs.qsave"
Jose.ctrl.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/Filtered_data/Jose_ctrl/"


# Source scripts
source(paste0(tool.dir, "transcriptome_tools.R"))


# Load Inserm and three COPD samples
integrated <- qs::qread(paste0(R.dir, integrated.file))
ncol(integrated)
colnames(integrated@meta.data)
rownames(integrated@meta.data)
integrated@meta.data$orig.ident
names(integrated@assays)
unique(integrated$Project)
project.v <- integrated@meta.data[, "Project"]
head(project.v)
names(project.v) <- rownames(integrated@meta.data)
head(project.v)


# Load the three ctrl samples from Jose
dir.list <- list.dirs(path = Jose.ctrl.dir)
dir.list <- dir.list[-1]
dir.list
names <- strsplit(dir.list, split = "\\/") %>% sapply(., "[[", (11)) %>% paste0("Jose_", .)
names
obj.list <- pblapply(seq_along(dir.list), function(i) {
  CreateSeuratObject(counts = Read10X(data.dir = dir.list[[i]]),
                     project = names[[i]], min.cells = 3, min.features = 200)
})
qs::qsave(obj.list, paste0(R.dir, "Jose_ctrl_objects.qsave"))


# Integrate the three Jose ctrl samples
names <- strsplit(dir.list, split = "\\/") %>% sapply(., "[[", (11)) %>% paste0("Jose_ctrl_", .)
names
Jose.integrated <- integrate_scRNA(ifnb.list = obj.list, name.list = names)
dim(Jose.integrated)
names(Jose.integrated@assays)
colnames(Jose.integrated@meta.data)
Jose.projects <- Jose.integrated$orig.ident
qs::qsave(Jose.integrated, paste0(R.dir, "Jose_three_ctrl.qsave"))


# Integrate Inserm and Jose's three ctrl data
# Inserm.integrated <- integrated
Inserm.integrated <- qs::qread(paste0(R.dir, integrated.file))
Jose.integrated <- qs::qread(paste0(R.dir, "Jose_three_ctrl.qsave"))
rm(integrated)
dim(Inserm.integrated)
DefaultAssay(Inserm.integrated) <- "RNA"
length(unique(Inserm.integrated$Project))
DefaultAssay(Jose.integrated) <- "RNA"
Jose.integrated@project.name
integrated <- integrate_scRNA(ifnb.list = list(Inserm = Inserm.integrated, 
                                               Jose = Jose.integrated), 
                              clustering = T, umap = T)
qs::qsave(integrated, paste0(R.dir, "Inserm_Jose_COPD_ctrl_integrated_cluster_0.5.qsave"))


# Check data
integrated <- qs::qread(paste0(R.dir, "Inserm_Jose_COPD_ctrl_integrated_cluster_0.5.qsave"))
ncol(integrated)
colnames(integrated)
colnames(integrated@meta.data)
unique(integrated$Project)
project.v <- integrated$Project
head(project.v)
which(is.na(integrated$Project))
length(project.v)
Jose.integrated <- qs::qread(paste0(R.dir, "Jose_three_ctrl.qsave"))
ncol(Jose.integrated)
colnames(Jose.integrated@meta.data)


# Rename the samples
# ctrl.list <- qs::qread(paste0(R.dir, "Jose_ctrl_objects.qsave"))
# length(ctrl.list)
# names(ctrl.list)
# ctrl.v <- do.call("c", lapply(seq_along(ctrl.list), function(i) {
#   x <- ctrl.list[[i]]
#   x.names <- names(x$orig.ident)
#   x.names <- paste0(x.names, "_", i)
#   orig <- x$orig.ident
#   names(orig) <- x.names
#   orig
# }))
# head(ctrl.v)
# tail(ctrl.v)
# length(ctrl.v) == sum(sapply(ctrl.list, ncol))
length(project.v)
project.copd <- project.v[which(!is.na(project.v))]
length(project.copd)
unique(project.copd)
ctrl.v <- integrated$orig.ident
head(ctrl.v)
length(ctrl.v)
ctrl.v <- ctrl.v[which(ctrl.v != "QUX")]
head(ctrl.v)
length(ctrl.v)
project.renamed <- c(project.copd, ctrl.v)
length(intersect(names(project.renamed), colnames(integrated))) == ncol(integrated)
integrated.renamed <- AddMetaData(integrated, metadata = project.renamed, 
                                  col.name = "Project")
unique(integrated.renamed$Project)
qs::qsave(integrated.renamed, paste0(R.dir, "Obj_Inserm_Jose_2batches.qsave"))
