#############################################################
#                                                           #
#                 Generate UMAP plot                        #
#                                                           #
#############################################################


# Load data
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/QUX/Tables/"


# Generate UMAP
source("/fs/ess/PCON0022/liyang/r_utilities/functions/visual_tools.R")
library(Seurat)
Idents(ifnb.integrated) <- ifnb.integrated$Project
DefaultAssay(ifnb.integrated) <- "integrated"
integrated.umap <- get_UMAP(object = ifnb.integrated, reduction.method = "umap", 
         txt = "Project", legend.position = "right")
png(filename = paste0(image.dir, "UMAP.png"), width = 5000, height = 4000, 
    res = 300)
print(integrated.umap)
dev.off()


# Print basic information
project.dt <- table(ifnb.integrated$Project)
project.dt <- data.frame(Project = names(project.dt), 
                         Cells = project.dt)
project.dt <- project.dt[, c(1, 3)]
colnames(project.dt) <- c("Project", "Cells")
write.csv(project.dt, paste0(table.dir, "Project_info.csv"), row.names = F, 
          quote = F)
