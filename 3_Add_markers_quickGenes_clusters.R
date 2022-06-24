######################################################
#                                                    #
#         Generate the following three files:        #
#         1. meta.tsv                                #
#         2. markers.tsv                             #
#         3. quickGenes.csv                          #
#                                                    #
######################################################


# Libraries and parameters
library(Seurat)
library(dplyr)
out.dir <- "/fs/ess/scratch/PCON0022/liyang/AstraZeneca/QUX/Rfiles/"
html.dir <- "/fs/ess/scratch/PCON0022/liyang/AstraZeneca/QUX/html_files/RDagher_CB/"
marker.inFile <- "Cell_cluster_markers_parallel.qsave"
obj.file <- "Cell_clusters_res_0.5.qsave"


# Generate the meta file composed of barcodes, project, and cluster
obj <- qs::qread(paste0(out.dir, obj.file))
identical(names(obj$Project), names(obj$seurat_clusters))
meta.df <- data.frame(cell = names(obj$Project), project = obj$Project, 
                      cluster = obj$seurat_clusters)
dim(meta.df)
head(meta.df)
write.table(meta.df, paste0(html.dir, "meta.tsv"), sep = "\t", quote = F, row.names = F)
message ("meta.tsv")


# Generate markers.tsv file composed of cluster, gene,	avg_diff, and	p_val
marker.df <- qs::qread(paste0(out.dir, marker.inFile))
dim(marker.df)
head(marker.df)
marker.df <- marker.df[marker.df$p_val_adj < 0.05 & (marker.df$avg_log2FC > 0.25 |
                                                   marker.df$avg_log2FC < -0.25),]
dim(marker.df)
head(marker.df)
markerOut.df <- data.frame(cluster = marker.df$cluster, gene = marker.df$gene, 
                           avg_diff = marker.df$avg_log2FC, 
                           p_val = marker.df$p_val_adj)
dim(markerOut.df)
head(markerOut.df)
write.table(markerOut.df, paste0(html.dir, "markers.tsv"), 
            sep = "\t", row.names = F, 
            quote = F)
message ("markers.tsv")


# Generate quickGenes.csv composed of gene symbols
df.list <- split(markerOut.df, f = markerOut.df$cluster)
length(df.list)
names(df.list)
sapply(df.list, nrow)
# Note that: clusters 24, 29, 32, 33, and 34 have no significant DEGs
quick.genes <- do.call("c", lapply(df.list, function(x) {
  x[order(-abs(x$avg_diff), x$p_val),] %>% pull("gene") %>% head(n = 50)
}))
length(quick.genes)
head(quick.genes)
names(quick.genes) <- NULL
head(quick.genes)
quickGenes.df <- data.frame(symbol = unique(quick.genes))
dim(quickGenes.df)
write.csv(quickGenes.df, paste0(html.dir, "quickGenes.csv"), quote = F, 
          row.names = F)
message ("quickGenes.csv")
