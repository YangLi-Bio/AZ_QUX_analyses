# AZ_QUX_analyses

This repository aims at all scRNA-Seq analysis codes for the QUX dataset from AstraZeneca.

# Pipeline:
Use Scanpy to integrate multiple scRNA-Seq datasets
Perform transformation and PCA, UMAP
Save the files using any one of the following two methods:
1. Saved as annData --> load annData in Python
2. Saved as triplets (matrix.mtx, barcodes.tsv, and genes.tsv) --> convert them to exprMatrix.tsv.gz using cbTools mtx2tsv
