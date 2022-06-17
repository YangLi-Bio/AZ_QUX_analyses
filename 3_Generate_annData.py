# Import modules
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)


# Load data
work_dir = "/fs/ess/PCON0022/liyang/astrazeneca/QUX/20220524_DITQUX_145_INSERM_Raw/csv_files/"
out_dir = "/fs/ess/PCON0022/liyang/astrazeneca/QUX/"
counts = pd.read_csv(f"{work_dir}/RNA_matrix.csv", header = None, engine = "pyarrow")
adata = ad.AnnData(counts)
adata.X
len(adata.X)
len(adata.X[1,:])
adata.to_pickle(f"{out_dir}/Matrix_annData.pkl.gz")


# Load observations
obs = pd.read_csv(f"{work_dir}/Barcodes_matrix.csv", header = 0)
adata.obs_names = obs["Barcodes"]
adata.obs["project"] = obs["Project"]


# Load variables
var = pd.read_csv(f"{work_dir}/Barcodes_matrix.csv", header = 0)
adata.var_names = var["Features"]


# Load UMAP embeddings
umap = pd.read_csv(f"{work_dir}/UMAP_matrix.csv", header = 0)
adata.obsm["X_umap"] = umap
adata.to_pickle(f"{out_dir}/Total_annData.pkl.gz")

