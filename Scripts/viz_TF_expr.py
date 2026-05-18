#============================================
# Visualize TF expression in scRNAseq data
#============================================

####################
### Introduction ###
####################

# Using homer analysis, sevaral TF motifs were find in wnt5a peaks sequence.
# The goul of this script is to visualize the TF factors in dermal fibroblasts,
# to get better understanding of wnt5a regulatory mechanism.


###################
### Preparation ###
###################

# Load modules
import scanpy as sc
from matplotlib import pyplot as plt

# Load the data
cncc = sc.read_h5ad("Data/scRNAseq_data/zebrafish_14dpf.h5ad")


##########################
### Data preprocessing ###
##########################

# Quality control calculation

## Find mitochondrial genes
cncc.var["mt"] = cncc.var_names.str.startswith("mt-")
## Find ribosomal genes
cncc.var["rb"] = cncc.var_names.str.startswith(("rps", "rpl"))
## Find hemoglobin genes
cncc.var["hb"] = cncc.var_names.str.contains("^hb[^(p)]")
## Calculate the qc metric
sc.pp.calculate_qc_metrics(cncc, qc_vars = ["mt", "rb", "hb"], inplace = True, log1p = True)

# Inspect violin plots
sc.pl.violin(cncc, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter = 0.4, multi_panel = True)

# Inspect the scatter plot
sc.pl.scatter(cncc, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# NOTE: No mt genes were identified (ensembl ids). However the data looks solid.

# Filter the cells
sc.pp.filter_cells(cncc, min_genes=100)
sc.pp.filter_genes(cncc, min_cells=3)

# Find doublets
sc.pp.scrublet(cncc, batch_key = "batch")

# Save count data
cncc.layers["counts"] = cncc.X.copy()

# Normalize to median total counts
sc.pp.normalize_total(cncc)

# Logarithmize the data
sc.pp.log1p(cncc)

# Feature selection
sc.pp.highly_variable_genes(cncc, n_top_genes = 2000, batch_key = "batch")

# Check the varible genes
sc.pl.highly_variable_genes(cncc)


###############################################
### Dimensionality reduction and clustering ###
###############################################

# Run the PCA
sc.tl.pca(cncc)

# Check the variance ratio
sc.pl.pca_variance_ratio(cncc, n_pcs = 50, log = True)

# NOTE: Use first 30 PCs

# Plot the PC1 and PC2
sc.pl.pca(cncc, color = ["batch"])

# Nearest neighbors graph construction
sc.pp.neighbors(cncc)

# Run umap dimensionality reguction
sc.tl.umap(cncc)

# Visulize umap
sc.pl.umap(cncc, color = "batch", size = 50)

# Cluster the cells
sc.tl.leiden(cncc, flavor = "igraph", n_iterations = 2, resolution = 0.3)

# Check the results of the clustering
sc.pl.umap(cncc, color = "leiden", size = 50)

# Re-asses the qulity control with doublets
sc.pl.umap(cncc, color = ["leiden", "predicted_doublet", "doublet_score"], wspace = 0.5, size = 3)

# NOTE: No doublets detected (probably filtered previously)


##############################
### Visualize possible TFs ###
##############################

# Select genes to be visualized
gene_list = {"wnt5a": ["ENSDARG00000104973"], "Selected TFs": ["ENSDARG00000042904", "ENSDARG00000059483", "ENSDARG00000042032", "ENSDARG00000095896", "ENSDARG00000062420"]}

# Plot the dotplot
sc.pl.dotplot(cncc, gene_list, groupby = "leiden", standard_scale = "var")
