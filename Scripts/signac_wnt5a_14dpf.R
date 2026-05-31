#==============================================================
# Analyze 14 dpf ATAC CNCCs dataset with focus on wnt5a locus
#==============================================================

####################
### Introduction ###
####################

# wnt5a locus has very specific chromatin activity based on
# snATACseq data. Aim of this script is to reanalyze the these
# data with special focus on this locus.

# NOTE: This specific script uses only 14 dpf stage.


#####################################
### Preparations and data loading ###
#####################################

# Load necessary packages
library(Signac)
library(Seurat)
library(GenomicRanges)
library(Matrix)
library(AnnotationHub)
library(ggplot2)
library(patchwork)
library(cicero)

# Load the ATAC data
## Load separate files
counts <- Matrix::readMM("Data/snATACseq_data/filtered_peak_bc_matrix/matrix.mtx")
barcodes <- readLines("Data/snATACseq_data/filtered_peak_bc_matrix/barcodes.tsv")
peaks <- read.csv("Data/snATACseq_data/filtered_peak_bc_matrix/peaks.bed", sep = "\t", header = FALSE)

## Extract peaks names
peaknames <- paste0(peaks[, 1], ":", peaks[, 2], "-", peaks[, 3])

## Add barcodes and peaknames to the count matrix
colnames(counts) <- barcodes
rownames(counts) <- peaknames

# Create chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "Data/snATACseq_data/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

# Create Seurat object
cncc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC"
)

# Use AnnotationHub to add Danio Rerio GRC11 rel. 98 genome
## Initialize AnnotationHub
ah <- AnnotationHub()

## Query the proper version
query_danre <- query(ah, c("Danio rerio", "98", "EnsDb")); query_danre
# NOTE: GRCz11 v98 AnnotationHub ID is AH74989

## Extract Ensembl genome
EnsDb_v98 <- ah[["AH74989"]]

## Extract annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb_v98)
genome(annotations) <- "danre98"

# Add gene information to the seurat object
Annotation(cncc) <- annotations

# Add pre-computed UMAP
## Read and reformat embedding coordinates
umap_data <- read.csv("Data/snATACseq_data/UMAP.csv")
umap_matrix <- as.matrix(umap_data[, 2:3])
colnames(umap_matrix) <- c("UMAP1", "UMAP2")
rownames(umap_matrix) <- umap_data[, 1]

## Filter only the cells that are in the precomputed UMAP
## NOTE: I suppose that the rest of the cells would be filtered
## during the QC step
cncc <- cncc[, rownames(umap_matrix)]

## Add the embedding to the Seurat object
cncc[["umap"]] <- CreateDimReducObject(
  embeddings = umap_matrix,
  assay = "ATAC",
  key = "UMAP"
)

# Load pre-computed clusters
## Read the clusters
LSI_snn <- read.csv("Data/snATACseq_data/LSI_snn.csv", row.names = "Barcode")

## Add the clusters to metadata
cncc <- AddMetaData(cncc, LSI_snn)

## NOTE: Use AddMetaData since it matches cell barcodes regardless the order

# Finally check the dimplots for different clustering resolutions
p1 <- DimPlot(cncc, group.by = "ATAC_snn_res.0.5")
p2 <- DimPlot(cncc, group.by = "ATAC_snn_res.1")
p3 <- DimPlot(cncc, group.by = "ATAC_snn_res.1.5")
p4 <- DimPlot(cncc, group.by = "ATAC_snn_res.2")
p1 + p2 + p3 + p4

# NOTE: Use resolution 1 (nice separation of wnt5a+ dermal fibroblasts)

# Create DF + Other cells "clusters"
cncc$DFvsRest <- ifelse(cncc$ATAC_snn_res.1 == "Cluster 19", "Dermal fibroblasts", "Other cells") 
DimPlot(cncc, group.by = "DFvsRest")

# Save the RDS object
saveRDS(cncc, paste0("Analysis/", Sys.Date(), "_cncc_ATAC_14dpf_precomputedUMAP.rds"))


#######################
### Explore dataset ###
#######################

# Set DF ns Rest as an active ident
Idents(cncc) <- "DFvsRest"

# Show wnt5a locus - specific wnt5a "promotor"
CoveragePlot(cncc, "wnt5a", group.by = "DFvsRest", extend.downstream = 2000)

# Show wnt5a locus - upstream "enhacers"
## NOTE: These peaks are really specific for DF cluster
CoveragePlot(cncc, "wnt5a", group.by = "DFvsRest", extend.downstream = 35000)


##################
### Run cicero ###
##################

# Extract count matrix from the cncc object
## NOTE: Dataset is filtered, so the original count matrix cannot be used
cncc_counts <- GetAssayData(
  cncc,
  assay = "ATAC",
  layer = "counts"
)

# Cicero prefers 1_1000_2000 format, but the data are now in 1-1000-2000 format
# Replace "-" with "_"
rownames(cncc_counts) <- gsub("-", "_", rownames(cncc_counts))
# NOTE: Check this. Kinda sus.

# Extract metadata for individual cells and peaks
cncc_metadata <- cncc@meta.data
peak_metadata <- data.frame(
  gene_short_name = rownames(cncc_counts), # Keep the colums naming
  row.names = rownames(cncc_counts)
)

# Create monocle CellDataSet (CDS)
pheno_data <- new("AnnotatedDataFrame", data = cncc_metadata)
feature_data <- new("AnnotatedDataFrame", data = peak_metadata)
cncc_cds <- newCellDataSet(
  cncc_counts,
  phenoData = pheno_data,
  featureData = feature_data,
  expressionFamily = binomialff()
)

# Add umap to CDS object
umap_coords <- Embeddings(cncc, "umap")
cncc_cds_cicero <- make_cicero_cds(
  cncc_cds,
  reduced_coordinates = umap_coords
)

# Use prefiously prepared annotationhub object to extract chromosome lengths
# NOTE: see line 67
genome_df <- seqlengths(EnsDb_v98)

# Keep only conventional chromosomes
genome_df <- genome_df[names(genome_df) %in% as.character(1:25)]

# Transform the vector into data.frame
genome_df <- data.frame(
  Chromosome = names(genome_df),
  Length = genome_df
)

# Run Cicero
# NOTE: This calculation takes cca 60 minutes
conns <- cicero::run_cicero(
  cncc_cds_cicero,
  genomic_coords = genome_df
)

# Save cicero results
saveRDS(conns, paste0("Analysis/", Sys.Date(), "_cncc_ATAC_14dpf_conns.rds"))

# Filter only strong links
conns_filter <- subset(conns, coaccess > 0.25)

# Calculate cis-coaccesibility networks (CCANs)
ccans <- generate_ccans(conns_filter)

# Fix "_" to "-"
conns_filter[, 1] <- gsub("_", "-", conns_filter[, 1])
conns_filter[, 2] <- gsub("_", "-", conns_filter[, 2])

# Convert connections to peak-to-peak network and add them to Seurat object
conns_links <- ConnectionsToLinks(conns = conns_filter, ccans = ccans)
Links(cncc) <- conns_links

# Save the seurat object with cicero links
saveRDS(cncc, paste0("Analysis/", Sys.Date(), "_cncc_ATAC_14dpf_SeuratObject_ciceroLinks.rds"))

pdf("Outs/wnt5a_distal.pdf", width = 8, height = 5)
CoveragePlot(cncc, "wnt5a", group.by = "DFvsRest", extend.upstream = 50000, extend.downstream = 50000)
dev.off()

