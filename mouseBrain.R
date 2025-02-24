library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)

# download the 1.3 Million Brain Cells from E18 Mice dataset (h5 output from 10x genomics).
# place the h5 file in your working directory.
brain.data <- open_matrix_10x_hdf5(
  path = "1M_neurons_filtered_gene_bc_matrices_h5.h5"
)

# Writing the matrix to a directory 'brain_counts'.
write_matrix_dir(
  mat = brain.data,
  dir = 'brain_counts')

# Matrix on disk; we can load it.
brain.mat <- open_matrix_dir(dir = "brain_counts")
# this converts Ensembl gene id to gene symbols as per mouse species. 
brain.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = brain.mat, species = "mouse")

# Creating Seurat Object.
brain <- CreateSeuratObject(counts = brain.mat)

#example analysis. We are working with a lot of data so this particular step is very memory intensive.
# VlnPlot(brain, features = c("Sox10", "Slc17a7", "Aif1"), ncol = 3, layer = "counts", alpha = 0.1)

# we can save this object for the future.
saveRDS(
  object = brain,
  file = "obj.Rds"
)
