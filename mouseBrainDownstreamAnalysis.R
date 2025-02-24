library(Seurat)
library(BPCells)
library(ggplot2)
# needs to be set for large dataset analysis
options(future.globals.maxSize = 6 * 1024^3) # setting max to 4 GB.

# We already have an on-disk representation of our dataset.

# Now, instead of loading our dataset into memory (1.3M samples), we simply create a connection to the data on-disk.
obj <- readRDS("obj.rds")

# obj has 27k features across 1.3 million samples within 1 assay.

# subsetting 60k cells out of the 1.3M cells based on gene co-variance matrix.
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- SketchData(
  object = obj,
  ncells = 40000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switching to analyzing the sketched dataset (in-memory) {this is actually already in place}.
# can switch to full dataset with DefaultAssay(obj) <- "RNA".
DefaultAssay(obj) <- "sketch"

# performing clustering on the sketched dataset.
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 2)

# running UMAP..
obj <- RunUMAP(obj, dims = 1:50, return.model = T)
DimPlot(obj, label = T, label.size = 3, reduction = "umap") + NoLegend()

# from literature, we know some genes are overexpressed in some cell types.
# endothelial cells (Igfbp7), Excitatory (Neurod6) and interneruon precursor (Dlx2) neurons, 
# Intermediate Progenitors (Eomes), Radial Glia (Vim), Cajal-Retzius cells (Reln),
# Oligodendroytes (Olig1), Macrophages (c1qa).
# feature plot for these genes & checking if they match with any particular cluster(s).
FeaturePlot(
  object = obj,
  features = c(
    "Igfbp7", "Neurod6", "Dlx2", "Gad2",
    "Eomes", "Vim", "Reln", "Olig1", "C1qa"
  ),
  ncol = 3
)

# extending results to full dataset.
obj <- ProjectData(
  object = obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(cluster_full = "seurat_clusters")
)
# now that we have projected the full dataset, switching back to analyzing all cells..
DefaultAssay(obj) <- "RNA"

DimPlot(obj, label = T, label.size = 3, reduction = "full.umap", 
        group.by = "cluster_full", alpha = 0.1) + NoLegend()

# iterative sub-clustering: zooming into Dlx2 gene producing cells.
DefaultAssay(obj) <- "sketch"
VlnPlot(obj, "Dlx2")

# from the violin plot, clusters 5, 15, 17, 21, 30, 37 are associated with Dlx2.
# subset cells in these clusters. Note that the data remains on-disk after subsetting
obj.sub <- subset(obj, subset = cluster_full %in% c(5, 15, 17, 21, 30, 37))
# over memory limit with this..so only using sketched dataset.
# Ideally: DefaultAssay(obj.sub) <- "RNA"

# now converting the RNA assay (previously on-disk) into an in-memory representation (sparse Matrix).
# we only convert the data layer, and keep the counts on-disk.
obj.sub[["RNA"]]$data <- as(obj.sub[["RNA"]]$data, Class = "dgCMatrix")

# reclustering the cells..
obj.sub <- FindVariableFeatures(obj.sub)
obj.sub <- ScaleData(obj.sub)
obj.sub <- RunPCA(obj.sub)
obj.sub <- RunUMAP(obj.sub, dims = 1:30)
obj.sub <- FindNeighbors(obj.sub, dims = 1:30)
obj.sub <- FindClusters(obj.sub)

DimPlot(obj.sub, label = T, label.size = 3) + NoLegend()


# subclusters within interneuron precursors.

# sub biomarkers: medial ganglionic eminence (Lhx6), caudal ganglionic eminence (Nr2f2).
# Sst (Sst) and Pvalb (Mef2c)-committed interneurons,
# CGE-derived (Meis2)-expressing progenitor population.
# (Gad2), (id2), (Dlx6a).
FeaturePlot(
  object = obj.sub,
  features = c(
    "Dlx2", "Gad2", "Lhx6", "Nr2f2", "Sst",
    "Mef2c", "Meis2", "Id2", "Dlx6os1"
  ),
  ncol = 3
)
