#read in raw data

Healthy <- read.table(("C:/My files/KB Lab/Data/GSE116222, human colon/Healthy_GSE116222_Expression_matrix.txt"), header=TRUE,row.names=1)
Non_inflamed <- read.table(("C:/My files/KB Lab/Data/GSE116222, human colon/Non_inflamed_GSE116222_Expression_matrix.txt"), header=TRUE,row.names=1)
Inflamed <- read.table(("C:/My files/KB Lab/Data/GSE116222, human colon/Inflamed_GSE116222_Expression_matrix.txt"), header=TRUE,row.names=1)


# Set up Healthy object
Healthy <- CreateSeuratObject(counts = Healthy, project = "h_colon", min.cells = 5)
Healthy$status <- "Healthy"
Healthy[["percent.mt"]] <- PercentageFeatureSet(Healthy, pattern = "MT-")
VlnPlot(Healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Healthy <- subset(Healthy, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
Healthy <- NormalizeData(Healthy, verbose = FALSE)
Healthy <- FindVariableFeatures(Healthy, selection.method = "vst", nfeatures = 3000)

# Set up Inflamed object
Inflamed <- CreateSeuratObject(counts = Inflamed, project = "h_atlas_f_intestine", min.cells = 5)
Inflamed$status <- "IBD"
Inflamed[["percent.mt"]] <- PercentageFeatureSet(Inflamed, pattern = "MT-")
VlnPlot(Inflamed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Inflamed <- subset(Inflamed, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 5)
Inflamed <- NormalizeData(Inflamed, verbose = FALSE)
Inflamed <- FindVariableFeatures(Inflamed, selection.method = "vst", nfeatures = 3000)

# Set up Non_inflamed object
Non_inflamed <- CreateSeuratObject(counts = Non_inflamed, project = "h_atlas_f_intestine", min.cells = 5)
Non_inflamed$status <- "Non_inflamed"
Non_inflamed[["percent.mt"]] <- PercentageFeatureSet(Non_inflamed, pattern = "MT-")
VlnPlot(Non_inflamed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Non_inflamed <- subset(Non_inflamed, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
Non_inflamed <- NormalizeData(Non_inflamed, verbose = FALSE)
Non_inflamed <- FindVariableFeatures(Non_inflamed, selection.method = "vst", nfeatures = 3000)






#Combining the GSE116222 datasets!
Col_healthy_non_IBD.anchors <- FindIntegrationAnchors(object.list = list(Healthy, Non_inflamed, Inflamed), dims = 1:30)
Col_healthy_non_IBD <- IntegrateData(anchorset = Col_healthy_non_IBD.anchors, dims = 1:30)





#########
Col.features <- SelectIntegrationFeatures(object.list = list(Healthy, Inflamed), nfeatures = 3000)
Col_list <- PrepSCTIntegration(object.list = list(Healthy, Inflamed), anchor.features = Col.features, 
                               verbose = FALSE)
Col.anchors <- FindIntegrationAnchors(object.list = Col_list, normalization.method = "SCT", 
                                      anchor.features = Col.features, verbose = FALSE)
Col_healthy_IBD <- IntegrateData(anchorset = Col.anchors, normalization.method = "SCT", 
                                 verbose = FALSE)
DefaultAssay(Col_healthy_non_IBD) <- "integrated"


# Run the standard workflow for visualization and clustering
Col_healthy_non_IBD <- ScaleData(Col_healthy_non_IBD, verbose = FALSE)

## continue!

Col_healthy_non_IBD <- RunPCA(Col_healthy_non_IBD, verbose = FALSE)

# Dimeonsionality reduction and Clustering

Col_healthy_non_IBD <- RunUMAP(Col_healthy_non_IBD, reduction = "pca", dims = 1:30)
Col_healthy_non_IBD <- FindNeighbors(Col_healthy_non_IBD, dims = 1:30)
Col_healthy_non_IBD <- FindClusters(Col_healthy_non_IBD, resolution = 0.3)


# Visualization
DimPlot(Col_healthy_non_IBD, reduction ="umap", label = TRUE, pt.size = 1) 

p1 <- DimPlot(Col_healthy_non_IBD, reduction = "umap", group.by = "status")
p2 <- DimPlot(Col_healthy_non_IBD, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1, p2), pt.size = 0.5)
plot_grid(p1, p2)


DimPlot(Col_healthy_non_IBD, reduction ="umap", split.by = "status", label = TRUE, pt.size = 1) + NoLegend()

Col_healthy_non_IBD$status <- factor(x = Col_healthy_non_IBD$status,
                                     levels = c("Healthy", "Non_inflamed", "IBD"))
