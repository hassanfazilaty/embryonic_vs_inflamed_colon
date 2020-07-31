library(Seurat)
library(Matrix)
library(dplyr)
library(cowplot)
library(ggplot2)

#read in raw data

E14_5_dev_std <- Read10X(data.dir = "~/E14_5/outs/filtered_feature_bc_matrix")
E15_5_dev_std <- Read10X(data.dir = "~/E15_5/outs/filtered_feature_bc_matrix")
E18_5_dev_std <- Read10X(data.dir = "~/E18_5/outs/filtered_feature_bc_matrix")
adult_std <- Read10X(data.dir = "~/adult/outs/filtered_feature_bc_matrix")

E14_5_dev_std
E15_5_dev_std

# Set up E14.5 object
E14_5_dev_std <- CreateSeuratObject(counts = E14_5_dev_std, project = "hindgut development", min.cells = 3)
E14_5_dev_std$stage <- "E14.5"
E14_5_dev_std[["percent.mt"]] <- PercentageFeatureSet(E14_5_dev_std, pattern = "mt-")
E14_5_dev_std <- subset(E14_5_dev_std, subset = nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 15)
E14_5_dev_std <- NormalizeData(E14_5_dev_std, verbose = FALSE)
E14_5_dev_std <- FindVariableFeatures(E14_5_dev_std, selection.method = "vst", nfeatures = 3000)

# Set up E15.5 object
E15_5_dev_std <- CreateSeuratObject(counts = E15_5_dev_std, project = "hindgut development", min.cells = 3)
E15_5_dev_std$stage <- "E15.5"
E15_5_dev_std[["percent.mt"]] <- PercentageFeatureSet(E15_5_dev_std, pattern = "mt-")
E15_5_dev_std <- subset(E15_5_dev_std, subset = nFeature_RNA > 1250 & nFeature_RNA < 7000 & percent.mt < 12)
E15_5_dev_std <- NormalizeData(E15_5_dev_std, verbose = FALSE)
E15_5_dev_std <- FindVariableFeatures(E15_5_dev_std, selection.method = "vst", nfeatures = 3000)


# Set up E18.5 object
E18_5_dev_std <- CreateSeuratObject(counts = E18_5_dev_std, project = "hindgut development", min.cells = 3)
E18_5_dev_std$stage <- "E18.5"
E18_5_dev_std[["percent.mt"]] <- PercentageFeatureSet(E18_5_dev_std, pattern = "mt-")
E18_5_dev_std <- subset(E18_5_dev_std, subset = nFeature_RNA > 1000 & nFeature_RNA < 7200 & percent.mt < 20)
E18_5_dev_std <- NormalizeData(E18_5_dev_std, verbose = FALSE)
E18_5_dev_std <- FindVariableFeatures(E18_5_dev_std, selection.method = "vst", nfeatures = 3000)
# Set up adult object
adult_std <- CreateSeuratObject(counts = adult_std, project = "hindgut development", min.cells = 3)
adult_std$stage <- "adult"
adult_std[["percent.mt"]] <- PercentageFeatureSet(adult_std, pattern = "mt-")
adult_std <- subset(adult_std, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 25)
adult_std <- NormalizeData(adult_std, verbose = FALSE)
adult_std <- FindVariableFeatures(adult_std, selection.method = "vst", nfeatures = 3000)



#Combining the datasets!
####only embryonic
hindgut_dev_std.anchors <- FindIntegrationAnchors(object.list = list(E14_5_dev_std, E15_5_dev_std, E18_5_dev_std), dims = 1:30)
hindgut_dev_std <- IntegrateData(anchorset = hindgut_dev_std.anchors, dims = 1:30)
DefaultAssay(hindgut_dev_std) <- "integrated"

####embryonic plus adult
Col_hindgut_std.anchors <- FindIntegrationAnchors(object.list = list(E14_5_dev_std, E15_5_dev_std, E18_5_dev_std, adult_std), dims = 1:30)
Col_hindgut_std <- IntegrateData(anchorset = Col_hindgut_std.anchors, dims = 1:30)


DefaultAssay(hindgut_dev_std) <- "integrated"

# standard workflow for visualization and clustering
hindgut_dev_std <- ScaleData(hindgut_dev_std, verbose = FALSE)


hindgut_dev_std <- RunPCA(hindgut_dev_std, npcs = 30, verbose = FALSE)

# Dimeonsionality reduction and Clustering

hindgut_dev_std <- RunUMAP(hindgut_dev_std, reduction = "pca", dims = 1:30)
hindgut_dev_std <- FindNeighbors(hindgut_dev_std, dims = 1:30)
hindgut_dev_std <- FindClusters(hindgut_dev_std, resolution = 0.5)


# Visualization

DimPlot(hindgut_dev_std, reduction ="umap", label = TRUE, pt.size = 0.5) 


p1 <- DimPlot(hindgut_dev_std, reduction = "umap", group.by = "stage")
p2 <- DimPlot(hindgut_dev_std, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1, p2), pt.size = 0.7)


DimPlot(hindgut_dev_std, reduction ="umap", split.by = "stage", pt.size = 1, label = TRUE) + NoLegend()



####GenesorteR

devtools::install_github("mahmoudibrahim/genesorteR", force = TRUE) 
#if "seuratObject" is the Seurat object that contains your data, I think this should work:
library(genesorteR)

gs = sortGenes(hindgut_dev_std@assays$RNA@data, Idents(hindgut_dev_std))
head(gs$specScore) #specificity scores for each gene in each cluster
mm = getMarkers(gs, quant = 0.99)
pp = plotMarkerHeat(gs$inputMat, gs$inputClass, mm$markers, clusterGenes=TRUE, outs = TRUE)

pp$gene_class_info #gene clusters



#dotplot of marker genes
DefaultAssay(hindgut_dev_std) <- "RNA"

markers.to.plot <- c("Epcam", "Cdh1", "Krt8", "Cdx2", "Muc2", "Neurod1", "Chga",
                     "Col1a1", "Col14a1", "Pdgfra", "Acta2", "Myh11",
                     "Stmn3", "Sox10", "Ascl1", "Wt1", "Cd93", "Pecam1", "Alas2", "Ctla2a")

DotPlot(hindgut_dev_std, features = rev(markers.to.plot), 
        cols = c("gray", "blue"), dot.scale = 15) + RotatedAxis()



#feature plot
FeaturePlot(hindgut_dev_std, features = c('Epcam', 'Col1a1'), 
            pt.size = 1, blend = TRUE, order = TRUE, blend.threshold = 0.3)


FeaturePlot(hindgut_dev_std, features = c("Cnn1", "Lmod1", "Sox10", "Ascl1", "Stmn3", "Wt1"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)



