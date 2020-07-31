
#read in raw data

DSS <- Read10X(data.dir = "~/DSS/outs/filtered_feature_bc_matrix")

adult_std <- Read10X(data.dir = "~/adult/outs/filtered_feature_bc_matrix")



# Set up DSS3d object
DSS <- CreateSeuratObject(counts = DSS, project = "DSS_treatment", min.cells = 5)
DSS$treatment <- "DSS"
DSS[["percent.mt"]] <- PercentageFeatureSet(DSS, pattern = "mt-")
DSS <- subset(DSS, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)
DSS <- NormalizeData(DSS, verbose = FALSE)
DSS <- FindVariableFeatures(DSS, selection.method = "vst", nfeatures = 3000)

# Set up adult object
adult_std <- CreateSeuratObject(counts = adult_std, project = "DSS_treatment", min.cells = 5)
adult_std$treatment <- "untreated"
adult_std[["percent.mt"]] <- PercentageFeatureSet(adult_std, pattern = "mt-")
adult_std <- subset(adult_std, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 25)
adult_std <- NormalizeData(adult_std, verbose = FALSE)
adult_std <- FindVariableFeatures(adult_std, selection.method = "vst", nfeatures = 3000)




#Combining the datasets!

DSS.anchors <- FindIntegrationAnchors(object.list = list(adult_std, DSS), dims = 1:30)
DSS.combined <- IntegrateData(anchorset = DSS.anchors, dims = 1:30)


DefaultAssay(DSS.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
DSS.combined <- ScaleData(DSS.combined, verbose = FALSE)

## Regress out variables!!!!

## continue!

DSS.combined <- RunPCA(DSS.combined, npcs = 30, verbose = FALSE)

# Dimeonsionality reduction and Clustering
DSS.combined <- RunUMAP(DSS.combined, reduction = "pca", dims = 1:30)
DSS.combined <- FindNeighbors(DSS.combined, dims = 1:30)
DSS.combined <- FindClusters(DSS.combined, resolution = 0.55)

DSS.combined$treatment <- factor(x = DSS.combined$treatment,
                                 levels = c("untreated", "DSS"))
# Visualization
DimPlot(DSS.combined, reduction ="umap", label = TRUE, pt.size = 1)



p1 <- DimPlot(DSS.combined, reduction = "umap", group.by = "treatment")
p2 <- DimPlot(DSS.combined, reduction = "umap", label = TRUE)

DimPlot(DSS.combined, reduction ="umap", split.by = "treatment", label = TRUE, pt.size = 1) + NoLegend()

CombinePlots(plots = list(p1, p2), pt.size = 1)

DimPlot(DSS.combined, reduction ="umap", label = TRUE) + NoLegend()
DimPlot(DSS.combined, reduction ="umap", label = TRUE, pt.size = 1)






###### Subseting the DSS.combined epithelial datasets!


epi_DSS3d.combined <- subset(DSS.combined, ident = c("11", "5", "2", "6", "7", "10",
                                                     "12", "9", "24", "18", "14", "23"))
epi_DSS3d.combined <- subset(epi_DSS3d.combined, ident = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "14", "12", "13"))
epi_DSS3d.combined <- subset(epi_DSS3d.combined, ident = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))

DefaultAssay(epi_DSS3d.combined) <- "RNA"
epithelium.list <- SplitObject(epi_DSS3d.combined, split.by = "treatment")
for (i in 1:2){
  epithelium.list[[i]] <- NormalizeData(epithelium.list[[i]], verbose = FALSE)}
for (i in 1:2){
  epithelium.list[[i]] <- FindVariableFeatures(epithelium.list[[i]], selection.method = "vst", nfeatures = 3000)}

#De novo combining the previously subsetted datasets!
epi_DSS3d.combined.anchors <- FindIntegrationAnchors(object.list = epithelium.list , dims = 1:30)
epi_DSS3d.combined <- IntegrateData(anchorset = epi_DSS3d.combined.anchors, dims = 1:30)
DefaultAssay(epi_DSS3d.combined) <- "integrated"
epi_DSS3d.combined <- ScaleData(epi_DSS3d.combined, verbose = FALSE)
epi_DSS3d.combined <- RunPCA(epi_DSS3d.combined, npcs = 30, verbose = TRUE)


#Dimensionality reduction
epi_DSS3d.combined <- RunUMAP(epi_DSS3d.combined, dims = 1:20)
epi_DSS3d.combined <- FindNeighbors(epi_DSS3d.combined, dims = 1:20)
epi_DSS3d.combined <- FindClusters(epi_DSS3d.combined, resolution = 0.35)
DimPlot(epi_DSS3d.combined, label = TRUE)


DefaultAssay(epi_DSS3d.combined) <- "RNA"
FeaturePlot(epi_DSS3d.combined, features = c("Epcam", "Chga", "Acta2", "Pdgfra"), 
            reduction = "umap", pt.size = 0.5, order = TRUE, min.cutoff = 0)

epi_DSS3d.combined$treatment <- factor(x = epi_DSS3d.combined$treatment,
                                       levels = c("untreated", "DSS"))
save(epi_DSS3d.combined, file = "epi_DSS3d.combined.Rdata")
# Visualization
p1 <- DimPlot(epi_DSS3d.combined, reduction = "umap", group.by = "treatment")
p2 <- DimPlot(epi_DSS3d.combined, reduction = "umap", label = TRUE)
CombinePlots(plots = list(p1, p2), pt.size = 1)



DimPlot(epi_DSS3d.combined, reduction ="umap", split.by = "treatment", label = TRUE, pt.size = 1) + NoLegend()


DimPlot(epi_DSS.combined, reduction ="umap", label = FALSE) + NoLegend()
DimPlot(epi_DSS3d.combined, reduction ="umap", label = TRUE, pt.size = 1)


DefaultAssay(epi_DSS3d.combined) <- "RNA"

#Plot features in dimensional reduction

FeaturePlot(epi_DSS3d.combined, features = c("Reg3b", "Clca4b"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0, split.by = "treatment")

FeaturePlot(epi_Col_hindgut_std, features = c("Reg3b"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0, split.by = "treatment")
FeaturePlot(epi_Col_hindgut_std, features = c("Fabp6", "Maf"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)


VlnPlot(epi_DSS3d.combined, features = c("Reg3a", "Reg3b", "Apoa1", "Maf",
                                         "Fabp6", "Clca4b"), split.by = "treatment", pt=FALSE)

#matchscore2
library(matchSCore2)
library(nnet)
library(Matrix)

#epithelial cells
epi_hindgut_markers <-FindAllMarkers(epi_Col_hindgut_std, only.pos = TRUE)
epi_DSS_markers <- FindAllMarkers(epi_DSS3d.combined, only.pos = TRUE)



gene_cl.ref <- cut_markers(levels(epi_hindgut_markers$cluster),epi_hindgut_markers,ntop=200)

gene_cl.obs <- cut_markers(levels(epi_DSS_markers$cluster),epi_DSS_markers,ntop=200)

ms <- matchSCore2(gene_cl.ref = gene_cl.ref,gene_cl.obs = gene_cl.obs,ylab = "hindgut_epi",xlab = "epi_DSS3d")

ms$ggplot


#epithelial human cells
mes_hindgut_markers <-FindAllMarkers(hindgut_Mes_dev, only.pos = TRUE)
mes_DSS_markers <- FindAllMarkers(mes_DSS3d.combined, only.pos = TRUE)



gene_cl.ref <- cut_markers(levels(mes_hindgut_markers$cluster),mes_hindgut_markers,ntop=100)

gene_cl.obs <- cut_markers(levels(mes_DSS_markers$cluster),mes_DSS_markers,ntop=100)

ms <- matchSCore2(gene_cl.ref = gene_cl.ref,gene_cl.obs = gene_cl.obs,ylab = "hindgut_mes",xlab = "mes_IBD")

ms$ggplot

##Violin plot

Idents(epi_DSS.combined) <- "treatment"

VlnPlot(epi_DSS3d.combined, features = c('Myo15b',
                                         'S100a11',
                                         'Cdv3',
                                         'Ktn1',
                                         'Mif',
                                         'Snrpd2',
                                         'Polr2f',
                                         'Gpx4',
                                         'Uba52'), split.by = "treatment", pt=FALSE)

Idents(epi_Col_hindgut_std) <- "stage"
VlnPlot(epi_Col_hindgut_std, features = c('Myo15b',
                                         'S100a11',
                                         'Cdv3',
                                         'Ktn1',
                                         'Mif',
                                         'Snrpd2',
                                         'Polr2f',
                                         'Gpx4',
                                         'Uba52'), split.by = "treatment", pt=FALSE)



###### Subseting the DSS.combined epithelial stem datasets!


IESC_DSS3d.combined <- subset(epi_DSS3d.combined, ident = c("0", "2"))

DefaultAssay(IESC_DSS3d.combined) <- "RNA"
epithelium.list <- SplitObject(IESC_DSS3d.combined, split.by = "treatment")
for (i in 1:3){
  epithelium.list[[i]] <- NormalizeData(epithelium.list[[i]], verbose = FALSE)}
for (i in 1:3){
  epithelium.list[[i]] <- FindVariableFeatures(epithelium.list[[i]], selection.method = "vst", nfeatures = 3000)}

#De novo combining the previously subsetted datasets!
IESC_DSS3d.combined.anchors <- FindIntegrationAnchors(object.list = epithelium.list , dims = 1:30)
IESC_DSS3d.combined <- IntegrateData(anchorset = IESC_DSS3d.combined.anchors, dims = 1:30)
DefaultAssay(IESC_DSS3d.combined) <- "integrated"
IESC_DSS3d.combined <- ScaleData(IESC_DSS3d.combined, verbose = FALSE)
IESC_DSS3d.combined <- RunPCA(IESC_DSS3d.combined, npcs = 10, verbose = TRUE)


#Dimensionality reduction
IESC_DSS3d.combined <- RunUMAP(IESC_DSS3d.combined, dims = 1:20)
IESC_DSS3d.combined <- FindNeighbors(IESC_DSS3d.combined, dims = 1:20)
IESC_DSS3d.combined <- FindClusters(IESC_DSS3d.combined, resolution = 0.2)
DimPlot(IESC_DSS3d.combined, label = TRUE)


DefaultAssay(IESC_DSS3d.combined) <- "RNA"


IESC_DSS3d.combined$treatment <- factor(x = IESC_DSS3d.combined$treatment,
                                        levels = c("untreated", "DSS"))
save(IESC_DSS3d.combined, file = "IESC_DSS3d.combined.Rdata")
# Visualization
p1 <- DimPlot(IESC_DSS3d.combined, reduction = "umap", group.by = "treatment")
p2 <- DimPlot(IESC_DSS3d.combined, reduction = "umap", label = TRUE)

DimPlot(IESC_DSS3d.combined, reduction ="umap", split.by = "treatment", label = TRUE, pt.size = 2) + NoLegend()

CombinePlots(plots = list(p1, p2), pt.size = 2)

DimPlot(IESC_DSS.combined, reduction ="umap", label = FALSE) + NoLegend()
DimPlot(IESC_DSS3d.combined, reduction ="umap", label = TRUE, pt.size = 2)


DefaultAssay(IESC_DSS3d.combined) <- "RNA"

#Plot features in dimensional reduction

FeaturePlot(IESC_DSS3d.combined, features = c("Mki67", "Lgr5", "Myo10", "Sox4", "Lars2", "Uba52"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)

##DEGs
#Identify differential expressed genes across conditions

DefaultAssay(mes_DSS3d.combined) <- "RNA"

mes_DSS3d.combined <- FindVariableFeatures(mes_DSS3d.combined, selection.method = "vst", nfeatures = 3000)

top20000 <- head(VariableFeatures(mes_DSS3d.combined), 20000)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Col_healthy_non_IBD)
plot2 <- LabelPoints(plot = plot1, points = top500, repel = TRUE)


CombinePlots(plots = list(plot1, plot2), ncol =1)


theme_set(theme_cowplot())

Idents(mes_DSS3d.combined) <- "treatment"


avg.IBDsix.cells <- log1p(AverageExpression(mes_DSS3d.combined, verbose = FALSE)$RNA)
avg.IBDsix.cells$gene <- rownames(avg.IBDsix.cells)


genes.to.label = c("Reg3b", "Reg3g", "Mki67", "Lrg1", "Apoa1",
                   "Lars2", "Lgr5", "Muc2", "Aqp8", "Mptx1",
                  "Prpf4b", "Sval1", "Myo15b",
                   "Wfdc2", "Fabp6")


p1 <- ggplot(avg.IBDsix.cells, aes(untreated, DSS)) + geom_point() + ggtitle("DSS")


p1 <- LabelPoints(plot = p1, points = top20000, repel = TRUE, xnudge = 0, ynudge = 0)

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


plot_grid(p1)



