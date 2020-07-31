#read in raw data

Healthy <- read.table(("~/Healthy_GSE116222_Expression_matrix.txt"), header=TRUE,row.names=1)
Non_inflamed <- read.table(("~/Non_inflamed_GSE116222_Expression_matrix.txt"), header=TRUE,row.names=1)
Inflamed <- read.table(("~/Inflamed_GSE116222_Expression_matrix.txt"), header=TRUE,row.names=1)


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




#read in raw data

Healthy_Kinchen <- read.table(("~/GSE114374_Human_HC_expression_matrix.txt"), header=TRUE,row.names=1)
IBD_Kinchen <- read.table(("~/GSE114374_Human_UC_expression_matrix.txt"), header=TRUE,row.names=1)



# Set up Healthy_Kinchen object
Healthy_Kinchen <- CreateSeuratObject(counts = Healthy_Kinchen, project = "h_colon", min.cells = 3)
Healthy_Kinchen$status <- "Healthy"
Healthy_Kinchen[["percent.mt"]] <- PercentageFeatureSet(Healthy_Kinchen, pattern = "MT-")
VlnPlot(Healthy_Kinchen, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Healthy_Kinchen <- subset(Healthy_Kinchen, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 3)
Healthy_Kinchen <- NormalizeData(Healthy_Kinchen, verbose = FALSE)
Healthy_Kinchen <- FindVariableFeatures(Healthy_Kinchen, selection.method = "vst", nfeatures = 3000)

# Set up IBD_Kinchen object
IBD_Kinchen <- CreateSeuratObject(counts = IBD_Kinchen, project = "h_atlas_f_intestine", min.cells = 3)
IBD_Kinchen$status <- "IBD"
IBD_Kinchen[["percent.mt"]] <- PercentageFeatureSet(IBD_Kinchen, pattern = "MT-")
VlnPlot(IBD_Kinchen, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
IBD_Kinchen <- subset(IBD_Kinchen, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 3)
IBD_Kinchen <- NormalizeData(IBD_Kinchen, verbose = FALSE)
IBD_Kinchen <- FindVariableFeatures(IBD_Kinchen, selection.method = "vst", nfeatures = 3000)


#Combining the Kinchen datasets!

Col_K.features <- SelectIntegrationFeatures(object.list = list(Healthy_Kinchen, IBD_Kinchen), nfeatures = 3000)
Col_K_list <- PrepSCTIntegration(object.list = list(Healthy_Kinchen, IBD_Kinchen), anchor.features = Col_K.features, 
                                 verbose = FALSE)
Col_K.anchors <- FindIntegrationAnchors(object.list = Col_K_list, normalization.method = "SCT", 
                                        anchor.features = Col_K.features, verbose = FALSE)
Col_K_healthy_IBD <- IntegrateData(anchorset = Col_K.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)

#Combining the Kinchen datasets STANDARD!

Col_K_healthy_IBD.anchors <- FindIntegrationAnchors(object.list = list(Healthy_Kinchen, IBD_Kinchen), dims = 1:30)
Col_K_healthy_IBD <- IntegrateData(anchorset = Col_K_healthy_IBD.anchors, dims = 1:30)



DefaultAssay(Col_K_healthy_IBD) <- "integrated"



# Run the standard workflow for visualization and clustering
Healthy_Kinchen <- ScaleData(Healthy_Kinchen, verbose = FALSE)

## continue!

Healthy_Kinchen <- RunPCA(Healthy_Kinchen, verbose = FALSE)

# Dimeonsionality reduction and Clustering

Healthy_Kinchen <- RunUMAP(Healthy_Kinchen, reduction = "pca", dims = 1:20)
Healthy_Kinchen <- FindNeighbors(Healthy_Kinchen, dims = 1:20)
Healthy_Kinchen <- FindClusters(Healthy_Kinchen, resolution = 0.3)


# Visualization
DimPlot(Healthy_Kinchen, reduction ="umap", label = TRUE, pt.size = 1) 

p1 <- DimPlot(Col_K_healthy_IBD, reduction = "umap", group.by = "status")
p2 <- DimPlot(Col_K_healthy_IBD, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1, p2), pt.size = 0.5)


DimPlot(Col_K_healthy_IBD, reduction ="umap", split.by = "status", label = TRUE, pt.size = 1) + NoLegend()



#### subset epithelial
################################################
epi_Col_healthy_non_IBD <- subset(Col_healthy_non_IBD,
                                  ident = c("0", "1", "2", "3", "4", "6", "7", "9"))
DefaultAssay(epi_Col_healthy_non_IBD) <- "RNA"
epi.list <- SplitObject(epi_Col_healthy_non_IBD, split.by = "status")
for (i in 1:3){
  epi.list[[i]] <- NormalizeData(epi.list[[i]], verbose = FALSE)}
for (i in 1:3){
  epi.list[[i]] <- FindVariableFeatures(epi.list[[i]], selection.method = "vst", nfeatures = 3000)}

#De novo combining the previously subsetted datasets!
epi_Col_healthy_non_IBD.anchors <- FindIntegrationAnchors(object.list = epi.list , dims = 1:30)
epi_Col_healthy_non_IBD <- IntegrateData(anchorset = epi_Col_healthy_non_IBD.anchors, dims = 1:30)
DefaultAssay(epi_Col_healthy_non_IBD) <- "integrated"
epi_Col_healthy_non_IBD <- ScaleData(epi_Col_healthy_non_IBD, verbose = FALSE)
epi_Col_healthy_non_IBD <- RunPCA(epi_Col_healthy_non_IBD, verbose = TRUE)

#Dimensionality reduction
epi_Col_healthy_non_IBD <- RunUMAP(epi_Col_healthy_non_IBD, dims = 1:10, verbose = TRUE)
epi_Col_healthy_non_IBD <- FindNeighbors(epi_Col_healthy_non_IBD, dims = 1:10, verbose = TRUE)
epi_Col_healthy_non_IBD <- FindClusters(epi_Col_healthy_non_IBD, verbose = TRUE, resolution = 0.2)


#visualization
DimPlot(epi_Col_healthy_non_IBD, label = TRUE)

epi_Col_healthy_non_IBD$status <- factor(x = epi_Col_healthy_non_IBD$status,
                                         levels = c("Healthy", "Non_inflamed", "IBD"))


p1.2 <- DimPlot(epi_Col_healthy_non_IBD, reduction = "umap", group.by = "status")
p2.2 <- DimPlot(epi_Col_healthy_non_IBD, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1.2, p2.2), pt.size = 0.5)

DimPlot(epi_Col_healthy_non_IBD, reduction ="umap", split.by = "status", label = TRUE, pt.size = 1) + NoLegend()

DimPlot(epi_Col_healthy_non_IBD, reduction ="umap", label = TRUE, pt.size = 1) 

# feature plot
FeaturePlot(epi_Col_healthy_non_IBD, features = c("ADH1C", "AQP8", "MDK", "MKI67"), 
            reduction = "umap", pt.size = 1, order = FALSE, min.cutoff = 0)
FeaturePlot(epi_Col_healthy_non_IBD, features = c("CA1", "LGR5", "MUC2", "SCGN"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)

FeaturePlot(epi_Col_hindgut_std, features = c("Car1"), 
            reduction = "umap", pt.size = 1, order = TRUE,split.by = "stage", min.cutoff = 0)

#violin plot
VlnPlot(epi_Col_healthy_non_IBD, features = c("PI3"), split.by = "status", pt=FALSE)
VlnPlot(epi_Col_hindgut_std, features = c("PI3"), split.by = "stage", pt=FALSE)


##matchscore2
##### compared to human
# epithelial


epi_DSS_markers <-FindAllMarkers(epi_DSS3d.combined, only.pos = TRUE)
human_epi_markers <- FindAllMarkers(epi_Col_healthy_non_IBD, only.pos = TRUE)


row.names(epi_DSS_markers) = toupper((row.names(epi_DSS_markers)))
epi_DSS_markers[,7] = toupper((epi_DSS_markers[,7]))



gene_cl.ref <- cut_markers(levels(human_epi_markers$cluster),human_epi_markers,ntop=200)
gene_cl.obs <- cut_markers(levels(epi_DSS_markers$cluster),epi_DSS_markers,ntop=200)
ms <- matchSCore2(gene_cl.ref = gene_cl.ref,gene_cl.obs = gene_cl.obs,ylab = "Human epi",xlab = "epi DSS")
ms$ggplot



##psupertime
#Convert Seurat object to SingleCellExperiment

epi_Col_healthy_non_IBD.sce <- as.SingleCellExperiment(epi_Col_healthy_non_IBD)

#Create vector with sequential cell labels to superimpose over pseudotime
treatment_labels <- 
  epi_Col_healthy_non_IBD@meta.data[["status"]]

#create psupertime dataset
epi_Col_healthy_non_IBD_psuper_best = psupertime(epi_Col_healthy_non_IBD.sce, factor(treatment_labels), sel_genes='all', penalization='best')

Col_healthy_non_IBD.psuper <- psupertime(Col_healthy_non_IBD.sce, factor(treatment_labels), sel_genes = 'all')

plot_train_results(epi_Col_healthy_non_IBD_psuper_best)

save(epi_Col_healthy_non_IBD_psuper_best, file = "epi_Col_healthy_non_IBD_psuper_best.Rdata")



plot_labels_over_psupertime(epi_Col_healthy_non_IBD_psuper_best, label_name='Status',palette="Set1")
plot_labels_over_psupertime(epi_Col_healthy_non_IBD_psuper_best, label_name='Status')


plot_identified_gene_coefficients(epi_Col_healthy_non_IBD_psuper_best)

plot_predictions_against_classes(epi_Col_healthy_non_IBD_psuper_best)


plot_identified_genes_over_psupertime(epi_Col_healthy_non_IBD_psuper_best, label_name='Status',palette="Set1")

plot_specified_genes_over_psupertime(epi_Col_healthy_non_IBD_psuper_best,
                                     c("MDK", "CCL20", "OLFM4", "PI3", "UBA7", "LCN2", "REG1A", "UBA6",
                                       "UBA1", "RBCK1", "AQP8", "LCN2", "XIST", "RNF19A", "REG1A", "CA1", "S100A11",
                                       "RNF19B", "RPLP1", "ECH1"), label_name ='Stages',palette="Set1")



#Identify differential expressed genes across conditions


epi_Col_healthy_non_IBD <- FindVariableFeatures(epi_Col_healthy_non_IBD, selection.method = "vst", nfeatures = 3000)

top20000 <- head(VariableFeatures(epi_Col_healthy_non_IBD), 20000)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Col_healthy_non_IBD)
plot2 <- LabelPoints(plot = plot1, points = top500, repel = TRUE)

CombinePlots(plots = list(plot1, plot2), ncol =1)


theme_set(theme_cowplot())

six.cells <- subset(Col_healthy_non_IBD, idents = "6")

Idents(epi_Col_healthy_non_IBD) <- "status"


avg.IBDsix.cells <- log1p(AverageExpression(epi_Col_healthy_non_IBD, verbose = FALSE)$RNA)
avg.IBDsix.cells$gene <- rownames(avg.IBDsix.cells)

genes.to.label = c("MDK", "CCL20", "OLFM4", "PI3", "LCN2", "REG1A",
                   "AQP8","CA1",
                   "CD74", "PDZK1IP1")


p1 <- ggplot(avg.IBDsix.cells, aes(Healthy, IBD)) + geom_point() + ggtitle("IBD")


p1 <- LabelPoints(plot = p1, points = top20000, repel = TRUE, xnudge = 0, ynudge = 0)

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


plot_grid(p1)


##violin plot
Idents(epi_Col_healthy_non_IBD) <- "status"
VlnPlot(epi_Col_healthy_non_IBD, features = c('MDK', 
                                              'OLFM4', 
                                              'CA1', 
                                              'AQP8',
                                              'MYO15B', 
                                              'S100A11', 
                                              'CDV3', 
                                              'FUT2', 
                                              'FKBP11', 
                                              'SPINK1', 
                                              'ITGA3', 
                                              'TFF1'), split.by = "status", pt=FALSE)






#########################
#subset fibroblasts
fb_Col_K_healthy_IBD <- subset(Col_K_healthy_IBD, ident = c("1", "2", "3", "7", "9"))
DefaultAssay(fb_Col_K_healthy_IBD) <- "RNA"
mes.list <- SplitObject(fb_Col_K_healthy_IBD, split.by = "status")
for (i in 1:2){
  mes.list[[i]] <- NormalizeData(mes.list[[i]], verbose = FALSE)}
for (i in 1:2){
  mes.list[[i]] <- FindVariableFeatures(mes.list[[i]], selection.method = "vst", nfeatures = 3000)}

#De novo combining the previously subsetted datasets!
fb_Col_K_healthy_IBD.anchors <- FindIntegrationAnchors(object.list = mes.list , dims = 1:30)
fb_Col_K_healthy_IBD <- IntegrateData(anchorset = fb_Col_K_healthy_IBD.anchors, dims = 1:30)
DefaultAssay(fb_Col_K_healthy_IBD) <- "integrated"
fb_Col_K_healthy_IBD <- ScaleData(fb_Col_K_healthy_IBD, verbose = FALSE)
fb_Col_K_healthy_IBD <- RunPCA(fb_Col_K_healthy_IBD, verbose = TRUE)

#Dimensionality reduction
fb_Col_K_healthy_IBD <- RunUMAP(fb_Col_K_healthy_IBD, dims = 1:20)
fb_Col_K_healthy_IBD <- FindNeighbors(fb_Col_K_healthy_IBD, dims = 1:20)
fb_Col_K_healthy_IBD <- FindClusters(fb_Col_K_healthy_IBD, resolution = 0.2)


fb_Col_K_healthy_IBD$treatment <- factor(x = fb_Col_K_healthy_IBD$status,
                                         levels = c("Healthy", "IBD"))


#visualization
DimPlot(fb_Col_K_healthy_IBD, label = TRUE)

p1.2 <- DimPlot(fb_Col_K_healthy_IBD, reduction = "umap", group.by = "status")
p2.2 <- DimPlot(fb_Col_K_healthy_IBD, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1.2, p2.2), pt.size = 0.5)

DimPlot(fb_Col_K_healthy_IBD, reduction ="umap", split.by = "status", label = TRUE, pt.size = 1) + NoLegend()

DimPlot(fb_Col_K_healthy_IBD, reduction ="umap", label = TRUE, pt.size = 1) 


#Plot features in dimensional reduction

FeaturePlot(fb_Col_K_healthy_IBD, features = c("PDGFRA", "ACTA2", "WNT5A", "WNT2B"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)

FeaturePlot(Col_hindgut_Mes_std, features = c("Grem1", "Sostdc1", "Wnt2b"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)
FeaturePlot(Col_hindgut_Mes_std, features = c("Col23a1", "Col15a1"), 
            reduction = "umap", pt.size = 1, order = TRUE, split.by = "status", min.cutoff = 0)

#violin plot
lnPlot(fb_Col_K_healthy_IBD, features = c("WNT2B", "SOSTDC1", "GREM1",
                                          "RSPO3"),  pt=FALSE)




#Identify differential expressed genes across conditions


fb_Col_K_healthy_IBD <- FindVariableFeatures(fb_Col_K_healthy_IBD, selection.method = "vst", nfeatures = 3000)

FeaturePlot(fb_Col_K_healthy_IBD, features = c("PDGFRA", "ACTA2", "WNT5A", "WNT2B"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)
top50000 <- head(VariableFeatures(fb_Col_K_healthy_IBD), 50000)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Col_K_healthy_IBD)
plot2 <- LabelPoints(plot = plot1, points = top500, repel = TRUE)

CombinePlots(plots = list(plot1, plot2), ncol =1)



theme_set(theme_cowplot())

Idents(fb_Col_K_healthy_IBD) <- "status"


fb_Col_K_healthy_IBD.cells <- log1p(AverageExpression(fb_Col_K_healthy_IBD, verbose = FALSE)$RNA)
fb_Col_K_healthy_IBD.cells$gene <- rownames(fb_Col_K_healthy_IBD.cells)

genes.to.label = c("MMP3", "BMP4", "IGFBP5", "PROCR", "FBLN1", "GSN",
                   "SMOC2", "EGFL7",
                   "ADAMDEC1", "COL7A1",
                   "COL15A1", "COL23A1", "BMP5")



p1 <- ggplot(fb_Col_K_healthy_IBD.cells, aes(Healthy, IBD)) + geom_point() + ggtitle("IBD_Kinchen")


p1 <- LabelPoints(plot = p1, points = top50000, repel = TRUE, xnudge = 0, ynudge = 0)

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


plot_grid(p1)


##violin plot

VlnPlot(fb_Col_K_healthy_IBD, features = c("PDGFRA", "BMP5", "WNT5A",
                                           "SOX6", "ADAMDEC1",
                                           "GSN", "MMP3", "IGFBP5","SMOC2"), split.by = "status", pt=FALSE)




Idents(fb_Col_K_healthy_IBD) <- "status"

VlnPlot(fb_Col_K_healthy_IBD, features = c('CXCL1',
                                              'HIF1A',
                                              'ITGA5',
                                              'NME1',
                                              'DUSP5',
                                              'SVIL',
                                              'GRAMD1A',
                                              'ARL5B'), split.by = "status", pt=FALSE)


