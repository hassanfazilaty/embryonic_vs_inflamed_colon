###### Subseting the DSS.combined mesenchymal datasets!


mes_DSS3d.combined <- subset(DSS.combined, ident = c("4", "1", "8", "16", "3"))
mes_DSS3d.combined <- subset(mes_DSS3d.combined, ident = c("0", "1", "2", "3"))
mes_DSS3d.combined <- subset(mes_DSS3d.combined, ident = c("0", "1", "2"))


DefaultAssay(mes_DSS3d.combined) <- "RNA"
mesenchyme.list <- SplitObject(mes_DSS3d.combined, split.by = "treatment")
for (i in 1:2){
  mesenchyme.list[[i]] <- NormalizeData(mesenchyme.list[[i]], verbose = FALSE)}
for (i in 1:2){
  mesenchyme.list[[i]] <- FindVariableFeatures(mesenchyme.list[[i]], selection.method = "vst", nfeatures = 3000)}

#De novo combining the previously subsetted datasets!
mes_DSS3d.combined.anchors <- FindIntegrationAnchors(object.list = mesenchyme.list, dims = 1:30)
mes_DSS3d.combined <- IntegrateData(anchorset = mes_DSS3d.combined.anchors, dims = 1:30)
DefaultAssay(mes_DSS3d.combined) <- "integrated"
mes_DSS3d.combined <- ScaleData(mes_DSS3d.combined, verbose = FALSE)
mes_DSS3d.combined <- RunPCA(mes_DSS3d.combined, npcs = 30, verbose = TRUE)


#Dimensionality reduction
mes_DSS3d.combined <- RunUMAP(mes_DSS3d.combined, dims = 1:10)
mes_DSS3d.combined <- FindNeighbors(mes_DSS3d.combined, dims = 1:10)
mes_DSS3d.combined <- FindClusters(mes_DSS3d.combined, resolution = 0.1)
DimPlot(mes_DSS3d.combined, label = TRUE)


DefaultAssay(mes_DSS3d.combined) <- "RNA"
FeaturePlot(mes_DSS3d.combined, features = c("Epcam", "Cdh1", "Acta2", "Pdgfra"), 
            reduction = "umap", pt.size = 0.5, order = TRUE, min.cutoff = 0)

mes_DSS3d.combined$treatment <- factor(x = mes_DSS3d.combined$treatment,
                                       levels = c("untreated", "DSS"))
save(mes_DSS3d.combined, file = "labeled_mes_DSS3d.combined.Rdata")
# Visualization
p1 <- DimPlot(mes_DSS3d.combined, reduction = "umap", group.by = "treatment")
p2 <- DimPlot(mes_DSS3d.combined, reduction = "umap", label = TRUE)

DimPlot(mes_DSS3d.combined, reduction ="umap", split.by = "treatment", label = TRUE, pt.size = 2) + NoLegend()

CombinePlots(plots = list(p1, p2), pt.size = 1)

DimPlot(mes_DSS.combined, reduction ="umap", label = FALSE) + NoLegend()
DimPlot(mes_DSS3d.combined, reduction ="umap", label = TRUE, pt.size = 2)

#violin plot
VlnPlot(mes_DSS3d.combined, features = c( "Bmp5", "Bmp3", "Sfrp1", "Grem1",
                                         "Bmp2", "Rspo3", "Uba52", "Adamdec1",
                                         "Snrpg"), split.by = "treatment", pt=FALSE)
VlnPlot(mes_DSS3d.combined, features = c('Cxcl1',
                                         'Cxcl10',
                                         'Hif1a',
                                         'Susd6',
                                         'Itga5'), split.by = "treatment", pt=FALSE)
VlnPlot(Col_hindgut_Mes_std, features = c('Cxcl1',
                                         'Cxcl10',
                                         'Hif1a',
                                         'Susd6',
                                         'Itga5'), split.by = "stage", pt=FALSE)


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



genes.to.label = c("Mmp3","Serpina3n", "Igfbp5", "Reg3b",
                   "Snrpg", "Adamdec1", "Lcn2", "Saa3",
                   "Bmp2", "Bmp5", "Col15a1", "Wnt4", "Gpx2",
                   "Rspo3", "Smoc2",  "Postn", "Gsn", "Ctgf", "Fbln1 ",
                   "Tnc")


p1 <- ggplot(avg.IBDsix.cells, aes(untreated, DSS)) + geom_point() + ggtitle("DSS")


p1 <- LabelPoints(plot = p1, points = top20000, repel = TRUE, xnudge = 0, ynudge = 0)

p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)


plot_grid(p1)




