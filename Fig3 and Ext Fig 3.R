

#mesenchymal subset
#Mesenchymal embryonic only
hindgut_Mes_std <- subset(hindgut_dev_std, ident = c("0", "14", "3", "2", "1", "10", "13", "16"))
DefaultAssay(hindgut_Mes_std) <- "RNA"
mes.list <- SplitObject(hindgut_Mes_std, split.by = "stage")
for (i in 1:3){
  mes.list[[i]] <- NormalizeData(mes.list[[i]], verbose = FALSE)}
for (i in 1:3){
  mes.list[[i]] <- FindVariableFeatures(mes.list[[i]], selection.method = "vst", nfeatures = 3000)}

#De novo combining the previously subsetted datasets!
hindgut_Mes_std.anchors <- FindIntegrationAnchors(object.list = mes.list , dims = 1:30)
hindgut_Mes_dev <- IntegrateData(anchorset = hindgut_Mes_std.anchors, dims = 1:30)
DefaultAssay(hindgut_Mes_dev) <- "integrated"
hindgut_Mes_dev <- ScaleData(hindgut_Mes_dev, verbose = FALSE)
hindgut_Mes_dev <- RunPCA(hindgut_Mes_dev, verbose = TRUE)

#Dimensionality reduction
hindgut_Mes_dev <- RunUMAP(hindgut_Mes_dev, dims = 1:10, verbose = TRUE)
hindgut_Mes_dev <- FindNeighbors(hindgut_Mes_dev, dims = 1:10, verbose = TRUE)
hindgut_Mes_dev <- FindClusters(hindgut_Mes_dev, verbose = TRUE, resolution = 0.4)


#visualization
DimPlot(hindgut_Mes_dev, label = TRUE)

p1.2 <- DimPlot(hindgut_Mes_dev, reduction = "umap", group.by = "stage")
p2.2 <- DimPlot(hindgut_Mes_dev, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1.2, p2.2), pt.size = 0.5)

DimPlot(hindgut_Mes_dev, reduction ="umap", split.by = "stage", label = TRUE, pt.size = 1) + NoLegend()

DimPlot(hindgut_Mes_dev, reduction ="umap", label = TRUE, pt.size = 1) 



#dotplot

markers.to.plot <- c("Acta2", "Myh11", "Smoc2", "Wnt4", "Slit2", "Dkk2", "Kit", "Mki67", "Cdk1", "Dcn", "Col14a1",  
                     "Dek", "Hmmr", "Aurkb", "Adamdec1", "Snai2", "Lum", 
                     "Pitx1", "Wnt5a", "Pdgfra",  "Bmp2",
                     "Bmp5")
DotPlot(hindgut_Mes_dev, features = rev(markers.to.plot), 
        cols = c("gray", "blue"), dot.scale = 15) + RotatedAxis()
#Cell cycle

# segregate this list into markers of G2/M phase and markers of S phase

CC_hindgut_Mes_dev <- CellCycleScoring(hindgut_Mes_dev, s.features, g2m.features, set.ident = FALSE)

s.features <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Hells"
                ,"Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45",
                "Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")

g2m.features <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67",
                  "Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
                  "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2",
                  "Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe",
                  "Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")


DimPlot(CC_hindgut_Mes_dev, label = FALSE, pt.size = 1, group.by = "Phase", cols = c("#9CC1E1", "#7A0000", "#BF9000"))


#feature plot
FeaturePlot(CC_hindgut_Mes_dev, features = c("Pdgfra", "Bmp2"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)

##################
#Mesenchymal embryonic and adult combined
plot <- FeaturePlot(Col_hindgut_Mes_std, features = c("Epcam"),
                    reduction = "umap", pt.size = 0.5, order = TRUE, min.cutoff = 0, cols = c('black','green'))



select.cells <- CellSelector(plot = plot)
Idents(Col_hindgut_Mes_std, cells = select.cells) <- "Mesenchyme"
DimPlot(Col_hindgut_Mes_std, reduction = "umap", label = TRUE)

#subset dataset
Col_hindgut_Mes_std <- subset(Col_hindgut_Mes_std, idents = c("Mesenchyme"))
DefaultAssay(Col_hindgut_Mes_std) <- "RNA"

Col_hindgut_Mes_std <- subset(Col_hindgut_std, ident = c("0", "2", "3", "5", "11", "9", "15", "17"))
DefaultAssay(Col_hindgut_Mes_std) <- "RNA"

#Combining the datasets!
Mes.list <- SplitObject(Col_hindgut_Mes_std, split.by = "stage")
for (i in 1:4){
  Mes.list[[i]] <- NormalizeData(Mes.list[[i]], verbose = FALSE)}
for (i in 1:4){
  mMs.list[[i]] <- FindVariableFeatures(Mes.list[[i]], selection.method = "vst", nfeatures = 3000)}

#De novo combining the previously subsetted datasets!
Mes.anchors <- FindIntegrationAnchors(object.list = Mes.list, dims = 1:30)
Col_hindgut_Mes_std <- IntegrateData(anchorset = Mes.anchors, dims = 1:30)



DefaultAssay(Col_hindgut_Mes_std) <- "integrated"



Col_hindgut_Mes_std <- ScaleData(Col_hindgut_Mes_std, verbose = FALSE)
Col_hindgut_Mes_std <- RunPCA(Col_hindgut_Mes_std, verbose = TRUE)

#Dimensionality reduction
Col_hindgut_Mes_std <- RunUMAP(Col_hindgut_Mes_std, dims = 1:30, verbose = TRUE)
Col_hindgut_Mes_std <- FindNeighbors(Col_hindgut_Mes_std, dims = 1:30, verbose = TRUE)
Col_hindgut_Mes_std <- FindClusters(Col_hindgut_Mes_std, verbose = TRUE, resolution = 0.6)


#visualization
DimPlot(Col_hindgut_Mes_std, label = TRUE)


DefaultAssay(Col_hindgut_Mes_std) <- "RNA"
FeaturePlot(Col_hindgut_Mes_std, features = c("Pdgfra"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)




p1.2 <- DimPlot(Col_hindgut_Mes_std, reduction = "umap", group.by = "stage")
p2.2 <- DimPlot(Col_hindgut_Mes_std, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1.2, p2.2), pt.size = 0.5)

DimPlot(Col_hindgut_Mes_std, reduction ="umap", split.by = "stage", label = TRUE, pt.size = 1) + NoLegend()

DimPlot(Col_hindgut_Mes_std, reduction ="umap", label = TRUE, pt.size = 1) 


DefaultAssay(Col_hindgut_Mes_std) <- "RNA"

Col_hindgut_Mes_std <- NormalizeData(Col_hindgut_Mes_std, verbose = FALSE)





Col_hindgut_Mes_std$stage <- factor(x = Col_hindgut_Mes_std$stage,
                                    levels = c("E14.5", "E15.5", "E18.5", "adult"))

#Violin plot
VlnPlot(Col_hindgut_Mes_std, features = c('Ly6a',
                                          'Cd9',
                                          'Postn',
                                          'Pdgfra',
                                          'Wnt5a',
                                          'Acta2'), pt=FALSE)


#psupertime
#Convert Seurat object to SingleCellExperiment

Col_hindgut_Mes_std.sce <- as.SingleCellExperiment(
  Col_hindgut_Mes_std)

#Create vector with sequential cell labels to superimpose over pseudotime
treatment_labels <- 
  Col_hindgut_Mes_std@meta.data[["stage"]]

#create psupertime dataset

Col_hindgut_Mes_std_psuper_best = psupertime(Col_hindgut_Mes_std.sce, factor(treatment_labels), sel_genes='all', penalization='best')


Col_hindgut_Mes_std.psuper <- psupertime(
  Col_hindgut_Mes_std.sce, factor(treatment_labels), sel_genes = 'all')

plot_train_results(Col_hindgut_Mes_std_psuper_best)

save(Col_hindgut_Mes_std_psuper_best, file = "Col_hindgut_Mes_std_psuper_best.Rdata")


plot_labels_over_psupertime(Col_hindgut_Mes_std_psuper_best, label_name='Stages',palette="Set1")
plot_labels_over_psupertime(Col_hindgut_Mes_std_psuper_best, label_name='Stages')


plot_identified_gene_coefficients(Col_hindgut_Mes_std_psuper_best)

plot_predictions_against_classes(Col_hindgut_Mes_std_psuper_best)


plot_identified_genes_over_psupertime(Col_hindgut_Mes_std_psuper_best, label_name='Stages',palette="Set1")


plot_specified_genes_over_psupertime(Col_hindgut_Mes_std_psuper_best,
                                     c("Tuba1b", "Dbp", "Rhoj", "Mptx1", "Uba52", "Egfl7", "Adamdec1", "Bmp4",
                                       "Igfbp5", "Ube2c", "Procr", "Postn", "Fbln1", "Snai2", "Wnt2b", "Wt5a", "Wnt4",
                                       "Wnt2", "Pdgfra"), label_name ='Stages',palette="Set1")

