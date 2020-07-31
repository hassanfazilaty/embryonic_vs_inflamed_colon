

#####subset 
#Epithelial


hindgut_epi_dev <- subset(hindgut_dev_std, ident = c("4", "5", "6", "9", "19"))
hindgut_epi_dev <- subset(hindgut_epi_dev, ident = c("0", "1", "2", "3", "4", "5", "6", "7", "9", "10", "12"))
hindgut_epi_dev <- subset(hindgut_epi_dev, ident = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))

DefaultAssay(hindgut_epi_dev) <- "RNA"
epithelium.list <- SplitObject(hindgut_epi_dev, split.by = "stage")
for (i in 1:3){
  epithelium.list[[i]] <- NormalizeData(epithelium.list[[i]], verbose = FALSE)}
for (i in 1:3){
  epithelium.list[[i]] <- FindVariableFeatures(epithelium.list[[i]], selection.method = "vst", nfeatures = 3000)}

#De novo combining the previously subsetted datasets!
hindgut_epi_dev.anchors <- FindIntegrationAnchors(object.list = epithelium.list , dims = 1:30)
hindgut_epi_dev <- IntegrateData(anchorset = hindgut_epi_dev.anchors, dims = 1:30)
DefaultAssay(hindgut_epi_dev) <- "integrated"
hindgut_epi_dev <- ScaleData(hindgut_epi_dev, verbose = FALSE)
hindgut_epi_dev <- RunPCA(hindgut_epi_dev, npcs = 30, verbose = TRUE)

#Dimensionality reduction
hindgut_epi_dev <- RunUMAP(hindgut_epi_dev, dims = 1:30, verbose = TRUE, reduction = 'pca')
hindgut_epi_dev <- FindNeighbors(hindgut_epi_dev, dims = 1:30, verbose = TRUE,reduction = 'pca')
hindgut_epi_dev <- FindClusters(hindgut_epi_dev, verbose = TRUE, resolution = 0.4, algorithm = 1)


#visualization
DimPlot(hindgut_epi_dev, label = TRUE)

p1.1 <- DimPlot(hindgut_epi_dev, reduction = "umap", group.by = "stage")
p2.1 <- DimPlot(hindgut_epi_dev, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1.1, p2.1), pt.size = 1)

DimPlot(hindgut_epi_dev, reduction ="umap", split.by = "stage", label = TRUE, pt.size = 1.5) + NoLegend()

DimPlot(hindgut_epi_dev, reduction ="umap", label = TRUE, pt.size = 1) 


#Cell cycle

# segregate this list into markers of G2/M phase and markers of S phase

CC_hindgut_epi_dev <- CellCycleScoring(hindgut_epi_dev, s.features, g2m.features, set.ident = FALSE)

s.features <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Hells"
                ,"Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45",
                "Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")

g2m.features <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67",
                  "Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
                  "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2",
                  "Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe",
                  "Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")


DimPlot(CC_hindgut_epi_dev, label = FALSE, pt.size = 1, group.by = "Phase", cols = c("#9CC1E1", "#7A0000", "#BF9000"))


#psupertime
library('psupertime')
library('SingleCellExperiment')

#Convert Seurat object to SingleCellExperiment

epi_Col_hindgut_std.sce <- as.SingleCellExperiment(
  epi_Col_hindgut_std)

#Create vector with sequential cell labels to superimpose over pseudotime
treatment_labels <- 
  epi_Col_hindgut_std@meta.data[["stage"]]

#create psupertime dataset

epi_Col_hindgut_std_psuper_best = psupertime(epi_Col_hindgut_std.sce, factor(treatment_labels), sel_genes='all', penalization='best')


plot_train_results(epi_Col_hindgut_std_psuper_best)



plot_labels_over_psupertime(epi_Col_hindgut_std_psuper_best, label_name='Stages',palette="Set1")
plot_labels_over_psupertime(epi_Col_hindgut_std_psuper_best, label_name='Stages')


plot_identified_gene_coefficients(epi_Col_hindgut_std_psuper_best)

plot_predictions_against_classes(epi_Col_hindgut_std_psuper_best)


plot_identified_genes_over_psupertime(epi_Col_hindgut_std.psuper, label_name='Stages',palette="Set1")

plot_specified_genes_over_psupertime(epi_Col_hindgut_std_psuper_best,
                                     c("Aqp8", "Muc2", "Uba52","Eno1", "A100a11", "Acin1", 
                                       "Ccl20", "Marcksl1", "Chga", "klk1", "Sox11", "Sox4", "Rpl36", "Hnrnpu",
                                       "Hnrnpl", "Lgr4", "Tacstd2", "Chgb", "Neurod1", "Lars2"), label_name ='Stages',palette="Set1")


plot_specified_genes_over_psupertime(epi_Col_hindgut_std_psuper_best,
                                     c("Tff3"), label_name ='Stages',palette="Set1")


#Plot features in dimensional reduction


FeaturePlot(hindgut_epi_dev, features = c("Mdk", "Apoa1", "Lgr5", "Ly6a", "Spdef", "Tacstd2"), 
            reduction = "umap", pt.size = 0.5, split.by = "stage", order = TRUE, min.cutoff = 0)

###### Subseting the embryonic and adult epithelial datasets!


plot <- FeaturePlot(Col_hindgut_std, features = c("Epcam"),
                    reduction = "umap", pt.size = 0.5, order = TRUE, min.cutoff = 0, cols = c('black','green'))

plot <- FeaturePlot(epi_Col_hindgut_std, features = c("Epcam", "Chga"),
                    reduction = "umap", pt.size = 0.5, order = TRUE, min.cutoff = 0, cols = c('black','green'))

select.cells <- CellSelector(plot = plot)
Idents(epi_Col_hindgut_std, cells = select.cells) <- "Epithelium"
DimPlot(epi_Col_hindgut_std, reduction = "umap", label = TRUE)

#subset dataset
epi_Col_hindgut_std <- subset(epi_Col_hindgut_std, idents = c("Epithelium"))
DefaultAssay(epi_Col_hindgut_std) <- "RNA"
epithelium.list <- SplitObject(epi_Col_hindgut_std, split.by = "stage")
for (i in 1:4){
  epithelium.list[[i]] <- NormalizeData(epithelium.list[[i]], verbose = FALSE)}
for (i in 1:4){
  epithelium.list[[i]] <- FindVariableFeatures(epithelium.list[[i]], selection.method = "vst", nfeatures = 3000)}
#De novo combining the previously subsetted datasets!
epi_Col_hindgut_std.anchors <- FindIntegrationAnchors(object.list = epithelium.list , dims = 1:30)
epi_Col_hindgut_std <- IntegrateData(anchorset = epi_Col_hindgut_std.anchors, dims = 1:30)
DefaultAssay(epi_Col_hindgut_std) <- "integrated"
epi_Col_hindgut_std <- ScaleData(epi_Col_hindgut_std, verbose = FALSE)
epi_Col_hindgut_std <- RunPCA(epi_Col_hindgut_std, verbose = TRUE)

#Dimensionality reduction
epi_Col_hindgut_std <- RunUMAP(epi_Col_hindgut_std, dims = 1:20, verbose = TRUE, reduction = 'pca')
epi_Col_hindgut_std <- FindNeighbors(epi_Col_hindgut_std, dims = 1:20, verbose = TRUE,reduction = 'pca')
epi_Col_hindgut_std <- FindClusters(epi_Col_hindgut_std, verbose = TRUE, resolution = 0., algorithm = 1)
DimPlot(epi_Col_hindgut_std, label = TRUE)

#feature plot
FeaturePlot(hindgut_epi_dev, features = c("Lgr5", "Muc2", "Aqp8", "Chga"), 
            reduction = "umap", pt.size = 1, order = TRUE, min.cutoff = 0)

#dotplot

markers.to.plot <- c("Epcam",  "Cdk1", "Mki67", "Dmbt1", "Apoa1",
                     "Lgr5",  "Sox4", "Tacstd2", "Sox11", "Alpi", "Axin2","Wfdc2",
                     "Muc2", "Reg4", "Aqp8", "Atp12a", "Ccl25", "Krt20", "Krt19", "Chga", "Scgn", "Fras1", 
                     "Spdef", "Atoh1", "Fabp2", "Klk1", "Muc3")
DotPlot(hindgut_epi_dev, features = rev(markers.to.plot), 
        cols = c("gray", "blue"), dot.scale = 15) + RotatedAxis()