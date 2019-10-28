library(Seurat)
library(monocle)
library(dplyr)
setwd("./Figure 1 and Extended Data Figure 1/")
data <- read.csv("./ALL_exprMatrix.tsv", sep = "\t", header = TRUE, row.names = 1)
metadata <- read.csv("./ALL_meta.tsv", sep = "\t", header = TRUE, row.names = 1)
gse <- CreateSeuratObject(counts = data, project = "gse", min.cells = 3, min.features = 200)
gse[['nGene']]=metadata$nGene
gse[['nUMI']]=metadata$nUMI
gse[['Cluster']]=metadata$Cluster
gse[['Time.Point']] = metadata$Time.Point
# gse <- subset(gse, subset = nFeature_RNA < 8000 & nFeature_RNA > 1500)

Endocardial <- subset(gse, subset = Cluster == 'Endocardial')
MP <- subset(gse, subset = Cluster == 'Multipot_progenitors')
Myo <- subset(gse, subset = Cluster == 'Myocardium')

gse <- NormalizeData(gse, normalization.method = "LogNormalize", scale.factor = 10000)
Endocardial <- NormalizeData(Endocardial, normalization.method = "LogNormalize", scale.factor = 10000)
MP <- NormalizeData(MP, normalization.method = "LogNormalize", scale.factor = 10000)
Myo <- NormalizeData(Myo, normalization.method = "LogNormalize", scale.factor = 10000)


all.genes <- rownames(gse)
gse <- ScaleData(gse, features = all.genes, vars.to.regress = "nUMI")
Endocardial <- ScaleData(Endocardial, features = all.genes, vars.to.regress = "nUMI")
MP <- ScaleData(MP, features = all.genes, vars.to.regress = "nUMI")
Myo <- ScaleData(Myo, features = all.genes, vars.to.regress = "nUMI")



gse <- FindVariableFeatures(gse, selection.method = "vst", nfeatures = 2000, mean.function = ExpMean, dispersion.function = LogVMR)
Endocardial <- FindVariableFeatures(Endocardial, selection.method = "vst", nfeatures = 2000, mean.function = ExpMean, dispersion.function = LogVMR)
MP <- FindVariableFeatures(MP, selection.method = "vst", nfeatures = 2000, mean.function = ExpMean, dispersion.function = LogVMR)
Myo <- FindVariableFeatures(Myo, selection.method = "vst", nfeatures = 2000, mean.function = ExpMean, dispersion.function = LogVMR)

gse <- ScaleData(gse, features = all.genes)
Endocardial <- ScaleData(Endocardial, features = all.genes)
MP <- ScaleData(MP, features = all.genes)
Myo <- ScaleData(Myo, features = all.genes)


gse <- RunPCA(gse, features = VariableFeatures(object = gse))
Endocardial <- RunPCA(Endocardial, features = VariableFeatures(object = Endocardial))
MP <- RunPCA(MP, features = VariableFeatures(object = MP))
Myo <- RunPCA(Myo, features = VariableFeatures(object = Myo))


DimHeatmap(gse, dims = 1:10)
DimHeatmap(Endocardial, dims = 11:20)
DimHeatmap(MP, dims = 1:10)
DimHeatmap(Myo, dims = 1:10)


ElbowPlot(gse,ndims = 50)
# gse <- FindNeighbors(gse, dims = 1:30)
# gse <- FindClusters(gse, resolution = 0.4)
pdf("./Figure 1 and Extended Data Figure 1/output.pdf")


gse <- RunUMAP(gse, dims = 1:10)
Endocardial <- RunUMAP(Endocardial, dims = 1:10)
MP <- RunUMAP(MP, dims = 1:10)
Myo <- RunUMAP(Myo, dims = 1:10)


DimPlot(gse, reduction = "umap", label = TRUE, group.by = "Cluster")

DimPlot(gse, reduction = "umap", label = TRUE, group.by = "Time.Point")
DimPlot(Endocardial, reduction = "umap", label = TRUE, group.by = "Time.Point")
DimPlot(MP, reduction = "umap", label = TRUE, group.by = "Time.Point")
DimPlot(Myo, reduction = "umap", label = TRUE, group.by = "Time.Point")


gse <- RunTSNE(gse, dims = 1:10)
DimPlot(gse, label = TRUE)
DimPlot(gse, label = TRUE, reduction = "tsne", group.by = "Cluster")
DimPlot(gse, label = TRUE, reduction = "tsne", group.by = "Time.Point")
FeaturePlot(gse, features = c("Flt4|Flt4", "Upp1|Upp1", "Rrad|Rrad", "Ank3|Ank3", "Prkaa2|Prkaa2"), reduction = "umap")
FeaturePlot(gse, features = c("Flt4|Flt4", "Upp1|Upp1", "Rrad|Rrad", "Ank3|Ank3", "Prkaa2|Prkaa2"), reduction = "tsne")
dev.off()


DimPlot(MP, reduction = "umap", label = TRUE, group.by = "Time.Point")
FindMarkers(gse, ident.1 = "E7.75", group.by = "Time.Point")



MP <- FindNeighbors(MP, dims = 1:30)
MP <- FindClusters(MP, resolution = 0.4)
MP <- RunTSNE(MP, dims = 1:10)
DimPlot(MP, label = TRUE, reduction = "umap")
DimPlot(MP, label = TRUE, reduction = "tsne")
MP.markers <- FindAllMarkers(MP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- MP.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MP, features = top10$gene) + NoLegend()


Endocardial <- FindNeighbors(Endocardial, dims = 1:30)
Endocardial <- FindClusters(Endocardial, resolution = 0.4)
Endocardial <- RunTSNE(Endocardial, dims = 1:10)
DimPlot(Endocardial, label = TRUE, reduction = "umap")
DimPlot(Endocardial, label = TRUE, reduction = "tsne")
Endocardial.markers <- FindAllMarkers(Endocardial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- Endocardial.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Endocardial, features = top10$gene) + NoLegend()

Myo <- FindNeighbors(Myo, dims = 1:30)
Myo <- FindClusters(Myo, resolution = 0.4)
Myo <- RunTSNE(Myo, dims = 1:10)
DimPlot(Myo, label = TRUE, reduction = "umap")
DimPlot(Myo, label = TRUE, reduction = "tsne")
Myo.markers <- FindAllMarkers(Myo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- Myo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Myo, features = top10$gene) + NoLegend()
FeaturePlot(Myo, features = c("Itga6|Itga6", "Acta1|Acta1"), reduction = "umap")










gse_cds = importCDS(gse, import_all = T)
Endocardial_cds = importCDS(Endocardial)
Myo_cds = importCDS(Myo)
MP_cds = importCDS(MP)
