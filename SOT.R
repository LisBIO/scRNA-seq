library(SOT)
library(DT)
library(ggplot2)
library(tidyr)
library(SingleCellExperiment)
library(clusterProfiler)
library(ReactomePA)
library(destiny)
library(Matrix)
library(dplyr)
file_dir = system.file("extdata", "c1_osk.csv.gz", package = "SOT")
exprs = read.csv(file_dir, row.names = 1, check.names = F)
exprs[1:5, 1:5]
mysample = exprs[,sample(1:912, 900, replace = F)]
exprs = mysample
dim(exprs)
sce = BuildSCE(normcounts = exprs, trim = 5.1/ncol(exprs))
# sce$Day
day = c("mef", "osk_d0", "osk_d1", "osk_d2", "osk_d3", "osk_d4", "osk_d5", "osk_d7", "osk_d8", "ips", "esc")
sce$Day = factor(gsub("(mef|osk_d0|osk_d1|osk_d2|osk_d3|osk_d4|osk_d5|osk_d7|osk_d8|ips|esc).*", "\\1", colnames(sce)), levels = day)
metadata(sce)$'Day color' = c("mef" = "#BEAED4",
                              "osk_d0" = "#fdac86",
                              "osk_d1" = "#FFFF99",
                              "osk_d2" = "#386CB0",
                              "osk_d3" = "#F0027F",
                              "osk_d5" = "#BF5B17",
                              "osk_d7" = "#666666",
                              "osk_d8" = "#ffd56f",
                              "ips" = "#7FBC41",
                              "esc" = "#A6CEE3")                 




sce = FilterLow(sce, minexp = 10, mincells = 10, datatype = "normcounts")
assay(sce, "logcounts") = Matrix(log2(assay(sce, "normcounts")+1), sparse = TRUE)
pdf("./variance_trend.pdf")
sce = FindHVGs(sce, datatype = "logcounts", thr.bio = 0, thr.FDR = 0.1)
dev.off()
sum(rowData(sce)$genes.use)
sce = ADtest(sce, 
             datatype = "logcounts",
             condition = "Day", 
             genes.use = "find.hvgs",
             useLevels = c("osk_d1","osk_d2","osk_d3","osk_d5","osk_d7","osk_d8"), 
             ncore = 2,
             thr.padj = 1e-3,
             r.implement = T)
sum(rowData(sce)$genes.use)
sce = ap.cluster(sce, datatype = "logcounts", genes.use = "genes.use", ncore = 2)


metadata(sce)$p.cluster.size
plot_embed(sce, usedEmbed = "lv1.score", anno_col = "Day", color_scale = metadata(sce)$'Day color', fontsize = 5)

sce = reduce.cluster(sce)

metadata(sce)$p.group.size


plot_embed(sce, usedEmbed = "lv2.score", anno_col = "Day", cluster_row = FALSE, color_scale = metadata(sce)$`Day color`, fontsize=8)





# GO 分析
gene_anno = rowData(sce)
gene_anno = gene_anno[gene_anno$genes.use, ]
grlabs = split(gene_anno$symbol, gene_anno$lv2.labels)
gcSample = lapply(grlabs, function(gr) as.numeric(bitr(gr, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID))
xx.mus.go <- compareCluster(gcSample, OrgDb='org.Mm.eg.db', fun='enrichGO', pvalueCutoff = 0.1, qvalueCutoff = 0.1, ont = "BP", readable=T)
EnrichResGO = xx.mus.go@compareClusterResult

dotplot(xx.mus.go, title = paste0("Mouse Gene Ontology"))
datatable(EnrichResGO, filter = "top", options = list(pageLength = 10))




lv2.score = as.data.frame(reducedDim(sce, "lv2.score"))
group.std = lv2.score %>% 
  cbind(Day = sce$Day) %>%
  group_by(Day) %>% 
  summarise_all(funs(sd)) %>%
  gather(Group, Std, -Day)

ggplot(data=group.std, aes(x=Day, y=Std, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point() + facet_wrap(~Group) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10)) + 
  labs(title = "Divergence of gene groups in each day")




dm = DiffusionMap(reducedDim(sce, "lv2.score")[,-c(1)], distance ="cosine")
plot(eigenvalues(dm), ylim = c(0.8,1), pch = 20,main="Scree plot", xlab = "Diffusion component (DC)", ylab = "Eigenvalue")


top_dcs = dm[[paste0("DC", 1:7)]]
reducedDim(sce, "dcs") = top_dcs



sce = layout_2d(sce, k = 35, seed = 13, usedEmbed = "dcs", layout_method = "FR", colour_condition = "Day", color_scale = metadata(sce)$`Day color`)
metadata(sce)$p.dcs_FR.Day



sce = layout_2d(sce, usedEmbed = "dcs_FR", colour_embed = "lv2.score")
metadata(sce)$p.dcs_FR.lv2.score


genes = c("Dppa5a", "Nanog", "Sall4", "Klk10", "Fxyd5", "Cd34","Gbp2", "Gbp5", "Ifng")
sce = layout_2d(sce, datatype = "logcounts", usedEmbed = "dcs_FR", colour_gene = genes)
metadata(sce)$p.dcs_FR.gene

