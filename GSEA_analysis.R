###GSEA analysis of NA and Nras
setwd("./gsea2")
Ctr_NA<-read.csv("Ctr_NA.csv",header=T,row.names = 1)

###### Ctr_NA  GSEA analysis
########GSEA analysis
Ctr_NA<-Ctr_NA[apply(Ctr_NA,1,mean)>10,]  ##---Set the threshhold ---##
row.names(Ctr_NA)<-toupper(row.names(Ctr_NA))
source(file="../GSEA_Function scripts.R",local=T,verbose=T)

####Analysis Nras_Ctr
condition<-c("Ctr","NA")  # input the sample_groups
groups<-c(4,4) #input the replicate numbers for each group

dir<-getwd()
new_dir<-paste(dir,"NA_gsea",sep = "/")
dir.create(new_dir,mode="0777") 
setwd(new_dir)

###--Running the GSEA
##--Part1 creat the gsea.cls----#####
i=1
j=2
sp_n<-sum(groups[i],groups[j])
line1<-paste(sp_n, "2 1")
sg<-paste(condition[i],condition[j],sep = " ")
line2<-paste("#",sg,sep = " ")
gp1<-condition[i]
for(g1 in 2:groups[i]){gp1<-paste(gp1,condition[i],sep=" ")}
gp2<-condition[j]
for(g2 in 2:groups[j]){gp2<-paste(gp2,condition[j],sep=" ")}
line3<-paste(gp1,gp2,sep=" ")
gsea_cls<-c(line1,line2,line3)
write(gsea_cls,file="gsea.cls")  ##Creat .cls files
cls_na_file<-paste(getwd(),"gsea.cls",sep="/")

gsea_dir<-getwd()

setwd(gsea_dir)
gsea_out<-GSEA(
  # Input/Output Files :------------------------------------------------
  
  # Input gene expression Affy dataset file in RES or GCT format
  input.ds=Ctr_NA,
  
  # Input class vector (phenotype) file in CLS format
  input.cls=cls_na_file,
  # Gene set database in GMT format
  #gs.db = mm_pathway,
  #gs.db="E:/Bioinformatic_learning/Annotation/GSEA/c2.all.v5.2.symbols.gmt",
  gs.db="../Cell_cycle_mcm_apoptosis_genesets.gmt",
  # Directory where to store output and results (default: "")
  output.directory = gsea_dir,
  
  # Program parameters :-----------------------------------------------
  doc.string = "mouse",
  non.interactive.run = T,
  reshuffling.type = "sample.labels",
  nperm = 1000,
  weighted.score.type = 1,
  nom.p.val.threshold = -1,
  fwer.p.val.threshold = -1,
  fdr.q.val.threshold = 0.25,
  topgs = 10,
  adjust.FDR.q.val = F,
  gs.size.threshold.min = 15,
  gs.size.threshold.max = 500,
  reverse.sign = F,
  preproc.type = 0,
  random.seed = 3338,
  perm.type = 0,
  fraction = 1.0,
  replace = F,
  save.intermediate.results = F,
  OLD.GSEA = F,
  use.fast.enrichment.routine = T
)



