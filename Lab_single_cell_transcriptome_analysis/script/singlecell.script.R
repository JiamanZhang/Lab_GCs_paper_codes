library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(cowplot)
library(limma)
library(DESeq2)
library(scater)
library(SingleR)
library(scRNAseq)
library(dplyr)
library(pheatmap)
library(SeuratData)


## filt low quality cells and merge all single cell data
dir=c('./D4_heart/','./D7_LV_heart/','./D7_RV_heart/','./D10_LV_heart/','./D10_RV_heart/','./D14_LV_heart/','./D14_RV_heart/','./SWF/outs/filtered_feature_bc_matrix/','./F1/outs/filtered_feature_bc_matrix/','./POF/outs/filtered_feature_bc_matrix/')

sample_list<-c('D4heart','D7LVheart','D7RVheart','D10LVheart','D10RVheart','D14LVheart','D14RVheart','SWF','F1','POF')
names(dir) = c('D4heart','D7LVheart','D7RVheart','D10LVheart','D10RVheart','D14LVheart','D14RVheart','SWF','F1','POF')
scRNAlist <- list()
for (i in 1:7){
    counts <- Read10X(data.dir = dir[i])
    scRNAlist[[i]] <- CreateSeuratObject(counts, project =sample_list[i], min.cells = 3, min.features = 200)
    head(scRNAlist[[i]]@meta.data, 2)
    scRNAlist[[i]][["ND2"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND2")
    scRNAlist[[i]][["ND1"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND1")
    scRNAlist[[i]][["COX1"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "COX1")
    scRNAlist[[i]][["COII"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "COII")
    scRNAlist[[i]][["ATP8"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ATP8")
    scRNAlist[[i]][["ATP6"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ATP6")
    scRNAlist[[i]][["COX3"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "COX3")
    scRNAlist[[i]][["ND3"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND3")
    scRNAlist[[i]][["ND4L"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND4L")
    scRNAlist[[i]][["ND4"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND4")
    scRNAlist[[i]][["ND5"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND5")
    scRNAlist[[i]][["CYTB"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "CYTB")
    scRNAlist[[i]][["ND6"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND6")
    scRNAlist[[i]][["percent.mt"]] <- apply(scRNAlist[[i]]@meta.data[,4:16], 1, sum)
    scRNAlist[[i]]@meta.data <- scRNAlist[[i]]@meta.data[,-c(4:16)]
    head(scRNAlist[[i]]@meta.data, 2)

    scRNAlist[[i]]<-subset(scRNAlist[[i]],subset = nFeature_RNA > 200 & percent.mt < 20)
}

for (i in 8:10){
    counts <- Read10X(data.dir = dir[i])
    scRNAlist[[i]] <- CreateSeuratObject(counts, project =sample_list[i], min.cells = 3, min.features = 200)
    head(scRNAlist[[i]]@meta.data, 2)
    scRNAlist[[i]][["ND2"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND2")
    scRNAlist[[i]][["ND1"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND1")
    scRNAlist[[i]][["COX1"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "COX1")
    scRNAlist[[i]][["COII"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "COII")
    scRNAlist[[i]][["ATP8"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ATP8")
    scRNAlist[[i]][["ATP6"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ATP6")
    scRNAlist[[i]][["COX3"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "COX3")
    scRNAlist[[i]][["ND3"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND3")
    scRNAlist[[i]][["ND4L"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND4L")
    scRNAlist[[i]][["ND4"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND4")
    scRNAlist[[i]][["ND5"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND5")
    scRNAlist[[i]][["CYTB"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "CYTB")
    scRNAlist[[i]][["ND6"]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = "ND6")
    scRNAlist[[i]][["percent.mt"]] <- apply(scRNAlist[[i]]@meta.data[,4:16], 1, sum)
    scRNAlist[[i]]@meta.data <- scRNAlist[[i]]@meta.data[,-c(4:16)]
    head(scRNAlist[[i]]@meta.data, 2)

    scRNAlist[[i]]<-subset(scRNAlist[[i]],subset = nFeature_RNA > 1200 & nFeature_RNA < 5000 & nCount_RNA > 3000 & percent.mt < 20)
}

scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]],scRNAlist[[7]],scRNAlist[[8]],scRNAlist[[9]],scRNAlist[[10]]),add.cell.ids = sample_list,project='pbmc3k')
pbmc<-scRNA2

## filt the contaminating genes in the blood
filter_genes_data<-read.table(file='filter_genes.chicken.xls',sep='',header=F)
filter_genes<-filter_genes_data[,1]
removed_filter<-rownames(pbmcnew@assays$RNA)[!(rownames(pbmcnew@assays$RNA)%in%filter_genes)]
pbmcnewer<-subset(pbmcnew,features = removed_filter)

## normalised the all merged data by Seurat
bm280k<-pbmcnew
bm280k.list <- SplitObject(bm280k, split.by = "orig.ident")
bm280k.list <- lapply(X = bm280k.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE,normalization.method = "LogNormalize",scale.factor = 10000)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = bm280k.list)
bm280k.list <- lapply(X = bm280k.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = bm280k.list, reduction = "rpca",dims = 1:20)
bm280k.integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
bm280k.integrated <- ScaleData(bm280k.integrated, verbose = FALSE)
bm280k.integrated <- RunPCA(bm280k.integrated, verbose = FALSE)
saveRDS(bm280k.integrated, file = "merge_times_all_heart_muscle_ollicle_processed.no.muscle.rds")

## cluster by Umap
pbmcnewer<-readRDS('merge_times_all_heart_muscle_ollicle_processed.no.muscle.rds')
table(pbmcnew@meta.data$orig.ident)
pbmcnewer <- FindNeighbors(pbmcnewer, dims = 1:20)
for (num in c(0.2)){
    pbmcnewer <- FindClusters(pbmcnewer, resolution = num)
    myresult<-table(pbmcnewer$orig.ident,pbmcnewer@meta.data$seurat_clusters)
    write.table(myresult,file=paste('mergedata_clusters.r',num,'.no.muscle.xls',sep=''),sep='\t',quote=F,row.names=T,col.names=T)
    p1 <- DimPlot(pbmcnewer, reduction = "umap", label = TRUE, repel = TRUE)
    oup<-paste('mergedata-UMAP-cluster-r',num,'.cluster.no.muscle.pdf',sep='')
    ggsave(oup, plot=p1, width = 8, height =8)

}

## find DEgenes
pbmcnewer<-readRDS('merge_times_all_heart_muscle_ollicle_processed.no.muscle.rds')
pbmcnewer <- FindNeighbors(pbmcnewer, dims = 1:20)
pbmcnewer <- FindClusters(pbmcnewer, resolution = 0.2)
cluster2.markers <- FindAllMarkers(pbmcnewer)
write.table(cluster2.markers,file=paste('resolution0.2.all.cluster.DEgenes.new.txt',sep=''),sep='\t',quote=F,row.names=T,col.names=T)

## makergenes violin and bubble plot
pbmcnewer<-readRDS('merge_times_all_heart_muscle_ollicle_processed.no.muscle.rds')
pbmcnewer <- FindNeighbors(pbmcnewer, dims = 1:20)
pbmcnewer <- FindClusters(pbmcnewer, resolution = 0.2)
# example
mymarkers_to_plot<-c('FSHR','CYP11A1','CHST8','TSPAN6','DSP','CCK','NOV','RLN3','EDN2','FGL2','RGS16','MSX1','HBZ','FABP5','LCP1','ALDH1A2')
p1<-VlnPlot(pbmcnewer,assay='RNA',ncol=1, features = mymarkers_to_plot, log = TRUE,pt.size=0,slot = "data")
ggsave(paste("example.violin.scale.pdf",sep=""), plot=p1, width = 6, height = 46)
p1<-DotPlot(pbmcnewer,scale=F,features = mymarkers_to_plot,group.by='seurat_clusters',assay='RNA',cols = c("lightgrey", "red")) +
    RotatedAxis()+scale_size(range=c(1,15))
ggsave(paste("example.bubble.new.pdf",sep=""), plot=p1, width = 20, height = 8)

## Umaps of makergenes expression
pbmcnewer<-readRDS('merge_times_all_heart_muscle_ollicle_processed.no.muscle.rds')
pbmcnewer <- FindNeighbors(pbmcnewer, dims = 1:20)
pbmcnewer <- FindClusters(pbmcnewer, resolution = 0.2)
# example
mymarkers<-c('POSTN','COL3A1','DCN','TNNC1','TNNT2','IRX4','PITX2','NKX2-5')
for (i in 1:length(mymarkers)){
    g1 <- FeaturePlot(object = pbmcnewer, features = mymarkers[i], cols = c("grey", "blue"), reduction = "umap",pt.size=0.1,min.cutoff=0)
    ggsave(paste(mymarkers[i],"_featurePlot.modify.no.muscle.new.pdf",sep=""), plot=g1, width = 8, height = 8)
}

## Get top50 vst genes for cluster
pbmcnewer<-readRDS('merge_times_all_heart_muscle_ollicle_processed.no.muscle.rds')
pbmcnewer <- FindNeighbors(pbmcnewer, dims = 1:20)
pbmcnewer <- FindClusters(pbmcnewer, resolution = 0.2)
pbmcnewer <- ScaleData(pbmcnewer, verbose = FALSE)
clusternew<-pbmcnewer
tpmmatrix<-as.matrix(clusternew[["RNA"]]@data)
colnames(tpmmatrix)<-clusternew@meta.data$seurat_clusters
top50 <- head(VariableFeatures(pbmcnewer), 50)
newmatrix<-tpmmatrix[rownames(tpmmatrix)%in%top50,]
oup<-'whole.top50.vst.genes.matrix.txt'
write.table(newmatrix,oup,quote=F,col.names=T,row.names=T,sep='\t')


