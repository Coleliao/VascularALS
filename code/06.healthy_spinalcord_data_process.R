# snRNA-seq data from healthy human spinal cord samples (6 individuals)
# 2023 Neuron
# # Cole 2025/11/26


# 1. process the raw data
library(Seurat)

counts <- data.table::fread('./data/human_spinalcord_healthsample_Neuron2023/snRNA-seq/GSE190442_aggregated_counts_postqc.csv.gz',data.table = F)
rownames(counts) <- counts$V1; counts$V1 <- NULL
meta <- data.table::fread('./data/human_spinalcord_healthsample_Neuron2023/snRNA-seq/GSE190442_aggregated_metadata_postqc.csv.gz')
rownames(meta) <- meta$V1; meta$V1 <- NULL

data <- CreateSeuratObject(counts = counts, meta.data = meta)
table(data$sample)
table(data$top_level_annotation)


data[["RNA"]] <- split(data[["RNA"]], f = data$sample)
data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)
data <- IntegrateLayers(
  object = data, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
data <- JoinLayers(data,layers=c('counts','data'))
data <- FindNeighbors(data, dims = 1:30, reduction = "harmony")
data <- RunUMAP(data, reduction = "harmony", dims = 1:30, reduction.name = "umap") #
#n.neighbors = 100,min.dist = 0.3,spread=0.8)
data <- FindClusters(data, resolution = 0.5)
DimPlot(data, reduction = "umap", group.by = "sample")
DimPlot(data, reduction = "umap", label = T)
DimPlot(data, reduction = "umap", label = T,group.by='top_level_annotation')
saveRDS(data,'./data/human_spinalcord_healthsample_Neuron2023/snRNA-seq/human_spinalcord_healthsample_Neuron2023_snRNAseq_integrated.rds')

# 2 subset the data for vacularture
Idents(data) <- 'top_level_annotation'
vas <- subset(data, idents = c('Endothelial', 'Meninges', 'Pericytes'))
vas[["RNA"]] <- split(vas[["RNA"]], f = vas$sample)

vas <- DoIntegration(vas,split='sample')
vas <- FindClusters(vas, resolution = 0.8)
DimPlot(vas, reduction = 'umap', label = T)
DimPlot(vas, reduction = 'umap', group.by = 'sample')
DimPlot(vas, reduction = 'umap', group.by = 'top_level_annotation', label = T)
DotPlot(vas,features = vas_markers)+RotatedAxis()
FeaturePlot(vas,features = c("LAMA1",'TIAM1','KCNMA1','FBLN5'),ncol = 2,order = T)
FeaturePlot(vas,features = c("SIX1",'SLC2A1','SNCAIP','MAPK4'),order = T)
FeaturePlot(vas,features = c("CEMIP",'NAV2','ABCA10','TRPM3'),ncol = 2,order = T)
FeaturePlot(vas,features = c('THBS2'),ncol = 1,order = T)

FeaturePlot(vas,features = endo_marker_sub,ncol = 3,order = T)

anno <- data.frame(
  cluster=c(7,2,5,12,8,9,4,10,0,3,11,1,6),
  celltype_class=c(rep('Endothelial cell',4), 'Mural cell','Mural cell',rep('Fibroblast',7)),
  celltype_type=c('Arterial',rep('Capillary',2),'Venous','Pericyte','Smooth muscle cell',
                  rep('M. Fibroblast',2),rep('P. Fibroblast type1',3),rep('P. Fibroblast type2',2)),
  celltype_sub=c('Arterial',rep('Capillary',2),'Venous','Pericyte','Smooth muscle cell',
                 rep('M. Fibroblast',2),rep('P. Fibroblast c1',3),'P. Fibroblast c2','P. Fibroblast c3')
)
vas$celltype_class <- anno$celltype_class[match(vas$seurat_clusters,anno$cluster)]
vas$celltype_type <- anno$celltype_type[match(vas$seurat_clusters,anno$cluster)]
vas$celltype_sub <- anno$celltype_sub[match(vas$seurat_clusters,anno$cluster)]

# fig 3C
col <- c(Arterial='#DC143C',
         Venous='#FF7F00',
         Capillary='#FB9A99',
         `M. Fibroblast`='#276419',
         `P. Fibroblast type1`='#B8E186',
         `P. Fibroblast type2`='#7FBC41',
         `Smooth muscle cell`='#B6B0FF',
         Pericyte='#5D69B1')
pdf('./vascular_celltype_health_0511.pdf',width=6,height=4)
DimPlot(health,group.by='celltype_type',label=T,label.size=4,cols = col)
dev.off()

saveRDS(vas,'human_spinalcord_healthsample_vascular_cells_1126.rds') # 1715cells


# 3 cell type correspondence between the vascular cells in D1 and D2
sc <- readRDS('./human_als_vascular_final_1111.rds') # D1
health <- vas # D2

sc$dataset <- 'D1'
health$dataset <- 'D2'
data <- merge(sc,health)
data <- DietSeurat(data,layers = c('counts','data'))
data <- JoinLayers(data,layer = c('counts','data'))
library(SingleCellExperiment)
library(MetaNeighbor)
sce <- as.SingleCellExperiment(data)
data <- FindVariableFeatures(data,selection.method = 'vst',nfeatures = 2000)
var_genes <- VariableFeatures(data)
aurocs <- MetaNeighborUS(var_genes = var_genes,
                         dat = sce,
                         study_id = sce$dataset,
                         cell_type = sce$celltype_type,
                         fast_version = TRUE)
# fig 3E
pdf('MetaNeighbor_heatmap.pdf',width=10,height=10)
plotHeatmap(aurocs,cex=1,margins = c(15, 15),cluster_rows = T,cluster_cols = T) # note, 是绝对标准化的AUC,就没必要scale了
dev.off()


