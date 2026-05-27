# process ALS single cell data -- spinal cord
# subset vascular cells (particularly for fibroblast)
# Related to Fig1 and FigS1
# Kuo 2025-04-30

library(Seurat)
library(hdf5r)
library(ggplot2)
source('~/functions/process.R')


#### Spinal cord data of human, Immunity 2025 ####
setwd('./01.SpinalCord_2025')
path = "./data/Spinaldata"
filelist = list.files(
  path = path,
  pattern = "h5")
datalist = list()
sampleInfo = sapply(strsplit(filelist,'_'),'[',3)
individualInfo = sapply(strsplit(filelist,'_'),'[',1)

pb <- txtProgressBar(min = 0, max = length(filelist), style = 3) # show progress
for (i in 1:length(filelist)) {
  tmp = Read10X_h5(paste0(path, '/', filelist[i]))
  tmp = CreateSeuratObject(counts = tmp, project = sampleInfo[i], min.cells = 3, min.features = 200)
  tmp$species = 'human'
  tmp$sample = sampleInfo[i] # pathological state
  tmp$individual = individualInfo[i] # individual
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  tmp <- subset(tmp, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &
                  nCount_RNA > 1000 & percent.mt < 15) 
  datalist[[i]] = tmp
  setTxtProgressBar(pb, i) 
}
close(pb)

# integrate the data using Seurat's integration workflow (all cells)
spinal = DoIntegration(datalist,method='log',nfeatures=3000,dims=1:50,resolution=0.5,gene_filter=T)
saveRDS(spinal, file = 'spinal_integrated_allcell_sc1_0430.rds')

DefaultAssay(spinal) = "RNA"
pdf('DimPlot_Spinal_allcell_sc1_0430.pdf',width=10,height=8)
DimPlot(spinal, reduction = "umap", group.by = "individual", label = TRUE, pt.size = 0.5)
DimPlot(spinal, label = TRUE, pt.size = 0.5)
FeaturePlot(spinal, features = c('COL3A1','BICC1','COL1A2','VIM'),order =T) # Fibroblast
FeaturePlot(spinal, features = c('PDGFRB','RGS5'),order =T) # Pericyte and Endothelial cell
dev.off()



# subset vascular cells
spinal <- readRDS('spinal_integrated_allcell_sc1_0430.rds')
vas <- subset(spinal, idents=c(16,13,11,22))
DefaultAssay(vas) = "RNA" # 4958
anno <- data.frame(suerat_clusters=c(16,13,11,22),
                   celltype=c('mural cell','Endothelial cell','Fibroblast','Fibroblast')) 
vas$celltype <- anno[match(vas$seurat_clusters, anno$suerat_clusters),'celltype']
pdf('DimPlot_Vascular_allcell_sc1_0501.pdf',width=10,height=8)
DimPlot(vas, reduction = "umap", group.by = 'celltype',label = TRUE, pt.size = 0.5)
dev.off()
# remove discrete cells in UMAP (regarded as low-quality) using DBSCAN
df <- Embeddings(vas, reduction = "umap") %>% as.data.frame()
library(dbscan)
kNNdistplot(df, k = 5) 
abline(h = 0.2, col = "red", lty = 2)
res <- dbscan(df, eps = 0.2, minPts = 5) # 
plot(x=df[,1], y=df[,2], col=res$cluster + 1L, pch=19, cex=0.5,
     xlab="UMAP_1", ylab="UMAP_2", main="DBSCAN Clustering") 
df$cluster <- res$cluster
df <- df[which(df$cluster %in% c(1,2)),]
vas <- vas[,rownames(df)] 

saveRDS(vas, file = 'spinal_integrated_VascularCell_sc1_0501.rds') # precise vascular cells for downstream analysis


# harmony integration in vacular cells
library(harmony)
DefaultAssay(vas) <- 'RNA'
vas_ha <- DietSeurat(vas,assays = 'RNA')
vas_ha <- vas_ha %>% SCTransform() %>% RunPCA() 
vas_ha <- RunHarmony(vas_ha, group.by.vars = 'donor_id', plot_convergence = F, max.iter.harmony = 15)
vas_ha <- vas_ha %>%  RunUMAP(reduction='harmony',dims = 1:50) %>% 
  FindNeighbors(reduction='harmony',dims = 1:50) %>% FindClusters(resolution = 0.25)
# cell type annotation
pdf('FeaurePlot_Vascular_dataset1_recluster_Harmony_0512.pdf',width=8,height=6)
FeaturePlot(vas_ha, features = c('COL3A1','BICC1','COL1A2','COL1A1'),order =T) # Fibroblast
FeaturePlot(vas_ha, features = c('PDGFRB','RGS5','CLDN5','ACTA2'),order =T) # pericyte + endo + smooth muscle cells
dev.off()
genes <- c('COL3A1','BICC1','COL1A2','COL1A1','PDGFRB','RGS5','CLDN5','ACTA2')
pdf('Dotplot_Vascular_dataset1_MajorMarkers_Harmony_0513.pdf',width=10,height=6)
DotPlot(vas_ha, features = genes, group.by = 'seurat_clusters', 
        cols = c('lightgrey', 'blue'), dot.scale = 6) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
genes <- c('VEGFC','ALPL','MFSD2A','SLC7A5','IL1R1','NR2F2')
pdf('FeaurePlot_Endo_subMarkers_0513.pdf',width=8,height=10)
FeaturePlot(vas_ha, features = genes,order =T) # Arterial - Capillary - Venous
dev.off()

anno <- data.frame(suerat_clusters=0:13,
                   celltype=c('Fibroblast','Endothelial cell','Fibroblast','Fibroblast',
                              'Endothelial cell','Smooth muscle cell','Pericyte','Fibroblast',
                              'MIX','Endothelial cell','MIX','Fibroblast','MIX','MIX'),
                   sub=c('P. Fibroblast c1','Capillary','P. Fibroblast c2','P. Fibroblast c3',
                         'Venous','Smooth muscle cell','Pericyte','M. Fibroblast',
                         'MIX','Arterial','MIX','P. Fibroblast c4','MIX','MIX'))
vas_ha$celltype_vas_new <- anno[match(vas_ha$seurat_clusters, anno$suerat_clusters),'celltype']
vas_ha$celltype_vas_sub <- anno[match(vas_ha$seurat_clusters, anno$suerat_clusters),'sub']
vas_ha$celltype_vas_sub <- factor(vas_ha$celltype_vas_sub, levels = c('Arterial','Capillary','Venous','P. Fibroblast c1','P. Fibroblast c2','P. Fibroblast c3',
                                                                      'P. Fibroblast c4','M. Fibroblast','Smooth muscle cell','Pericyte','MIX'))
saveRDS(vas_ha, file = 'ALS_vascular_spinalcord_Harmony_annotated_0513.rds')


# curation and final annotation
vas <- readRDS('ALS_vascular_spinalcord_Harmony_annotated_0513.rds')
Idents(vas)='celltype_vas_new'
vas <- subset(vas,idents='MIX',invert=T) # 4622 cells
vas$celltype_vas_final <- ifelse(vas$celltype_vas_sub%in%c('P. Fibroblast c1','P. Fibroblast c4'),'P. Fibroblast type1',
                                ifelse(sc$celltype_vas_sub%in%c('P. Fibroblast c2','P. Fibroblast c3'),'P. Fibroblast type2',sc$celltype_vas_sub))

saveRDS(vas,'human_als_vascular_final_1111.rds')



