# process ALS single cell data -- cortex
# subset vascular cells (particularly for fibroblast)
# Related to Fig1 and FigS1
# Kuo 2025-05-07


library(Seurat)
library(hdf5r)
library(ggplot2)
source('~/functions/process.R')
source('~/functions/tools.R')

setwd('~/MC_PFC_cell2024')
rm(list = ls());gc()

fpath <- '~/data/als_CTX_cell2024'
als_files <- list.files(fpath, pattern = 'ALS', full.names = T)
pn_files <- list.files(fpath, pattern = 'PN',full.names = T)

# 1. ALS
als <- list()
pb <- txtProgressBar(min = 0, max = length(als_files), style = 3)
for(i in 1:length(als_files)){
  setwd(als_files[i])
  tmp <- ReadMtx(mtx = 'counts_fil.mtx',cells = 'col_metadata.tsv',features = 'row_metadata.tsv',
                 cell.column = 1,feature.column = 1,skip.feature = 1, skip.cell = 1) # ENSEMBL-ID 
  tmp <- CreateSeuratObject(counts = tmp, project = 'ALS', assay = 'RNA',
                            min.cells = 3, min.features = 200)
  meta <- read.table('col_metadata.tsv', sep = '\t', header = T, row.names = 1)
  meta <- cbind(tmp@meta.data,meta)
  if(all(rownames(meta)==colnames(tmp))){tmp@meta.data <- meta}else{cat(als_files[i],'meta data not match!\n')}
  als[i] <- tmp
  setTxtProgressBar(pb, i) 
}
close(pb)
als <- merge(als[[1]], als[-1])
gc()
saveRDS(als, file = 'ALS_CTX+MC_dataset2_allcell_Noprocess_0507.rds') # 333211 cells
Idents(als) <- 'CellClass'
vas <- subset(als, idents = c('Vasc'))
saveRDS(vas, file = 'ALS_CTX+MC_dataset2_Vasc_0507.rds') # 8624


# 2. PN 
pn <- list()
pb <- txtProgressBar(min = 0, max = length(pn_files), style = 3)
for(i in 1:length(pn_files)){
  setwd(pn_files[i])
  tmp <- ReadMtx(mtx = 'counts_fil.mtx',cells = 'col_metadata.tsv',features = 'row_metadata.tsv',
                 cell.column = 1,feature.column = 1,skip.feature = 1, skip.cell = 1) # 
  tmp <- CreateSeuratObject(counts = tmp, project = 'ALS', assay = 'RNA',
                            min.cells = 3, min.features = 200)
  meta <- read.table('col_metadata.tsv', sep = '\t', header = T, row.names = 1)
  meta <- cbind(tmp@meta.data,meta)
  if(all(rownames(meta)==colnames(tmp))){tmp@meta.data <- meta}else{cat(pn_files[i],'meta data not match!\n')}
  pn[i] <- tmp
  setTxtProgressBar(pb, i)
}
close(pb)
pn <- merge(pn[[1]], pn[-1])
gc()
saveRDS(pn, file = 'CTRL_CTX+MC_dataset2_allcell_Noprocess_0507.rds') # 123196 cells
Idents(pn) <- 'CellClass'
vas <- subset(pn, idents = c('Vasc')) # 2225
saveRDS(vas, file = 'PN_CTX+MC_dataset2_Vasc_0507.rds') 


# 3. merge vascular cells
vas <- readRDS('ALS_CTX+MC_dataset2_Vasc_0507.rds')
vas2 <- readRDS('PN_CTX+MC_dataset2_Vasc_0507.rds')
vas <- merge(vas, vas2)
# transfer ENSEMBL-ID to gene name by the row_metadata.tsv 
en2gene <- read.table('~/data/als_CTX_cell2024/210526_PN_301_snRNA-D7/row_metadata.tsv', 
                      sep = '\t', header = T, row.names = 1)
gene <- ifelse(en2gene[rownames(vas),'Gene']!='N/A',
               en2gene[rownames(vas),'Gene'], rownames(vas))
rownames(vas@assays$RNA@counts) <- gene
rownames(vas@assays$RNA@data) <- gene



# Harmony integration
m=table(vas$Donor)
Idents(vas) <- 'Donor'
vas <- subset(vas, idents = names(m[m>100]) ) # samples with too few cells will be removed to avoid over-correction by harmony
library(harmony)
cp <- vas
vas <- vas %>%
  subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000) %>% 
  SCTransform() %>% RunPCA() 
vas <- RunHarmony(vas, group.by.vars = 'Donor', plot_convergence = F, max.iter.harmony = 20)
vas <- vas %>%  RunUMAP(reduction='harmony',dims = 1:30) %>% 
  FindNeighbors(reduction='harmony',dims = 1:30) %>% FindClusters(resulution = 0.5)


pdf('Vasc_dataset2_CTX_Harmony_reculsterd_0507.pdf', width = 10, height = 8)
DimPlot(vas, reduction ="umap", label=T,pt.size = 1)
DimPlot(vas, reduction ="umap", group.by = "full_label",label=T,pt.size = 1)
DimPlot(vas, reduction ="umap", group.by = "Donor",label=F,pt.size = 1)
DimPlot(vas, reduction ="umap", group.by = "Condition",label=F,pt.size = 1)
DimPlot(vas, reduction ="umap", group.by = "Region",label=F,pt.size = 1)
dev.off()

res <- consist(vas@meta.data,'full_label','Condition',plot=T)

DefaultAssay(vas) <- 'RNA';vas <- NormalizeData(vas)
pdf('FeaurePlot_Vascular_dataset2_recluster_0508.pdf',width=8,height=6)
FeaturePlot(vas, features = c('COL3A1','BICC1','COL1A2','COL1A1'),order =T) # Fibroblast
FeaturePlot(vas, features = c('PDGFRB','RGS5','CLDN5','ACTA2'),order =T) # pericyte + endo + smooth muscle cells
dev.off()


saveRDS(vas, file = 'Vasc_dataset2_CTX_Harmony_reculsterd_0507.rds')




# data clearing and annotation
vas <- readRDS('./Vasc_dataset2_CTX_Harmony_reculsterd_0507.rds')
vas <- subset(vas,idents=c(8,11,10,16),invert=T) # t cells and sample-specific clusters removed
vas <- FindClusters(vas,resolution=0.5)
DimPlot(vas,reduction='umap',group.by='seurat_clusters',label=T)
DotPlot(vas,features=c('CLDN5','FLT1','PECAM1','VWF','PDGFRB','RGS5','ACTA2','COL1A1','COL1A2','COL3A1','DCN','LUM'),group.by='seurat_clusters') + RotatedAxis()

FeaturePlot(vas,features=endo_marker_sub,ncol=3,order = T)

# reclustering after removing contaminating clusters
Idents(vas) <- 'Donor'
m <- table(vas$Donor)
vas <- subset(vas,idents=names(m[m>=100])) # 8835 cells
vas <- DietSeurat(vas,assay='RNA',layers=c('counts','data'))
vas[["RNA"]] <- split(vas[["RNA"]], f = vas$Donor) 
vas <- NormalizeData(vas)
vas <- FindVariableFeatures(vas)
vas <- ScaleData(vas)
vas <- RunPCA(vas)
vas <- IntegrateLayers(
  object = vas, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
vas <- JoinLayers(vas,layers=c('counts','data'))
vas <- FindNeighbors(vas, dims = 1:30, reduction = "harmony")
vas <- RunUMAP(vas, reduction = "harmony", dims = 1:30, reduction.name = "umap",
                spread = 0.8, min.dist = 0.3)
vas <- FindClusters(vas, resolution = 0.5)
vas <- subset(vas,idents='11',invert=T) # removing clusters with cells < 50, 8787 cells

DimPlot(vas,reduction='umap',group.by='seurat_clusters',label=T)
FeaturePlot(vas,features=c('LAMA1','LAMA2','KCNMA1','SLC4A4'),order = T)
FeaturePlot(vas,features=c('TIAM1','CDH19','NAV2','CEMIP'),order = T)
FeaturePlot(vas,features=c('SLC6A13','LAMA1'),order = T)

anno <- data.frame(seurat_clusters=c(0,8,10,4,3,7,9,6,1,2,5),
                   celltype_vas=c(rep('Pericyte',3),'Smooth muscle cell',
                                  rep('Fibroblast',3),rep('Endothelial cell',4)),
                   celltype=c(rep('Pericyte',3),'Smooth muscle cell',
                              'P. Fibroblast','P. Fibroblast','M. Fibroblast',
                              'Arterial',rep('Capillary',2),'Venous'))
vas$celltype_vas_final <- anno$celltype[match(vas$seurat_clusters,anno$seurat_clusters)]
vas$celltype_vas <- anno$celltype_vas[match(vas$seurat_clusters,anno$seurat_clusters)]
ctx <- vas
ctx$celltype_vas_final <- factor(ctx$celltype_vas_final, 
                                 levels = c('Arterial','Capillary','Venous','Pericyte',
                                            'Smooth muscle cell','P. Fibroblast', 'M. Fibroblast'))


DimPlot(ctx,reduction='umap',group.by='celltype_vas_final',label=T,
              repel = T,cols=col)+xlab('UMAP_1')+ylab('UMAP_2')

saveRDS(vas,'ALS_Cortex_Vasculature_processed_1116.rds')

