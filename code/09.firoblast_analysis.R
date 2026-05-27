# 2026/05/13本地跑一下
# type2 的细分和trajectory

pacman::p_load(Seurat, dplyr, ggplot2, monocle3, patchwork)

col2 <- c(Arterial='#DC143C',
          Capillary='#FB9A99',
          Venous='#FF7F00',
          `P. Fibroblast c1`='#B8E186',
          `P. Fibroblast c4`='#7FBC41',
          `P. Fibroblast c2`='#5ecc74',
          `P. Fibroblast c3`='#4D9221',
          `M. Fibroblast`='#276419',
          `Smooth muscle cell`='#B6B0FF',
          Pericyte='#5D69B1')

vas <- readRDS('./human_als_vascular_final_1126.rds')
Idents(vas) <- vas$celltype_type


# cell proportion
meta <- vas@meta.data;meta$donor_id <- as.character(meta$donor_id)
df <- meta %>% group_by(donor_id,celltype_type) %>% summarise(count=n())%>% 
  group_by(donor_id) %>% mutate(ration=count/sum(count)*100)
df$disease_status <- ifelse(df$donor_id%in%c('13085','13122','13394','13446'),'Health','ALS')
library(ggalluvial)
# fig 4A
pdf('vascular_celltype_alluvial_donor_0512.pdf',width = 6,height =4)
ggplot(df, aes(x = donor_id,y=ration,fill = celltype_type, 
               stratum = celltype_type, alluvium = celltype_type)) +
  geom_col(width = 0.7,color=NA)+
  #geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0)  
  geom_alluvium(width = 0.5,alpha = 0.3)+  
  scale_fill_manual(values = col)+
  theme_bw()+
  theme(axis.text.x=element_text(size=10,vjust = 0.5,angle = 90,hjust = 0.5,colour = 'black'),
        legend.position = 'right',axis.title.x=element_text(colour = 'black'))
dev.off()


# GO analysis - fig S3
Idents(vas) <- 'celltype_type'
DimPlot(vas,group.by='celltype_type',label=F,label.size=4,cols = col)
fibro <- subset(vas,idents=c('P. Fibroblast type1','P. Fibroblast type2','M. Fibroblast'))

library(clusterProfiler)
Idents(fibro)='celltype_type' # celltype_vas_major
degs <- FindAllMarkers(fibro,only.pos = T,min.pct = 0.2,logfc.threshold = 0.25)
top <- degs %>% filter(p_val_adj<0.05&avg_log2FC>1.5) %>% as.data.frame()
table(top$cluster)
genelist <- list()
for(i in unique(top$cluster)){
  gene <- top[which(top$cluster==i),'gene']
  ids=bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ids <- na.omit(ids)
  gene <- ids$ENTREZID
  genelist[[i]] <- gene
}
sam.go.BP <- compareCluster(genelist,
                            fun = "enrichGO",
                            OrgDb = "org.Hs.eg.db",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05
)
p <- dotplot(sam.go.BP,font.size=5,title="GO of stromal cells",showCategory = 5) 
df <- sam.go.BP@compareClusterResult
pdf('./GO_deg_fibroblast_subtypes_Harmony_0715.pdf',9,14)
ggplot(df,aes(x=Cluster,y = Description,size = GeneRatio, color = p.adjust))+
  geom_point() +
  scale_y_discrete(position = "right",labels=function(x) stringr::str_wrap(x, width=80)) + 
  scale_color_gradientn(colours = rev(viridis::viridis(20))) + 
  cowplot::theme_cowplot() +
  ylab("") + xlab("") + theme_bw() + 
  theme(
    axis.text.x = element_text(size=10, angle=90, hjust=0.5, color="black"),
    axis.text.y = element_text(size=10, color="black",margin = c(5,1)),
    axis.title = element_text(size=14)
  )
dev.off()



# focus on P. Fibroblast type2
t2 <- subset(vas,idents='P. Fibroblast type2')
DimPlot(t2,group.by = 'disease_status',shuffle = T)
umap <- Embeddings(t2, 'umap') %>% as.data.frame()
cells <- umap %>% filter(UMAP_2 < 0) %>% rownames() # there is an outlier
t2 <- t2[,cells]
cellprop(t2@meta.data,idents = 'disease_status',sample = 'celltype_sub',show = T)
Idents(t2) <- 'celltype_sub'

# rename
t2$celltype_sub <- plyr::revalue(t2$celltype_sub,replace = c('P. Fibroblast c2'='c2','P. Fibroblast c3'='c1')) # 因为我不需要c1/c4了

colt2 <- c('c1'='#4D9221','c2'='#5ecc74')
col3 <- c('control' = '#8DD3C7', 'ALS' = '#FDB462')
p1 <- DimPlot(t2,shuffle = T,cols = colt2)
p2 <- DimPlot(t2,group.by = 'disease_status',shuffle = T,cols = col3)

# fig 4E
pdf('./t2_umap.pdf',width = 8,height = 4)
p1+p2
dev.off()

FeaturePlot(t2,features = c('ACTA2','COL6A2','KLF3'),order = T,label = T,split.by = 'disease_status')&scale_color_viridis_c()
VlnPlot(t2,features = c('ACTA2','ACTB'),pt.size = 0,group.by = 'celltype_sub',stack = T)


# major fibroblast state module score
# gene list from https://www.nature.com/articles/s41588-025-02284-1 

library(readxl)
genes <- read_excel('./fibroblast-genes.xlsx') %>% as.data.frame() 

myof <- genes[which(genes$`transcriptional states`=='myofibroblast'),'gene']
t2 <- AddModuleScore(t2,features = list(myof),name = 'myof_score')
FeaturePlot(t2,features = 'myof_score1',order = T,label = T,split.by = 'disease_status')&scale_color_viridis_c()
VlnPlot(t2,features = 'myof_score1',pt.size = 0,group.by = 'celltype_sub')

inf <- genes[which(genes$`transcriptional states`=='Inflammatory'),'gene']
t2 <- AddModuleScore(t2,features = list(inf),name = 'inf_score')
FeaturePlot(t2,features = 'inf_score1',order = T,label = T,split.by = 'disease_status')&scale_color_viridis_c()
VlnPlot(t2,features = 'inf_score1',pt.size = 0,group.by = 'celltype_sub')
FeaturePlot(t2,features= c('LOXL2','COL5A2', 'COL8A1','SPARC'),
            order = T,label = T)&scale_color_viridis_c()

pan <- genes[which(genes$`transcriptional states`=='Universal'),'gene']
t2 <- AddModuleScore(t2,features = list(pan),name = 'pan_score')
FeaturePlot(t2,features = 'pan_score1',order = T,label = T,split.by = 'disease_status')&scale_color_viridis_c()
VlnPlot(t2,features = 'pan_score1',pt.size = 0,group.by = 'celltype_sub')

antigen <- genes[which(genes$`transcriptional states`=='antigen presentation'),'gene']
t2 <- AddModuleScore(t2,features = list(antigen),name = 'antigen_score')
FeaturePlot(t2,features = 'antigen_score1',order = T,label = T,split.by = 'disease_status')&scale_color_viridis_c()
VlnPlot(t2,features = 'antigen_score1',pt.size = 0,group.by = 'celltype_sub')

df <- t2@meta.data[,c('celltype_sub','disease_status','myof_score1','inf_score1','pan_score1','antigen_score1')]
library(ggpubr)
library(reshape2)
df2 <- melt(df,id.vars = c('celltype_sub','disease_status'),variable.name = 'module',value.name = 'score')
df2$module <- plyr::revalue(df2$module,replace = c('myof_score1'='Myofibroblast','inf_score1'='Inflammatory','pan_score1'='Universal','antigen_score1'='Antigen presentation'))
df2$module <- factor(df2$module,levels = c('Universal','Myofibroblast','Inflammatory','Antigen presentation'))
# fig 4F
pdf('./t2_fibro_module_score.pdf',width = 6,height = 2)
ggplot(df2,aes(x=module,y=score,fill=celltype_sub))+
  geom_violin(position=position_dodge(width = 0.8))+
  geom_boxplot(outlier.size = 0,width=0.2,position=position_dodge(width = 0.8))+
  #facet_wrap(~module,scales = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_rect(color = 'black',size=1))+
  scale_fill_manual(values = colt2)+ # c('#B6B0FF','#5D69B1')
  stat_compare_means(aes(group = celltype_sub),method = 'wilcox.test',label = 'p.signif')
dev.off()




# trajectory analysis using monocle3
cds <- monocle3::new_cell_data_set(
  expression_data = GetAssayData(t2, assay = "RNA", layer = "counts"),
  cell_metadata = t2@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(t2), row.names = rownames(t2))
)
# preprocess
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "donor_id") 

# reduce dimension (use UMAP from Seurat -- optional)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(t2, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed

p1 <- plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "celltype_sub",
                 group_label_size = 5,cell_size = 2,trajectory_graph_color = 'red',
                 label_cell_groups = F)+
  scale_color_manual( values = col2)+
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_rect(color = 'black',size=1))
p1

# cluster cells
cds <- cluster_cells(cds)
plot_cells(cds, label_groups_by_cluster = TRUE, color_cells_by = "partition") 

# learn graph
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype_sub", 
           label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE)

# order cells
cds <- order_cells(cds) 
p2 <- plot_cells(cds, color_cells_by = "pseudotime", 
                 label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
                 cell_size = 2,trajectory_graph_color = 'red')+theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_rect(color = 'black',size=1))
p2

save_monocle_objects(cds, '.')
cds <- load_monocle_objects(directory_path='.')

# fig 4G
pdf('./monocle3_pseudotime.pdf',width = 10,height = 4)
p1+p2
dev.off()

# DEG along pseudotime
deg_pseudotime <- graph_test(cds, neighbor_graph="principal_graph", cores=4,reduction_method = "UMAP")
top20 <- deg_pseudotime %>%
  top_n(n = 20, morans_I) %>%
  pull(gene_short_name) %>%
  as.character()

pdf('./monocle3_pseudotime_top_genes.pdf',width = 7,height = 16)
plot_genes_in_pseudotime(cds[top20,],
                         color_cells_by="celltype_sub",
                         min_expr=0.5)+scale_color_manual( values = col2)
dev.off()

genes <- c('APOD','COL28A1','CDH19','TNFRSF19','THSD4','SORBS1','GFRA1','SLC2A1')
genes2 <- c('CACNA2D3','SLC7A11','CSMD3','GFRA1','PDZRN4','ABCA8')

p1 <- plot_genes_in_pseudotime(cds[genes[1:4],],
                               color_cells_by="celltype_sub",
                               min_expr=0.5)+scale_color_manual(values = col2)
p2 <- plot_genes_in_pseudotime(cds[genes[5:8],],
                               color_cells_by="celltype_sub",
                               min_expr=0.5)+scale_color_manual(values = col2)
pdf('./monocle3_pseudotime_top_genes_selected.pdf',width = 14,height = 4)
p1+p2
dev.off()

gene1 <- c('AEBP1','APOD','CDKN1A','CDKN2A','CIITA','HLA-DOA','COL6A2')
gene2 <- c('FGL2','HGF','ICAM1','ITIH5','LAMA1','LAMB1','NNMT','PCOLCE')
metaolism <- c('HK2','PFKP','SLC2A1','SLC2A3','SDHD','NDUFA9','PDK4','PDK1','SLC1A1')
plot_genes_in_pseudotime(cds[gene2,],
                         color_cells_by="celltype_sub",
                         min_expr=0.5)+scale_color_manual( values = col2)
plot_genes_in_pseudotime(cds[metaolism,],
                         color_cells_by="celltype_sub",
                         min_expr=0.5)+scale_color_manual( values = col2)






