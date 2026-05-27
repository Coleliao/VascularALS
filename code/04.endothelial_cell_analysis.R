# Endothelial cell analysis
# Related to Fig2
# 2026.03.06

setwd('./combine_analysis/endo')
pacman::p_load(dplyr, Seurat, ggplot2, patchwork, stringr, cowplot, ggpubr,ggrepel)
source("process.R")
source("Visualization.R")
source("tools.R")
col <- c(Arterial='#DC143C',
         Venous='#FF7F00',
         Capillary='#FB9A99',
         `M. Fibroblast`='#276419',
         `P. Fibroblast type1`='#B8E186',
         `P. Fibroblast type2`='#7FBC41',
         `Smooth muscle cell`='#B6B0FF',
         Pericyte='#5D69B1')


sc <- readRDS('../human_als_vascular_final_1116.rds')
ctx <- readRDS('../ALS_Cortex_Vasculature_final_1116.rds')

endo_sc <- subset(sc, celltype_vas_new %in% c('Endothelial cell'))
p <- DimPlot(endo_sc)
cells <- CellSelector(p)
endo_sc <- endo_sc[,cells]

endo_ctx <- subset(ctx, celltype_vas %in% c('Endothelial cell'))
p <- DimPlot(endo_ctx)
cells <- CellSelector(p)
endo_ctx <- endo_ctx[,cells]

p1 <- DimPlot(endo_sc, group.by = 'celltype_vas_final',cols = col,pt.size = 1)+theme_blank()
p2 <- DimPlot(endo_ctx, group.by = 'celltype_vas_final',cols = col,pt.size = 1)+theme_blank()

# Fig2A
pdf('endo_subcluster.pdf',width = 8,height = 3.5)
cowplot::plot_grid(p1,p2,ncol = 2,rel_widths = c(1,1.25)) 
dev.off()


# proportion significant quantitation
p3 <- cellprop(endo_sc@meta.data,sample = 'disease_status',idents = 'celltype_vas_final',col = col,show = T)
p4 <- cellprop(endo_ctx@meta.data,sample = 'disease_status',idents = 'celltype_vas_final',col = col,show = T)
pdf('endo_cell_proportion.pdf',width = 7,height = 4) # Fig2B
p3|p4
dev.off()



# BBB and periphery module
genes <- c(
  "SPOCK2","SLCO1A2","SLC2A1","SLCO1C1","SLC22A8","ABCB1","ITM2A",
  "LEF1","SLC7A5","ITIH5","BSG","MFSD2A","LRP8","CXCL12","PGLYRP1",
  "ZIC3","SLC38A3","ABHD2","TFRC","APCDD1","SLC6A6","PROM1","APOD",
  "CCDC141","SLC19A3","FLT1","PTN","FOXF2","SLC39A10","SPARCL1","ZIC2",
  "SLC1A1","ATP10A","TNFRSF19","FOXQ1","FOSB","STRA6",
  "SLC7A1","VWA1","PLTP","EFR3B","GPCPD1","SGPP2","HIAT1",
  "IGF1R","OCLN","SLC16A4","ABLIM1"
)
endo <- AddModuleScore(endo,features = list(genes),name = 'endo_score')
VlnPlot(endo,features = 'endo_score1',group.by = 'region_celltype',pt.size = 0)
endo$region_celltype_condition <- paste(endo$region_celltype,endo$disease_status,sep = '_')
VlnPlot(endo,features = 'endo_score1',group.by = 'region_celltype_condition',pt.size = 0)

meta <- endo@meta.data
group_cols <- c(control='#8DD3C7',ALS='#FDB462')
group_cols <- c('#F1BB72','#53A85F')
p5 <- ggplot(data = meta,aes(x=region_celltype,y=endo_score1,fill=disease_status))+
  geom_violin(position = position_dodge(0.8),scale = 'width')+
  geom_boxplot(outlier.size = 0,width=0.15,position = position_dodge(0.8))+
  theme_classic()+ theme(axis.text.x = element_text(size=12,color = 'black',angle = 90,vjust = 0.5),
                         axis.text.y = element_text(size=12,color = 'black'),
                         axis.title.y=element_text(size=12,colour = 'black',face = 'bold'),
                         plot.title=element_text(hjust = 0.5))+
  #geom_jitter(aes(x=region_celltype,y=periphery_score1,fill=disease_status), size=0.8, alpha=0.9)+
  labs(x='',title='',y='BBB score')+
  #theme(panel.grid = element_blank())+
  scale_fill_manual(values = group_cols)+
  scale_y_continuous(expand = c(0,.3))+
  stat_compare_means(method = 'wilcox.test',
                     label = 'p.signif',hide.ns = F,show.legend = F,size=6)
p5


genes <- c(
  "AQP1","IGFBP5","CD36","FABP4","GPIHBP1","MRC1","CYP4B1","MEIS2",
  "EDNRA","SLC43A3","LYVE1","TSPAN7","BTNL9","MMRN1","SELP","HPGD","F8",
  "PLVAP","LYZ","F2R","SERPINA3","PLPP3","NEURL3",
  "CXCL9","IGFBP4","ESM1","PTGS1","CYTL1","HLA-DRA","CD300LG","GPR182",
  "CCDC80","RELN","CD74","MLXIP","GIMAP4","ITGA9",
  "LBP","PLPP1","LDB2","FSTL1","BACE2","NEURL3"
)

endo <- AddModuleScore(endo,features = list(genes),name = 'periphery_score')
VlnPlot(endo,features = 'periphery_score1',group.by = 'region_celltype_condition',pt.size = 0)
# 有显著性检验
meta <- endo@meta.data
p6 <- ggplot(data = meta,aes(x=region_celltype,y=periphery_score1,fill=disease_status))+
  geom_violin(position = position_dodge(0.8),scale = 'width')+
  geom_boxplot(outlier.size = 0,width=0.15,position = position_dodge(0.8))+
  theme_classic()+ theme(axis.text.x = element_text(size=12,color = 'black',angle = 90,vjust = 0.5),
                         axis.text.y = element_text(size=12,color = 'black'),
                         axis.title.y=element_text(size=12,colour = 'black',face = 'bold'),
                         plot.title=element_text(hjust = 0.5))+
  #geom_jitter(aes(x=region_celltype,y=periphery_score1,fill=disease_status), size=0.8, alpha=0.9)+
  labs(x='',title='',y='Periphery score')+
  #theme(panel.grid = element_blank())+
  scale_fill_manual(values = group_cols)+
  scale_y_continuous(expand = c(0,.1))+
  stat_compare_means(method = 'wilcox.test',
                     label = 'p.signif',hide.ns = F,show.legend = F,size=6)


# Fig2C
pdf('endo_BBB_periphery_score_col2.pdf',width = 7,height = 5)
p5+theme_scp()
p6+theme_scp()
dev.off()

pdf('endo_BBB_periphery_score_v2.pdf',width = 6,height = 6)
p5+coord_flip()
p6+coord_flip()
dev.off()


FeaturePlot(endo_sc,features = c('CLDN5','OCLN','PECAM1'),split.by = 'disease_status')
VlnPlot(endo_sc,features = c('CLDN5','OCLN','PECAM1'),group.by = 'disease_status',pt.size = 0)

# cell proportion and BBB score will lead to a foucs on capillary endo, which is the most interesting part of this analysis.
Idents(endo_sc) <- 'celltype_vas_final'
cap <- subset(endo_sc,idents= 'Capillary')
Idents(cap) <- 'disease_status'
deg1 <- FindMarkers(cap,ident.1 = 'ALS',only.pos = F,min.pct = 0.2)
deg1$gene <- rownames(deg1)
palette <- c("#228B22", "#8FBC8F", "#FFFFE0", "#FFA07A", "#DC143C")
# fig2D
pdf('capillary_endo_DEG.pdf',width = 6,height = 6)
sVp1(df=deg1,logFC_col = "avg_log2FC", FDR_col = "p_val_adj", Symbol_col = "gene", 
     title = 'DEGs of capillary endothelial cells',
     y_increased=5,colours = palette) + xlab('log2FC')
dev.off()

top5 <- deg1 %>% filter(p_val_adj < 0.05) %>% slice_max(avg_log2FC,n = 5) %>% rownames()
btm5 <- deg1 %>% filter(p_val_adj < 0.05) %>% slice_max(-avg_log2FC,n = 5) %>% rownames()
# fig2E
pdf('capillary_endo_UP_DEG_top.pdf',width = 5,height = 10)
FeaturePlot(endo_sc,features = c("TNFRSF6B","PDLIM1","MPZL2","SOCS3"),
            cols = c("lightgrey", "red"),ncol = 2,split.by = 'disease_status',order = T,label = F)
dev.off()
FeaturePlot(endo_sc,features = c("CTNNA2","CSMD1","ACTN2","LRRC17"),
            cols = c("lightgrey", "red"),ncol = 2,split.by = 'disease_status',order = T,label = F)

# Go/Kegg for deg1
library(clusterProfiler)
library(org.Hs.eg.db)
deg1 <- deg1 %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>% arrange(avg_log2FC)
genelist <- list(
  up=deg1[which(deg1$avg_log2FC>1),'gene'],
  down=deg1[which(deg1$avg_log2FC<(-1)),'gene']
)
genelist <- lapply(genelist, function(x){
  ids=bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ids <- na.omit(ids)
  gene <- ids$ENTREZID
  return(gene)
})

go.bp.up = enrichGO(gene = genelist$up,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

p <- dotplot(go.bp.up,font.size=5,showCategory = 10) 
p


allgg <- as.data.frame(go.bp.up)
allgg <- allgg %>% filter(p.adjust < 0.05) %>% arrange(p.adjust)
items <- c(
  'alpha-beta T cell activation',
  'alpha-beta T cell differentiation',
  'T-helper cell differentiation',
  'lipid transport',
  'response to steroid hormone',
  'stress response to copper ion',
  'cellular response to zinc ion',
  'stress response to metal ion',
  'artery morphogenesis',
  'vascular process in circulatory system')
df <- allgg[which(allgg$Description%in%items),]
df$geneID <- sapply(df$geneID,function(x){
  genes <- str_split(x,'/')[[1]]
  symbols <- bitr(genes, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
  return(paste(symbols,collapse = '/'))
})
df$Description <- factor(df$Description,levels = rev(items))
#df <- df %>% arrange(desc(p.adjust))

options(repr.plot.width=7, repr.plot.height=6.5)
cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")
bar_cols <- viridis(option=cmap[7], direction=-1,n=20)[1:10]
library(viridis)
p <- ggplot(data = df, aes(x = Count, y = Description,fill = p.adjust)) +
  geom_bar(width = 0.5,stat = 'identity') +
  theme_classic() +
  scale_x_continuous(expand = c(0,0.5)) +
  scale_fill_gradientn(colors = bar_cols)
p <- p + theme(axis.text.y = element_blank(),axis.ticks = element_blank()) + 
  geom_text(data = df, aes(x = 0.1, y = Description, label = Description), size = 4.8,hjust = 0) 
p <- p + geom_text(data = df,aes(x = 0.1, y = Description, label = geneID),
                   size = 3.5, fontface = 'italic',hjust = 0,vjust = 2.7) 
# scale_colour_viridis(option=cmap[7], direction=1) 
# fig2F
pdf('capillary_endo_DEG_GO.pdf',width = 7,height = 6.5)
p
dev.off()



# save data
saveRDS(cap,file = 'capillary_endo.rds') # spinal cord capillary endo



