# spinal cord and coretex combined analysis
# Related to fig1 and figS1
# Kuo

library(Seurat)
setwd('./combine_analysis')
source('Visualization.R')
source('tools.R')


ctx <- readRDS('ALS_Cortex_Vasculature_processed_1116.rds')
col <- c(Arterial='#DC143C',
         Venous='#FF7F00',
         Capillary='#FB9A99',
         `M. Fibroblast`='#276419',
         `P. Fibroblast`='#B8E186',
         `Smooth muscle cell`='#B6B0FF',
         Pericyte='#5D69B1')
p1 <- DimPlot(ctx,reduction='umap',group.by='celltype_vas_final',label=T,
              repel = T,cols=col)+xlab('UMAP_1')+ylab('UMAP_2')



# spinal cord
sc <- readRDS('human_als_vascular_final_1111.rds')
col <- c(Arterial='#DC143C',
         Venous='#FF7F00',
         Capillary='#FB9A99',
         `M. Fibroblast`='#276419',
         `P. Fibroblast type1`='#B8E186',
         `P. Fibroblast type2`='#7FBC41',
         `Smooth muscle cell`='#B6B0FF',
         Pericyte='#5D69B1')
sc$celltype_vas_final <- factor(sc$celltype_vas_final, 
                                levels = c('Arterial','Capillary','Venous','Pericyte',
                                           'Smooth muscle cell','P. Fibroblast type1',
                                           'P. Fibroblast type2', 'M. Fibroblast'))

p2 <- DimPlot(sc,reduction='umap',group.by='celltype_vas_final',label=T,
              repel = T,cols=col)

pdf('Fig1_vasculature_UMAP_CTX_SC_0128.pdf',width=12,height=4)
p2 + p1
dev.off()




# deg
Idents(sc) <- sc$disease_status
table(sc$disease_status)
deg_sc <- data.frame()
for(i in unique(sc$celltype_vas)){
  print(i)
  sub_data <- subset(sc,celltype_vas==i)
  deg <- FindMarkers(sub_data,ident.1='ALS',ident.2='control',logfc.threshold = 0.05,min.pct = 0.1)
  deg$gene <- rownames(deg)
  deg$celltype_vas <- i
  deg_sc <- rbind(deg_sc,deg)
}

df <- subset(deg_sc,deg_sc$p_val_adj<0.05) # & deg_sc$pct.1 > 0.25
col4 <- c(Pericyte='#5D69B1',
          `Smooth muscle cell`='#B6B0FF',
          `Fibroblast`='#4D9221',
          `Endothelial cell`='#d73027')
p1 <- multiVolcano(df,group = 'celltype_vas',title.col = col4,log2FC.cutoff = 0.5,
                   title.size = 4,gene.size = 4)+NoLegend()
p1

Idents(ctx) <- ctx$Condition
table(ctx$Condition)
deg_ctx <- data.frame()
for(i in unique(ctx$celltype_vas)){
  print(i)
  sub_data <- subset(ctx,celltype_vas==i)
  deg <- FindMarkers(sub_data,ident.1='ALS',ident.2='PN',logfc.threshold = 0.05,min.pct = 0.1)
  deg$gene <- rownames(deg)
  deg$celltype_vas <- i
  deg_ctx <- rbind(deg_ctx,deg)
}
df2 <- subset(deg_ctx,deg_ctx$p_val_adj<0.05)
p2 <- multiVolcano(df2,group = 'celltype_vas',title.col = col4,log2FC.cutoff = 0.5,
                   celltypeSize=4,gene.size = 4)+NoLegend()

pdf('Fig2_vasculature_DEG_volcano_CTX_SC_1116.pdf',width=9,height=5)
print(p1);print(p2);
dev.off()

write.csv(deg_sc,'ALS_SC_vasculature_DEG_by_celltype_1116.csv',row.names = F)
write.csv(deg_ctx,'ALS_CTX_vasculature_DEG_by_celltype_1116.csv',row.names = F)

saveRDS(ctx,'ALS_Cortex_Vasculature_final_1116.rds')
saveRDS(sc,'human_als_vascular_final_1116.rds')

#deg_sc <- read.csv('ALS_SC_vasculature_DEG_by_celltype_1116.csv')
#deg_ctx <- read.csv('ALS_CTX_vasculature_DEG_by_celltype_1116.csv')

deg_sc <- split(deg_sc,deg_sc$celltype_vas)
names( deg_sc ) <- paste0( 'SC_', names(deg_sc))
deg_ctx <- split(deg_ctx,deg_ctx$celltype_vas)
names( deg_ctx ) <- paste0( 'CTX_', names(deg_ctx))
deg <- c(deg_sc,deg_ctx)
up_list <- lapply(deg,function(x) {
  sub <- subset(x,x$avg_log2FC>0.5 & x$p_val_adj<0.05)
  return(sub$gene)
})
library(UpSetR)
p3 <- upset(fromList(up_list),  
            nsets = 100,     
            nintersects = 40, 
            order.by = "freq", 
            keep.order = T, 
            mb.ratio = c(0.6,0.4),  
            text.scale = 2 
)

inter <- VennDiagram::get.venn.partitions(up_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ',')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
inter <- inter[order(inter$..count..,decreasing = T),]


down_list <- lapply(deg,function(x) {
  sub <- subset(x,x$avg_log2FC< -0.5 & x$p_val_adj<0.05)
  return(sub$gene)
})
p4 <- upset(fromList(down_list),  
            nsets = 100,     
            nintersects = 40, 
            order.by = "freq", 
            keep.order = T, 
            mb.ratio = c(0.6,0.4),   
            text.scale = 2 
)

intersect(up_list$`SC_Endothelial cell`,up_list$`CTX_Endothelial cell`)
# "FKBP5"   "AKAP12"  "ERO1A"   "PDE10A"  "NEDD9"   "CCNY"    "HILPDA"  "ACSL5"   "NXN"     "GALNT18" "SAMD4A"

# Fig S1G and S1F
pdf('Upset_vasculature_DEG_up_down_overlap_CTX_SC_1116.pdf',width=12,height=6.5)
print(p3)
print(p4)
dev.off()





