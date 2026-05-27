#  GRN by pyscenic on capillary endothelial cells
# 
# library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(ggplot2)
library(BiocParallel)
cap <- readRDS('capillary_endo.rds')
auc <- read.table('./results/aucell_output.csv',sep=',', header=TRUE, row.names=1,stringsAsFactors=F,check.names=F)
auc2 <- CreateAssay5Object(data=t(auc))
cap[["AUC"]] <- auc2
DefaultAssay(cap) <- 'AUC'
FeaturePlot(cap,features = 'AHDC1(+)')

auc <- t(auc) 
auc <- auc[,colnames(cap)] 
all(colnames(auc) == colnames(cap)) 
rss<-calcRSS(AUC=auc,cellAnnotation=cap@meta.data$disease_status)
rssPlot <- plotRSS(rss)
pdf('capillary_endo_Scenic_RSS.pdf',width = 4,height = 7)
rssPlot$plot
dev.off()
ctl_top10 <- rss %>% as.data.frame() %>% rownames_to_column('regulon') %>% arrange(desc(control)) %>% slice_head(n=30) %>% pull(regulon)
# clt_top5 <- c("EBF1(+)","ZNF70(+)","LCORL(+)","NFIA(+)","ZNF518A(+)")
als_top10 <- rss %>% as.data.frame() %>% rownames_to_column('regulon') %>% arrange(desc(ALS)) %>% slice_head(n=30) %>% pull(regulon)
top_regulons <- unique(c(ctl_top10, als_top10))
rss_top <- rss[top_regulons, ]

rssPlot <- plotRSS(rss_top)
pdf('capillary_endo_Scenic_RSS_top10.pdf',width = 4,height = 7)
rssPlot$plot
dev.off()


#  auc heatmap for each cell
df <- auc[c(ctl_top5,als_top5),]
df %<>% t() %>% scale() %>% t() %>% as.data.frame()
pdf('heatmap_capillary_endo_Scenic_top5.pdf',width = 8,height = 4)
col_anno <- HeatmapAnnotation(
  group = cap@meta.data$disease_status,
  col = list(group = c('control' = '#8DD3C7', 'ALS' = '#FDB462'))
)
Heatmap(df, name = "AUC", show_row_names = T, show_column_names = F,
        column_title = "Top 5 regulons in control and ALS",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        cluster_rows = F, cluster_columns = F,
        top_annotation = col_anno,
        col = colorRamp2(c(0, 0.5, 1), c("lightblue", "white", "red")))
dev.off()

regulon <- c('CEBPD(+)','BHLHE40(+)','EBF1(+)','ZNF70(+)')
FeaturePlot(cap,features = c('CEBPD(+)','BHLHE40(+)'),split.by = 'disease_status')
FeaturePlot(cap,features = c('EBF1(+)','ZNF70(+)'),split.by = 'disease_status')

p <- RidgePlot(cap, features = regulon, fill.by = 'ident',stack = F,ncol = 1,cols = c('#8DD3C7','#FDB462'))
pdf('capillary_endo_RegulonTop2_Scenic_RidgePlot.pdf',width = 5,height = 7)
p
dev.off()

wilcox.test(auc['CEBPD(+)',cap$disease_status=='ALS'],auc['CEBPD(+)',cap$disease_status=='control'])
wilcox.test(auc['BHLHE40(+)',cap$disease_status=='ALS'],auc['BHLHE40(+)',cap$disease_status=='control'])
wilcox.test(auc['EBF1(+)',cap$disease_status=='ALS'],auc['EBF1(+)',cap$disease_status=='control'])
wilcox.test(auc['ZNF70(+)',cap$disease_status=='ALS'],auc['ZNF70(+)',cap$disease_status=='control'])
RidgePlot(cap,features = 'CEBPD(+)',group.by = 'disease_status',cols = c('#8DD3C7','#FDB462'))


# regulon: TF and its target genes
ctx <- fread("./results/ctx_output.csv")
head(ctx)

cebpd_targets <- ctx[which(ctx$V1=='CEBPD'), 'V9']




