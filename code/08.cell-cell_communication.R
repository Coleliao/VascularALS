# cell-cell communication analysis using CellChat
# spinal cord vascular data


library(CellChat)
vas <- readRDS('human_als_vascular_final_1126.rds')
Idents(vas) <- 'celltype_type'
vas <- subset(vas,idents='M. Fibroblast',invert=T) # too few in the ctl group
Idents(vas)='disease_status'
als <- subset(vas,idents='ALS')
ctr <- subset(vas,idents='control')

# run cellchat separately
future::plan("multisession", workers = 1) ## do parallel
for(i in list(als,ctr)){
  DefaultAssay(i)='RNA'
  i <- DietSeurat(i,assays = 'RNA')
  i <- NormalizeData(i)
  i$samples <- i$disease_status
  mtr <- LayerData(i, assay = "RNA", layer = "data")
  cellchat <- createCellChat(object = mtr, meta = i@meta.data, 
                             group.by = "celltype_type")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10) # M.F. in ctrl = 10
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  cat('Process Done!\n')
  try(saveRDS(cellchat,glue('./cellchat.{unique(i$samples)}.rds')))
  gc()
}
rm(i,als,augur,ctr);gc()


als=readRDS('./cellchat.ALS.rds');ctr <- readRDS('./cellchat.control.rds')
object.list <- list(ctr = ctr,als = als)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T) 
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# 2D scatter plot of signaling role analysis
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "P. Fibroblast type1")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Arterial")
patchwork::wrap_plots(plots = list(gg1,gg2))

pdf('cellchat_als_vs_ctrl_rankNet_0513.pdf',7,8)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,cutoff.pvalue = 0.01)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,cutoff.pvalue = 0.01)
gg1 + gg2
dev.off()
pdf('cellchat_als_vs_ctrl_PF2Capillary_0513.pdf',5,8)
netVisual_bubble(cellchat, sources.use = 6:7, targets.use = 3,  
                 comparison = c(1, 2), angle.x = 45,
                 signaling = c('CXCL',"LAMININ",'FN1','COLLAGEN','AGRN'))
dev.off()


pathways.show <- c("LAMININ") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

