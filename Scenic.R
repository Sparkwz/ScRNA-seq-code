rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(readr)
source("~/Spark/R/custom_seurat_functions.R")

####提取out_SCENIC.loom 信息
loom <- open_loom('~/Scenic/out_SCENIC.loom') 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
embeddings <- get_embeddings(loom)  
close_loom(loom)
rownames(regulonAUC)

###加载SeuratData
seurat.data = read_rds(file = "~/Spark/Step1.pySCENIC_test_seurat.rds")
seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = F)

n.pcs = 30
seurat.data <- seurat.data %>% 
  RunUMAP(reduction = "pca", dims = 1:n.pcs, verbose = F) %>% 
  FindNeighbors(reduction = "pca", k.param = 10, dims = 1:n.pcs)

seurat.data$seurat_annotations[is.na(seurat.data$seurat_annotations)] = "B"
Idents(seurat.data) <- "seurat_annotations"

options(repr.plot.width = 6, repr.plot.height = 4.5)
DimPlot(seurat.data, reduction = "umap", label=T) 

##把SCENIC的结果与seurat数据合并
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(seurat.data))
dim(sub_regulonAUC)

cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$seurat_annotations)
head(cellTypes)
sub_regulonAUC[1:4,1:4] 

scenic_res = assay(sub_regulonAUC) %>% as.matrix()
seurat.data[["scenic"]] <- SeuratObject::CreateAssayObject(counts = scenic_res)
seurat.data <- SeuratObject::SetAssayData(seurat.data, slot = "scale.data",
                                          new.data = scenic_res, assay = "scenic")
seurat.data
#保存数据
save(sub_regulonAUC,cellTypes,seurat.data,
     file = '~/Spark/Step2.for_rss_and_visual.Rdata')


##常规转录因子可视化(挑选转录因子)
regulonsToPlot = c('RFX1(+)','EOMES(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
# Vis
p1 = DotPlot(seurat.data, features = unique(regulonsToPlot), assay = "scenic") + coord_flip() + RotatedAxis()
p2 = VlnPlot(seurat.data, features = regulonsToPlot, assay = "scenic",pt.size = 0)&labs(y="TF activity")
p3 = RidgePlot(seurat.data, features = regulonsToPlot, assay = "scenic" , ncol = 2)&labs(x="TF activity")
p4 = FeaturePlot(seurat.data,features = regulonsToPlot)

options(repr.plot.width = 15, repr.plot.height = 7)
#wrap_plots(p1,p2,p3,p4)
wrap_plots(p1,p2,p3)

#亚群特异性转录因子分析及可视化
###TF活性均值
# 看看不同单细胞亚群的转录因子活性平均值
# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
options(repr.plot.width = 4.5, repr.plot.height = 15)
Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

###rss查看特异TF
#参考文章：The RSS was first used by Suo et al. in: Revealing the Critical Regulators of Cell Identity in the Mouse Cell Atlas. Cell Reports (2018). doi: 10.1016/j.celrep.2018.10.045
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot&coord_flip())