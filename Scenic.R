library(Seurat)
library(dplyr)
library(readr)
###加载示例数据
install.packages('~/Rpackages/pbmc3k.SeuratData_3.1.4.tar.gz', repos = NULL, type = "source")
library(pbmc3k.SeuratData)	
# 加载该数据集
data("pbmc3k")
# 查看数据
pbmc3k = UpdateSeuratObject(pbmc3k)
pbmc3k
#注意矩阵一定要转置，不然会报错
write.csv(t(as.matrix(pbmc3k@assays$RNA@counts)),file = "~/Spark/for.pyscenic.csv")
write_rds(pbmc3k, file = "~/Spark/Step1.pySCENIC_test_seurat.rds")

################可视化################
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

#####把SCENIC的结果与seurat数据合并####
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

#热图查看TF分布：
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

### 4.3 其他查看TF方式
rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat = data.frame(
                 path = rownames(rss),
                 cluster = colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)

options(repr.plot.width = 4.5, repr.plot.height = 8)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)

#Rank图
library(parallel)
library(philentropy)
library(ggrepel)
library(latex2exp)
### 5.1 读入RAS矩阵
rasMat <- seurat.data@assays$scenic@counts %>% t() %>% as.data.frame()

### 5.2 读入细胞类型矩阵
cell.info <- seurat.data@meta.data
cell.info$celltype <- seurat.data@meta.data$seurat_annotations

cell.types <- names(table(cell.info$celltype))
ctMat <- lapply(cell.types, function(i) {
  as.numeric(cell.info$celltype == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- cell.types
rownames(ctMat) <- rownames(cell.info)
head(ctMat)
### 5.3 计算RSS矩阵(Regulon Specificity Score)
options(ggrepel.max.overlaps = Inf)
rssMat <- mclapply(colnames(rasMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
}, mc.cores = 10)
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(rasMat)
colnames(rssMat) <- colnames(ctMat)
p.rank = lapply(colnames(rssMat), function(ct){
  PlotRegulonRank(rssMat, ct, topn = 5, front.size = 10, point.size = 1)
})
names(p.rank) = colnames(rssMat)

options(repr.plot.width = 9, repr.plot.height = 9)
wrap_plots(p.rank)
