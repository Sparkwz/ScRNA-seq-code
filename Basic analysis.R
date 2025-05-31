library(Seurat) 
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)
library(harmony)
#check the current active plan
plan()
# change the current plan to access parallelization
plan("multisession", workers =4)
plan()

######################数据读入######################
## 1.1 lappy循环依次读入数据
fileID = list.files("~/Spark/matrix/")
path = "~/Spark/matrix/"
fileID

# 运行lapply函数读入数据
seurat.list = lapply(fileID, function(file){
  seurat_data <- Read10X(data.dir = paste0(path, file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.cells = 0,
                                   min.features = 0,
                                   project = file)
  return(seurat_obj)
})
names(seurat.list) = fileID
seurat.list

##merge多个seurat对象
# Create a merged Seurat object
merged_seurat <- merge(x = seurat.list[[1]],
                       y = seurat.list[-1],
                       add.cell.id = fileID)
merged_seurat
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

##表型数据载入
merged_seurat$sampleID = merged_seurat$orig.ident

#根据样本信息重命名编组
merged_seurat$group <- recode(merged_seurat$sampleID,
                              "SRR7722939" = "PBMC_Pre",
                              "SRR7722940" = "PBMC_Disc_Early",
                              "SRR7722941" = "PBMC_Disc_Resp",
                              "SRR7722942" = "PBMC_Disc_AR",
                              "SRR7722937" = "Tumor_Disc_Pre",
                              "SRR7722938" = "Tumor_Disc_AR")
head(merged_seurat@meta.data)

#保存
saveRDS(merged_seurat,file = "~/Spark/Step1.RawCount_merged_seurat.rds")

######################质量控制######################
##读入数据##
seurat.data = read_rds(file = "~/Spark/Step1.RawCount_merged_seurat.rds")
seurat.data
##查看数据##
head(seurat.data@meta.data)
table(seurat.data$group)
##筛选样本##
pbmc = subset(seurat.data, group %in% c("PBMC_Disc_AR", "PBMC_Disc_Early", 
                                        "PBMC_Disc_Resp", "PBMC_Pre"))
table(pbmc$group)

###计算QC指标，如线粒体、核糖体和血红细胞###
##线粒体百分比###
#注意人的是 ^MT-；小鼠的则是 ^mt-
mito_genes=rownames(pbmc)[grep("^MT-", rownames(pbmc))]
mito_genes

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)

#可视化
#options(repr.plot.width=10, repr.plot.height=5)
VlnPlot(pbmc,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        group.by = "group")
dev.off()

##计算核糖体基因比例##
#注意人的是^RP[SL]；小鼠的则是 ^Rp[sl]
ribo_genes=rownames(pbmc)[grep("^RP[SL]", rownames(pbmc),ignore.case = T)]
ribo_genes

pbmc=PercentageFeatureSet(pbmc, "^RP[SL]",col.name = "percent.ribo")

##计算红血细胞基因比例##
#注意人的是 ^HB[^(P)]；小鼠的 ^Hb[^(p)]
hb_genes <- rownames(pbmc)[grep("^HB[^(P)]", rownames(pbmc),ignore.case = T)]
hb_genes
pbmc=PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent.hb")

##可视化
#options(repr.plot.width=10, repr.plot.height=10)
VlnPlot(pbmc,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
        ncol = 3,
        group.by = "group")

##过滤
pbmc.qc <- subset(pbmc, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & 
                    nCount_RNA> 500 & percent.mt < 15)

#如果需要根据线粒体、核糖体、血红蛋白过滤
#pbmc.qc <- subset(pbmc, subset = percent.mt < 15 & percent.ribo> 3 & percent.hb < 0.1)

#保存
saveRDS(pbmc.qc,file = "~/Spark/Step2.PBMC_afterQC.rds")

######################降维注释######################
#设置可用内存(!!!)
options(future.globals.maxSize = 10 * 1024^3)
###读入数据
seurat.data = read_rds(file = "~/Spark/Step2.PBMC_afterQC.rds")
seurat.data

###三步走：标准化、特征选择和归一化分析
seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F)
seurat.data

##降维和聚类，检查批次
###降维聚类(PCA+UMAP)
seurat.data = seurat.data %>% 
  RunPCA(npcs = 30, verbose = F) %>% 
  #RunTSNE(reduction = "pca", dims = 1:30, verbose = F) %>% 
  RunUMAP(reduction = "pca", dims = 1:30, verbose = F)
seurat.data
##检查批次
options(repr.plot.width = 10, repr.plot.height = 4.5)
p1.compare=wrap_plots(ncol = 2,
                      DimPlot(seurat.data, reduction = "pca", group.by = "sampleID")+NoAxes()+ggtitle("Before_PCA"),
                      DimPlot(seurat.data, reduction = "umap", group.by = "sampleID")+NoAxes()+ggtitle("Before_UMAP"),
                      guides = "collect"
)
p1.compare #存在批次效应
#保存去批次前图像
ggsave(plot=p1.compare, filename="~/Spark/Step3.Before_inter_sum.pdf", width = 10 ,height = 4.5)

##去批次
seurat.data <- seurat.data %>% RunHarmony("sampleID", plot_convergence = T)
seurat.data #出现pca,umap,harmony

#Check the generated embeddings:
harmony_embeddings <- Embeddings(seurat.data, 'harmony')
harmony_embeddings[1:5,1:5]

###RunUMAP及聚类(去批次后)
n.pcs = 20
seurat.data <- seurat.data %>% 
  RunUMAP(reduction = "harmony", dims = 1:n.pcs, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:n.pcs)
##检查批次
p2.compare=wrap_plots(ncol = 2,
                      DimPlot(seurat.data, reduction = "harmony", group.by = "sampleID")+NoAxes()+ggtitle("After_PCA (harmony)"),
                      DimPlot(seurat.data, reduction = "umap", group.by = "sampleID")+NoAxes()+ggtitle("After_UMAP"),
                      guides = "collect"
)
p2.compare #无明显批次效应
#展示前后图像
options(repr.plot.width = 10, repr.plot.height = 9)
wrap_plots(p1.compare, p2.compare, ncol = 1)
ggsave(plot=p2.compare, filename="~/Spark/Step3.After_inter_Harmony.pdf", width = 10 ,height = 4.5)

#绘制Umap图
#设置批量分辨率
for (res in c(0.05,0.1,0.3,0.5,0.8,1,1.2,1.4,1.5,2)){
  print(res)
  seurat.data <- FindClusters(seurat.data, resolution = res, algorithm = 1)%>% 
    identity()
}
options(repr.plot.width = 20, repr.plot.height = 8)
#umap可视化
cluster_umap <- wrap_plots(ncol = 5,
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.05", label = T) & NoAxes(),  
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T)& NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T) & NoAxes(), 
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.4", label = T)& NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.5", label = T)& NoAxes()
)
cluster_umap
ggsave(cluster_umap,filename = "~/Spark/Step3.After_inter.cluster_umap_Harmony.pdf",
       width = 25, height = 9)
dev.off()

#选择合适的分辨率
Idents(object = seurat.data) <- "RNA_snn_res.0.3"
options(repr.plot.width = 6, repr.plot.height = 5)
DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T)& NoAxes()

#########################亚群注释########################
#####COSG方法#####
library(COSG)
marker_cosg <- COSG::cosg(
  seurat.data,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  expressed_pct=0.1,
  remove_lowly_expressed = T,
  n_genes_user=200)
markers = as.data.frame(marker_cosg$names)
head(markers, 10)
write.csv(markers, file = "~/Spark/Step3.COSG_res.csv")

options(repr.plot.width = 7.5, repr.plot.height = 4)

##绘制标志基因气泡图
options(repr.plot.width = 7.5, repr.plot.height = 7)
check_genes = c("RGS1","PTPRC",'TYROBP', #Leukocytes (Leu) 
                'CD3D','CD3E',"CD3G","CD2",'TRAC','IL32', #T cells
                'SELL',"CCR7","LEF1","TCF7",'IL7R', #Naive
                "GZMA", "GZMB","IFNG","PRF1","GNLY", #Effect/cytotoxic
                'CD4','CD40LG','CD8A','CD8B','IL2RA','FOXP3',
                "NKG7","KLRD1","IFIT1", #NK cells
                "IGHM","CD22","CD79A","CD19","MS4A1","SDC1", #B cells
                "JCHAIN","MZB1","PRDM1","IGJ" #Plasma cells
)
DotPlot(object = seurat.data, features = check_genes,assay = "RNA",scale = T) + 
  coord_flip()

###分配细胞名称
celltype=data.frame(ClusterID=0:14,celltype='NA')

## Others
celltype[celltype$ClusterID %in% c(7),2]='Platelets'
celltype[celltype$ClusterID %in% c(11),2]='Erythroid cells'

## NK/T
celltype[celltype$ClusterID %in% c(0),2]='NK'
celltype[celltype$ClusterID %in% c(3),2]='NKT'
celltype[celltype$ClusterID %in% c(4),2]='CD4T'
celltype[celltype$ClusterID %in% c(5),2]='CD8T'
celltype[celltype$ClusterID %in% c(10),2]='T cycling'

## B细胞
celltype[celltype$ClusterID %in% c(2),2]='B'
celltype[celltype$ClusterID %in% c(8),2]='B plasma'

## 髓系
celltype[celltype$ClusterID %in% c(1, 6, 12),2]='Mono/Mac'
celltype[celltype$ClusterID %in% c(9),2]='DCs'
celltype[celltype$ClusterID %in% c(13),2]='pDCs'
celltype[celltype$ClusterID %in% c(14),2]='Mast cells'

colnames(celltype) = c("ClusterID","celltype_main")
seurat.data@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  seurat.data@meta.data[which(seurat.data@active.ident == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(seurat.data@meta.data$celltype)
options(repr.plot.width = 6, repr.plot.height = 5)
DimPlot(seurat.data, reduction = "umap", group.by = "celltype", label = T)& NoAxes()
#将celltype设置为默认插槽
Idents(object = seurat.data) <- "celltype"
# 保存
saveRDS(seurat.data,file = "~/Spark/Step3.PBMC_annotation.rds")

#####FindAllMarkers方法#####
marker.Find <- FindAllMarkers(object = seurat.data, 
                              only.pos = T,
                              min.pct = 0.2,
                              logfc.threshold = 0.25,
                              test.use="wilcox")
head(marker.Find, n=5)
###寻找marker基因（每个cluster与其它所有总和的cluster做差异基因）
all.markers =marker.Find %>% dplyr::select(gene, everything())%>%subset(p_val<0.05)
top10= all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
###手动注释后读取
celltype <- read.delim("celltype.txt")
###将注释结果添加到Seurat对象的meta.data中
seurat.data@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){seurat.data@meta.data[which(seurat.data@meta.data$seurat_clusters == celltype$cluster[i]),'celltype'] <- celltype$cell_type[i]}
dev.off()