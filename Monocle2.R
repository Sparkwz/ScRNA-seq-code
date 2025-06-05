library(Seurat)
library(dplyr)
library(monocle)
##安装旧版monocle
#install.packages("Rawdata/monocle_2.26.0.tar.gz", repos = NULL, type = "source")
library(patchwork)
library(ggpubr)
library(pbmc3k.SeuratData)	#加载seurat数据集 

#InstallData("pbmc3k") 
data("pbmc3k")
#提取注释后seurat对象
seurat.data <- pbmc3k.final
seurat.data
options(repr.plot.width = 6, repr.plot.height = 4.5)
DimPlot(seurat.data, reduction = "umap", label=T) 
#提取纳入拟时序分析细胞类型
seurat.data$celltype = Idents(seurat.data)
#需要严格选择具有分化关系的细胞进行分析，不随意All in
sce = subset(seurat.data, celltype %in% c('CD14+ Mono','FCGR3A+ Mono'))
table(Idents(sce))
sce <- NormalizeData(sce) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

#1 标准Seurat转为Monocle2 cds数据
HSMM <- as.CellDataSet(sce)
HSMM

#2 如果是其他格式的数据，需要自行构建cds格式
{#2.1 表型数据（sample id/celltype)
sample_ann <- sce@meta.data  
head(sample_ann)

#2.2 基因信息
gene_ann <- data.frame(
  gene_short_name = rownames(sce@assays$RNA), 
  row.names = rownames(sce@assays$RNA) 
)
head(gene_ann)

#2.3 表达矩阵
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
ct=as.data.frame(sce@assays$RNA@counts)
ct[1:4,1:4]

#2.4 构建cds对象
HSMM2 <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
HSMM2
}
#########运行Monocle2########
##01 Estimate size factor##
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
saveRDS(HSMM,file = "~/Spark/Step1.HSMM.rds")

##02 数据质控##
HSMM <- detectGenes(HSMM, min_expr = 1)
print(head(fData(HSMM)))
#获取细胞表达量阈值的基因(在10个细胞以上有表达）
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

head(pData(HSMM))
length(expressed_genes)

##03 选择输入的基因用于降维聚类##
#（1）选择clusters差异表达基因（算法预测）；

#（2）选择离散程度高的基因（例如Seurat的高变基因）；

#（3）选择发育差异表达基因（需要结合背景知识）；

#（4）自定义发育marker基因（需要结合背景知识）

#1.选择clusters差异表达基因
##选择clusters差异表达基因
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~celltype",
                                      cores = 10)
# 挑选差异最显著的基因可视化
ordering_genes <- subset(diff_test_res, qval < 0.05)
ordering_genes=ordering_genes[order(ordering_genes$pval),]
head(ordering_genes[,c("gene_short_name", "pval", "qval")] )
cg=as.character(head(ordering_genes$gene_short_name))
ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))
length(ordering_genes)
HSMM <- setOrderingFilter(HSMM, ordering_genes[1:3000])
plot_ordering_genes(HSMM)
#红线表示单片基于这种关系对色散的期望。用于聚类的基因用黑点表示，其他带过滤的基因用灰点表示

#2.选择高变基因（推荐）
length(VariableFeatures(sce))
HSMM <- setOrderingFilter(HSMM, VariableFeatures(sce))
plot_ordering_genes(HSMM)

#3.基于表达量进行过滤
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
length(unsup_clustering_genes$gene_id)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

##04 降维 & 排序##
HSMM <- reduceDimension(HSMM,
                        max_components = 2,
                        num_dim = 20,
                        #residualModelFormulaStr = "~SampleID", #如果存在批次则指定批次
                        method = 'DDRTree') # DDRTree方式
HSMM <- orderCells(HSMM)
pData(HSMM) %>% head()

##05 构建拟时序轨迹（识别细胞轨迹的起点和终点）##
GM_state <- function(cds, starting_point, cluster){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds@phenoData@data[,cluster])[,starting_point]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
root_start = GM_state(cds = HSMM,starting_point = "CD14+ Mono",cluster = "celltype")
root_start
HSMM <- monocle::orderCells(HSMM, root_state = root_start)

##06 可视化##
#配色
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

#主题
if(T){
  text.size = 12
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
}
#配色
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

#主题
if(T){
  text.size = 12
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
}

#分别为根据 seurat cluster ，State ，Pseudotime 和 singleR注释后的cell type 着色。
a1 <- plot_cell_trajectory(HSMM, color_by = "celltype") + scale_color_manual(values = colour)
a2 <- plot_cell_trajectory(HSMM, color_by = "State") + scale_color_manual(values = colour)
a3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime") + ggsci::scale_color_gsea()
options(repr.plot.width = 14, repr.plot.height = 4.5)
wrap_plots(a1, a2, a3)

#分面展示
options(repr.plot.width = 4, repr.plot.height = 6)
plot_cell_trajectory(HSMM, color_by = "Pseudotime") +
  facet_wrap(~celltype, nrow = 2) + ggsci::scale_color_gsea()

###添加“树形图”
#plot_complex_cell_trajectory函数添加“树形图”
p1 <- plot_cell_trajectory(HSMM, x = 1, y = 2, color_by = "celltype") + 
  theme(legend.position='none',panel.border = element_blank()) + 
  scale_color_manual(values = colour) 
p2 <- plot_complex_cell_trajectory(HSMM, x = 1, y = 2,
                                   color_by = "celltype")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 

options(repr.plot.width = 4.5, repr.plot.height = 9)
wrap_plots(p1, p2, ncol = 1, heights = c(2,1.5))

input.data = data.frame(celltype = HSMM$celltype,
                        Pseudotime = HSMM$Pseudotime)

###添加箱式图
input.data = data.frame(celltype = HSMM$celltype,
                        Pseudotime = HSMM$Pseudotime)

options(repr.plot.width = 4.5, repr.plot.height = 4)
ggboxplot(data = input.data,
          fill = "celltype", 
          x = "celltype",  
          y = "Pseudotime")+mytheme

#特定基因表达
input.gene <- c("CD14","S100A8", "FCGR3A")
cds_subset <- HSMM[input.gene,]

options(repr.plot.width = 12, repr.plot.height = 3)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype", ncol = 3)
table(row.names(pData(HSMM)) == colnames(sce))
pData(HSMM)$FCGR3A = as.numeric(GetAssayData(object = sce, assay = "RNA",slot = "data")["FCGR3A",])
pData(HSMM)$CD14 = as.numeric(GetAssayData(object = sce, assay = "RNA",slot = "data")["CD14",])
p1=plot_cell_trajectory(HSMM, color_by = "FCGR3A") + ggsci::scale_color_gsea()
p2=plot_cell_trajectory(HSMM, color_by = "CD14") + ggsci::scale_color_gsea()
options(repr.plot.width = 12, repr.plot.height = 4.5)
wrap_plots(p1, p2, ncol = 2)


#####gene/通路的相关性折线图
input.data = data.frame(celltype = HSMM$celltype,
                        Pseudotime = HSMM$Pseudotime)
## 基因
input.data$FCGR3A = as.numeric(GetAssayData(object = sce, assay = "RNA",slot = "data")["FCGR3A",])
options(repr.plot.width = 6, repr.plot.height = 4)
## Vis
ggplot(input.data, aes(x = Pseudotime, y = FCGR3A)) +
  labs(x="Pseudotime",y = "Expression level (log2)")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 3,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "FCGR3A")+
  mytheme + theme(legend.position = "none")

## 通路
sce = AddModuleScore(object = sce, features = 
                       list(test = c("LGALS2", "FCGR3A", "S100A9", "CCL3", "CD14", "LYZ", "FOLR3", "ECHDC1", "S100A8")),
                     name = "test_pathway")
input.data$test_pathway = sce$test_pathway1
ggplot(input.data, aes(x = Pseudotime, y = test_pathway)) +
  labs(x="Pseudotime",y = "Pathway activity")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 0,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "Test pathway")+
  mytheme + theme(legend.position = "none")

##感兴趣基因热图
input.gene <- c("LGALS2", "FCGR3A", "S100A9", "CCL3", "CD14", "LYZ", "FOLR3", "ECHDC1", "S100A8")
options(repr.plot.width = 6, repr.plot.height = 4)
plot_pseudotime_heatmap(HSMM[input.gene,], 
                        num_cluster = 4, 
                        show_rownames = T, 
                        return_heatmap = T)

##########branch分析#########
#确定了分化起点后，Monocle可以模拟出每个细胞所处的分化时间，并寻找随着分化时间逐渐升高或降低的基因，即Beam分析
BEAM_res <- BEAM(HSMM, branch_point = 1, 
                 cores = 4, 
                 progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
dim(BEAM_res)
write.csv(BEAM_res,file = "~/Spark/Step1.BEAM_res.csv")
input.gene = row.names(subset(BEAM_res, qval < 0.05))
input.gene
length(input.gene)
options(repr.plot.width = 4.5, repr.plot.height = 6)
plot_genes_branched_heatmap(HSMM[input.gene,],
                            branch_point = 1,
                            num_clusters = 6,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)
#monocle2 拟时间分支点分析结果解读：https://www.jianshu.com/p/9995cd707002
