library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(future)
library(stringr)
library(enrichplot)
library(scRNAtoolVis)
# check the current active plan
plan()
# change the current plan to access parallelization
plan("multisession", workers =4)
plan()

#设置可用的内存
options(future.globals.maxSize = 10 * 1024^3)
#加载绘图函数
source("~/Spark/R/custom_seurat_functions.R")
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

###读入数据
seurat.data = read_rds(file = "~/Spark/Step3.PBMC_annotation.rds")
seurat.data

#UMAP降维可视化
options(repr.plot.width = 6, repr.plot.height = 5)
p1 = DimPlot(seurat.data, 
             reduction = "umap",
             group.by = "celltype",
             label = T);p1

Idents(seurat.data) = "celltype"
ct.marker <- FindAllMarkers(object = seurat.data, 
                            only.pos = T,#只纳入正表达基因
                            logfc.threshold = 0.25,#保留logFC>0.25的基因
                            test.use="wilcox")
head(ct.marker)
table(ct.marker$cluster)
save(ct.marker,
     file = "~/Spark/Step6.celltype.markers.0.25.Rdata")
table(ct.marker$cluster)
#设定差异基因阈值
deg.df = filter(ct.marker, avg_log2FC>=0.8, p_val_adj<0.05)
table(deg.df$cluster)
#产生差异基因序列
gene.list = lapply(split(deg.df, deg.df$cluster), function(x){x$gene})
str(gene.list)

############批量GO分析############
#GOBP运行比较慢(Myenrich为自定义函数)
SigGOBP <- lapply(gene.list, Myenrich, 
                  category = "go", geneid = "SYMBOL")
##可视化
SigGOBPplot <- lapply(names(SigGOBP), function(z)barplot(SigGOBP[[z]],
                                                         font.size = 12,
                                                         showCategory = 5) +
                        ggtitle(z))
options(repr.plot.width = 22, repr.plot.height = 8)
gobpcol <- which(sapply(SigGOBPplot, function(z)nrow(z$data)) != 0)
if (length(gobpcol) != 0){
  SigGOBPplot <- SigGOBPplot[gobpcol]
  p.gobp = wrap_plots(plotlist = SigGOBPplot, ncol = 5) 
}
p.gobp

############批量KEGG分析############
Sigkegg <- lapply(gene.list, Myenrich, 
                  category = "kegg", geneid = "SYMBOL")
#柱状图可视化
Sigkeggplot <- lapply(names(Sigkegg), function(z)barplot(Sigkegg[[z]],
                                                         font.size = 12,
                                                         showCategory = 5) + ggtitle(z))

keppcol <- which(sapply(Sigkeggplot, function(z)nrow(z$data)) != 0)
if (length(keppcol) != 0){
  Sigkeggplot <- Sigkeggplot[keppcol]
  p.kepp = cowplot::plot_grid(plotlist = Sigkeggplot, ncol = 5) 
}
options(repr.plot.width = 22, repr.plot.height = 8)
p.kepp

##整合式气泡图代替柱状图##
deg.df = merge(deg.df,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(deg.df$ENTREZID, deg.df$cluster)
head(gcSample$`B plasma`) # entrez id , compareCluster
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
options(repr.plot.width = 12, repr.plot.height = 6)
p.kegg = dotplot(xx,font.size = 10, showCategory = 3, label_format = 100) + 
  mytheme + ggtitle("KEGG")
p.kegg

############GSEA############
#GSEA需要所有的基因，然后基于foldchange排序
##计算所有的差异基因
Idents(seurat.data) = "celltype"
#纳入所有基因，only.pos只保留正表达选F，Logfc选最小
ct.marker.all <- FindAllMarkers(object = seurat.data, 
                                only.pos = F,
                                min.pct = 0.01,
                                logfc.threshold = 0.01,
                                test.use="wilcox")
table(ct.marker.all$cluster)
str(ct.marker.all)
save(ct.marker.all,
     file = "~/Spark/Step6.celltype.markers.for_GSEA.Rdata")
#load(file = "~/Spark/Step6.celltype.markers.for_GSEA.Rdata")
##根据avg_log2FC排序
sce.list = lapply(split(ct.marker.all,ct.marker.all$cluster), function(x){
  tpm = x
  geneList = tpm$avg_log2FC
  names(geneList) = tpm$gene
  geneList = sort(geneList, decreasing = TRUE)
  return(geneList)
})
str(sce.list[1:5])
##基于Hallmark基因集的GSEA富集分析
gene.sets <- clusterProfiler::read.gmt("~/Spark/R/h.all.v7.2.symbols.gmt")
gene.sets$term = str_replace_all(string = gene.sets$term,
                                 pattern = "HALLMARK_",replacement = "")
unique(gene.sets$term)
##批量GSEA分析
GSEA_analy <- lapply(sce.list, function(x){
  clusterProfiler::GSEA(x,TERM2GENE = gene.sets,
                        pvalueCutoff = 0.05,
                        eps = 1e-100)
})
##批量可视化
GSEAplot <- lapply(names(GSEA_analy), 
                   function(z)enrichplot::dotplot(GSEA_analy[[z]], 
                                                  showCategory = 3,
                                                  x = "NES") +
                     ggtitle(z) + 
                     theme(plot.title = element_text(color="black",hjust = 0.5)))

options(repr.plot.width = 30, repr.plot.height = 8)
p.gsea = wrap_plots(plotlist=GSEAplot, ncol= 5)
p.gsea

##单个细胞前20通路可视化
options(repr.plot.width = 6, repr.plot.height = 4.5)
enrichplot::dotplot(GSEA_analy[["T cycling"]], 
                    showCategory = 20,
                    x = "NES") +
  ggtitle("GSEA (T cycling)") + 
  theme(plot.title = element_text(color="black",hjust = 0.5))

##单个细胞单通路可视化
options(repr.plot.width = 5, repr.plot.height = 5)
gseaplot2(GSEA_analy[["T cycling"]], 
          geneSetID = 'G2M_CHECKPOINT',
          pvalue_table=T,
          title = "G2M_CHECKPOINT")

##汇总气泡图
gseaTab = NULL
for (x in names(GSEA_analy)) {
  gsea.df = as.data.frame(GSEA_analy[[x]])
  gseaTab <- rbind.data.frame(gseaTab,
                              data.frame(term = gsea.df$ID,
                                         NES = gsea.df$NES,
                                         FDR = gsea.df$p.adjust,
                                         celltype = x,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
  
}
gseaTab$NES = ifelse(gseaTab$NES>2.5,2.5,
                     ifelse(gseaTab$NES< -2.5,-2.5,gseaTab$NES))
options(repr.plot.width = 11.5, repr.plot.height = 5.5)
p.gsea = ggplot(gseaTab, aes(x=celltype,y=term)) +
  geom_point(shape=21,aes(size=-log10(FDR),fill=NES),position =position_dodge(0)) +
  labs(title = "Hallmark GSEA analysis",y=NULL,x=NULL) + 
  scale_color_manual(values = c("white","black"))+
  scale_fill_gradient2(
    low = "#0072B5FF",
    mid = "white",
    high = "#FF7F00",
    midpoint = 0
  )+
  theme_bw() + mytheme + coord_flip()+
  scale_size(name = "-log10 (FDR)",range = c(1,7))+
  theme(legend.key.size = unit(0.4, "cm"))
p.gsea
head(gseaTab)

#####AUCell#####
#assay.names与assay需要统一
seurat.data <- sc.Pathway.Seurat(obj = seurat.data, 
                                 method = "AUCell", #可选方法："AUCell", "VISION", "ssGSEA","gsva"，单细胞推荐使用AUCell
                                 ncores = 4,
                                 assay.names = "pathway",
                                 geneList = "~/Spark/R/h.all.v7.2.symbols.gmt")
seurat.data
##小提琴图指定通路可视化
options(repr.plot.width = 6, repr.plot.height = 5)
VlnPlot_2(seurat.data, 
          features = c("HALLMARK-G2M-CHECKPOINT"),
          assay = "pathway", 
          y.lab = "Pathway enrichemnt score (AUCell)")

##featureplot通路UMAP图可视化
options(repr.plot.width = 5, repr.plot.height = 4.5)
FeaturePlot(seurat.data, features = c("HALLMARK-P53-PATHWAY"))&
  scale_color_distiller(palette = 'RdBu')

##气泡图可视化
options(repr.plot.width = 7, repr.plot.height = 2.5)
check_terms = c('HALLMARK-G2M-CHECKPOINT','HALLMARK-INTERFERON-GAMMA-RESPONSE',
                'HALLMARK-P53-PATHWAY',"HALLMARK-KRAS-SIGNALING-UP","HALLMARK-GLYCOLYSIS",
                "HALLMARK-FATTY-ACID-METABOLISM")
p.dot = DotPlot_2(object = seurat.data,
                  assay = "pathway",
                  text.size = 10,
                  Combine = F, dot.range.max = 3,
                  dot.range.min = 0, label.size = 3,
                  features = check_terms, legend.key.size = 0.4)&
  scale_color_distiller(palette = 'RdYlBu')&coord_flip();p.dot

##热图
#input.path = c('HALLMARK-G2M-CHECKPOINT','HALLMARK-INTERFERON-GAMMA-RESPONSE',
#                'HALLMARK-P53-PATHWAY',"HALLMARK-KRAS-SIGNALING-UP","HALLMARK-GLYCOLYSIS",
#               "HALLMARK-FATTY-ACID-METABOLISM")
input.path = row.names(seurat.data@assays$pathway@data)
#pdf(file = "./Outplot/Step6.AUCell_Pathway_heatmap.pdf", height = 4.3, width = 5)
options(repr.plot.width = 7.5, repr.plot.height = 8)
AverageHeatmap(object = seurat.data, 
               assays = "pathway",
               htRange = c(-2,0,2),
               markerGene = input.path,
               cluster_rows = T,
               cluster_columns = F,
               show_column_dend = F,
               show_row_dend = F,
               row_title = NULL)
#dev.off()
##保存数据
write_rds(seurat.data, file = "~/Spark/Step6.PBMC_annotation_AUCell_Hallmark.rds")
