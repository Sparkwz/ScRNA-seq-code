rm(list = ls())
library(infercnv) #安装inferCNV需要先install.packages("rjags")
library(Seurat)
library(dplyr)
library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(scales)
library(pheatmap)

outdir = "~/Spark/inferCNV/"
dir.create(outdir)
### 1. 数据读入
seurat.data = read_rds("~/Spark/test_for_inverCNV.rds")
table(seurat.data$celltype)
seurat.data

### 2. 数据准备(三个文件)
## 2.1 细胞表型数据
#Create slim cell annotation file (malignant vs non-malignant) from Seurat object 
inferCNV.anno = data.frame(cell.id = rownames(seurat.data@meta.data),
                           group = seurat.data$celltype)
table(inferCNV.anno$group)

## 2.2 count数据
count.data <- GetAssayData(seurat.data, slot='counts',assay='RNA') %>% as.data.frame()
count.data = count.data[,inferCNV.anno$cell.id]
dim(count.data)

## 2.3 基因组文件(外用)
geneInfor = read.table("~/Spark/gene_pos.txt")
comm.gene = intersect(geneInfor$V1, rownames(count.data) )
head(comm.gene)

# Write table为txt格式
write.table(inferCNV.anno,file = paste0(outdir,'groupFiles.txt'),
            sep = '\t',quote = F,col.names = F,row.names = F)

write.table(count.data,file = paste0(outdir,'expFile.txt'),
            sep = '\t',quote = F)

write.table(geneInfor,
            file = paste0(outdir,"geneFile.txt"),
            row.names = F, col.names = F, quote=F, sep="\t")

## 2.4 构建inferCNV对象
expFile = paste0(outdir,'expFile.txt')
groupFiles=paste0(outdir,'groupFiles.txt')
geneFile=paste0(outdir,'geneFile.txt')

group.data = read.table(groupFiles, sep="\t")
dim(group.data)
table(group.data$V2)
##查看分组信息
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = expFile,
                                    annotations_file = groupFiles,
                                    gene_order_file = geneFile,
                                    ref_group_names = c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"),
                                    delim="\t")  ## 这个取决于自己的分组信息里面的
saveRDS(infercnv_obj,file = "~/Spark/inferCNV/Step1.infercnv_obj_inputfor_run.rds")

## 2.5 运行inferCNV
setwd("~/Spark/inferCNV/") #设置结果输出目录
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='Infercnv-results/', 
                             cluster_by_groups = T, # 默认False; 先区分细胞来源，再做层次聚类
                             analysis_mode = 'subclusters',
                             tumor_subcluster_partition_method = 'random_trees',
                             denoise=TRUE,
                             HMM=TRUE,
                             HMM_type = 'i6',
                             write_expr_matrix = T, # 默认False  
                             num_threads=parallel::detectCores())
saveRDS(infercnv_obj, file = "~/Spark/inferCNV/Step1.inferCNV_res.rds")

###CNV heatmap可视化
rm(list = ls())
dat <- read.table("~/Spark/inferCNV/Infercnv-results/infercnv.observations.txt",header=T,row.names=1,sep=" ",stringsAsFactors=F)
colnames(dat) <- str_replace(colnames(dat),"^X","")
ann <- read.table("~/Spark/inferCNV/test.cell_groupings",header = T,sep = "\t",stringsAsFactors = F)
ann$cell_group_name <- str_replace(ann$cell_group_name,"malignant_93.malignant_93.","")
ann$cell_group_name <- str_replace(ann$cell_group_name,'malignant_97.malignant_97.','')
ann$cell_group_name <- str_replace(ann$cell_group_name,'malignant_MGH36.malignant_MGH36.','')
ann$cell_group_name <- str_replace(ann$cell_group_name,'malignant_MGH53.malignant_MGH53.','')
rownames(ann) <- ann$cell
ann$cell <- NULL
ann$cell_group_name <- factor(ann$cell_group_name,levels = sort(unique(ann$cell_group_name)))
dat <- dat[,rownames(ann)]

gn <- rownames(dat)
geneFile <- read.table("~/Spark/gene_pos.txt",header = F,sep = "\t",stringsAsFactors = F)
dat <- dat[intersect(geneFile$V1,gn),]
sub_geneFile <- geneFile[match(rownames(dat),geneFile$V1),]

##热图注释
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), 
                                               labels = 1:22,labels_gp = gpar(cex = 1.5)))
len_c <- length(names(table(ann$cell_group_name)))
color_c <- RColorBrewer::brewer.pal(12, "Paired")[1:(len_c)]
names(color_c) <- names(table(ann$cell_group_name))
left_anno <- rowAnnotation(df = ann,col=list(cell_group_name=color_c),name='Groups',border = F)

#自定义哪些基因显示
cosmic=read.table("~/Spark/inferCNV/Census_allMon.small.tsv",header = F,sep = "\t",stringsAsFactors = F)

#发生了CNV的基因
pre_res <- read.table("~/Spark/inferCNV/Infercnv-results/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                      header = T,sep = "\t",stringsAsFactors = F)
pre_res <- pre_res %>% filter(state != 3) # State 3: 1x--neutral 
pre_res <- pre_res[!str_detect(pre_res$cell_group_name,"Microglia|Oligodendrocytes"),]
pre_res$cell_group_name <- str_replace(pre_res$cell_group_name,"malignant_93.malignant_93.","")
pre_res$cell_group_name <- str_replace(pre_res$cell_group_name,'malignant_97.malignant_97.','')
pre_res$cell_group_name <- str_replace(pre_res$cell_group_name,'malignant_MGH36.malignant_MGH36.','')
pre_res$cell_group_name <- str_replace(pre_res$cell_group_name,'malignant_MGH53.malignant_MGH53.','')

####同一基因 不同CNV事件 要区别开来
pre_res$cnv_type=""
chr_pq <- read.table("~/Spark/inferCNV/chr_pq.txt",header = F,sep = "\t",stringsAsFactors = F) 
colnames(chr_pq) <- c("chr","arm","cutoff")
chr_pq$chr <- paste("chr",chr_pq$chr,sep = "")

for (i in 1:nrow(pre_res)) {
  tmp_chr_pq <- chr_pq%>%filter(chr==pre_res[i,"chr"])
  if(pre_res[i,"start"] <= tmp_chr_pq[1,"cutoff"] & pre_res[i,"end"] <= tmp_chr_pq[1,"cutoff"]) {
    pre_res[i,"cnv_type"] <- paste(tmp_chr_pq[1,"chr"],tmp_chr_pq[1,"arm"],sep = "")
  }else if (pre_res[i,"start"] >= tmp_chr_pq[2,"cutoff"] & pre_res[i,"end"] >= tmp_chr_pq[2,"cutoff"]) {
    pre_res[i,"cnv_type"] <- paste(tmp_chr_pq[2,"chr"],tmp_chr_pq[2,"arm"],sep = "")
  }else {
    pre_res[i,"cnv_type"] <- paste(pre_res[i,"chr"],"p,",pre_res[i,"chr"],"q",sep = "")
  }
  if (pre_res[i,"state"] < 3) {
    pre_res[i,"cnv_type"] <- paste(pre_res[i,"cnv_type"],"_loss",sep = "")
  }
  if (pre_res[i,"state"] > 3) {
    pre_res[i,"cnv_type"] <- paste(pre_res[i,"cnv_type"],"_gain",sep = "")
  }
  if (str_detect(pre_res[i,"cnv_type"],",")) {
    tmp1 <- strsplit(pre_res[i,"cnv_type"],",")[[1]][1]
    tmp2 <- strsplit(strsplit(pre_res[i,"cnv_type"],",")[[1]][2],"_")[[1]][1]
    tmp3 <- strsplit(strsplit(pre_res[i,"cnv_type"],",")[[1]][2],"_")[[1]][2]
    pre_res[i,"cnv_type"] <- paste(tmp1,"_",tmp3,",",tmp2,"_",tmp3,sep = "")
    rm(list = c("tmp1","tmp2","tmp3"))
  }
}
pre_res$gene_cnv_type <- paste(pre_res$gene,pre_res$cnv_type,sep = "_")

##发生在极少数细胞的CNV事件——不考虑
tmpdf <- as.data.frame(table(ann$cell_group_name))
colnames(tmpdf) <- c("cell_group_name","cellnum")
pre_res <- pre_res%>%inner_join(tmpdf,by="cell_group_name")
tmpdf2 <- aggregate(pre_res$cellnum,by=list(pre_res$gene_cnv_type),FUN=sum)
colnames(tmpdf2) <- c("gene_cnv_type","num")
tmpdf2 <- tmpdf2%>%filter(num >= dim(dat)[2] * 0.05)
pre_res <- pre_res[pre_res$gene_cnv_type %in% tmpdf2$gene_cnv_type,]

#交集
key_gene <- intersect(rownames(dat),intersect(pre_res$gene,cosmic$V1))
ha <- columnAnnotation(foo=anno_mark(at=which(rownames(dat) %in% key_gene),
                                     labels = rownames(dat)[which(rownames(dat) %in% key_gene)],
                                     which = "column",side = "bottom"))
pdf("~/Spark/inferCNV/Step4.CNV_heatmap.clone.pdf",width = 15,height = 11)
ht_tp = Heatmap(t(dat),
                col = colorRamp2(c(0.6,1,1.4), c("#377EB8","#F0F0F0","#E41A1C")),
                cluster_rows = F,cluster_columns = F,
                show_column_names = F,show_row_names = F,
                column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")),
                column_gap = unit(2, "mm"),
                row_split = ann$cell_group_name,
                heatmap_legend_param = list(title = "Score",
                                            at=c(0.4,1,1.6),
                                            legend_height = unit(3,"cm")),
                top_annotation = top_anno,
                left_annotation = left_anno,
                bottom_annotation = ha,
                row_title = NULL,
                column_title = "Chromosome",
                column_title_side = "top",
                column_title_gp = gpar(fontsize = 20),
)

draw(ht_tp, heatmap_legend_side = "right")
dev.off()

options(repr.plot.width =15, repr.plot.height = 11)
ht_tp

##可视化2
rm(list = ls())
cell_group <- read.table("~/Spark/inferCNV/Infercnv-results/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",
                         header = T,sep = "\t",stringsAsFactors = F)
cell_group <- cell_group[!str_detect(cell_group$cell_group_name,"Microglia|Oligodendrocytes"),]
cell_group$cell_group_name <- str_replace(cell_group$cell_group_name,"malignant_93.malignant_93.","")
cell_group$cell_group_name <- str_replace(cell_group$cell_group_name,'malignant_97.malignant_97.','')
cell_group$cell_group_name <- str_replace(cell_group$cell_group_name,'malignant_MGH36.malignant_MGH36.','')
cell_group$cell_group_name <- str_replace(cell_group$cell_group_name,'malignant_MGH53.malignant_MGH53.','')

group_cellcount <- as.data.frame(table(cell_group$cell_group_name))
colnames(group_cellcount) <- c("cell_group_name","cellcount")
group_cellcount$cellratio <- group_cellcount$cellcount / sum(group_cellcount$cellcount)

group_cnvtype <- read.table("CNVgroup_and_CNVtype_in_sampleA.txt",header = T,sep = "\t",stringsAsFactors = F)
group_cnvtype <- group_cnvtype%>%inner_join(group_cellcount,by="cell_group_name")

alltype=c()
for (i in 1:22) {
  for (j in c("p","q")) {
    for (k in c("gain","loss")) {
      alltype <- append(alltype,paste("chr",i,j,"_",k,sep = ""))
    }
  }
}

cellpercent <- c()
for (i in alltype) {
  if(i %in% unique(group_cnvtype$cnv_type)){
    tmp <- sum(group_cnvtype[group_cnvtype$cnv_type == i,"cellratio"])
    cellpercent <- append(cellpercent,tmp)
  }else{
    cellpercent <- append(cellpercent,0)
  }
}
names(cellpercent) <- alltype
one.sample.stat <- as.data.frame(cellpercent)

#### 当有多个样本时，直接合并就可以了。这里我为了演示，人为多加了几列 ####
some.sample.stat=one.sample.stat
colnames(some.sample.stat)="sampleA"
some.sample.stat$sampleB=sample(some.sample.stat$sampleA,88,replace = T)
some.sample.stat$sampleC=sample(some.sample.stat$sampleA,88,replace = T)
some.sample.stat$sampleD=sample(some.sample.stat$sampleA,88,replace = T)
some.sample.stat$sampleE=sample(some.sample.stat$sampleA,88,replace = T)
some.sample.stat$sampleF=sample(some.sample.stat$sampleA,88,replace = T)
some.sample.stat$sampleG=sample(some.sample.stat$sampleA,88,replace = T)
some.sample.stat$sampleH=sample(some.sample.stat$sampleA,88,replace = T)

pheatmap(some.sample.stat,
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(brewer.pal(9,"PuRd"))(100),
         #rev(brewer.pal(n = 7, name ="RdYlBu"))
         #brewer.pal(9,"PuRd")
         #brewer.pal(9,"RdPu"),
         border_color = "grey",
         filename = "~/Spark/inferCNV/Step5.test.CNV.heatmap.pdf",width = 4,height = 11
)

options(repr.plot.width = 4, repr.plot.height = 13)
pheatmap(some.sample.stat,
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(brewer.pal(9,"PuRd"))(100),
         #rev(brewer.pal(n = 7, name ="RdYlBu"))
         #brewer.pal(9,"PuRd")
         #brewer.pal(9,"RdPu"),
         border_color = "grey"
)

#可视化3
### 1.读入数据
rm(list = ls())
infercnv_obj = readRDS("~/Spark/inferCNV/Infercnv-results/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices%>%unlist()
tumor_loc <- infercnv_obj@observation_grouped_cell_indices%>%unlist()

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[tumor_loc]),
  class=c(rep("normal",length(normal_loc)),rep("tumor",length(tumor_loc)))
)
head(anno.df)

gn <- rownames(expr)
geneFile <- read.table("~/Spark/gene_pos.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile) <- geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr <- expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)

expr[1:4,1:4]

#聚成5类
set.seed(1234)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(Cluster=paste0('C',kmeans.result$cluster))
kmeans_df$CB <- colnames(expr)
kmeans_df <- kmeans_df%>%inner_join(anno.df,by="CB") 
kmeans_df_s <- arrange(kmeans_df,Cluster)
rownames(kmeans_df_s) <- kmeans_df_s$CB
kmeans_df_s$CB <- NULL
kmeans_df_s$Cluster <- as.factor(kmeans_df_s$Cluster) 
head(kmeans_df_s)

#定义热图的注释，及配色
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:5]
names(color_v)=paste0('C',1:5)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("tumor"="red","normal" = "blue"),Cluster=color_v))

pdf("~/Spark/inferCNV/Step6.CNV-kmeans-heatmap.pdf",width = 20,height = 20)
##报错未解决
ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
             column_gap = unit(2, "mm"),
             heatmap_legend_param = list(title = "Scores",
                                         at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

options(repr.plot.width = 15, repr.plot.height = 10)
ht

cnvScore <- function(data){
  require(tidyverse)
  require(scales)
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    scales::rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- data.frame(ID=colnames(data),Score=colSums(data * data),row.names = NULL)
  return(cnv_score)
}
cnv_score <- cnvScore(expr)
data <- merge(kmeans_df_s,cnv_score,by.x=0,by.y=1)
head(data)

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

options(repr.plot.width = 5, repr.plot.height = 5)
ggplot(data,aes(Cluster,Score))+
  geom_boxplot(aes(fill=class),outlier.colour = 'grey30',outlier.size = 0.3)+
  labs(x=NULL,y='CNV score',fill=NULL)+
  ggsci::scale_fill_jco()+theme_bw()+mytheme
ggsave('~/Spark/inferCNV/Step6.CNV_score_type.pdf',width = 5,height = 3)

ggplot(data,aes(Cluster,Score))+
  geom_boxplot(aes(fill=Cluster),outlier.colour = 'grey30',outlier.size = 0.3)+
  labs(x=NULL,y='CNV score',fill=NULL)+
  ggsci::scale_fill_jco()+theme_bw()+mytheme
ggsave('~/Spark/inferCNV/Step6.CNV_score_cluster.pdf',width = 5,height = 3)
