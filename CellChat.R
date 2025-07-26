library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
#InstallData("pbmc3k") 
library(dplyr)
library(aplot)
library(ggplotify)
library(readr)
devtools::install_github('immunogenomics/presto')
library(presto)
##############Step1. 构建cellchat对象##############
setwd("~/Spark/Cellchat")
## 1.1 读入数据
data("pbmc3k")  
pbmc3k= UpdateSeuratObject(object = pbmc3k)
## 1.2 构建cellchat对象
#pbmc3k里的seurat_annotations有一些NA注释，过滤掉
data.input = pbmc3k@assays$RNA@data  #NormalizeData()归一化后表达矩阵
meta.data =  pbmc3k@meta.data  #存储每个细胞的元数据，用于质控，注释及可视化
meta.data = meta.data[!is.na(meta.data$seurat_annotations),]
data.input = data.input[,row.names(meta.data)]
table(meta.data$seurat_annotations)
#设置因子水平
#levels输入所有的seurat_annotations细胞类型
meta.data$seurat_annotations = factor(meta.data$seurat_annotations,
                                      levels = c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", 
                                                 "FCGR3A+ Mono", "NK", "DC", "Platelet"))
class(meta.data$seurat_annotations)
### 1.3 Create a CellChat object Step1. 构建CellChat对象
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "seurat_annotations")

###可在cellchat对象的meta插槽中添加表型信息
# 添加meta.data信息
cellchat <- addMeta(cellchat, meta = meta.data)

# 设置默认的labels
levels(cellchat@idents) # show factor levels of the cell labels
#cellchat <- setIdent(cellchat, ident.use = "new.labels") 
#groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


##############Step2. 加载CellChatDB数据库##############
### 1.4 加载CellChat受配体数据库
CellChatDB <- CellChatDB.human #保证受配体数据库正确，此处为human
showDatabaseCategory(CellChatDB)
dev.off()
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB  # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

##############Step3. 对表达数据进行预处理##############
### 1.5 对表达数据进行预处理，用于细胞间通讯分析
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 2) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)        #01 识别高表达基因
cellchat <- identifyOverExpressedInteractions(cellchat) #02 识别高表达通路

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

##############Step4. 计算通讯概率，推断细胞通讯网络##############
cellchat <- computeCommunProb(cellchat,population.size = F)
#报错解决 #因为 RSubunitsV 变量是 NULL。这种情况发生在 cellchat@LR 槽为空时。换句话说，CellChat 无法在你提供的基因列表中识别任何相关的配体-受体对。这个槽会在你运行 identifyOverExpressedInteractions()时被填充
#1. Use the correct CellChatDB--human
#cellchat@data #2. check the input data matrix 
#cellchat@data.signaling #3. check the subset data matrix
#unique(cellchat@idents) #4. check the cell group information is correct
#cellchat <- computeCommunProbPathway(cellchat)
# 过滤掉通信数量少的细胞-细胞通信
cellchat <- filterCommunication(cellchat, min.cells = 10)

##############Step5. 提取预测的细胞通讯网络为data frame##############
### cellchat取子集(报错)
#barcode.use = sample(row.names(cellchat@meta),100)
#cellchat.subset = subsetCellChat(cellchat,cells.use = barcode.use)
#获取所有的配受体对以及其通讯概率
df.net <- subsetCommunication(cellchat)
head(df.net)

#以通路为单位提取通讯信息
df.pathway = subsetCommunication(cellchat,slot.name = "netP")
head(df.pathway)
# 对感兴趣的细胞提取受配体信息
# 这里的 source.use = c(1) 指的是Naive CD4 T，2和3分别对应Memory CD4 T和CD14+ Mono：
levels(cellchat@idents)
df.net.sub <- subsetCommunication(cellchat, sources.use = c(1), targets.use = c(2,3)) #source.use顺序根据levels(cellchat@idents)顺序
head(df.net.sub)

# 对感兴趣的通路提取受配体信息
df.net.sub <- subsetCommunication(cellchat, signaling = c("MIF", "MHC-I"))
head(df.net.sub)

##############Step6. 在信号通路水平推断细胞通讯##############
cellchat <- computeCommunProbPathway(cellchat)
head(cellchat@net)
head(cellchat@netP)

##############Step7. 计算加和的cell-cell通讯网络##############
#可视化加和的细胞间通讯网络。例如，使用circle plot显示任意两个细胞亚群之间的通讯次数或总通讯强度(权重)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
p1 = netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                      weight.scale = T, label.edge= F,
                      title.name = "Number of interactions")
p2 = netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                      weight.scale = T, label.edge= F,
                      title.name = "Interaction weights/strength")

#由于细胞间通讯网络的复杂性，我们可以对每个细胞亚群发出的信号进行检测。这里我们还控制参数edge.weight.max，以便我们可以比较不同网络之间的边权值：
#单独绘制每个细胞亚群发出的信号
mat <- cellchat@net$weight
par(mar = c(1, 1, 2, 1))  # 下，左，上，右的边距
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  p1 = netVisual_circle(mat2, vertex.weight = groupSize,
                        weight.scale = T, edge.weight.max = max(mat),
                        title.name = rownames(mat)[i])
}
saveRDS(cellchat,file = "~/Spark/Cellchat/Step1.CellCha_Res.rds")

##############Step8. 可视化##############
rm(list = ls())
####读入数据####
cellchat = read_rds(file = "~/Spark/Cellchat/Step1.CellCha_Res.rds")
pathways.show <- c("MIF") 

#### 01 Hierarchy plot 层次图####
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = c(1,2,3,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
levels(cellchat@idents)
levels(cellchat@idents)[c(1,2,3,4)]

#### 02 Circle plot show pathway 圆圈图####
par(mfrow=c(1,2))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",label.edge= T)
# Circle plot show L-R pairs 
# 用extractEnrichedLR函数提取指定pathways内的所有受配体信号值
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
LR.show
#"MIF_CD74_CXCR4"
# Vis
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = "MIF_CD74_CXCR4", layout = "circle")

#### 03 Chord diagram 和弦图####
par(mfrow = c(1,2), xpd=TRUE)
png("~/Spark/chord_plot.png", width = 3000, height = 3000, res = 300)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord",title.name = "Chord diagram  1")
dev.off()
# Chord diagram 2 show L-R pairs 显示配受体层面的和弦图，指定slot.name为net
netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(5:8), lab.cex = 0.5,title.name = "Chord diagram  2: show gene",slot.name = "net")
# Chord diagram 3 show pathway 显示通路层面的和弦图，指定slot.name为netP
netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(5:8), lab.cex = 0.5,slot.name = "netP",title.name = "Chord diagram  2: show pathway")

#### 04 Heatmap 热图####
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

##############Step9. 计算每个配体-受体对L-R pairs对整个信号通路的贡献，并可视化单个配体-受体对介导的细胞-细胞通信##############
netAnalysis_contribution(cellchat, signaling = pathways.show)
# Chord diagram
png("~/Spark/Cellchat/chord_plot2.png", width = 3000, height = 3000, res = 300)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = "MIF_CD74_CXCR4", layout = "chord")
dev.off()

##############Step10. 自动保存所有推断网络##############
#加载绘图函数
source("~/Spark/Cellchat/custom_seurat_functions.R")
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
# Vis
gg.list = list()
for (i in 1:length(pathways.show.all)) {
  ## 可视化1：hierarchy plot 可视化配受体对 netVisual.V2是在作者的netVisual基础上补充了out.dir参数
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual.V2(cellchat, signaling = pathways.show.all[i],out.format = c("png"),
               vertex.receiver = vertex.receiver, layout = "hierarchy",out.dir = "~/Spark/")#修改存储目录
  
  ## 可视化2： 柱状图可视化配受体对 
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  gg.list[[pathways.show.all[i]]] = gg
  #ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
options(repr.plot.width = 14, repr.plot.height = 8)
gg.plot = wrap_plots(gg.list,ncol = 5)
gg.plot
ggsave(gg.plot,filename = "~/Spark/Step2.L_R_contribution.pdf",width = 14, height = 10)

##############Step11. 观察多种配体受体或信号通路介导的细胞-细胞通信##############
### Bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
options(repr.plot.width = 4, repr.plot.height = 5)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:7), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:5), 
                 signaling = c("MIF","MHC-I"), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MIF","MHC-I"))
pairLR.use  = pairLR.use[c(1,3),,drop=F]
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(1:5), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)

##############Step12. 使用小提琴/气泡图绘制信号基因表达分布##############
options(repr.plot.width = 8, repr.plot.height = 4.5)
plotGeneExpression(cellchat, signaling = "MIF")
dev.off()
##############Step13. 识别细胞亚群的信号作用（例如主要的发送者，接收者）以及主要的贡献信号##############
##计算并可视化网络中心性得分##
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
cellchat@netP$pathways
pathways.show = "MHC-I"
netAnalysis_signalingRole_network(cellchat,
                                  signaling = pathways.show,
                                  width = 8, height = 2.5,
                                  font.size = 10)
# 识别对某些细胞亚群的输出或输入信号中贡献最大的信号
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1+ht2
##多分组代码
