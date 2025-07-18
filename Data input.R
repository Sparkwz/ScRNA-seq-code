#######10x数据读取######
#barcodes.tsv.gz
#features.tsv.gz
#matrix.mtx.gz
#####1.单样本10x数据读取#####
### 加载R包
library(Seurat)
### 数据读取
counts <- Read10X(data.dir = "文件夹名称")#文件夹名称即存放上述三个标准文件
### 创建Seurat对象
Seu_obj <- CreateSeuratObject(counts,
                              min.cells = 3,#gene至少在3个细胞中表达，否则去除该gene
                              min.features = 300,#细胞中至少有300gene存在，否则去除该细胞
                              project="name" #Seurat对象的Project名称，对于单个样本影响不大
                              )

#####2.多样本10x数据读取#####
###加载R包
library(Seurat)
###列出想要加载数据的文件夹
dir <- list.files()#将想要读取的10x数据文件夹放置同一目录下，并指定该目录为工作目录
###for循环读取数据并创建Seurat对象
Seu_obj_list <- list()
for(i in dir){
  counts <- Read10X(data.dir = i)  
  Project <- i  
  Seu_obj_list[[Project]] <- CreateSeuratObject(counts,min.cells = 3,min.features =1000, project=Project)
}

#######Txt格式读取######
#####1.单样本Txt格式读取#####
##加载包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(data.table)
library(limma)
###读取并处理数据
counts <- read.table("GSM4808947.txt.gz",sep="\t",header=T,check.names=F)
rt=as.matrix(counts)
rownames(rt)=rt[,1]
#指定基因名
exp=rt[,4:ncol(rt)]#表达矩阵是从第4列开始
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
Project <- strsplit("GSM4808947.txt.gz",split = "[.]")[[1]][1]#对GSM4808947.txt.gz进行按照.分割取GSM4808947
###创建Seurat对象
Seu_obj <- CreateSeuratObject(data,
                              min.cells = 3,
                              min.features = 300,
                              project=Project
                              )
#####2.多样本Txt格式读取#####
###加载包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(data.table)
library(limma)
###读取文件名称
dir <- list.files()#将文件放入同一文件夹下并指定该文件夹为工作路径
###for循环读取并创建seurat对象
Seu_obj_list <- list()
for(i in 1:length(dir)){  
  counts <- read.table(dir[i],sep="\t",header=T,check.names=F)  
  rt=as.matrix(counts)  
  rownames(rt)=rt[,1]###这里根据自己数据进行修改
  exp=rt[,4:ncol(rt)]###这里根据自己数据进行修改
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  Project <- strsplit(dir[i],split = "[.]")[[1]][1]
  Seu_obj_list[[i]] <- CreateSeuratObject(data,
                                          min.cells = 3,
                                          min.features = 300,
                                          project=Project)
}
#######CSV格式读取######
#CSV文件与txt文件的不同就在于csv文件是用“,”进行分割每列的，txt文件是用"\t"分隔每列
#####1.单样本CSV格式读取#####
###加载包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(data.table)
library(limma)
###读取并处理数据
counts <- read.table("GSM3348304.csv.gz",sep=",",header=T,check.names=F)
rt <- t(counts) #转置为行为基因名列为样本名
rt=as.matrix(rt)
colnames(rt) <- rt[1,]
exp=rt[2:nrow(rt),]#表达矩阵是从第2行开始
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
Project <- strsplit("GSM3348304.csv.gz",split = "[.]")[[1]][1]#对GSM3348304.csv.gz按照"."分割取GSM4808947
###创建Seurat对象
Seu_obj <- CreateSeuratObject(data,
                              min.cells = 3,
                              min.features = 300,
                              project=Project
                              )
#####2.多样本CSV格式读取#####
###加载包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(data.table)
library(limma)
###读取文件名称
dir <- list.files()#将文件放入同一文件夹下并指定该文件夹为工作路径
###for循环读取并创建seurat对象
Seu_obj_list <- list()
for(i in (dir)){
  counts <- read.table(i,sep=",",header=T,check.names=F)  
  rt <- t(counts)  
  rt=as.matrix(rt)  
  colnames(rt) <- rt[1,]  
  exp=rt[2:nrow(rt),]#表达矩阵是从第2行开始
  dimnames=list(rownames(exp),colnames(exp))  
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)  
  data=avereps(data)  
  Project <- strsplit(i,split = "[.]")[[1]][1]  
  Seu_obj_list[[i]] <- CreateSeuratObject(data,
                                          min.cells = 3,
                                          min.features = 300,
                                          project=Project)
}

#######RDS或Rdata格式读取######
#####1.RDS格式读取#####
Seu_obj=readRDS("文件名")
#####2.Rdata格式读取#####
Seu_obj=readRDS("文件名")