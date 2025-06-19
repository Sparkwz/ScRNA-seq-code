devtools::install_github('sunduanchen/Scissor')
library(Seurat)
library(preprocessCore)
library(scAB) #用于获取示例数据
library(Scissor)
library(Matrix)  
library(matrixStats)
library(MatrixGenerics)
library(EpiDISH) 
library(ComplexHeatmap)
library(sctransform)
library(gplots)
#########1. 数据加载#########
data("data_survival")

# 单细胞RNA-seq数据Seurat对象
dim(sc_dataset)
sc_dataset

# Bulk RNA-seq表达矩阵
head(bulk_dataset[,1:10])

#表型数据
head(phenotype)

#样本名称需保持一致
table(colnames(bulk_dataset) == row.names(phenotype))

#########2. 处理单细胞数据########
#降维聚类等预处理
sc_dataset <- run_seurat(sc_dataset,verbose = FALSE)
sc_dataset

UMAP_celltype <- DimPlot(sc_dataset, reduction ="umap",
                         group.by="celltype",label = T)

options(repr.plot.width = 6, repr.plot.height = 4.5)
UMAP_celltype
save(bulk_dataset, sc_dataset, phenotype,file = "~/Spark/scissor.RData")##下载至本地
#########本地运行normalize.quantiles########
load("C:\\Users\\56426\\Desktop\\Doing\\scissor.RData")
family = "cox"
tag = NULL
alpha = NULL
cutoff = 0.2
Save_file = "Scissor_inputs.RData"
common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
if (length(common) == 0) {    
  stop("There is no common genes between the given single-cell and bulk samples.")  
}
sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)
network <- as.matrix(sc_dataset@graphs$RNA_snn)


library(Seurat)
library(Matrix)
library(preprocessCore)

diag(network) <- 0
network[which(network != 0)] <- 1
dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ])
dataset1 <- normalize.quantiles(as.matrix(dataset0))
saveRDS(dataset1,file = "dataset1.rds")#上传至本地

#########定义并运行Scissor########
Scissor <- function (bulk_dataset, sc_dataset, phenotype, tag = NULL, alpha = NULL,                      
                     cutoff = 0.2, family = c("gaussian", "binomial", "cox"),                      
                     Save_file = "Scissor_inputs.RData", Load_file = NULL) 
{  
  library(Seurat)  
  library(Matrix)  
  library(preprocessCore)  
  if (is.null(Load_file)) {    
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))    
    if (length(common) == 0) {      
      stop("There is no common genes between the given single-cell and bulk samples.")    
    }    
    if (class(sc_dataset) == "Seurat") {      
      sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)      
      network <- as.matrix(sc_dataset@graphs$RNA_snn)    
    }    
    else {      
      sc_exprs <- as.matrix(sc_dataset)      
      Seurat_tmp <- CreateSeuratObject(sc_dataset)      
      Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst",                                          
                                         verbose = F)      
      Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)      
      Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp),                            
                           verbose = F)      
      Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10,                                   
                                  verbose = F)      
      network <- as.matrix(Seurat_tmp@graphs$RNA_snn)    
    }    
    diag(network) <- 0    
    network[which(network != 0)] <- 1    
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common,     ])
    dataset1 <- readRDS("~/Spark/dataset1.rds")    
    #dataset1 <- normalize.quantiles(as.matrix(dataset0))    
    rownames(dataset1) <- rownames(dataset0)    
    colnames(dataset1) <- colnames(dataset0)    
    Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]    
    Expression_cell <- dataset1[, (ncol(bulk_dataset) +                                      
                                     1):ncol(dataset1)]    
    X <- cor(Expression_bulk, Expression_cell)    
    quality_check <- quantile(X)    
    print("|**************************************************|")    
    print("Performing quality-check for the correlations")    
    print("The five-number summary of correlations:")    
    print(quality_check)    
    print("|**************************************************|")    
    if (quality_check[3] < 0.01) {      
      warning("The median correlation between the single-cell and bulk samples is relatively low.")    
    }    
    if (family == "binomial") {      
      Y <- as.numeric(phenotype)      
      z <- table(Y)      
      if (length(z) != length(tag)) {        
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")      
      }      
      else {        
        print(sprintf("Current phenotype contains %d %s and %d %s samples.",                       
                      z[1], tag[1], z[2], tag[2]))        
        print("Perform logistic regression on the given phenotypes:")      
      }    
    }    
    if (family == "gaussian") {      
      Y <- as.numeric(phenotype)      
      z <- table(Y)      
      if (length(z) != length(tag)) {        
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")      
      }      
      else {        
        tmp <- paste(z, tag)        
        print(paste0("Current phenotype contains ",                      
                     paste(tmp[1:(length(z) - 1)], collapse = ", "),                      
                     ", and ", tmp[length(z)], " samples."))        
        print("Perform linear regression on the given phenotypes:")      
      }    
    }    
    if (family == "cox") {      
      Y <- as.matrix(phenotype)      
      if (ncol(Y) != 2) {        
        stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")      
      }      
      else {        
        print("Perform cox regression on the given clinical outcomes:")      
      }    
    }    
    save(X, Y, network, Expression_bulk, Expression_cell,          
         file = Save_file)  
  }  
  else {   
    load(Load_file)  
  }  
  if (is.null(alpha)) {    
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,                
               0.6, 0.7, 0.8, 0.9)  }  
  for (i in 1:length(alpha)) {    
    set.seed(123)    
    fit0 <- APML1(X, Y, family = family, penalty = "Net",                   
                  alpha = alpha[i], Omega = network, nlambda = 100,                   
                  nfolds = min(10, nrow(X)))    
    fit1 <- APML1(X, Y, family = family, penalty = "Net",                   
                  alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)    
    if (family == "binomial") {      
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])    
    }    
    else {      Coefs <- as.numeric(fit1$Beta)    
    }    
    Cell1 <- colnames(X)[which(Coefs > 0)]    
    Cell2 <- colnames(X)[which(Coefs < 0)]    
    percentage <- (length(Cell1) + length(Cell2))/ncol(X)    
    print(sprintf("alpha = %s", alpha[i]))    
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.",                   
                  length(Cell1), length(Cell2)))    
    print(sprintf("The percentage of selected cell is: %s%%",                   
                  formatC(percentage * 100, format = "f", digits = 3)))    
    if (percentage < cutoff) 
    {      
      break    
    }    
    cat("\n")  
  }  
  print("|**************************************************|")  
  return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min,                           
                          family = family), Coefs = Coefs, Scissor_pos = Cell1,               
              Scissor_neg = Cell2))
}
infos1 <- Scissor(bulk_dataset, 
                  sc_dataset, 
                  phenotype, 
                  alpha = NULL, 
                  cutoff = 0.2,
                  #alpha = 0.05, 
                  family = "cox", #二分类变量选binomial，连续性变量选gaussian，生存表型选cox
                  Save_file = '~/Spark/Step2.Scissor_survival.RData')
names(infos1)
length(infos1$Scissor_pos)
infos1$Scissor_pos[1:4]
length(infos1$Scissor_neg)

#########可视化########
Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor+"
Scissor_select[infos1$Scissor_neg] <- "Scissor-"
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
UMAP_scissor <- DimPlot(sc_dataset, reduction = 'umap', 
                        group.by = 'scissor',
                        cols = c('grey','royalblue','indianred1'), 
                        pt.size = 0.001, order = c("Scissor+","Scissor-"))

options(repr.plot.width = 12, repr.plot.height = 4.5)
patchwork::wrap_plots(plots = list(UMAP_celltype,UMAP_scissor), ncol = 2)
table(sc_dataset$scissor,sc_dataset$celltype)
balloonplot(table(sc_dataset$scissor,sc_dataset$celltype))
