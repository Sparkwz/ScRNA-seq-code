library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
devtools::install_github("navinlabcode/CellTrek")
remotes::install_local("~/Rpackages/CellTrek-main",upgrade = F,dependencies = T)
library(CellTrek)
library(viridis)
library(ConsensusClusterPlus)
library(SeuratObject)
library(magrittr)

#####Step1 加载单细胞和空转数据#####
brain_sc <- read_rds("~/Spark/brain_sc.rds")
brain_st_cortex <- read_rds("~/Spark/brain_st_cortex.rds")

## Rename the cells/spots with syntactically valid names
brain_st_cortex <- RenameCells(brain_st_cortex, new.names=make.names(Cells(brain_st_cortex)))
brain_sc <- RenameCells(brain_sc, new.names=make.names(Cells(brain_sc)))

##单细胞数据：
brain_sc

##空转数据
brain_st_cortex

#####Step2 使用CellTrek进行细胞映射#####
# We first co-embed ST and scRNA-seq datasets using traint
brain_traint <- CellTrek::traint(st_data=brain_st_cortex, 
                                 sc_data=brain_sc, 
                                 sc_assay='RNA', 
                                 cell_names='cell_type')

# 在这里，我们使用非线性插值（intp = T，intp_lin=F）方法来增强ST的spots
brain_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=brain_sc, sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                     dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek

#####Step3 CellTrek细胞共定位分析#####
glut_cell <- c('L2/3 IT', 'L4', 'L5 IT', 'L5 PT', 'NP', 'L6 IT', 'L6 CT',  'L6b')
names(glut_cell) <- make.names(glut_cell)
brain_celltrek_glut <- subset(brain_celltrek, subset=cell_type %in% glut_cell)

#然后使用scoloc进行共定位分析：
brain_celltrek_glut$cell_type <- factor(brain_celltrek_glut$cell_type,
                                        levels=glut_cell)

## 我们从图中提取最小生成树（MST）的结果。
brain_sgraph_KL_mst_cons <- brain_sgraph_KL$mst_cons
rownames(brain_sgraph_KL_mst_cons) <- colnames(brain_sgraph_KL_mst_cons) <- glut_cell[colnames(brain_sgraph_KL_mst_cons)]

## 然后，我们提取meta.data数据（包括细胞类型及其频率信息）。
brain_cell_class <- brain_celltrek@meta.data %>% dplyr::select(id=cell_type) %>% unique
brain_celltrek_count <- data.frame(freq = table(brain_celltrek$cell_type))
brain_cell_class_new <- merge(brain_cell_class, brain_celltrek_count, by.x ="id", by.y = "freq.Var1")

brain_sgraph_KL <- CellTrek::scoloc(brain_celltrek_glut, 
                                    col_cell='cell_type',
                                    use_method='KL', eps=1e-50)
CellTrek::scoloc_vis(brain_sgraph_KL_mst_cons,
                     meta_data=brain_cell_class_new )


#####Step4 目标细胞类型的空间加权共表达分析#####
brain_celltrek_l5 <- FindVariableFeatures(brain_celltrek_l5)
vst_df <- brain_celltrek_l5@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
nz_test <- apply(as.matrix(brain_celltrek_l5[['RNA']]@data), 1, function(x) mean(x!=0)*100)
hz_gene <- names(nz_test)[nz_test<20]
mt_gene <- grep('^Mt-', rownames(brain_celltrek_l5), value=T)
rp_gene <- grep('^Rpl|^Rps', rownames(brain_celltrek_l5), value=T)
vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(.,-vst.variance.standardized)
feature_temp <- vst_df$id[1:2000]

# 我们使用 scoexp 进行空间加权基因共表达分析。
brain_celltrek_l5_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=brain_celltrek_l5, 
                                                    assay='RNA', 
                                                    approach='cc', 
                                                    gene_select = feature_temp, 
                                                    sigm=140, 
                                                    avg_cor_min=.4, 
                                                    zero_cutoff=3, 
                                                    min_gen=40, max_gen=400)
