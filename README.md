# 单细胞数据分析学习笔记
## 基础处理 Basic analysis
### 数据储存格式
- Seurat官网示例数据
- 10X Cellranger输出标准三件套（barcodes.tsv.gz/features.tsv.gz/matrix.mtx.gz)
- 文章中自带数据，如Rdata `read::read_rds`, RDS `load()` ,h5ad [link](https://mp.weixin.qq.com/s/7eUQ_yvJslizM3rKh0I1xA)
- GEO数据库各种格式数据集
- **数据处理准则**
 1. 标准三件套：第一步重命名文件，第二步`list.files`+`lappy`循环（`Read10X()`+`CreateSeuratObject()`）读取
 2. 非标准格式：`CreateSeuratObject()`读取
### 01 数据读入
- `Read10X()` 读入10X标准格式为稀疏矩阵
- `CreateSeuratObject()` 构建Seurat对象
- `merge()` 合并多个样本
### 02 数据质控
1. 双细胞预测及过滤（可选）
2. 质控指标(不同领域过滤标准可能不同）：`PercentageFeatureSet`
 - **nFeature_RNA**: 每个细胞中检测到的唯一基因数 200-2000
 - **nCount_RNA**: 每个细胞检测到的分析总数 <500
 - **percent.mt**：低质量或者死亡细胞含有很高的线粒体基因 >15%
3. 低质量细胞过滤 `subset`
4. 基因水平过滤
### 03 下游分析（三步走）
- 标准化 确保数据在不同样本之间具有可比性 `NormalizeData`
- 特征选择 识别细胞类型特征 `FindVariableFeatures`
- 归一化 减小样本间的技术差异 `ScaleData`
### 04 降维聚类去批次（两步不存在先后关系）
1. PCA线性降维： 减少数据维度、去除技术噪声、保留关键特征
2. UMAP非线性降维： 映射数据及可视化
3. 若存在批次效应，则去除批次效应： `Harmony`+`Dimplot`后再进行UMAP非线性降维
4. FindClusters亚群聚类 
### 05 分群注释
- 首先推荐COSG，通过`remotes::install_local("~/Rpackages/COSGR",upgrade = F,dependencies = T)`本地安装
- 外周组织可分为：上皮细胞、免疫细胞、基质细胞
- 前列腺癌与乳腺癌类似：即上皮细胞可分为管腔和基底细胞，管腔细胞一般为恶性肿瘤细胞，基底细胞为正常上皮细胞
- 上皮、髓质、T细胞、B细胞、肥大细胞和基质细胞谱系相似性低，一般UMAP图可见泾渭分明
- T细胞中，一般按照CD4T、CD8T、NKT和NK细胞依次排列，其中CD8T与NK/T细胞谱系更为接近
## 富集分析 Enrichment analysis
### 01 基于TOP差异基因的富集分析（GO/KEGG）
- 使用**TOP差异基因**作为输入数据，得到富集分析结果
### 02 基于秩次排序的富集分析（GSEA)
- 使用**任意基因集**进行GSEA富集分析，一般纳入差异分析的所有基因，按照Log FoldChang进行排序后进行富集分析
### 03 单样本/单细胞打分
- 评估单个样本/细胞中特定基因集的活性水平，可用于量化单个样本/细胞内的生物学特征或通路活性，使用基因表达量矩阵，结合已知基因集作用输入数据
 1. ssGSEA、GSVA和Z-score:常用于Bulk数据
 2. AUCell、VISION和AddModuleScore:常用于单细胞数据
## 转录因子分析 SCENIC
### 01 SCENIC三步骤
#### First step
- **GENIE3（随机森林)/GRNBoost (Gradient Boosting)** 推断转录因子与候选靶基因之间的共表达模块，每个模块包含一个转录因子及其靶基因，纯粹基于共表达
#### Second step 
- **RcisTatget** 分析每个共表达模块中的基因，以鉴定enriched motifs，仅保留TF motif富集的模块和targets，构建TF-targets网络，每个TF及其潜在的直接targets gene被称作一个调节因子（Regulons）；
#### Third step 
- **AUCell** 计算调节因子（Regulons）的活性，这将确定Regulon在哪些细胞中处于“打开”状态
### 02 SCENIC参考文件准备x3
#### 转录起始子信息：
- hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
- hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
#### 转录因子信息
- allTFs_hg38.txt
#### motif信息
- motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
### 03 SCENIC流程详细解释及报错解决方案
1. **输出表达矩阵为CSV格式**
```library(Seurat)
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
```
2. **pySCENIC**
- **创建文件目录，目录下包含（1.change.py 2.scenic_human.bash 3.SCENIC参考文件(转录起始子/转录因子/motif) 4.for.pyscenic.csv）**
- **change.py详细信息**
```
import os,sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("for.pyscenic.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs);
```
- **scenic_human.bash详细信息**
```
### Step1.运行change.py
python change.py
### Step2.设置路径
# 不同物种的数据库不一样，这里是人类是human 
dir=/home/iyun42/Scenic #改成自己的目录
tfs=$dir/allTFs_hg38.txt
feather=$dir/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
# 一定要保证上面的数据库文件完整无误哦 
input_loom=./sample.loom
ls $tfs  $feather  $tbl  
CORE=10

### Step3.运行pySCENIC
#3.1 grn
pyscenic grn \
--num_workers $CORE \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
$tfs #转录因子文件，human or mouse

#3.2 cistarget
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers $CORE  \
--mask_dropouts

#3.3 AUCell
pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers $CORE
```
- **报错解读**
```
报错1：TypeError: Must supply at least one delayed object
解决措施：**numpy版本降至1.23.5后运行** `pip install numpy==1.23.5` **python版本调整** `pip install dask-expr==0.5.3 distributed==2024.2.1`
报错2：AttributeError: 'Series' object has no attribute 'iteritems'
解决措施：**pandas版本调整** `pip install pandas==1.5.3`
```
- **运行结果解读**
- GRNBoost结果文件adj.sample.tsv: 基因调控网络邻接矩阵/描述了基因之间的调控关系，共三列，第1列为TF，第2列为target,第3列为importance,表示一个基因对的调控强度
- RcisTatget结果文件reg.csv: 包含基因调控网络的上下游关系信息，其中列出了每个基因及其预测的上游和下游调控基因
- AUCell打分结果文件out_SCENIC_loom：类似基于转录因子靶基因集的富集分析打分（细胞X转录因子表达矩阵）
3. **AUCell输出文件可视化**
