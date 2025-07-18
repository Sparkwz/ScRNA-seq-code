# 单细胞数据分析学习笔记
## 原始测序数据下载及Cellranger上游分析
### 数据下载 
1. SRA网站获取SRR_Acc_List.txt  [link](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA662018&o=acc_s%3Aa)
2. Linux环境安装Kingfisher
```
conda create -n kingfisher python=3.8
conda activate kingfisher
conda install -c bioconda kingfisher
```
3. Kingfisher多种形式下载(从ENA下载，意味着直接下载的是 FASTQ文件)
```
#下载整个Bioproject
kingfisher get -p PRJNA486534 -m ena-ascp ena-ftp prefetch aws-http 1>down_prjan486534.log 2>&1

##下载单个样本
kingfisher get -r SRR14615558 -m ena-ascp ena-ftp prefetch aws-http --download-threads 10  1>down.log 2>&1

##下载多个样本
kingfisher get --run-identifiers-list SRR_Acc_List.txt -m ena-ascp ena-ftp prefetch --download-threads 10 --check-md5sums 1>down_srr_list.log 2>&1
```
### 数据上游分析
1. 修改样本名称
```
#确认文件存在
ls *.fastq.gz       # 查看当前目录下的fastq.gz文件
cat SRR_ACC_List.txt  # 查看SRR列表文件内容
#执行重命名命令
while IFS= read -r i; do
    mv -v "${i}_1"*.fastq.gz "${i}_S1_L001_R1_001.fastq.gz"
    mv -v "${i}_2"*.fastq.gz "${i}_S1_L001_R2_001.fastq.gz"
done < SRR_ACC_List.txt
```
2. 下载cellranger
```
#下载cellranger
#https://www.10xgenomics.com/support/software/cell-ranger/downloads
curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1752260309&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=BdJZmx4rS~yKbvI5RmENSGIMNQeY5sJZSdhPVpIDH~vkDBA0130Ih~ph-Rq~VtNLd2nu41Aifh1v1RkBT-Snrb9BFWU~57oRf7jyn6vYJperQwUzKyQOSxCtvdEu3EfdKM~MMxuKCDfI-tsGOD0D4NpWamjtfSH5OZ3e6q2LB8n2FXTyuvzxj1ps~ueKpEfTK2UyVUIDIhnWQGZhfcLj-x29pYWlZ5T73N3RCTNmXEbqClg5PghYrNtFS2pCb9rS64kkSg5lA6wNXghINoRg~2ZQ82cxQGyUALOhHMX7zTlNPlXe4t7Kzai7VruClbHLXkI7f569b-WNxQOHj6B7dw__"
tar -zxvf cellranger-9.0.1.tar.gz
```
3. 下载人的参考基因组数据refdata文件
```
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -zxvf refdata-gex-GRCh38-2024-A.tar.gz
```
4. cellranger运行
```
# 临时添加到PATH
export PATH=/home/spark/cellranger-9.0.1/bin:$PATH
# 确认是否添加成功
which cellranger
# 创建新文件夹并移动SRR文件至文件夹
mkdir merged_SRR
#生成sample样本文件
ls *.fastq.gz | cut -d'_' -f1 | sort -u > samples.txt
# 设置参考路径
# 设定参考基因组
ref="/home/spark/refdata-gex-GRCh38-2024-A"
# 设定fastq文件路径
fastqs="/home/spark/merged_SRR"
# 运行cellranger
while read sample_id; do
    echo "Processing sample: $sample_id"
    
    nohup cellranger count \
        --id="$sample_id" \
        --transcriptome="$ref" \
        --fastqs="$fastqs" \
        --sample="$sample_id" \
        --create-bam=true \
        --nosecondary \
        --localcores=4 \
        --localmem=30 \
        > "${sample_id}.log" 2>&1 &
    
    echo "Submitted $sample_id, log: ${sample_id}.log"
done < samples.txt
# 查看运行进度
jobs -l
```

## 基础处理 Basic analysis
### 数据储存格式
- Seurat官网示例数据
- 10X Cellranger输出标准三件套（barcodes.tsv.gz/features.tsv.gz/matrix.mtx.gz)
- 文章中自带数据，如Rdata `read::read_rds`, RDS `load()` ,h5ad [link](https://mp.weixin.qq.com/s/7eUQ_yvJslizM3rKh0I1xA)
- GEO数据库各种格式数据集
- **数据处理准则**
 1. 标准三件套：第一步重命名文件，第二步`list.files`+`lappy`循环（`Read10X()`+`CreateSeuratObject()`）读取
 2. 非标准格式：`CreateSeuratObject()`读取
### 01 数据读入 Data input
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
### 06 认识Seurat对象数据结构
- **assays** 用于存储不同模态或不同分析阶段的数据
  - **RNA**
     - **layers**
        - **counts** 原始UMI计数矩阵（稀疏矩阵格式）`pbmc[["RNA"]]$counts`
        - **data** NormalizeData()归一化后表达矩阵 `pbmc[["RNA"]]$data`
          - **scale.data** ScaleData()标准化矩阵 `pbmc[["RNA"]]$scale.data[c(1:4),c(1:4)]`
        - **meta.data** FindVariableFeatures()储存高变基因
- **meta.data** 存储每个细胞的元数据，用于质控，注释及可视化
  - **nFeature_RNA**：每个细胞检测到的基因数
  - **nCount_RNA**：每个细胞总UMI数
  - **percent.mt**：线粒体基因百分比
  - **seurat_clusters**：细胞聚类的结果
  - **orig.ident**：样本来源标识（如不同样本或批次）
- **active.ident** 默认分组信息
- **reductions** 用于存储降维结果
  - **pca** 运行RunPCA后reductions存储PCA的降维结果
  - **umap** UMAP降维
  - **t-SNE** t-SNE降维

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
## 拟时序分析 Monocle2
### 基本概念及流程
- 用于分析单细胞转录组数据的方法，旨在推断细胞在发育或分化过程中的顺序
- 基于反向图嵌入学习单细胞轨迹，分为5个步骤
1. 构建monocle2数据对象（cds对象）
2. 数据过滤（QC）
3. 基于离散型基因进行降维聚类
4. 构建拟时序轨迹
5. 可视化
## 细胞通讯分析 CellChat
- 配体细胞表达配体，配体作用于受体和受体细胞，引起受体细胞的生物学变化
### 细胞通讯的三种通讯方式
1. 化学信号分子通讯：细胞通过**释放化学信号分子**与细胞进行通讯
2. 相邻细胞表面分子的黏着：细胞通过**直接接触相邻细胞表面分子**进行通讯
3. 细胞与细胞外基质的黏着：细胞通过**与细胞外基质相互作用**，通过细胞表面受体感知环境信号
### 细胞通讯所回答的三个层面的问题
1. **细胞层面**：哪些细胞亚群之间存在互作关系，及互作的数量和强度
2. **通路层面**：细胞亚群A和细胞亚群B存在哪些通讯相关的通路，而本质上，配受体对应的基因构成了通讯通路
3. **基因层面**：亚群A和亚群B存在哪些互作的配受体，其通讯概率多大
### Cellchat大概流程
1. 构建Cellchat对象
2. 载入CellchatDB数据库
3. 过滤表达数据
4. 基于受配体数据库计算细胞间的受配体通讯概率和P值
5. 可视化
- **Hierarchy plot 层次图**: USER should define vertex.receiver, which is a numeric vector giving the index of the cell groups as targets in the left part of hierarchy plot. This hierarchical plot consist of two components: the left portion shows autocrine and paracrine signaling to certain cell groups of interest (i.e, the defined vertex.receiver), and the right portion shows autocrine and paracrine signaling to the remaining cell groups in the dataset. Thus, hierarchy plot provides an informative and intuitive way to visualize autocrine and paracrine signaling communications between cell groups of interest. For example, when studying the cell-cell communication between fibroblasts and immune cells, USER can define vertex.receiver as all fibroblast cell groups.
- **Circle plot 圆圈图**：Visualization of cell-cell communication at different levels**: One can visualize the inferred communication network of signaling pathways using netVisual_aggregate, and visualize the inferred communication networks of individual L-R pairs associated with that signaling pathway using netVisual_individual.
- **Chord diagram 和弦图**: CellChat provides two functions netVisual_chord_cell and netVisual_chord_gene for visualizing cell-cell communication with different purposes and different levels. netVisual_chord_cell is used for visualizing the cell-cell communication between different cell groups (where each sector in the chord diagram is a cell group), and netVisual_chord_gene is used for visualizing the cell-cell communication mediated by mutiple ligand-receptors or signaling pathways (where each sector in the chord diagram is a ligand, receptor or signaling pathway.)
## 拷贝数变异分析 inferCNV
### inferCNV 三步骤
1. 基于单细胞数据构建inferCNV对象
2. `infercnv::run()`函数一键运行
3. 个性化分析和可视化
## 使用单细胞数据联合Bulk表型分析（ScAB/Scissor算法）
- 参考文献：PMID 37652986/PMID: 36368318
### Scissor
- step 1 输入数据包括1.单细胞数据（scATAC-seq/scRNA-seq),其次是2.对应组织的Bulk RNA-seq数据，以及3.表型数据
- step 2 计算每对细胞和bulk样本的Pearson相关系数构建相关系数矩阵，通过优化样本表型Y与相关矩阵S的回归模型
- step 3 由上述优化模型求解的非零系数β用于选择与目标表型相关的细胞亚群。其中Scissor+表示所选择的细胞与目标表型呈正相关，Scissor为负相关。
- step 4 可靠性检验+差异表达基因分析+功能富集分析+motif分析
- **报错解决方案** `Error in preprocessCore::normalize.quantiles(Y) :   ERROR; return code from pthread_create() is 22`: 本地运行`normalize.quantiles`后传输中间文件
### CellTrek
scRNA-seq测序有丰富的细胞类型和基因表达量信息，但是缺乏基因表达的空间位置信息，无法了解基因在组织中的位置，会遗漏细胞间相互作用等关键生物学信息
- step 1 准备单细胞数据及空间转录组数据
- step 2 使用CellTrek进行细胞映射
- step 3 可视化分析
## 问题合辑
### 如何判断多样本或多队列单细胞数据整合分析时是否要去除批次效应？
- **`Batch.Quant`函数**计算Z-score评分，**输入数据**为前期经过质控并批次校正整合后的Seurat对象，其中包括原始矩阵对象和校正后矩阵的对象。若Z-score大于1.96（p=0.05)表明细胞与同数据集的其他细胞显著重叠，即批次效应显著
### 细胞亚群被划分为多少类是相对合适的，即Resolution参数选择多少合适？
- 在运行FingClusters函数时，选取不同的resolution，**以迭代地增加已识别的簇的数量**，并用FindAllMarkers分别提取不同resolution下各簇的特征基因，统计最少特征基因数量，当某簇特征基因极少时，反映该簇转录特征并不显著，可能是被强行分出的一类细胞簇，提示Resolution值过大。
### 细胞注释后如何验证注释结果的可靠性？
1. 展示细胞特异性标记基因的表达情况
2. 构建细胞特异性Signature,通过AddModuleScore打分评估细胞群的表达模式是否接近所注释细胞类型
3. 根据已注释细胞的分布情况来判断


