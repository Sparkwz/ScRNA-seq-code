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
### 分群注释
- 外周组织可分为：上皮细胞、免疫细胞、基质细胞
- 前列腺癌与乳腺癌类似：即上皮细胞可分为管腔和基底细胞，管腔细胞一般为恶性肿瘤细胞，基底细胞为正常上皮细胞
- 上皮、髓质、T细胞、B细胞、肥大细胞和基质细胞谱系相似性低，一般UMAP图可见泾渭分明
- T细胞中，一般按照CD4T、CD8T、NKT和NK细胞依次排列，其中CD8T与NK/T细胞谱系更为接近
