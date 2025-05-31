# 单细胞数据分析学习笔记
## 基础处理
### 数据储存格式
- Seurat官网示例数据
- 10X Cellranger输出标准三件套（barcodes.tsv.gz/features.tsv.gz/matrix.mtx.gz)
- 文章中自带数据，如Rdata `read::read_rds`, RDS `load()` ,h5ad [link](https://mp.weixin.qq.com/s/7eUQ_yvJslizM3rKh0I1xA)
- GEO数据库各种格式数据集
- **数据处理准则**
 1. 标准三件套：第一步重命名文件，第二步`list.files`+`lappy`循环（`Read10X()`+`CreateSeuratObject()`）读取
 2. 非标准格式：`CreateSeuratObject()`读取
### 数据读入关键步骤解释(10X 标准三件套)
- **`Read10X()`** 读入10X标准格式为稀疏矩阵
- **`CreateSeuratObject()`** 构建Seurat对象
