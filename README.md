# 单细胞数据分析学习笔记
## 基础处理
### 数据储存格式
- Seurat官网示例数据
- 10X Cellranger输出标准三件套（barcodes.tsv.gz/features.tsv.gz/matrix.mtx.gz)
- 文章中自带数据，如Rdata, RDS,h5ad
- GEO数据库各种格式数据集
### 数据读入关键步骤解释(10X 标准三件套)
`Read10X()` 读入10X标准格式为稀疏举证
`CreateSeuratObject()` 构建Seurat对象
