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
