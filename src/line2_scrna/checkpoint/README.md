# checkpoint/ 占位说明
#
# 跨语言传文件（Python 出 .h5ad 给 R 读，或反过来）后，
# 必须先跑本目录的检查点脚本，四项任一不符即停：基因名、矩阵方向、细胞名、embedding。
#
# 命名：checkpoint_<步骤>_<内容>.<R|py>，例如 checkpoint_02_qc_to_seurat.R。
#
# 当前为空：等线内真实桥接步骤出现时再写，不预造。
