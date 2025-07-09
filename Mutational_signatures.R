library('NMF')
library('pheatmap')
library("barplot3d")
library("devtools")
library(BSgenome.Hsapiens.UCSC.hg38)
# 先构建三连核苷酸矩阵
maf.tnm = trinucleotideMatrix(maf = m1, 
                              #prefix = 'chr', 
                              #add = TRUE, 
                              ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
# 运行 NMF非负矩阵分解，并拟合
maf.sign = estimateSignatures(mat = maf.tnm, nTry = 6, pConstant=0.05)
# 确定最佳突变特征数量
plotCophenetic(res = maf.sign)
# 使用非负矩阵分解将矩阵分解为n签名
maf.sig = extractSignatures(mat = maf.tnm, n = 5,pConstant = 0.1)
# 与 COSMIC 的突变特征比较，计算余弦相似度
maf.v3.cosm = compareSignatures(nmfRes = maf.sig, sig_db = "SBS")
# 热图展示余弦相似度
pheatmap::pheatmap(mat = maf.v3.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
# 可视化突变特征
maftools::plotSignatures(nmfRes = maf.sig, title_size = 1.2, sig_db = "SBS")