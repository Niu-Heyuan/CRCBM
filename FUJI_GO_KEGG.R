library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(RColorBrewer)
#导入需要富集分析的列表列名称为Gene
diff <- read.csv("monogene.txt",header = T)
diff <- as.character(diff$Gene)
columns(org.Hs.eg.db)
diff_entrez <- bitr(diff,
                    fromType = "SYMBOL",#现有的ID类型
                    toType = "ENTREZID",#需转换的ID类型
                    OrgDb = "org.Hs.eg.db")
head(diff_entrez) #提示少量无法映射，不同数据库间ID转换存在少量缺失是正常现象
KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 10,
                        maxGSSize = 500)
KEGG_diff <- setReadable(KEGG_diff,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID")
View(KEGG_diff@result)
#计算Rich Factor（富集因子）：
KEGG_diff2 <- mutate(KEGG_diff,
                     RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
#计算Fold Enrichment（富集倍数）：
KEGG_diff2 <- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]
KEGG_result <- KEGG_diff2@result
#保存富集结果到本地：
save(KEGG_diff2, KEGG_result, file = c("KEGG_diff.Rdata"))
write.csv(KEGG_result, file = c('KEGG_CRCBM_Driver.csv'))
#富集气泡图绘制：
dotplot(
  KEGG_diff2,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 20,
  font.size = 12,
  title = "Top 20 of Pathway Enrichment",
  label_format = 30
)
#GO富集分析
library(dplyr)
library(stringr)
library(org.Hs.eg.db) #智人注释包
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
diff <- read.csv("laml_300.txt",header = T)
diff <- as.character(diff$Gene)
diff_entrez <- bitr(diff,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")
head(diff_entrez)
GO_all_diff <- enrichGO(gene = diff_entrez$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        minGSSize = 5,
                        maxGSSize = 500,
                        readable = T)
