library(maftools)
library(stringr)
library(tidyr)
library(data.table)
library(dplyr)

#clinical.data = read.table("clinical.txt", header = T)
# 读取 vep 注释的 maf 文件 ----
vep_maf_nopairedBM = data.table::fread("nopaired_BM.maf")
vep_maf_pairedBM = data.table::fread("paired_BM.maf")
vep_maf_CRC= data.table::fread("CRC.maf")

#提取配对样本中的原发灶
vep_maf_pairedP <- vep_maf_pairedBM[grep("P$", vep_maf_pairedBM$Tumor_Sample_Barcode), ]
#提取配对样本中的转移灶
vep_maf_pairedBM1 <- vep_maf_pairedBM[grep("BM$", vep_maf_pairedBM$Tumor_Sample_Barcode), ]
#合并转移灶样本
vep_maf_BM<- rbind(vep_maf_pairedBM1, vep_maf_nopairedBM)

#分别传递vep_maf进行过滤
vep_maf<-vep_maf_CRC
vep_maf<-vep_maf_pairedP 

vep_maf<-vep_maf_pairedBM1

vep_maf<-vep_maf_BM

# Add VAF field
vep_maf$VAF <- vep_maf$t_alt_count/(vep_maf$t_alt_count+vep_maf$t_ref_count)

#注释突变位点
vep_maf$mutation <- paste(vep_maf$Hugo_Symbol, 
                          vep_maf$Start_Position, 
                          vep_maf$Chromosome,
                          vep_maf$End_Position, 
                          sep = "_")
#设置过滤条件
Variants = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
             "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", 
             "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")

vep_maf <- vep_maf[vep_maf$t_depth >=8, ]
vep_maf <- vep_maf[vep_maf$VAF >= 0.08, ]
vep_maf <- vep_maf[vep_maf$t_alt_count >=2, ]

#高突变处理流程
# 读取 MAF 文件
colnames(vep_maf)
maf <- read.maf(vep_maf)
# 提取原始突变表
maf_df <- maf@data
# 仅保留 SNVs
snv_df <- maf_df[maf_df$Variant_Type == "SNP", ]
# 按样本统计 SNV 数量
sample_snv_counts <- table(snv_df$Tumor_Sample_Barcode)
# 查看中位数与 IQR
med <- median(sample_snv_counts)
iqr_val <- IQR(sample_snv_counts)
threshold <- med + 1.5 * iqr_val
# 找出高突变样本
hyper_samples <- names(sample_snv_counts[sample_snv_counts > threshold])
split_hyper_mutants <- function(maf_df) {
  repeat {
    # 统计突变数
    snv_counts <- table(maf_df$Tumor_Sample_Barcode)
    med <- median(snv_counts)
    iqr_val <- IQR(snv_counts)
    threshold <- med + 1.5 * iqr_val
    # 筛选高突变样本
    hyper_samples <- names(snv_counts[snv_counts > threshold])
    if (length(hyper_samples) == 0) break  # 停止条件
    for (sample in hyper_samples) {
      sample_mut <- maf_df[maf_df$Tumor_Sample_Barcode == sample, ]
      n <- nrow(sample_mut)
      # 打乱并平均分
      shuffled <- sample_mut[sample(1:n), ]
      half <- floor(n / 2)
      part1 <- shuffled[1:half, ]
      part2 <- shuffled[(half + 1):n, ]   
      # 修改名称
      part1$Tumor_Sample_Barcode <- paste0(sample, "_1")
      part2$Tumor_Sample_Barcode <- paste0(sample, "_2")
      # 替换原始样本
      maf_df <- maf_df[maf_df$Tumor_Sample_Barcode != sample, ]
      maf_df <- rbind(maf_df, part1, part2)
    }
  }
  return(maf_df)
}

##vep_maf_CRC1、m1为CRC组
vep_maf_CRC1<-vep_maf
m1<-read.maf(vep_maf_CRC1)

##vep_maf_pairedP1、m2为CRCBM_P组
vep_maf_pairedP1<-vep_maf 
m2<-read.maf(vep_maf_pairedP1)

##vep_maf_BM1、m3为CRCBM_BM组
vep_maf_BM1<-vep_maf 
m3<-read.maf(vep_maf_BM1)
unique(vep_maf_BM1$Tumor_Sample_Barcode)

##vep_maf_pairedBM1、m4为CRCBM_Paired组
vep_maf_pairedBM1<-vep_maf 
m4<-read.maf(vep_maf_pairedBM1)
unique(vep_maf_pairedBM1$Tumor_Sample_Barcode)

#输出maf文件汇总内容
laml <-m1
laml <-m2
laml <-m3
laml <-m4
#Shows sample summry.
getSampleSummary(laml)
#Shows gene summary.
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
#MAF文件中的224列都有哪些内容
getFields(laml)
#将MAF文件写出来
write.mafSummary(maf = laml, basename = 'CRC')
# 定义自定义颜色
col = c(
  "Missense_Mutation" = "#146fa7", 
  "Nonsense_Mutation" = "#91ce8c", 
  "Frame_Shift_Del" = "#b6aacb",   
  "Frame_Shift_Ins" = "#fe3730",    
  "Splice_Site" = "#ff792f",        
  "In_Frame_Del" = "#389c46",       
  "Translation_Start_Site" = "#fdbf6f",
  "In_Frame_Ins" = "#fdb16a",     
  "Multi_Hit" = "#88cee5"
)
plotmafSummary(
  maf = maf,
  rmOutlier = TRUE,
  dashboard = TRUE,
  titvRaw = FALSE,
  log_scale = FALSE,
  addStat = NULL,
  showBarcodes = TRUE,
  fs = 1,
  textSize = 0.5,
  color = col,
  titleSize = c(1, 0.6),
  titvColor = NULL,
  top = 20
)
VAF<-vep_maf[,c("Hugo_Symbol","VAF")]
colnames(VAF)
table(VAF$Hugo_Symbol)
table(VAF$VAF)
VAF_summary <- VAF %>%
  group_by(Hugo_Symbol) %>%
  summarise(VAF = if_else(n() > 1, mean(VAF), VAF[1])) %>%
  ungroup()
aml_genes_vaf<-VAF_summary
# 绘制瀑布图
oncoplot(maf = maf, 
         genes = gene_top,  # 选择要绘制的基因
         colors = col, 
         leftBarData = aml_genes_vaf,
         drawRowBar = TRUE, 
         fontSize = 0.6,
         barcodeSrt = 45,
         drawColBar = TRUE, 
         sampleOrder= sample_order,
         removeNonMutated =F, 
         showTumorSampleBarcodes = F,
         clinicalFeatures = c("Site","Paired_State"),
         annotationColor = list(Site=assign_Site,
                                Paired_State=assign_aired_State),
         sortByAnnotation = T,
         titleText = "Top30 Mutation Genes in 68 (100%) of 68 CRCBM samples")
