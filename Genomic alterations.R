#提取CRC突变特征(Variant_Type)
library(dplyr)
table(vep_maf$Variant_Type)
table(vep_maf$Tumor_Sample_Barcode)
Sum_Variant_Type<-vep_maf %>%
  group_by(Tumor_Sample_Barcode, Variant_Type) %>%
  summarise(count = n()) %>%
  spread(key = Variant_Type, value = count, fill = 0)
write.csv(Sum_Variant_Type,"Sum_Variant_Type.csv")
library(dplyr)
library(tidyr)
library(ggplot2)
# 读取数据并进行汇总
Sum_Variant_Type <- vep_maf %>%
  group_by(Tumor_Sample_Barcode, Variant_Type) %>%
  summarise(count = n(), .groups = "drop") %>%
  spread(key = Variant_Type, value = count, fill = 0)
# 计算每个样本的总变异数，只考虑数值型列
Sum_Variant_Type <- Sum_Variant_Type %>%
  mutate(total = rowSums(select(., -Tumor_Sample_Barcode)))
# 计算百分比
Sum_Variant_Type_percent <- Sum_Variant_Type %>%
  pivot_longer(cols = -c(Tumor_Sample_Barcode, total), names_to = "Variant_Type", values_to = "count") %>%
  mutate(percent = count / total * 100)
# 绘制堆积柱状图
ggplot(Sum_Variant_Type_percent, aes(x = Tumor_Sample_Barcode, y = percent, fill = Variant_Type)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample Barcode", y = "Percentage of Variant Type (%)", 
       title = "Percentage of Variant Types per Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") # 使用不同的颜色

# 按照 Tumor_Sample_Barcode 和 Variant_Classification 统计变异数
Sum_Variant_Classification <- vep_maf %>%
  group_by(Tumor_Sample_Barcode, Variant_Classification) %>%
  summarise(count = n(), .groups = "drop") %>%
  spread(key = Variant_Classification, value = count, fill = 0)

# 计算每个样本的总变异数
Sum_Variant_Classification <- Sum_Variant_Classification %>%
  mutate(total = rowSums(select(., -Tumor_Sample_Barcode), na.rm = TRUE))
# 转换为长格式，方便绘图
Sum_Variant_Classification_long <- Sum_Variant_Classification %>%
  pivot_longer(cols = -c(Tumor_Sample_Barcode, total), 
               names_to = "Variant_Classification", 
               values_to = "count") %>%
  mutate(percent = count / total * 100)
# 定义自定义排序顺序 (按照你提供的顺序，从下至上)
variant_order <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", 
                   "In_Frame_Del", "Frame_Shift_Ins", "Splice_Site", 
                   "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site")
# 定义颜色方案
color_palette <- c(
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
# 确保Variant_Classification按照自定义顺序排序
Sum_Variant_Classification_long <- Sum_Variant_Classification_long %>%
  mutate(Variant_Classification = factor(Variant_Classification, levels = rev(variant_order)))
# 绘制堆积柱状图
ggplot(Sum_Variant_Classification_long, aes(x = Tumor_Sample_Barcode, y = percent, fill = Variant_Classification)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample Barcode", y = "Percentage of Variant Classification (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = color_palette) +  # 应用自定义颜色方案
  theme(panel.background = element_blank(),  # 去除灰色底纹
        plot.background = element_blank())  # 去除整个图形背景

#提取CRC突变特征(VARIANT_CLASS)
library(dplyr)
table(vep_maf$VARIANT_CLASS)
table(vep_maf$Tumor_Sample_Barcode)
Sum_VARIANT_CLASS <- vep_maf %>%
  group_by(Tumor_Sample_Barcode, VARIANT_CLASS) %>%
  summarise(count = n()) %>%
  spread(key = VARIANT_CLASS, value = count, fill = 0)
library(dplyr)
library(tidyr)
library(ggplot2)
# 按照 Tumor_Sample_Barcode 和 VARIANT_CLASS 统计变异数
Sum_VARIANT_CLASS <- vep_maf %>%
  group_by(Tumor_Sample_Barcode, VARIANT_CLASS) %>%
  summarise(count = n(), .groups = "drop") %>%
  spread(key = VARIANT_CLASS, value = count, fill = 0)
# 计算每个样本的总变异数
Sum_VARIANT_CLASS <- Sum_VARIANT_CLASS %>%
  mutate(total = rowSums(select(., -Tumor_Sample_Barcode), na.rm = TRUE))
# 转换为长格式，方便绘图
Sum_VARIANT_CLASS_long <- Sum_VARIANT_CLASS %>%
  pivot_longer(cols = -c(Tumor_Sample_Barcode, total), 
               names_to = "VARIANT_CLASS", 
               values_to = "count") %>%
  mutate(percent = count / total * 100)
# 定义排序顺序 (按照要求，反转顺序)
variant_order <- c("SNV", "deletion", "insertion", "substitution")
# 确保VARIANT_CLASS按照自定义顺序排序
Sum_VARIANT_CLASS_long <- Sum_VARIANT_CLASS_long %>%
  mutate(VARIANT_CLASS = factor(VARIANT_CLASS, levels = rev(variant_order)))  # 反转排序顺序
# 定义配色方案
color_scheme <- c("SNV" = "#62b7e7", 
                  "deletion" = "#f8b02e", 
                  "insertion" = "#e2292a", 
                  "substitution" = "#6db628")
# 绘制堆积柱状图
ggplot(Sum_VARIANT_CLASS_long, aes(x = Tumor_Sample_Barcode, y = percent, fill = VARIANT_CLASS)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample Barcode", y = "Percentage of VARIANT_CLASS (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = color_scheme) +  # 使用自定义颜色
  theme(panel.background = element_blank(),  # 去除灰色底纹
        plot.background = element_blank())  # 去除整个图形背景

# 1. 过滤出 VARIANT_CLASS 为 SNV 的数据
snv_data <- subset(vep_maf, VARIANT_CLASS == "SNV")
# 2. 创建 SNV 变异形式（如 T>G, T>C）
snv_data$SNV_Form <- paste(snv_data$Reference_Allele, ">", snv_data$Tumor_Seq_Allele2, sep = "")
# 3. 根据映射规则合并 SNV 类型
snv_data$SNV_Type_Merged <- snv_data$SNV_Form
# 定义映射规则，合并相关的变异类型
snv_data$SNV_Type_Merged <- recode(snv_data$SNV_Type_Merged,
                                   "C>A" = "C>A/G>T", 
                                   "G>T" = "C>A/G>T", 
                                   "C>G" = "C>G/G>C", 
                                   "G>C" = "C>G/G>C", 
                                   "C>T" = "C>T/G>A", 
                                   "G>A" = "C>T/G>A", 
                                   "T>A" = "T>A/A>T", 
                                   "A>T" = "T>A/A>T", 
                                   "T>C" = "T>C/A>G", 
                                   "A>G" = "T>C/A>G", 
                                   "T>G" = "T>G/A>C", 
                                   "A>C" = "T>G/A>C", 
                                   .default = snv_data$SNV_Type_Merged)
# 按样本名称和合并后的 SNV 变异形式统计频数
snv_summary <- as.data.frame(table(snv_data$Tumor_Sample_Barcode, snv_data$SNV_Type_Merged))
colnames(snv_summary) <- c("Sample", "SNV_Type", "Count")
# 计算每个样本中不同 SNV 类型的比例
snv_summary <- snv_summary %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count))
# 配置颜色（根据您的图像中的配色方案）
snv_colors <- c("C>A/G>T" = "#1fbef0",  # 红色
                "C>G/G>C" = "#231f20",  # 蓝色
                "C>T/G>A" = "#e72725",  # 绿色
                "T>A/A>T" = "#b7b6b6",  # 黄色
                "T>C/A>G" = "#a1ce63",  # 紫色
                "T>G/A>C" = "#eec8c4")  # 深粉色
# 使用 ggplot2 绘制堆叠柱状图
ggplot(snv_summary, aes(x = Sample, y = Proportion, fill = SNV_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = snv_colors) +  # 使用指定的配色方案
  scale_y_continuous(labels = scales::percent) +  # 转换为百分比
  theme_minimal() +
  labs(x = "Samples",
       y = "Proportion of Total SNVs",
       fill = "SNV Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # 旋转样本标签
