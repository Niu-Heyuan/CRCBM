# 读取数据
data <- read.csv("filtered_CRC.csv", stringsAsFactors = FALSE)
# 计算每个基因突变位点在不同分组的突变数量，并计算总数
mutation_counts <- data %>%
  group_by(Hugo_Symbol, HGVSc, Group) %>%
  summarise(Mutation_Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = Mutation_Count, values_fill = 0) %>%
  rowwise() %>%
  mutate(Total_Count = sum(c_across(where(is.numeric)))) %>%
  ungroup()
# 按总数降序排列，并将顺序反转
mutation_counts <- mutation_counts %>%
  arrange(desc(Total_Count)) %>%
  mutate(Gene_Mutation = paste(Hugo_Symbol, HGVSc, sep = ":"),
         Gene_Mutation = factor(Gene_Mutation, levels = rev(paste(Hugo_Symbol, HGVSc, sep = ":"))))
# 转换为长格式数据用于绘图
plot_data <- mutation_counts %>%
  pivot_longer(cols = -c(Hugo_Symbol, HGVSc, Total_Count, Gene_Mutation),
               names_to = "Group", values_to = "Count") %>%
  filter(Count > 0)
# 根据突变总数设置Gene_Mutation顺序（从大到小）
plot_data <- plot_data %>%
  mutate(Gene_Mutation = factor(Gene_Mutation, levels = rev(levels(mutation_counts$Gene_Mutation))))

ggplot(plot_data, aes(x = Count, y = Gene_Mutation, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("CRC" = "#53867a", "CRCBM_P" = "#7d518f", "CRCBM_BM" = "#b04835")) +
  labs(x = "Mutation Count",
       y = "Mutation Site") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # 原y轴现在是x轴
    axis.title = element_text(size = 12)
  ) +
  coord_flip()
