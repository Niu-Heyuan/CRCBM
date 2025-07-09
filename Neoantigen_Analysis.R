config = read.table("neo.config",header = F)
colnames(config) = c("Tumor_Sample_Barcode","group")
head(config)
# 临床信息
clincal.all = read.table("clinical_all.txt",header = T,sep = "\t")
clincal.all = as.data.frame(clincal.all)
head(clincal.all)
# 合并新抗原结果
merge.neo  = data.frame()
for (i in config$Tumor_Sample_Barcode) {
  # i="B16_BM"
  neo.raw = read.table(paste0("./pvacseq/",i,".filtered.fc1.tsv"),header = T,sep = "\t")
  # 判断是否为空，否添加 id 
  if (nrow(neo.raw) != 0 ) {
    neo.raw$Tumor_Sample_Barcode = i
    # 合并所有结果
    merge.neo = rbind(merge.neo,neo.raw)
  }
}
# 对于同一个样本同一个基因有多个候选新抗原，取 Transcript.Length 较长的候选新抗原，
rmdup.neo = merge.neo %>% group_by(Gene.Name, Tumor_Sample_Barcode) %>%
  slice(which.max(Transcript.Length)) %>%
  ungroup()
# 长变宽
neo.df = rmdup.neo[,c("Gene.Name","Tumor_Sample_Barcode","Best.MT.IC50.Score")]
neo.mtx = spread(neo.df, key = "Tumor_Sample_Barcode", value = "Best.MT.IC50.Score")
neo.mtx = as.data.frame(neo.mtx)
rownames(neo.mtx) = neo.mtx[,1]
neo.mtx = neo.mtx[,-1]
# 添加上没检测到新抗原的样本
neo.mtx$XQC_P = NA
neo.mtx$LQ_P = NA
neo.mtx$B23_P = NA
# 样本排序：先分组 CRC  CRCBM_P CRCBM_BM，然后组内按新抗原数排序
neo.sum = data.frame(Tumor_Sample_Barcode = names(colSums(neo.mtx>0,na.rm = T)),
                     counts = colSums(neo.mtx>0,na.rm = T))
sample.df = merge(config,neo.sum,by="Tumor_Sample_Barcode")
sample.df = sample.df[order(sample.df$group,sample.df$counts,decreasing = T),] 
sample.df = rbind(sample.df[sample.df$group=="CRC",],sample.df[sample.df$group!="CRC",])
sample.df
write.table(neo.mtx[,sample.df$Tumor_Sample_Barcode],file = "neo_ic50.txt",col.names = T,sep = "\t",quote = F)
# 挑选前30个基因来可视化
topgene = sort(rowSums(neo.mtx>0,na.rm = T),decreasing = T)[1:30]
topgene = names(topgene)
topgene
# 添加注释信息
# 添加注释 neo counts 柱状图
neoantigens = sample.df$counts
names(neoantigens) = as.vector(sample.df$Tumor_Sample_Barcode)
write.table(sample.df,file = "sample_group_neo.txt",col.names = T,sep = "\t",quote = F)
column_top = HeatmapAnnotation(neoantigens = anno_barplot(neoantigens,
                                                          gp = gpar(fill = "#0E6CA3"),
                                                          border = F),
                               annotation_name_side = "left",
                               height = unit(3,"cm"))
# 肿瘤分组信息group
group = sample.df$group
column_bot = HeatmapAnnotation(group = group,annotation_name_side ="left")
# 基因的 neo 数
gene.sum = data.frame(gene = names(rowSums(neo.mtx>0,na.rm = T)),
                      counts = rowSums(neo.mtx>0,na.rm = T))
row_right = rowAnnotation("gene_neo" = anno_barplot(gene.sum[topgene,"counts"],
                                                    gp = gpar(fill = c("#0E6CA3")),
                                                    add_numbers = TRUE,
                                                    border = F),
                          width = unit(3,"cm"))
neo.ht2 = Heatmap(
  neo.mtx[topgene,sample.df$Tumor_Sample_Barcode], 
  name = "Best.MT.IC50.Score",   
  col = colorRampPalette(c( "#0E6CA3","white"))(100),  
  show_row_names = TRUE,   
  show_column_names = FALSE,  
  # row_title = "symbol",    
  column_title = "CRC and CRCBM_Paired Neo Best.MT.IC50.Score",  
  na_col = "white", 
  row_names_side = "left",
  use_raster = F,
  rect_gp = gpar(col = "grey"),
  cluster_rows = F,
  cluster_columns = F,
  top_annotation = column_top,
  bottom_annotation = column_bot,
  right_annotation = row_right
  
)
neo.ht2