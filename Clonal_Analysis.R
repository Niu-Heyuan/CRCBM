config = read.table("CRCBM_Paired.config",header = F)
colnames(config) = c("normal","tumor")
head(config)
# 临床信息
clincal.CRCBM = readxl::read_xlsx("CRCBM_Clinical.xlsx")
clincal.CRCBM = as.data.frame(clincal.CRCBM)
colnames(clincal.CRCBM) = gsub(" ","_",colnames(clincal.CRCBM))
head(clincal.CRCBM)
# 克隆结果
merge.clone  = data.frame()
for (i in config$tumor) {
  # i="B16_BM"
  clone.raw = read.table(paste0("./9.pyclone/",i,".h5.tsv"),header = T)
  # 根据 CCF = CP / 纯度 ，只能作为参考，大于1的赋值为1
  ccf = clone.raw$cellular_prevalence / clone.raw$tumour_content
  ccf[ccf>1] = 1
  clone.raw$ccf = ccf 
  # 区分 clone 和 subclone，每个样本最大的 CCF 为克隆，其余为 subclone
  clone.raw$clone = ifelse(clone.raw$ccf >= max(clone.raw$ccf),"clone","subclone")
  # 合并所有结果
  merge.clone = rbind(merge.clone,clone.raw)
}
# 添加 gene symbol
merge.clone$symbol = str_split(merge.clone$mutation_id,":",simplify = T)[,3]
# 因为后续以基因为单位，可能存在一个基因多个突变，即重复
table(duplicated(merge.clone[,c(2,10)]))
# 去重复:
# 对于 multi_hit 基因，取 cluster_assignment_prob 较大的突变位点
clone.rmdup = merge.clone %>% group_by(symbol, sample_id) %>%
  slice(which.max(cluster_assignment_prob)) %>%
  ungroup()
# 取关键列画热图
clone.df = clone.rmdup[,c("sample_id","ccf","symbol")]
# 长变宽
clone.mtx = spread(clone.df, key = "sample_id", value = "ccf")
clone.mtx = as.data.frame(clone.mtx)
rownames(clone.mtx) = clone.mtx[,1]
clone.mtx = clone.mtx[,-1]
write.table(clone.mtx,file = "gene_ccf.txt",col.names = T,sep = "\t",quote = F)
#从 maftools 获取 top30 基因
maf.df = fread("CRCBM_Paired.maf.gz")
maf = read.maf(maf.df)
topgene = getGeneSummary(maf)[1:30,1]
topgene = topgene$Hugo_Symbol
# 从 maf 获取样本排序,TMB
sample.sum = getSampleSummary(maf)
sample.sum = as.data.frame(sample.sum)
# 样本排序，分原发和转移，再加进行组内TMB排序
sample.od = as.vector(sample.sum$Tumor_Sample_Barcode)
sample.od = c(sample.od[grepl("P",sample.od)],sample.od[grepl("BM",sample.od)])
clone.topgene = clone.mtx[topgene,]
clone.topgene = clone.topgene[,sample.od]
write.table(clone.topgene,file = "topgene_ccf.txt",col.names = T,sep = "\t",quote = F)
# 添加注释 TMB 柱状图
tmb = sample.sum$total
names(tmb) = as.vector(sample.sum$Tumor_Sample_Barcode)
tmb = tmb[sample.od]
write.table(tmb,file = "sample_tmb.txt",col.names = T,sep = "\t",quote = F)
column_top = HeatmapAnnotation(TMB = anno_barplot(tmb,
                                                  gp = gpar(fill = "#0E6CA3"),
                                                  border = F),
                               annotation_name_side = "left",
                               height = unit(3,"cm"))
# 肿瘤原发和转移status、临床信息Sex、Age、Smoking、Location、Pathologic、Neoadjuvant therapy
status = ifelse(grepl("P",sample.od),"Primary","Brain")
sample.p = str_split(sample.od,"_",simplify = T)[,1]
rownames(clincal.CRCBM) = clincal.CRCBM$ID
column_bot = HeatmapAnnotation(status = status,
                               Sex = clincal.CRCBM[sample.p,"Sex"],
                               Age = clincal.CRCBM[sample.p,"Age_at__Primary_Dx"],
                               Smoking = clincal.CRCBM[sample.p,"History_of_smoking"],
                               Location = clincal.CRCBM[sample.p,"Location_of_the_primary_lesion"],
                               Pathologic = clincal.CRCBM[sample.p,"Pathologic_type"],
                               Neoadjuvant = clincal.CRCBM[sample.p,"Neoadjuvant_therapy"],
                               annotation_name_side ="left")
# 基因的 clone 和 subclone 数
clone.df2 = clone.rmdup[,c("sample_id","clone","symbol")]
clone.mtx2 = spread(clone.df2, key = "sample_id", value = "clone")
clone.mtx2 = as.data.frame(clone.mtx2)
rownames(clone.mtx2) = clone.mtx2[,1]
clone.mtx2 = clone.mtx2[,-1]
clone.mtx2$clone = rowSums(clone.mtx2=="clone",na.rm = T)
clone.mtx2$subclone = rowSums(clone.mtx2=="subclone",na.rm = T)
clone.mtx2[topgene,c("clone","subclone")]
clone.ht3 = Heatmap(
  clone.topgene[rownames(clone.p),], 
  name = "CCF",   
  col = colorRampPalette(c("white", "#0E6CA3"))(100),  
  show_row_names = TRUE,   
  show_column_names = FALSE,  
  #row_title = "symbol",    
  column_title = "CRCBM Paired Clone Event and CCF",  
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
clone.ht3


library(timescape)
library(plotly)
library(htmlwidgets)
library(webshot)
library(tidyr)
config = read.table("citup/p.config")
config = config[-c(1,3,13),]
for (patient in config[,1]) {
  # patient = "B28"
  #使用 timescape 进行可视化
  tree = read.table(paste0("citup/",patient,"_tree.txt"), header = T )
  tree_edges = tree
  clonefreq = read.table(paste0("citup/",patient,"_clonefreq.txt"), header = T )
  clonefreq = clonefreq[,c(2,1)]
  clonefreq$clone_id = rownames(clonefreq)
  clonal_prev = gather(clonefreq, key="timepoint", value = "clonal_prev", -clone_id)
  clonal_prev = clonal_prev[,c(2,1,3)]
  #clonal_prev = clonal_prev[order(clonal_prev$timepoint),]
  # targeted mutations
  # color
  clone_colours <- data.frame(clone_id = clonefreq$clone_id,
                              colour = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")[1:length(clonefreq$clone_id)]
  )
  p = timescape(clonal_prev = clonal_prev, 
                tree_edges = tree_edges, 
                height=260,
                genotype_position = "space",
                clone_colours = clone_colours)
  timescape(clonal_prev = clonal_prev, 
            tree_edges = tree_edges, 
            height=260,
            genotype_position = "space",
            clone_colours = clone_colours)
  saveWidget(p, paste0( patient ,"_timescape", ".html"))
  
}

